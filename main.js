(() => {
  const canvas = document.getElementById('sim');
  const ctx = canvas.getContext('2d');

  const ui = {
    grid: document.getElementById('grid'),
    rain: document.getElementById('rain'),
    rainLabel: document.getElementById('rainLabel'),
    ramp: document.getElementById('ramp'),
    strength: document.getElementById('strength'),
    strengthLabel: document.getElementById('strengthLabel'),
    drag: document.getElementById('drag'),
    dragLabel: document.getElementById('dragLabel'),
    visc: document.getElementById('visc'),
    viscLabel: document.getElementById('viscLabel'),
    reset: document.getElementById('reset'),
  };

  const g = 9.81;        // m/s^2
  const rho = 1000.0;    // kg/m^3
  let N = parseInt(ui.grid.value, 10);
  let W = canvas.width, H = canvas.height;
  let dx = 1.0 / N;

  // State arrays
  let h, u, v, strength, roof, drain;

  // Params
  let k = parseFloat(ui.drag.value);
  let nu = parseFloat(ui.visc.value);
  let rain_mmhr = parseFloat(ui.rain.value);
  let rainRamp = ui.ramp.checked;

  function alloc() {
    const size = N * N;
    h = new Float32Array(size);
    u = new Float32Array(size);
    v = new Float32Array(size);
    strength = new Float32Array(size);
    roof = new Uint8Array(size);
    drain = new Uint8Array(size);
    const base = parseFloat(ui.strength.value);
    for (let i = 0; i < size; i++) {
      h[i] = 0.0;
      u[i] = 0.0;
      v[i] = 0.0;
      strength[i] = base;
      roof[i] = 0;
      drain[i] = 0;
    }
  }

  function idx(x, y) { return y * N + x; }
  function clamp(a, lo, hi) { return Math.max(lo, Math.min(hi, a)); }

  function reset() {
    N = parseInt(ui.grid.value, 10);
    dx = 1.0 / N;
    alloc();
  }

  ui.reset.onclick = reset;
  ui.grid.onchange = reset;

  ui.rain.oninput = () => { rain_mmhr = parseFloat(ui.rain.value); ui.rainLabel.textContent = ui.rain.value; };
  ui.ramp.oninput = () => { rainRamp = ui.ramp.checked; };
  ui.strength.oninput = () => { ui.strengthLabel.textContent = ui.strength.value; };
  ui.drag.oninput = () => { k = parseFloat(ui.drag.value); ui.dragLabel.textContent = (+k).toFixed(2); };
  ui.visc.oninput = () => { nu = parseFloat(ui.visc.value); ui.viscLabel.textContent = (+nu).toFixed(4); };

  let mouseDown = false, erase = false;
  canvas.addEventListener('mousedown', (e) => {
    mouseDown = true;
    erase = (e.button === 2);
    paint(e);
  });
  canvas.addEventListener('mousemove', (e) => { if (mouseDown) paint(e); });
  canvas.addEventListener('mouseup', () => { mouseDown = false; });
  canvas.addEventListener('contextmenu', (e) => e.preventDefault());

  function paint(e) {
    const rect = canvas.getBoundingClientRect();
    const cx = ((e.clientX - rect.left) / rect.width) * N | 0;
    const cy = ((e.clientY - rect.top) / rect.height) * N | 0;
    const R = Math.max(2, (N/128)|0);
    for (let y = cy - R; y <= cy + R; y++) {
      for (let x = cx - R; x <= cx + R; x++) {
        if (x < 0 || y < 0 || x >= N || y >= N) continue;
        const d = Math.hypot(x - cx, y - cy);
        if (d <= R) {
          const id = idx(x,y);
          if (!erase) { roof[id] = 1; strength[id] = parseFloat(ui.strength.value); drain[id] = 0; }
          else { roof[id] = 0; drain[id] = 0; }
        }
      }
    }
  }

  function laplace(A, x, y) {
    const xm = x>0? x-1:x, xp = x<N-1? x+1:x;
    const ym = y>0? y-1:y, yp = y<N-1? y+1:y;
    return (A[idx(xp,y)] + A[idx(xm,y)] + A[idx(x,yp)] + A[idx(x,ym)] - 4*A[idx(x,y)]) / (dx*dx);
  }

  function step(dt, t) {
    // Rain in m/s from mm/hr
    let rain = rain_mmhr / 1000.0 / 3600.0;
    if (rainRamp) {
      // Ramp up over 2 minutes to 60 mm/hr peak
      const phase = Math.min(1, t / 120.0);
      const peak = 60.0 / 1000.0 / 3600.0;
      rain = (1 - phase) * rain + phase * peak;
    }

    // Velocity update (semi-implicit drag/visc)
    const u_new = new Float32Array(u.length);
    const v_new = new Float32Array(v.length);
    let umax = 0, vmax = 0, hmax = 0;

    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const id = idx(x,y);
        // Height gradients (central)
        const hx = (h[idx(Math.min(x+1,N-1),y)] - h[idx(Math.max(x-1,0),y)]) / (2*dx);
        const hy = (h[idx(x,Math.min(y+1,N-1))] - h[idx(x,Math.max(y-1,0))]) / (2*dx);

        const uL = u[id] + dt * ( -g*hx - k*u[id] + nu * laplace(u, x, y) );
        const vL = v[id] + dt * ( -g*hy - k*v[id] + nu * laplace(v, x, y) );

        // No-slip at outer boundary
        const isEdge = (x===0||y===0||x===N-1||y===N-1);
        u_new[id] = isEdge ? 0 : uL;
        v_new[id] = isEdge ? 0 : vL;

        umax = Math.max(umax, Math.abs(u_new[id]));
        vmax = Math.max(vmax, Math.abs(v_new[id]));
        hmax = Math.max(hmax, h[id]);
      }
    }

    // Height update (divergence + advection + weak dissipation)
    const h_new = new Float32Array(h.length);
    const He = 0.0, eps = 0.001;
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const id = idx(x,y);
        const ux = (u_new[idx(Math.min(x+1,N-1),y)] - u_new[idx(Math.max(x-1,0),y)]) / (2*dx);
        const vy = (v_new[idx(x,Math.min(y+1,N-1))] - v_new[idx(x,Math.max(y-1,0))]) / (2*dx);
        const hx = (h[idx(Math.min(x+1,N-1),y)] - h[idx(Math.max(x-1,0),y)]) / (2*dx);
        const hy = (h[idx(x,Math.min(y+1,N-1))] - h[idx(x,Math.max(y-1,0))]) / (2*dx);

        let dh = - (ux + vy) * (h[id] + He) - (hx * u_new[id] + hy * v_new[id]) - eps * h[id];
        dh += rain;

        // Roof drains if failed
        const p = rho * g * h[id];
        if (roof[id]) {
          if (p > strength[id]) { roof[id] = 0; drain[id] = 1; }
        }
        if (drain[id]) {
          dh -= Math.min(h[id]/dt, 0.5 * h[id] / dt); // sink quickly but not instant
        }

        // Keep nonnegative height
        h_new[id] = Math.max(0, h[id] + dt * dh);
      }
    }

    h = h_new;
    u = u_new;
    v = v_new;

    // CFL recompute suggestion
    const wavespeed = Math.sqrt(g * Math.max(1e-6, hmax));
    const cfl = 0.4;
    const dtCfl = cfl * dx / (wavespeed + Math.max(umax, vmax) + 1e-6);
    return dtCfl;
  }

  function render() {
    const img = ctx.createImageData(W, H);
    const data = img.data;
    const sx = W / N, sy = H / N;
    // Simple block rendering
    for (let y = 0; y < H; y++) {
      const gy = Math.min(N-1, (y / sy) | 0);
      for (let x = 0; x < W; x++) {
        const gx = Math.min(N-1, (x / sx) | 0);
        const id = idx(gx, gy);
        const hh = h[id];
        const rr = roof[id] ? 255 : (drain[id] ? 255 : 0);
        const gg = roof[id] ? 200 : (drain[id] ? 80 : 0);
        // Water color intensity by depth
        const b = Math.max(0, Math.min(255, (Math.sqrt(hh) * 300)|0));
        const off = (y * W + x) * 4;
        data[off+0] = Math.max(rr, 0);
        data[off+1] = Math.max(gg, 0);
        data[off+2] = b;
        data[off+3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }

  let t = 0;
  let dt = 0.01;
  let last = performance.now();

  function loop(now) {
    const elapsed = Math.min(0.05, (now - last) / 1000);
    last = now;
    let acc = elapsed;
    while (acc > 1e-6) {
      const stept = Math.min(dt, acc);
      const dtSuggest = step(stept, t);
      t += stept;
      acc -= stept;
      // Smoothly adapt dt to CFL suggestion
      dt = 0.9 * dt + 0.1 * dtSuggest;
      dt = clamp(dt, 0.0015, 0.02);
    }
    render();
    requestAnimationFrame(loop);
  }

  reset();
  requestAnimationFrame(loop);
})();
