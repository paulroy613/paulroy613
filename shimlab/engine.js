// ShimLab valving engine — independent model built from textbook physics.
// Shim stack: axisymmetric Kirchhoff plate, clamped at the clamp radius, free at
// the outer edge, with bending rigidity D(r) = sum of the rigidities of every shim
// that reaches radius r (shims bend together, no separation or friction).
// Flow: turbulent orifice equation, Q = Cd * A * sqrt(2 dP / rho).
// All internal units are SI (m, Pa, N, m^3/s).
(function (root) {
  'use strict';

  function stackModel(cfg) {
    const E = cfg.E, nu = cfg.nu;
    const b = cfg.clampD / 2000;
    const shims = cfg.shims
      .filter(s => s.qty > 0 && s.od > cfg.clampD && s.t > 0)
      .map(s => ({ r: s.od / 2000, t: s.t / 1000, qty: s.qty }));
    if (!shims.length) return { ok: false, reason: 'No shim is larger than the clamp diameter.' };
    const a = Math.max(...shims.map(s => s.r));
    let rs = cfg.seatOD / 2000;
    const warnings = [];
    if (rs > a) { warnings.push('Seat is wider than the largest shim, so the stack cannot seal.'); rs = a; }
    if (rs <= b) return { ok: false, reason: 'Seat diameter must be larger than the clamp diameter.' };
    let rp = Math.min(Math.max(cfg.portID / 2000, b), rs);
    const k = cfg.portArc;
    const Dfac = E / (12 * (1 - nu * nu));
    const Dat = r => { let s = 0; for (const sh of shims) if (sh.r > r) s += sh.qty * sh.t ** 3; return Dfac * s; };
    // Load per unit circumference at radius r, for unit pressure (1 Pa).
    const P = r => r >= rs ? 0 : k * Math.PI * (rs * rs - Math.max(r, rp) ** 2);
    const Qs = r => P(r) / (2 * Math.PI * r);

    const N = 600, h = (a - b) / N;
    const r = i => b + i * h;
    const Dm = [], rm = [];
    for (let i = 0; i < N; i++) { rm[i] = b + (i + 0.5) * h; Dm[i] = Dat(rm[i]); }
    // Tridiagonal system for slope phi_1..phi_N.
    const A = new Float64Array(N + 1), B = new Float64Array(N + 1), C = new Float64Array(N + 1), R = new Float64Array(N + 1);
    for (let i = 1; i < N; i++) {
      const ri = r(i), Dn = 0.5 * (Dm[i - 1] + Dm[i]);
      A[i] = Dm[i - 1] * (rm[i - 1] / h - nu / 2) / h + Dn * nu / (2 * h);
      B[i] = (Dm[i] * (-rm[i] / h + nu / 2) - Dm[i - 1] * (rm[i - 1] / h + nu / 2)) / h - Dn / ri;
      C[i] = Dm[i] * (rm[i] / h + nu / 2) / h - Dn * nu / (2 * h);
      R[i] = -ri * Qs(ri);
    }
    {
      const i = N, ri = r(N), Dn = Dm[N - 1];
      A[i] = -2 / h * Dm[N - 1] * (-rm[N - 1] / h + nu / 2) + Dn * nu / h;
      B[i] = -2 / h * Dm[N - 1] * (rm[N - 1] / h + nu / 2) - Dn / ri - Dn * nu / h;
      C[i] = 0; R[i] = -ri * Qs(ri);
    }
    // Thomas algorithm (phi_0 = 0).
    const cp = new Float64Array(N + 1), dp = new Float64Array(N + 1), phi = new Float64Array(N + 1);
    cp[1] = C[1] / B[1]; dp[1] = R[1] / B[1];
    for (let i = 2; i <= N; i++) {
      const m = B[i] - A[i] * cp[i - 1];
      cp[i] = C[i] / m; dp[i] = (R[i] - A[i] * dp[i - 1]) / m;
    }
    phi[N] = dp[N];
    for (let i = N - 1; i >= 1; i--) phi[i] = dp[i] - cp[i] * phi[i + 1];
    // Deflection profile (per Pa), positive away from seat.
    const w = new Float64Array(N + 1);
    for (let i = 1; i <= N; i++) w[i] = w[i - 1] + 0.5 * (phi[i] + phi[i - 1]) * h;
    const sgn = w[N] < 0 ? -1 : 1;
    for (let i = 0; i <= N; i++) w[i] *= sgn;
    const interp = rr => { const x = (rr - b) / h, i = Math.min(N - 1, Math.max(0, Math.floor(x))); return w[i] + (w[i + 1] - w[i]) * (x - i); };
    const compliance = interp(rs); // m of seat lift per Pa
    // Peak bending stress per Pa, in the thickest shim present at each radius.
    let stressPerPa = 0, stressR = b;
    for (let i = 0; i < N; i++) {
      const rr = rm[i];
      const dphi = (phi[i + 1] - phi[i]) / h, ph = 0.5 * (phi[i + 1] + phi[i]);
      let tmax = 0; for (const sh of shims) if (sh.r > rr) tmax = Math.max(tmax, sh.t);
      const s = E * tmax / (2 * (1 - nu * nu)) * Math.abs(dphi + nu * ph / rr);
      if (s > stressPerPa) { stressPerPa = s; stressR = rr; }
    }
    const profile = [];
    for (let i = 0; i <= N; i += 10) profile.push([r(i), w[i]]);
    return { ok: true, warnings, a, b, rs, rp, compliance, stressPerPa, stressR, profile, shims };
  }

  function flowArea(mode, boreD, rodD) {
    const bore = Math.PI / 4 * (boreD / 1000) ** 2, rod = Math.PI / 4 * (rodD / 1000) ** 2;
    if (mode === 'bore') return bore;
    if (mode === 'rod') return rod;
    return bore - rod;
  }

  function simulate(cfg) {
    const st = stackModel(cfg);
    if (!st.ok) return st;
    const rho = cfg.rho;
    const Ad = flowArea(cfg.areaMode, cfg.boreD, cfg.rodD);
    const Ab = Math.PI / 4 * (cfg.bleedD / 1000) ** 2;
    const Ap = cfg.portArea * 1e-6;
    const pre = cfg.preload / 1000;
    const c = st.compliance;
    const crack = pre / c;
    const lift = dP => Math.max(0, c * dP - pre);
    const flowAt = dP => {
      const vel = Math.sqrt(2 * dP / rho);
      const Ac = 2 * Math.PI * st.rs * lift(dP);
      const Aeff = Ac > 0 ? Ap * Ac / Math.sqrt(Ap * Ap + Ac * Ac) : 0;
      return { bleed: cfg.cdBleed * Ab * vel, stack: cfg.cdStack * Aeff * vel };
    };
    const solve = Qt => {
      if (Qt <= 0) return 0;
      let lo = 0, hi = 1e5;
      while ((() => { const f = flowAt(hi); return f.bleed + f.stack; })() < Qt) { hi *= 2; if (hi > 1e11) return hi; }
      for (let i = 0; i < 80; i++) {
        const mid = 0.5 * (lo + hi), f = flowAt(mid);
        if (f.bleed + f.stack < Qt) lo = mid; else hi = mid;
      }
      return 0.5 * (lo + hi);
    };
    const n = 80, pts = [];
    for (let j = 0; j <= n; j++) {
      const v = cfg.vmax * j / n, Qt = Ad * v, dP = solve(Qt), f = flowAt(dP);
      const tot = f.bleed + f.stack;
      pts.push({ v, dP, F: dP * Ad, lift: lift(dP), bleedFrac: tot > 0 ? f.bleed / tot : 1, stress: st.stressPerPa * dP });
    }
    const at = v => { const Qt = Ad * v, dP = solve(Qt), f = flowAt(dP), tot = f.bleed + f.stack; return { v, dP, F: dP * Ad, lift: lift(dP), bleedFrac: tot > 0 ? f.bleed / tot : 1, stress: st.stressPerPa * dP }; };
    // Velocity where the stack starts to lift (bleed alone reaches crack pressure).
    const vCrack = Ab > 0 ? cfg.cdBleed * Ab * Math.sqrt(2 * crack / rho) / Ad : 0;
    return { ok: true, st, pts, at, Ad, crack, vCrack, warnings: st.warnings };
  }

  root.ShimLab = { stackModel, simulate, flowArea };
  if (typeof module !== 'undefined') module.exports = root.ShimLab;
})(typeof window !== 'undefined' ? window : globalThis);
