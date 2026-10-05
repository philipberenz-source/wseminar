/* fdtd.js — JavaScript-Port von fdtd_core.py (FDTD2D) + Szenenmodell aus
 * fdtd_studio.py. Die Update-Gleichungen, die CPML-Koeffizienten und die
 * Quellen entsprechen dem Python-Kern Zeile fuer Zeile; Arrays liegen
 * zeilenweise mit Index [i, j] -> i * ny + j (wie numpy, Form (nx, ny)).
 *
 * Wird als klassisches Skript geladen (kein ES-Modul), damit die Seiten auch
 * per Doppelklick ueber file:// funktionieren.
 */
(function (global) {
  "use strict";

  const EPS0 = 8.8541878128e-12;
  const MU0 = 4e-7 * Math.PI;
  const C0 = 1.0 / Math.sqrt(EPS0 * MU0);

  // ------------------------------------------------------------- Kern
  class FDTD2D {
    constructor(nx, ny, d, opts = {}) {
      const { courant = 0.99, eps_r = null, sigma = null,
              npml = 10, pml_m = 3, r0 = 1e-6 } = opts;
      this.nx = nx; this.ny = ny; this.d = d;
      this.courant = courant;
      this.dt = courant * d / (C0 * Math.SQRT2);

      const N = nx * ny;
      this.eps_r = eps_r ? Float64Array.from(eps_r) : new Float64Array(N).fill(1);
      this.sigma = sigma ? Float64Array.from(sigma) : new Float64Array(N);
      this.Ez = new Float64Array(N);
      this.Hx = new Float64Array(nx * (ny - 1));
      this.Hy = new Float64Array((nx - 1) * ny);
      this.Jz = new Float64Array(N);
      this.pec = new Uint8Array(N);
      this.pecIdx = new Int32Array(0);

      this._updateMaterialCoeffs();

      this.npml = npml; this.pml_m = pml_m; this.r0 = r0;
      this._setupPml();

      this.softSources = [];      // {i, j, fn}         -> Ez += fn(t)
      this.currentSources = [];   // {i, j, fn, amp}    -> Jz += amp*fn(t)
      this.lineSources = [];      // {idx, g, delay, fn, mode}
      this.probes = [];

      this.timeStep = 0;
      this.t = 0;
    }

    _updateMaterialCoeffs() {
      const N = this.nx * this.ny, dt = this.dt, d = this.d;
      this.Ca = new Float64Array(N);
      this.Cb = new Float64Array(N);
      for (let k = 0; k < N; k++) {
        const loss = this.sigma[k] * dt / (2 * EPS0 * this.eps_r[k]);
        this.Ca[k] = (1 - loss) / (1 + loss);
        this.Cb[k] = (dt / (EPS0 * this.eps_r[k] * d)) / (1 + loss);
      }
    }

    setMaterialBox(i0, i1, j0, j1, eps_r = 1, sigma = 0) {
      for (let i = i0; i < i1; i++) for (let j = j0; j < j1; j++) {
        const k = i * this.ny + j; this.eps_r[k] = eps_r; this.sigma[k] = sigma;
      }
      this._updateMaterialCoeffs();
    }

    setPecMask(mask) {
      const idx = [];
      for (let k = 0; k < mask.length; k++) if (mask[k]) { this.pec[k] = 1; idx.push(k); this.Ez[k] = 0; }
      this.pecIdx = Int32Array.from(idx);
    }

    _pmlCoeffs(pos, n) {
      const b = new Float64Array(pos.length), C = new Float64Array(pos.length);
      for (let k = 0; k < pos.length; k++) {
        if (this.npml === 0) { b[k] = 1; C[k] = 0; continue; }
        const dl = this.npml - pos[k];
        const dr = pos[k] - (n - 1 - this.npml);
        const tiefe = Math.max(dl, dr, 0);
        const rho = Math.min(tiefe / this.npml, 1);
        const sig = this.sigma_max * Math.pow(rho, this.pml_m);
        b[k] = Math.exp(-sig * this.dt / EPS0);
        C[k] = b[k] - 1;
      }
      return [b, C];
    }

    _setupPml() {
      if (this.npml > 0) {
        const eta0 = Math.sqrt(MU0 / EPS0);
        const D = this.npml * this.d;
        this.sigma_max = -(this.pml_m + 1) * Math.log(this.r0) / (2 * eta0 * D);
      } else this.sigma_max = 0;
      const { nx, ny } = this;
      const ar = (n, off) => Float64Array.from({ length: n }, (_, k) => k + off);
      [this.bx_i, this.Cx_i] = this._pmlCoeffs(ar(nx, 0), nx);
      [this.by_j, this.Cy_j] = this._pmlCoeffs(ar(ny, 0), ny);
      [this.bx_h, this.Cx_h] = this._pmlCoeffs(ar(nx - 1, 0.5), nx);
      [this.by_h, this.Cy_h] = this._pmlCoeffs(ar(ny - 1, 0.5), ny);
      this.psi_hx = new Float64Array(this.Hx.length);
      this.psi_hy = new Float64Array(this.Hy.length);
      this.psi_ez_x = new Float64Array((nx - 2) * (ny - 2));
      this.psi_ez_y = new Float64Array((nx - 2) * (ny - 2));
    }

    addSoftSource(i, j, fn) { this.softSources.push({ i, j, fn }); }
    addCurrentSource(i, j, fn, amp = 1) { this.currentSources.push({ i, j, fn, amp }); }
    addLineSource(src) { this.lineSources.push(src); }
    addProbe(i, j) { this.probes.push({ i, j, t: [], Ez: [] }); return this.probes.length - 1; }

    _updateH() {
      const { nx, ny, d, Ez, Hx, Hy, psi_hx, psi_hy, by_h, Cy_h, bx_h, Cx_h } = this;
      const f = this.dt / MU0, nyh = ny - 1;
      for (let i = 0; i < nx; i++) {
        const r = i * ny, rh = i * nyh;
        for (let j = 0; j < nyh; j++) {
          const g = (Ez[r + j + 1] - Ez[r + j]) / d;
          const k = rh + j;
          const p = by_h[j] * psi_hx[k] + Cy_h[j] * g;
          psi_hx[k] = p;
          Hx[k] -= f * (g + p);
        }
      }
      for (let i = 0; i < nx - 1; i++) {
        const r = i * ny, bx = bx_h[i], cx = Cx_h[i];
        for (let j = 0; j < ny; j++) {
          const k = r + j;
          const g = (Ez[k + ny] - Ez[k]) / d;
          const p = bx * psi_hy[k] + cx * g;
          psi_hy[k] = p;
          Hy[k] += f * (g + p);
        }
      }
    }

    _updateE() {
      const { nx, ny, d, Ez, Hx, Hy, Jz, Ca, Cb, psi_ez_x, psi_ez_y, bx_i, Cx_i, by_j, Cy_j } = this;
      const nyh = ny - 1, nyi = ny - 2;
      for (let i = 1; i < nx - 1; i++) {
        const bx = bx_i[i], cx = Cx_i[i];
        for (let j = 1; j < ny - 1; j++) {
          const k = i * ny + j;
          const dHy = Hy[k] - Hy[k - ny];
          const dHx = Hx[i * nyh + j] - Hx[i * nyh + j - 1];
          const q = (i - 1) * nyi + (j - 1);
          const px = bx * psi_ez_x[q] + cx * (dHy / d);
          const py = by_j[j] * psi_ez_y[q] + Cy_j[j] * (dHx / d);
          psi_ez_x[q] = px; psi_ez_y[q] = py;
          Ez[k] = Ca[k] * Ez[k] + Cb[k] * (dHy - dHx + d * (px - py) - Jz[k] * d);
        }
      }
      const pi = this.pecIdx;
      for (let m = 0; m < pi.length; m++) Ez[pi[m]] = 0;
    }

    _injectCurrent(t) {
      this.Jz.fill(0);
      for (const s of this.currentSources) this.Jz[s.i * this.ny + s.j] += s.amp * s.fn(t);
      for (const L of this.lineSources) if (L.mode !== "weich") L.apply(this.Jz, t);
    }

    _injectSoft(t) {
      for (const s of this.softSources) this.Ez[s.i * this.ny + s.j] += s.fn(t);
      for (const L of this.lineSources) if (L.mode === "weich") L.apply(this.Ez, t);
    }

    step() {
      const tHalb = this.t + 0.5 * this.dt, tNeu = this.t + this.dt;
      this._updateH();
      this._injectCurrent(tHalb);
      this._updateE();
      this._injectSoft(tNeu);
      this.timeStep += 1;
      this.t = this.timeStep * this.dt;
      for (const p of this.probes) { p.t.push(this.t); p.Ez.push(this.Ez[p.i * this.ny + p.j]); }
    }

    totalEnergy() {
      let ue = 0, uh = 0;
      for (let k = 0; k < this.Ez.length; k++) ue += 0.5 * EPS0 * this.eps_r[k] * this.Ez[k] ** 2;
      for (let k = 0; k < this.Hx.length; k++) uh += 0.5 * MU0 * this.Hx[k] ** 2;
      for (let k = 0; k < this.Hy.length; k++) uh += 0.5 * MU0 * this.Hy[k] ** 2;
      return (ue + uh) * this.d * this.d;
    }

    /* Verschraenkte Energie wie im Studio:
       W^n = d^2 sum( eps0 eps_r (Ez^n)^2/2 + mu0 H^(n-1/2).H^(n+1/2)/2 ) */
    interleavedEnergy(ezAlt, hxAlt, hyAlt) {
      let we = 0, wh = 0;
      for (let k = 0; k < ezAlt.length; k++) we += this.eps_r[k] * ezAlt[k] * ezAlt[k];
      for (let k = 0; k < hxAlt.length; k++) wh += hxAlt[k] * this.Hx[k];
      for (let k = 0; k < hyAlt.length; k++) wh += hyAlt[k] * this.Hy[k];
      return (0.5 * EPS0 * we + 0.5 * MU0 * wh) * this.d * this.d;
    }
  }

  // ------------------------------------------------------ Wellenformen
  const WELLENFORMEN = { sinus: "Sinus (Dauerstrich)", ricker: "Ricker-Puls", gauss: "Gauss-Puls" };
  const ALTNAMEN = { rampe: "sinus", moduliert: "ricker", nadel: "ricker" };

  function machWellenform(art, f0, amp = 1, zyklen = 3, verzoegerung = 0) {
    art = ALTNAMEN[art] || art;
    f0 = Math.max(+f0, 1);
    const T = 1 / f0, tv = +verzoegerung * T;
    if (art === "ricker") {
      const a = (Math.PI * f0) ** 2, t0 = 1.2 * T + tv;
      return t => { const x = a * (t - t0) ** 2; return amp * (1 - 2 * x) * Math.exp(-x); };
    }
    if (art === "gauss") {
      const a = (Math.PI * f0) ** 2, t0 = 1.2 * T + tv;
      return t => amp * Math.exp(-a * (t - t0) ** 2);
    }
    if (art === "sinus") {
      const tr = Math.max(+zyklen, 0) * T;
      return t => {
        if (t < tv) return 0;
        const s = Math.sin(2 * Math.PI * f0 * (t - tv));
        if (tr <= 0) return amp * s;
        const x = Math.min(Math.max((t - tv) / tr, 0), 1);
        return amp * 0.5 * (1 - Math.cos(Math.PI * x)) * s;
      };
    }
    throw new Error("unbekannte Wellenform: " + art);
  }

  // ------------------------------------------------------- Szenenmodell
  // Feldnamen identisch zu den Python-Dataclasses -> JSON-kompatibel.
  const Objekt = o => Object.assign({ form: "rechteck", art: "dielektrikum", i0: 0, j0: 0, i1: 10, j1: 10,
                                      eps_r: 2.25, sigma: 0, aktiv: true }, o);
  const Quelle = q => {
    const r = Object.assign({ typ: "ebene", art: "sinus", i: 30, j: 100, f0: 15e9, amp: 1, zyklen: 4,
                              theta: 0, modus: "strom", taper: 0, arm: true, richtung: false,
                              verzoegerung: 0, aktiv: true }, q);
    r.art = ALTNAMEN[r.art] || r.art; return r;
  };
  const Sonde = s => Object.assign({ i: 100, j: 100, name: "" }, s);
  const Szene = s => {
    const r = Object.assign({ nx: 300, ny: 240, d: 1e-3, courant: 0.99, npml: 12, pml_m: 3, r0: 1e-6,
                              schirm: 0, objekte: [], quellen: [], sonden: [] }, s);
    r.objekte = r.objekte.map(Objekt); r.quellen = r.quellen.map(Quelle); r.sonden = r.sonden.map(Sonde);
    return r;
  };

  function objBox(o, nx, ny) {
    let [i0, i1] = [Math.trunc(o.i0), Math.trunc(o.i1)].sort((a, b) => a - b);
    let [j0, j1] = [Math.trunc(o.j0), Math.trunc(o.j1)].sort((a, b) => a - b);
    i0 = Math.max(i0, 0); j0 = Math.max(j0, 0);
    return [i0, j0, Math.min(Math.max(i1, i0 + 1), nx), Math.min(Math.max(j1, j0 + 1), ny)];
  }

  function objEnthaelt(o, i, j, nx, ny) {
    const [i0, j0, i1, j1] = objBox(o, nx, ny);
    if (o.form === "ellipse") {
      const ci = 0.5 * (i0 + i1 - 1), cj = 0.5 * (j0 + j1 - 1);
      const ri = Math.max(0.5 * (i1 - i0), 0.5), rj = Math.max(0.5 * (j1 - j0), 0.5);
      return ((i - ci) / ri) ** 2 + ((j - cj) / rj) ** 2 <= 1;
    }
    return i >= i0 && i < i1 && j >= j0 && j < j1;
  }

  function karten(sz) {
    const { nx, ny } = sz, N = nx * ny;
    const eps = new Float64Array(N).fill(1), sig = new Float64Array(N), pec = new Uint8Array(N);
    for (const o of sz.objekte) {
      if (!o.aktiv) continue;
      const [i0, j0, i1, j1] = objBox(o, nx, ny);
      for (let i = i0; i < i1; i++) for (let j = j0; j < j1; j++) {
        if (!objEnthaelt(o, i, j, nx, ny)) continue;
        const k = i * ny + j;
        if (o.art === "pec") pec[k] = 1;
        else { eps[k] = Math.max(o.eps_r, 1e-3); sig[k] = Math.max(o.sigma, 0); }
      }
    }
    return { eps, sig, pec };
  }

  function linienPunkte(nx, ny, iSrc, thetaGrad, arm, richtung, f0, d, taper = 0) {
    const th = thetaGrad * Math.PI / 180, kx = Math.cos(th), ky = Math.sin(th);
    const versatz = [0];
    if (richtung) {
      const lam = C0 / Math.max(f0, 1);
      versatz.push(Math.max(1, Math.round(lam / (4 * Math.max(Math.abs(kx), 0.25) * d))));
    }
    const seen = new Map();
    const push = (i, j, g) => { const key = i * ny + j; if (!seen.has(key)) seen.set(key, [i, j, g]); };
    for (const v of versatz) {
      const i0 = Math.min(Math.max(Math.trunc(iSrc + v), 0), nx - 1);
      for (let j = 0; j < ny; j++) {
        let g = 1;
        if (taper > 0) { const r = Math.min(Math.max(Math.min(j, ny - 1 - j) / taper, 0), 1); g = Math.sin(0.5 * Math.PI * r) ** 2; }
        push(i0, j, g);
      }
      if (arm && Math.abs(ky) > 1e-9) {
        const jk = ky > 0 ? 0 : ny - 1;
        for (let i = i0; i < nx; i++) push(i, jk, 1);
      }
    }
    const pts = [...seen.values()];
    const idx = new Int32Array(pts.length), g = new Float64Array(pts.length), delay = new Float64Array(pts.length);
    let mn = Infinity;
    pts.forEach(([i, j, w], k) => { idx[k] = i * ny + j; g[k] = w; delay[k] = (i * kx + j * ky) * d / C0; mn = Math.min(mn, delay[k]); });
    for (let k = 0; k < delay.length; k++) delay[k] -= mn;
    return { idx, g, delay, pts };
  }

  function registriereQuelle(sim, q) {
    const puls = machWellenform(q.art, q.f0, q.amp, q.zyklen, q.verzoegerung);
    if (q.typ === "punkt") {
      const i = Math.min(Math.max(Math.trunc(q.i), 0), sim.nx - 1);
      const j = Math.min(Math.max(Math.trunc(q.j), 0), sim.ny - 1);
      if (q.modus === "weich") sim.addSoftSource(i, j, puls); else sim.addCurrentSource(i, j, puls);
      return 1;
    }
    const { idx, g, delay } = linienPunkte(sim.nx, sim.ny, q.i, q.theta, q.arm, q.richtung, q.f0, sim.d, q.taper);
    sim.addLineSource({
      mode: q.modus, idx, g, delay,
      apply(target, t) { for (let k = 0; k < idx.length; k++) target[idx[k]] += g[k] * puls(t - delay[k]); }
    });
    return idx.length;
  }

  function baue(sz) {
    const { eps, sig, pec } = karten(sz);
    const sim = new FDTD2D(sz.nx, sz.ny, sz.d, { courant: sz.courant, eps_r: eps, sigma: sig,
                                                 npml: sz.npml, pml_m: sz.pml_m, r0: sz.r0 });
    sim.setPecMask(pec);
    for (const q of sz.quellen) if (q.aktiv) registriereQuelle(sim, q);
    for (const s of sz.sonden)
      sim.addProbe(Math.min(Math.max(Math.trunc(s.i), 0), sz.nx - 1), Math.min(Math.max(Math.trunc(s.j), 0), sz.ny - 1));
    return { sim, eps, sig, pec };
  }

  function kennzahlen(sz) {
    const dt = sz.courant * sz.d / (C0 * Math.SQRT2);
    const f0 = sz.quellen.length ? sz.quellen[0].f0 : 0;
    return { dt, breite: sz.nx * sz.d, hoehe: sz.ny * sz.d, n_lambda: f0 > 0 ? C0 / (f0 * sz.d) : Infinity, f0 };
  }

  // ------------------------------------------------------------ Presets
  const wand = (i0, i1, ny, schlitze) => {
    const k = [0]; for (const [a, b] of schlitze) k.push(a, b); k.push(ny);
    const out = [];
    for (let m = 0; m < k.length; m += 2) if (k[m + 1] > k[m]) out.push(Objekt({ form: "rechteck", art: "pec", i0, j0: k[m], i1, j1: k[m + 1] }));
    return out;
  };

  const PRESETS = {
    "Freiraum: ebene Welle": () => Szene({ nx: 320, ny: 260,
      quellen: [{ typ: "ebene", art: "ricker", i: 30, f0: 15e9 }], sonden: [{ i: 120, j: 130 }, { i: 220, j: 130 }] }),
    "Schräge ebene Welle (30°)": () => Szene({ nx: 320, ny: 320,
      quellen: [{ typ: "ebene", art: "sinus", i: 25, f0: 15e9, theta: 30, arm: true, richtung: true, zyklen: 4 }], sonden: [{ i: 200, j: 200 }] }),
    "Doppelspalt (Beugung)": () => Szene({ nx: 340, ny: 300, schirm: 280,
      quellen: [{ typ: "ebene", art: "sinus", i: 30, f0: 15e9, zyklen: 4 }],
      objekte: wand(130, 134, 300, [[115, 125], [175, 185]]), sonden: [{ i: 280, j: 150 }, { i: 280, j: 202 }] }),
    "Einzelspalt": () => Szene({ nx: 340, ny: 300, schirm: 280,
      quellen: [{ typ: "ebene", art: "sinus", i: 30, f0: 15e9, zyklen: 4 }],
      objekte: wand(130, 134, 300, [[130, 170]]), sonden: [{ i: 280, j: 150 }] }),
    "Beugungsgitter (5 Spalte)": () => {
      const s = []; for (let j = 40; j <= 280; j += 60) s.push([j - 5, j + 5]);
      return Szene({ nx: 360, ny: 320, schirm: 300, quellen: [{ typ: "ebene", art: "sinus", i: 30, f0: 15e9, zyklen: 4 }],
        objekte: wand(140, 144, 320, s), sonden: [{ i: 300, j: 160 }] });
    },
    "Brechung an Grenzfläche": () => Szene({ nx: 360, ny: 300,
      quellen: [{ typ: "ebene", art: "sinus", i: 25, f0: 15e9, theta: 30, arm: true, richtung: true, zyklen: 4 }],
      objekte: [{ form: "rechteck", art: "dielektrikum", i0: 190, j0: 0, i1: 360, j1: 300, eps_r: 2.25 }],
      sonden: [{ i: 120, j: 150 }, { i: 280, j: 150 }] }),
    "Totalreflexion (50°)": () => Szene({ nx: 360, ny: 300,
      quellen: [{ typ: "ebene", art: "sinus", i: 25, f0: 15e9, theta: 50, arm: true, richtung: true, zyklen: 4 }],
      objekte: [{ form: "rechteck", art: "dielektrikum", i0: 0, j0: 0, i1: 190, j1: 300, eps_r: 2.25 }],
      sonden: [{ i: 210, j: 150 }] }),
    "Dielektrische Linse": () => Szene({ nx: 360, ny: 280,
      quellen: [{ typ: "ebene", art: "sinus", i: 25, f0: 15e9, zyklen: 4 }],
      objekte: [{ form: "ellipse", art: "dielektrikum", i0: 140, j0: 70, i1: 200, j1: 210, eps_r: 2.25 }],
      sonden: [{ i: 250, j: 140 }, { i: 300, j: 140 }] }),
    "Streuung am PEC-Zylinder": () => Szene({ nx: 340, ny: 300,
      quellen: [{ typ: "ebene", art: "ricker", i: 30, f0: 15e9 }],
      objekte: [{ form: "ellipse", art: "pec", i0: 150, j0: 120, i1: 200, j1: 180 }],
      sonden: [{ i: 280, j: 150 }, { i: 80, j: 150 }] }),
    "Parallelplatten-Hohlleiter": () => Szene({ nx: 360, ny: 260,
      quellen: [{ typ: "punkt", art: "sinus", i: 60, j: 130, f0: 15e9, zyklen: 3 }],
      objekte: [{ form: "rechteck", art: "pec", i0: 40, j0: 104, i1: 360, j1: 110 },
                { form: "rechteck", art: "pec", i0: 40, j0: 150, i1: 360, j1: 156 }],
      sonden: [{ i: 200, j: 130 }, { i: 320, j: 130 }] }),
    "Bragg-Spiegel (Schichtstapel)": () => {
      const o = []; let i = 140;
      for (let n = 0; n < 8; n++) { o.push({ form: "rechteck", art: "dielektrikum", i0: i, j0: 0, i1: i + 5, j1: 200, eps_r: 4 }); i += 15; }
      return Szene({ nx: 380, ny: 200, quellen: [{ typ: "ebene", art: "ricker", i: 30, f0: 7.5e9 }], objekte: o,
        sonden: [{ i: 90, j: 100 }, { i: 330, j: 100 }] });
    },
    "Punktquelle im Freiraum": () => Szene({ nx: 280, ny: 280,
      quellen: [{ typ: "punkt", art: "ricker", i: 140, j: 140, f0: 15e9 }], sonden: [{ i: 200, j: 140 }, { i: 140, j: 200 }] }),
    "PEC-Resonator (ohne PML)": () => Szene({ nx: 141, ny: 141, npml: 0,
      quellen: [{ typ: "punkt", art: "ricker", i: 70, j: 70, f0: 15e9 }], sonden: [{ i: 40, j: 95 }] }),
  };

  // ------------------------------------------------------------ Farben
  const hex = h => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
  const CMAP_STOPS = {
    "Nacht": ["#39d0ff", "#1466c7", "#0b1a3a", "#07090d", "#3a1208", "#e2531f", "#ffc36b"],
    "RdBu": ["#053061", "#2166ac", "#4393c3", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d", "#b2182b", "#67001f"],
    "Seismic": ["#00004c", "#0000ff", "#ffffff", "#ff0000", "#800000"],
    "Coolwarm": ["#3b4cc0", "#7b9ff9", "#c0d4f5", "#dddddd", "#f2cbb7", "#ee8468", "#b40426"],
    "Twilight": ["#2f1436", "#5e3f9c", "#7f99c9", "#e2d9e2", "#c99277", "#a14a3f", "#2f1436"],
    "Inferno": ["#000004", "#320a5e", "#781b6c", "#bc3754", "#ed6925", "#fbb61a", "#fcffa4"],
    "Magma": ["#000004", "#2c115f", "#721f81", "#b73779", "#f1605d", "#feb078", "#fcfdbf"],
    "Viridis": ["#440154", "#414487", "#2a788e", "#22a884", "#7ad151", "#fde725"],
    "Turbo": ["#30123b", "#4686fb", "#1ae4b6", "#a2fc3c", "#fabb3a", "#e4460a", "#7a0403"],
    "Grau": ["#000000", "#ffffff"],
  };
  function makeLut(name, n = 512) {
    const stops = (CMAP_STOPS[name] || CMAP_STOPS.Nacht).map(hex);
    const lut = new Uint8ClampedArray(n * 3);
    for (let k = 0; k < n; k++) {
      const x = k / (n - 1) * (stops.length - 1), a = Math.min(Math.floor(x), stops.length - 2), f = x - a;
      for (let c = 0; c < 3; c++) lut[k * 3 + c] = stops[a][c] + (stops[a + 1][c] - stops[a][c]) * f;
    }
    return lut;
  }

  /* Schreibt ein (nx, ny)-Feld in ImageData (Breite nx, Hoehe ny, j nach oben). */
  function renderField(img, data, nx, ny, lut, vmin, vmax, extra) {
    const px = img.data, n = lut.length / 3 - 1, inv = 1 / ((vmax - vmin) || 1);
    for (let j = 0; j < ny; j++) {
      const row = (ny - 1 - j) * nx;
      for (let i = 0; i < nx; i++) {
        const k = i * ny + j;
        let x = (data[k] - vmin) * inv; x = x < 0 ? 0 : x > 1 ? 1 : x;
        const l = ((x * n) | 0) * 3, p = (row + i) * 4;
        px[p] = lut[l]; px[p + 1] = lut[l + 1]; px[p + 2] = lut[l + 2]; px[p + 3] = 255;
        if (extra) extra(px, p, k);
      }
    }
  }

  // ---------------------------------------------------------------- FFT
  function fftMag(y, dt) {
    const n0 = y.length; if (n0 < 32) return null;
    let n = 1; while (n < n0) n <<= 1;
    const re = new Float64Array(n), im = new Float64Array(n);
    let wsum = 0;
    for (let k = 0; k < n0; k++) { const w = 0.5 - 0.5 * Math.cos(2 * Math.PI * k / (n0 - 1)); re[k] = y[k] * w; wsum += w; }
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
      if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const a = -2 * Math.PI / len, wr = Math.cos(a), wi = Math.sin(a);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const ur = re[i + k], ui = im[i + k];
          const vr = re[i + k + len / 2] * cr - im[i + k + len / 2] * ci;
          const vi = re[i + k + len / 2] * ci + im[i + k + len / 2] * cr;
          re[i + k] = ur + vr; im[i + k] = ui + vi; re[i + k + len / 2] = ur - vr; im[i + k + len / 2] = ui - vi;
          const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
    const m = n / 2 + 1, f = new Float64Array(m), mag = new Float64Array(m);
    for (let k = 0; k < m; k++) { f[k] = k / (n * dt); mag[k] = Math.hypot(re[k], im[k]) * 2 / wsum; }
    return { f, mag };
  }

  global.FDTD = { EPS0, MU0, C0, FDTD2D, machWellenform, WELLENFORMEN, Szene, Objekt, Quelle, Sonde,
                  objBox, objEnthaelt, karten, linienPunkte, registriereQuelle, baue, kennzahlen,
                  PRESETS, CMAP_STOPS, makeLut, renderField, fftMag };
})(typeof window !== "undefined" ? window : globalThis);
