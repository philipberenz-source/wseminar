/* main.js — Animationen und interaktive Demos der Startseite. */
(function () {
  "use strict";
  const F = window.FDTD;
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => [...r.querySelectorAll(s)];
  const reduced = matchMedia("(prefers-reduced-motion: reduce)").matches;
  const fmt = (x, n = 2) => x.toLocaleString("de-DE", { minimumFractionDigits: n, maximumFractionDigits: n });

  // ----------------------------------------------------------- Allgemein
  function renderTex() {
    if (!window.katex) return setTimeout(renderTex, 60);
    $$(".tex").forEach(el => {
      try { katex.render(el.textContent, el, { displayMode: el.hasAttribute("data-display"), throwOnError: false, strict: false }); }
      catch (e) { /* Rohtext bleibt stehen */ }
    });
    initYee();
  }
  renderTex();

  // Lightbox
  const lb = $("#lightbox");
  $$("[data-zoom]").forEach(p => p.addEventListener("click", () => { lb.querySelector("img").src = p.querySelector("img").src; lb.classList.add("open"); }));
  lb.addEventListener("click", () => lb.classList.remove("open"));
  addEventListener("keydown", e => { if (e.key === "Escape") lb.classList.remove("open"); });

  // Slider-Fuellstand
  const fill = r => r.style.setProperty("--p", ((r.value - r.min) / (r.max - r.min) * 100) + "%");
  $$("input[type=range]").forEach(r => { fill(r); r.addEventListener("input", () => fill(r)); });

  /* Sichtbarkeitsgesteuerte Animationsschleife: laeuft nur, wenn das
     Element im Bild ist und der Tab aktiv ist. */
  let motion = !reduced;
  try { const m = localStorage.getItem("fdtd-motion"); if (m) motion = m === "an"; } catch (e) { /* egal */ }
  const loops = [];
  function loop(el, frame) {
    let vis = false, raf = 0, first = true;
    const tick = () => {
      if (!vis) { raf = 0; return; }
      if (!motion && !first) { raf = 0; return; }
      first = false; frame(); raf = requestAnimationFrame(tick);
    };
    loops.push(() => { if (vis && !raf) raf = requestAnimationFrame(tick); });
    new IntersectionObserver(es => { vis = es[0].isIntersecting; if (vis && !raf) raf = requestAnimationFrame(tick); }, { threshold: 0.01 }).observe(el);
  }

  const mBtn = $("#motion-btn");
  function setMotion(on) {
    motion = on;
    try { localStorage.setItem("fdtd-motion", on ? "an" : "aus"); } catch (e) { /* egal */ }
    if (mBtn) { mBtn.setAttribute("aria-pressed", String(!on)); mBtn.querySelector("span").textContent = on ? "Animation aus" : "Animation an"; }
    loops.forEach(f => f());
  }
  if (mBtn) mBtn.addEventListener("click", () => setMotion(!motion));
  setMotion(motion);

  // Divergierende "Nacht"-Farben mit Wurzel-Kompression: schwache Felder
  // bleiben sichtbar, Nullfeld = Hintergrund.
  const LUT = F.makeLut("Nacht", 1024), LUT_ABS = F.makeLut("Inferno", 1024);
  /* Flackerschutz: Die Skala folgt dem Spitzenwert nur sehr traege (kein
     Helligkeitspumpen), und mit abs=true wird |Ez| gezeigt - dann gibt es
     keinen Rot/Blau-Wechsel pro Schwingung. */
  function drawField(ctx, img, sim, st, opts = {}) {
    const { nx, ny, Ez, pec } = sim, px = img.data, abs = !!opts.abs;
    const lut = abs ? LUT_ABS : LUT, n = lut.length / 3 - 1;
    let peak = 0;
    for (let k = 0; k < Ez.length; k++) { const a = Math.abs(Ez[k]); if (a > peak) peak = a; }
    if (!st.peak) st.peak = peak;
    // langsam anpassen: hoch schneller als runter, beides weich
    st.peak += (peak - st.peak) * (peak > st.peak ? 0.05 : 0.004);
    st.peak = Math.max(st.peak, 1e-30);
    const inv = 1 / st.peak, gain = opts.gain ?? 1;
    for (let j = 0; j < ny; j++) {
      const row = (ny - 1 - j) * nx;
      for (let i = 0; i < nx; i++) {
        const k = i * ny + j, p = (row + i) * 4;
        if (pec[k]) { px[p] = 44; px[p + 1] = 50; px[p + 2] = 62; px[p + 3] = 255; continue; }
        let v = Ez[k] * inv * gain; v = v > 1 ? 1 : v < -1 ? -1 : v;
        const l = abs ? ((Math.abs(v) * 0.85 * n) | 0) * 3 : (((v + 1) * 0.5 * n) | 0) * 3;
        px[p] = lut[l]; px[p + 1] = lut[l + 1]; px[p + 2] = lut[l + 2]; px[p + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }

  function blit(dst, src, smooth = true) {
    const c = dst.getContext("2d");
    c.imageSmoothingEnabled = smooth; c.imageSmoothingQuality = "high";
    c.drawImage(src, 0, 0, dst.width, dst.height);
  }

  function offscreen(nx, ny) {
    const c = document.createElement("canvas"); c.width = nx; c.height = ny;
    const ctx = c.getContext("2d");
    return { c, ctx, img: ctx.createImageData(nx, ny) };
  }


  // ------------------------------------------------ Ebene Welle (analytisch)
  /* Senkrechter Einfall aus Luft (n1 = 1) auf ein Dielektrikum (n2 = sqrt(eps_r)).
     Medium 1: E = cos(wt - k1 x) + G cos(wt + k1 x),  eta0 H = cos(wt - k1 x) - G cos(wt + k1 x)
     Medium 2: E = T cos(wt - k2 x),                   eta0 H = n2 T cos(wt - k2 x)
     Die Zeit ist stark verlangsamt (0,3 Schwingungen pro Sekunde). */
  (function wave() {
    const cv = $("#wave-canvas"); if (!cv) return;
    const c = cv.getContext("2d"), W = cv.width, H = cv.height;
    const sf = $("#w-f"), se = $("#w-eps"), ce = $("#w-e"), ch = $("#w-h"), cs = $("#w-sep");
    const L = 0.12, x0 = 70, x1 = W - 40, xi = (x0 + x1) / 2, ya = H * 0.52;   // 120 mm Ausschnitt
    const skew = { x: 0.62, y: 0.42 };                                         // Tiefenachse (H)
    let phase = 0, last = performance.now();
    let p = {};
    function params() {
      const f = +sf.value * 1e9, er = +se.value, n = Math.sqrt(er);
      const G = (1 - n) / (1 + n), T = 2 / (1 + n), lam1 = F.C0 / f;
      p = { f, er, n, G, T, k1: 2 * Math.PI / lam1, k2: 2 * Math.PI * n / lam1, lam1 };
      $("#w-f-out").textContent = fmt(+sf.value, 1) + " GHz";
      $("#w-eps-out").textContent = fmt(er);
      $("#w-n").textContent = fmt(n, 3);
      $("#w-lam").textContent = `${fmt(lam1 * 1e3, 1)} mm / ${fmt(lam1 / n * 1e3, 1)} mm`;
      $("#w-g").textContent = (G < 0 ? "−" : "") + fmt(Math.abs(G), 3);
      $("#w-t").textContent = fmt(T, 3);
      $("#w-p").textContent = `${fmt(G * G * 100, 1)} % / ${fmt(n * T * T * 100, 1)} %`;
    }
    [sf, se].forEach(el => el.addEventListener("input", () => { params(); if (!motion) draw(); }));
    [ce, ch, cs].forEach(el => el.addEventListener("change", () => { if (!motion) draw(); }));
    params();

    const amp = H * 0.16;                          // Pixel je Feldeinheit
    const X = xm => x0 + (xm / L + 0.5) * (x1 - x0);
    function curve(fn, color, depth, width = 2.4, alpha = 1, stems = true) {
      c.globalAlpha = alpha;
      c.strokeStyle = color; c.lineWidth = width; c.beginPath();
      const pts = [];
      for (let px = x0; px <= x1; px += 2) {
        const xm = ((px - x0) / (x1 - x0) - 0.5) * L, v = fn(xm);
        const dx = depth ? v * amp * skew.x : 0, dy = depth ? v * amp * skew.y : -v * amp;
        pts.push([px, px + dx, ya + dy]);
      }
      pts.forEach(([, x, y], k) => k ? c.lineTo(x, y) : c.moveTo(x, y)); c.stroke();
      if (stems) {
        c.lineWidth = 1.2; c.globalAlpha = alpha * 0.35; c.beginPath();
        pts.forEach(([px, x, y], k) => { if (k % 5 === 0) { c.moveTo(px, ya); c.lineTo(x, y); } }); c.stroke();
      }
      c.globalAlpha = 1;
    }
    function draw() {
      c.clearRect(0, 0, W, H);
      // Medium 2
      const g = c.createLinearGradient(xi, 0, x1, 0);
      const a = Math.min(0.05 + (p.er - 1) * 0.025, 0.22);
      g.addColorStop(0, `rgba(51,217,178,${a})`); g.addColorStop(1, `rgba(51,217,178,${a * 0.5})`);
      c.fillStyle = g; c.fillRect(xi, 16, x1 - xi + 30, H - 32);
      c.strokeStyle = "rgba(51,217,178,.5)"; c.setLineDash([6, 6]); c.lineWidth = 1.5;
      c.beginPath(); c.moveTo(xi, 16); c.lineTo(xi, H - 16); c.stroke(); c.setLineDash([]);
      c.font = "500 19px JetBrains Mono"; c.fillStyle = "#8d97a9";
      c.fillText("Luft · n₁ = 1", x0, 44);
      c.fillStyle = "#33d9b2"; c.fillText(`εr = ${fmt(p.er)} · n₂ = ${fmt(p.n, 2)}`, xi + 16, 44);
      // Achsen
      c.strokeStyle = "rgba(255,255,255,.25)"; c.lineWidth = 1.2;
      c.beginPath(); c.moveTo(x0 - 20, ya); c.lineTo(x1 + 20, ya); c.stroke();
      c.beginPath(); c.moveTo(x0, ya); c.lineTo(x0, ya - amp * 2.3); c.stroke();
      c.beginPath(); c.moveTo(x0, ya); c.lineTo(x0 + amp * 2.3 * skew.x, ya + amp * 2.3 * skew.y); c.stroke();
      c.fillStyle = "#5d6679"; c.font = "500 17px JetBrains Mono";
      c.fillText("x", x1 + 24, ya + 6); c.fillText("E", x0 - 6, ya - amp * 2.3 - 10);
      c.fillText("H", x0 + amp * 2.3 * skew.x + 6, ya + amp * 2.3 * skew.y + 14);

      const { G, T, k1, k2, n } = p, w = phase;
      const inc = x => x < 0 ? Math.cos(w - k1 * x) : 0;
      const ref = x => x < 0 ? G * Math.cos(w + k1 * x) : 0;
      const tr = x => x >= 0 ? T * Math.cos(w - k2 * x) : 0;
      const sep = cs.checked;
      if (ch.checked) {
        if (sep) {
          curve(inc, "#39c2ff", true, 1.6, 0.45, false); curve(x => -ref(x), "#39c2ff", true, 1.6, 0.45, false);
          curve(x => n * tr(x), "#39c2ff", true, 1.6, 0.45, false);
        } else curve(x => inc(x) - ref(x) + n * tr(x), "#39c2ff", true, 2.2, 0.85);
      }
      if (ce.checked) {
        if (sep) {
          curve(inc, "#ff7a3d", false, 2, 0.9, false); curve(ref, "#ffc36b", false, 2, 0.9, false); curve(tr, "#ff7a3d", false, 2, 0.9, false);
          c.font = "500 16px JetBrains Mono"; c.fillStyle = "#ffc36b"; c.fillText("reflektiert", x0 + 8, H - 30);
          c.fillStyle = "#ff7a3d"; c.fillText("einfallend", x0 + 150, H - 30); c.fillText("transmittiert", xi + 16, H - 30);
        } else curve(x => inc(x) + ref(x) + tr(x), "#ff7a3d", false, 2.6, 1);
      }
      // Ausbreitungsrichtung
      c.fillStyle = "rgba(255,255,255,.35)"; c.font = "500 16px JetBrains Mono";
      c.fillText("→ Ausbreitung", x0, H - 60);
    }
    loop(cv, () => {
      const now = performance.now();
      phase += 2 * Math.PI * 0.3 * Math.min((now - last) / 1000, 0.1);
      last = now; draw();
    });
    draw();
  })();

  // ------------------------------------------------------------------ Yee
  function initYee() {
    const svg = $("#yee-svg"); if (!svg || svg.dataset.ok) return; svg.dataset.ok = 1;
    const NS = "http://www.w3.org/2000/svg";
    const I = 6, J = 4, s = 80, ox = 60, oy = 370;
    const X = i => ox + i * s, Y = j => oy - j * s;
    const nodes = [];
    const mk = (tag, a, parent = svg) => { const e = document.createElementNS(NS, tag); for (const k in a) e.setAttribute(k, a[k]); parent.appendChild(e); return e; };

    // Zellen
    for (let i = 0; i < I; i++) for (let j = 0; j < J; j++)
      mk("rect", { x: X(i), y: Y(j + 1), width: s, height: s, fill: "none", stroke: "rgba(255,255,255,.07)" });
    mk("text", { x: X(I) + 14, y: Y(0) + 5, fill: "#5d6679", "font-family": "JetBrains Mono", "font-size": 13 }).textContent = "x";
    mk("text", { x: X(0) - 5, y: Y(J) - 16, fill: "#5d6679", "font-family": "JetBrains Mono", "font-size": 13 }).textContent = "y";

    const arrow = (g, dir, col) => {
      const L = 17;
      const d = dir === "x" ? `M${-L} 0 L${L - 7} 0 M${L - 12} -6 L${L} 0 L${L - 12} 6` : `M0 ${L} L0 ${-L + 7} M-6 ${-L + 12} L0 ${-L} L6 ${-L + 12}`;
      mk("path", { d, stroke: col, "stroke-width": 3, fill: "none", "stroke-linecap": "round", "stroke-linejoin": "round" }, g);
      mk("circle", { r: 20, fill: "transparent" }, g);
    };
    const add = (type, i, j) => {
      const outer = mk("g", { class: "node", transform: `translate(${X(i)} ${Y(j)})`, tabindex: 0 });
      const g = mk("g", { class: "node-in" }, outer);
      const col = type === "ez" ? "#ff7a3d" : type === "hx" ? "#39c2ff" : "#a58bff";
      g.style.color = col;
      if (type === "ez") { mk("circle", { r: 7.5, fill: col }, g); mk("circle", { r: 18, fill: "transparent" }, g); }
      else arrow(g, type === "hx" ? "x" : "y", col);
      outer.style.color = col;
      const n = { type, i, j, g: outer }; nodes.push(n);
      outer.addEventListener("pointerenter", () => { auto = false; focus(n); });
      outer.addEventListener("focus", () => { auto = false; focus(n); });
      return n;
    };
    for (let i = 0; i <= I; i++) for (let j = 0; j <= J; j++) add("ez", i, j);
    for (let i = 0; i <= I; i++) for (let j = 0; j < J; j++) add("hx", i, j + 0.5);
    for (let i = 0; i < I; i++) for (let j = 0; j <= J; j++) add("hy", i + 0.5, j);
    svg.addEventListener("pointerleave", () => { auto = true; });

    const find = (t, i, j) => nodes.find(n => n.type === t && n.i === i && n.j === j);
    const eqEl = $("#yee-eq"), titleEl = $("#yee-title"), txtEl = $("#yee-text");
    const TEX = {
      ez: String.raw`E_z\big|^{n+1}_{i,j} = C_a E_z\big|^{n}_{i,j} + C_b\Big[\big({\color{#a58bff}H_y|_{i+\frac12,j}} - {\color{#a58bff}H_y|_{i-\frac12,j}}\big) - \big({\color{#39c2ff}H_x|_{i,j+\frac12}} - {\color{#39c2ff}H_x|_{i,j-\frac12}}\big)\Big]^{n+\frac12}`,
      hx: String.raw`H_x\big|^{n+\frac12}_{i,j+\frac12} = H_x\big|^{n-\frac12}_{i,j+\frac12} - \frac{\Delta t}{\mu\Delta}\Big({\color{#ff7a3d}E_z|^{n}_{i,j+1}} - {\color{#ff7a3d}E_z|^{n}_{i,j}}\Big)`,
      hy: String.raw`H_y\big|^{n+\frac12}_{i+\frac12,j} = H_y\big|^{n-\frac12}_{i+\frac12,j} + \frac{\Delta t}{\mu\Delta}\Big({\color{#ff7a3d}E_z|^{n}_{i+1,j}} - {\color{#ff7a3d}E_z|^{n}_{i,j}}\Big)`,
    };
    const TXT = {
      ez: "Die Rotation von H um diesen Punkt: zwei H_y-Werte links und rechts, zwei H_x-Werte darunter und darüber. Sie bilden einen kleinen Umlauf um die Zelle.",
      hx: "H_x ändert sich mit dem Anstieg von E_z in y-Richtung. Dafür braucht es genau die beiden E_z-Punkte darunter und darüber.",
      hy: "H_y ändert sich mit dem Anstieg von E_z in x-Richtung. Dafür braucht es die beiden E_z-Punkte links und rechts.",
    };
    const TITLE = { ez: "UPDATE · E<sub>z</sub> (ganzzahliger Zeitschritt)", hx: "UPDATE · H<sub>x</sub> (halber Zeitschritt)", hy: "UPDATE · H<sub>y</sub> (halber Zeitschritt)" };
    let cur = null;
    function focus(n) {
      if (cur === n) return; cur = n;
      nodes.forEach(m => m.g.classList.remove("on", "src"));
      svg.classList.add("focus");
      n.g.classList.add("on", "src");
      const nb = n.type === "ez"
        ? [find("hy", n.i + 0.5, n.j), find("hy", n.i - 0.5, n.j), find("hx", n.i, n.j + 0.5), find("hx", n.i, n.j - 0.5)]
        : n.type === "hx" ? [find("ez", n.i, n.j + 0.5), find("ez", n.i, n.j - 0.5)]
        : [find("ez", n.i + 0.5, n.j), find("ez", n.i - 0.5, n.j)];
      nb.forEach(m => m && m.g.classList.add("on"));
      if (window.katex) katex.render(TEX[n.type], eqEl, { displayMode: true, throwOnError: false });
      titleEl.innerHTML = TITLE[n.type]; txtEl.innerHTML = TXT[n.type].replace(/([EH])_([xyz])/g, "$1<sub>$2</sub>");
    }
    let auto = true, k = 0;
    const seq = [find("ez", 2, 2), find("hx", 4, 1.5), find("hy", 3.5, 3), find("ez", 4, 1), find("hy", 1.5, 1), find("hx", 2, 2.5)];
    focus(seq[0]);
    setInterval(() => { if (auto && motion) focus(seq[++k % seq.length]); }, 2600);
  }

  // ------------------------------------------------------------- Leapfrog
  (function leap() {
    const svg = $("#leap"); if (!svg) return;
    const NS = "http://www.w3.org/2000/svg";
    const mk = (tag, a) => { const e = document.createElementNS(NS, tag); for (const k in a) e.setAttribute(k, a[k]); svg.appendChild(e); return e; };
    const N = 5, x0 = 70, dx = 100, yE = 60, yH = 140;
    mk("text", { x: 18, y: yE + 5, fill: "#ff7a3d", "font-family": "JetBrains Mono", "font-size": 15 }).textContent = "E";
    mk("text", { x: 18, y: yH + 5, fill: "#a58bff", "font-family": "JetBrains Mono", "font-size": 15 }).textContent = "H";
    mk("line", { x1: 50, x2: 545, y1: yE, y2: yE, stroke: "rgba(255,255,255,.1)" });
    mk("line", { x1: 50, x2: 545, y1: yH, y2: yH, stroke: "rgba(255,255,255,.1)" });
    mk("defs", {}).innerHTML = `<marker id="ah" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="#e9edf4"/></marker>`;
    const steps = [];
    for (let n = 0; n < N; n++) {
      const xe = x0 + n * dx, xh = xe + dx / 2;
      const e = mk("circle", { cx: xe, cy: yE, r: 7, fill: "#1a202b", stroke: "#ff7a3d", "stroke-width": 2, class: "t-node" });
      mk("text", { x: xe, y: yE - 18, fill: "#5d6679", "text-anchor": "middle", "font-family": "JetBrains Mono", "font-size": 11.5 }).textContent = n === 0 ? "n" : `n+${n}`;
      steps.push({ el: e, col: "#ff7a3d", label: n === 0 ? "E bei n" : `E bei n+${n}` });
      if (n < N - 1) {
        const a1 = mk("line", { x1: xe + 6, y1: yE + 7, x2: xh - 6, y2: yH - 9, stroke: "#e9edf4", "stroke-width": 1.5, "marker-end": "url(#ah)", class: "t-arrow" });
        const h = mk("circle", { cx: xh, cy: yH, r: 7, fill: "#1a202b", stroke: "#a58bff", "stroke-width": 2, class: "t-node" });
        mk("text", { x: xh, y: yH + 28, fill: "#5d6679", "text-anchor": "middle", "font-family": "JetBrains Mono", "font-size": 11.5 }).textContent = n === 0 ? "n+½" : `n+${n}½`;
        const a2 = mk("line", { x1: xh + 6, y1: yH - 7, x2: xe + dx - 6, y2: yE + 9, stroke: "#e9edf4", "stroke-width": 1.5, "marker-end": "url(#ah)", class: "t-arrow" });
        steps.push({ el: h, col: "#a58bff", arrow: a1, label: n === 0 ? "H bei n+½" : `H bei n+${n}½` });
        h._next = a2;
      }
    }
    // Pfeil von H nach E gehoert zum folgenden E-Schritt
    for (let k = 2; k < steps.length; k += 2) steps[k].arrow = steps[k - 1].el._next;
    const label = $("#leap-label");
    let k = 0;
    function show() {
      steps.forEach((s, m) => {
        s.el.setAttribute("fill", m <= k ? s.col : "#1a202b");
        s.el.setAttribute("r", m === k ? 9.5 : 7);
        if (s.arrow) s.arrow.classList.toggle("on", m <= k && m > k - 2);
      });
      label.textContent = steps[k].label;
    }
    show();
    let vis = false;
    new IntersectionObserver(es => { vis = es[0].isIntersecting; }).observe(svg);
    setInterval(() => { if (!vis || !motion) return; k = (k + 1) % steps.length; show(); }, 850);
  })();

  // ------------------------------------------------------------------ CFL
  (function cfl() {
    const cv = $("#cfl-canvas"), ecv = $("#cfl-energy"); if (!cv) return;
    const slider = $("#cfl-s"), out = $("#cfl-out"), badge = $("#cfl-badge");
    const n = 96;
    let sim, off = offscreen(n, n), st = { peak: 0 }, E = [], blown = 0, S = 0.99;
    cv.width = 480; cv.height = 480;
    function build() {
      sim = new F.FDTD2D(n, n, 1e-3, { courant: S, npml: 0 });
      const pec = new Uint8Array(n * n);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (i === 0 || j === 0 || i === n - 1 || j === n - 1) pec[i * n + j] = 1;
      sim.setPecMask(pec);
      sim.addCurrentSource(Math.round(n * 0.4), Math.round(n * 0.55), F.machWellenform("ricker", 15e9, 1));
      E = []; st.peak = 0; blown = 0;
      const ok = S <= 1;
      badge.className = "badge " + (ok ? "good" : "bad");
      badge.textContent = ok ? "stabil" : "S > 1 …";
    }
    function set(v) {
      S = +v; slider.value = S; slider.dispatchEvent(new Event("input"));
    }
    slider.addEventListener("input", () => { S = +slider.value; out.textContent = fmt(S); build(); });
    $$("#cfl-demo [data-s]").forEach(b => b.addEventListener("click", () => set(b.dataset.s)));
    $("#cfl-reset").addEventListener("click", build);
    build();

    const ectx = ecv.getContext("2d");
    function drawEnergy() {
      const w = ecv.width, h = ecv.height; ectx.clearRect(0, 0, w, h);
      if (E.length < 2) return;
      const lg = E.map(v => Math.log10(Math.max(v, 1e-40)));
      const ref = Math.max(...lg.slice(0, Math.min(lg.length, 200)));
      const lo = ref - 2, hi = Math.max(ref + 2, ...lg);
      ectx.strokeStyle = "rgba(255,255,255,.08)"; ectx.lineWidth = 1;
      for (let d = Math.ceil(lo); d <= hi; d += Math.max(1, Math.round((hi - lo) / 5))) {
        const y = h - 16 - (d - lo) / (hi - lo) * (h - 30); ectx.beginPath(); ectx.moveTo(0, y); ectx.lineTo(w, y); ectx.stroke();
      }
      const bad = S > 1 && lg[lg.length - 1] > ref + 0.5;
      ectx.strokeStyle = bad ? "#ff5d6c" : "#4cd38f"; ectx.lineWidth = 3; ectx.beginPath();
      const N = 900;  // Schritte je Durchlauf
      lg.forEach((v, m) => { const x = m / N * w, y = h - 16 - (Math.min(v, hi) - lo) / (hi - lo) * (h - 30); m ? ectx.lineTo(x, y) : ectx.moveTo(x, y); });
      ectx.stroke();
      ectx.fillStyle = "#8d97a9"; ectx.font = "500 22px JetBrains Mono"; ectx.fillText("log W(t)", 14, 30);
    }

    loop(cv, () => {
      if (blown) { if (performance.now() - blown > 4000) build(); return; }
      for (let s = 0; s < 1; s++) {
        const ez = sim.Ez.slice(), hx = sim.Hx.slice(), hy = sim.Hy.slice();
        sim.step();
        E.push(sim.interleavedEnergy(ez, hx, hy));
      }
      const last = E[E.length - 1];
      if (!isFinite(last) || last > 1e10 * Math.max(...E.slice(0, 150))) {
        blown = performance.now(); badge.className = "badge bad"; badge.textContent = "instabil: Energie explodiert";
      } else if (S > 1 && E.length > 150 && last > 3 * E[149]) {
        badge.className = "badge bad"; badge.textContent = "instabil: Energie wächst";
      }
      if (E.length >= 900) build();
      // |Ez| statt Ez: die instabile Mode wechselt jeden Schritt das
      // Vorzeichen und wuerde sonst mit Bildwiederholrate flackern.
      drawField(off.ctx, off.img, sim, st, { abs: true });
      blit(cv, off.c);
      drawEnergy();
    });
  })();

  // ----------------------------------------------------------- Dispersion
  (function disp() {
    const cv = $("#disp-canvas"); if (!cv) return;
    const ctx = cv.getContext("2d");
    const sN = $("#disp-n"), sS = $("#disp-s");
    // Numerische Dispersionsrelation nach k aufloesen (Bisektion).
    // Mit dt = S*d/(c*sqrt2):  sin^2(a cos th) + sin^2(a sin th) = (2/S^2) sin^2(pi S /(sqrt2 N)),  a = k d / 2
    function vp(N, S, th) {
      const rhs = 2 / (S * S) * Math.sin(Math.PI * S / (Math.SQRT2 * N)) ** 2;
      let lo = 0, hi = Math.PI / 2;
      const g = a => Math.sin(a * Math.cos(th)) ** 2 + Math.sin(a * Math.sin(th)) ** 2 - rhs;
      for (let it = 0; it < 60; it++) { const m = (lo + hi) / 2; g(m) > 0 ? hi = m : lo = m; }
      const a = (lo + hi) / 2;
      return Math.PI / (N * a);
    }
    function draw() {
      const N = +sN.value, S = +sS.value;
      $("#disp-n-out").textContent = N; $("#disp-s-out").textContent = fmt(S);
      const W = cv.width, cx = W / 2, cy = W / 2, R = W * 0.36;
      ctx.clearRect(0, 0, W, W);
      const pts = []; let vmin = 1;
      for (let d = 0; d <= 360; d += 1) { const v = vp(N, S, d * Math.PI / 180); pts.push(v); vmin = Math.min(vmin, v); }
      // Ueberhoehung: 1 - vp/c wird auf einen festen Bereich gestreckt
      const span = 0.05, rad = v => R * (1 - (1 - v) / span * 0.6);
      // Ringe
      ctx.lineWidth = 1.5;
      [0.98, 0.99, 1].forEach(v => {
        ctx.strokeStyle = v === 1 ? "rgba(255,255,255,.35)" : "rgba(255,255,255,.08)";
        ctx.setLineDash(v === 1 ? [8, 8] : []);
        ctx.beginPath(); ctx.arc(cx, cy, rad(v), 0, 2 * Math.PI); ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = "#5d6679"; ctx.font = "500 20px JetBrains Mono";
        ctx.fillText(v === 1 ? "c" : fmt(v), cx + rad(v) * 0.707 + 6, cy - rad(v) * 0.707 - 6);
      });
      for (let a = 0; a < 360; a += 45) {
        ctx.strokeStyle = "rgba(255,255,255,.06)"; ctx.beginPath(); ctx.moveTo(cx, cy);
        ctx.lineTo(cx + R * 1.15 * Math.cos(a * Math.PI / 180), cy - R * 1.15 * Math.sin(a * Math.PI / 180)); ctx.stroke();
        ctx.fillStyle = "#5d6679"; ctx.font = "500 19px JetBrains Mono"; ctx.textAlign = "center";
        ctx.fillText(a + "°", cx + R * 1.25 * Math.cos(a * Math.PI / 180), cy - R * 1.25 * Math.sin(a * Math.PI / 180) + 6);
        ctx.textAlign = "left";
      }
      // Kurve
      const grad = ctx.createLinearGradient(0, 0, W, W);
      grad.addColorStop(0, "#39c2ff"); grad.addColorStop(1, "#ff7a3d");
      ctx.strokeStyle = grad; ctx.lineWidth = 4; ctx.shadowColor = "rgba(255,122,61,.5)"; ctx.shadowBlur = 18;
      ctx.beginPath();
      pts.forEach((v, d) => { const r = Math.max(rad(v), 4), a = d * Math.PI / 180; const x = cx + r * Math.cos(a), y = cy - r * Math.sin(a); d ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
      ctx.closePath(); ctx.stroke(); ctx.shadowBlur = 0;
      const v0 = pts[0], v45 = pts[45];
      $("#disp-v0").textContent = v0.toFixed(5).replace(".", ",");
      $("#disp-v45").textContent = v45.toFixed(5).replace(".", ",");
      const lag = 100 * (1 / v0 - 1);            // in Wellenlaengen
      $("#disp-lag").textContent = `${fmt(lag * 360, 0)}° ≙ ${fmt(lag * N, 1)} Zellen`;
    }
    sN.addEventListener("input", draw); sS.addEventListener("input", draw);
    draw();
  })();

  // ------------------------------------------------------------ PEC / PML
  (function rand() {
    const make = (id, npml) => {
      const cv = $(id); if (!cv) return;
      const n = 110; cv.width = 440; cv.height = 440;
      const sim = new F.FDTD2D(n, n, 1e-3, { courant: 0.99, npml });
      if (!npml) {
        const pec = new Uint8Array(n * n);
        for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (i === 0 || j === 0 || i === n - 1 || j === n - 1) pec[i * n + j] = 1;
        sim.setPecMask(pec);
      }
      const f = F.machWellenform("ricker", 12e9, 1);
      let t0 = 0;
      sim.addCurrentSource(55, 55, t => f(t - t0));
      const off = offscreen(n, n), st = { peak: 0 };
      let k = 0;
      loop(cv, () => {
        sim.step(); k++;
        if (k > 420) { k = 0; t0 = sim.t; }
        drawField(off.ctx, off.img, sim, st, { gain: 1.4 });
        blit(cv, off.c);
        if (npml) {
          const c = cv.getContext("2d"), p = npml / n * cv.width;
          c.strokeStyle = "rgba(255,255,255,.18)"; c.setLineDash([6, 6]); c.lineWidth = 1.5;
          c.strokeRect(p, p, cv.width - 2 * p, cv.height - 2 * p); c.setLineDash([]);
        }
      });
    };
    make("#pec-canvas", 0);
    make("#pml-canvas", 12);
  })();

})();
