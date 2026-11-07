/* studio.js — Oberflaeche des Web-Studios.
 *
 * Gerechnet wird ausschliesslich in Python: py-worker.js laedt Pyodide und
 * fuehrt fdtd_core.py (unveraendert aus der Arbeit) plus bridge.py aus. Diese
 * Datei baut nur die Szene (JSON im Format des Desktop-Studios), schickt
 * Befehle an den Worker und zeichnet Feld, Overlays und Diagramme.
 */
(function () {
  "use strict";
  const F = window.FDTD;
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => [...r.querySelectorAll(s)];
  const EPS0 = F.EPS0, C0 = F.C0;
  const PROBE_COLORS = ["#39c2ff", "#ff7a3d", "#4cd38f", "#e5c07b", "#a58bff", "#ff6fae", "#9be7ff", "#ffd0a8"];
  const DPR = () => Math.min(devicePixelRatio || 1, 2);
  const de = (x, n = 2) => (+x).toLocaleString("de-DE", { maximumFractionDigits: n });

  // ------------------------------------------------------------ Zustand
  const S = {
    scene: null, info: null, gen: 0,
    ez: null, intens: null, eps: null, pec: null,
    running: false, pending: false, ready: false, autoplay: true,
    probes: [], energy: { t: [], w: [] }, spectra: [], schirmProf: null, nMittel: 0,
    tool: "info", spf: 4, peak: 0,
    view: { mode: "Ez", cmap: "Nacht", scale: "auto", vmax: 1, pml: true, mat: true, marks: true, schirm: true, smooth: true, antiflicker: true },
    material: { art: "dielektrikum", eps_r: 2.25, sigma: 0 },
    open: { quellen: new Set([0]), objekte: new Set(), sonden: new Set() },
    rate: 0, pyMs: 0, lastFrame: 0, drag: null, hover: null,
  };
  try { const v = JSON.parse(localStorage.getItem("fdtd-studio-view") || "null"); if (v) Object.assign(S.view, v); } catch (e) { /* egal */ }
  // Gespeicherte Werte pruefen: eine unbekannte Darstellung oder Farbskala
  // wuerde sonst ein leeres Feld zeichnen.
  const MODUS_WERTE = ["Ez", "|Ez|", "Ez in dB", "Intensität", "Intensität in dB", "el. Energiedichte"];
  if (!MODUS_WERTE.includes(S.view.mode)) S.view.mode = "Ez";
  if (!F.CMAP_STOPS[S.view.cmap]) S.view.cmap = "Nacht";
  const saveView = () => { try { localStorage.setItem("fdtd-studio-view", JSON.stringify(S.view)); } catch (e) { /* egal */ } };

  // ------------------------------------------------------------ Worker
  let worker = null;
  const post = (msg, transfer) => { msg.gen = S.gen; worker.postMessage(msg, transfer || []); };

  function startWorker() {
    if (location.protocol === "file:") {
      loaderError("Das Studio braucht einen lokalen Webserver, weil Python in einem Web Worker läuft.<br><br>" +
        "Im Ordner <code>website</code> ausführen:<br><code>python -m http.server 8000</code><br>und dann <code>http://localhost:8000/studio.html</code> öffnen.");
      return;
    }
    try { worker = new Worker("assets/js/py-worker.js"); }
    catch (e) { loaderError("Web Worker konnte nicht gestartet werden: " + e.message); return; }
    worker.onmessage = onWorker;
    worker.onerror = e => loaderError("Fehler im Python-Worker: " + (e.message || e));
    worker.postMessage({ type: "init" });
  }

  function onWorker(e) {
    const m = e.data;
    switch (m.type) {
      case "status":
        $("#loader-text").textContent = m.text; $("#engine-text").textContent = m.text; break;
      case "ready":
        S.ready = true; S.coreSource = m.coreSource;
        $("#engine").classList.add("ok");
        $("#engine-text").textContent = m.version;
        $("#engine").title = "Rechenkern: fdtd_core.py in Pyodide (WebAssembly)";
        $("#src-pre").textContent = m.coreSource;
        rebuild(true);
        break;
      case "built":
        if (m.gen !== S.gen) return;
        S.info = m.info; S.eps = m.eps; S.pec = m.pec; S.ez = m.ez; S.intens = null;
        S.probes = S.scene.sonden.map(() => []); S.spectra = []; S.energy = { t: [], w: [] };
        S.peak = 0; S.pending = false; S.step = 0; S.t = 0; S.schirmProf = null; S.nMittel = 0;
        $("#schirm-i").max = S.info.nx - 1; $("#schirm-i").value = schirmI();
        buildMaterialLayer(); fit(); render(); drawCharts(); updateStatus(); kpi();
        $("#loader").classList.add("gone");
        ["#btn-play", "#btn-step", "#btn-reset"].forEach(s => $(s).disabled = false);
        if (S.autoplay) { S.autoplay = false; setRunning(true); } else if (S.running) pump();
        break;
      case "frame": {
        if (m.gen !== S.gen) return;
        S.pending = false;
        S.ez = m.ez; if (m.intens) S.intens = m.intens;
        if (m.profil) { S.schirmProf = m.profil; S.nMittel = m.nMittel; }
        S.step = m.step; S.t = m.t; S.pyMs = m.ms;
        m.probes.forEach((y, k) => { const a = S.probes[k]; if (a) for (let q = 0; q < y.length; q++) a.push(y[q]); });
        S.energy.t.push(m.tw); S.energy.w.push(m.w);
        const now = performance.now();
        if (S.lastFrame) S.rate = S.rate * 0.85 + 0.15 * (m.n / ((now - S.lastFrame) / 1000));
        S.lastFrame = now;
        render(); updateStatus();
        if (now - (S.lastCharts || 0) > 120) { drawCharts(); S.lastCharts = now; }
        if (now - (S.lastSpec || 0) > 900) { requestSpectra(); S.lastSpec = now; }
        if (!m.ok) { setRunning(false); toast("Die Rechnung divergiert und wurde gestoppt. Courant-Zahl S ≤ 1 wählen.", "bad"); break; }
        if (S.running) requestAnimationFrame(pump);
        break;
      }
      case "spectrum":
        if (m.gen !== S.gen) return;
        S.spectra[m.k] = { f: m.f, mag: m.mag }; drawSpectrum(); break;
      case "schirm":
        if (m.gen !== S.gen) return;
        S.schirmProf = m.profil; S.nMittel = m.nMittel; drawSchirm(); break;
      case "exec": conOut(m.out || "", m.out && /Error|Traceback/.test(m.out) ? "er" : ""); break;
      case "error":
        if (m.during === "init") loaderError("Python konnte nicht geladen werden: " + m.text + "<br><br>Besteht eine Internetverbindung? Pyodide wird beim ersten Start von cdn.jsdelivr.net geladen.");
        else { toast("Python-Fehler: " + m.text.split("\n").slice(-2).join(" "), "bad"); S.pending = false; setRunning(false); }
        break;
    }
  }

  function loaderError(html) {
    const l = $("#loader"); l.classList.add("err"); l.classList.remove("gone");
    $(".st-load-title", l).textContent = "Python nicht verfügbar";
    $("#loader-text").innerHTML = html; $(".st-load-sub", l).textContent = "";
    $("#engine").classList.add("err"); $("#engine-text").textContent = "kein Python";
  }

  // --------------------------------------------------- Simulation steuern
  function rebuild(autoplay) {
    if (!S.ready) return;
    if (autoplay !== undefined) S.autoplay = autoplay || S.running;
    else S.autoplay = S.running;
    S.gen++; S.pending = true; S.running = false; S.lastFrame = 0;
    const b = $("#btn-play"); b.classList.remove("on"); $(".lbl", b).textContent = "Start";
    post({ type: "build", scene: JSON.stringify(S.scene) });
    renderLists(); overlay();
  }
  let rebuildT = 0;
  const rebuildSoon = () => { clearTimeout(rebuildT); rebuildT = setTimeout(() => rebuild(), 250); };

  /* Flackerschutz: Laeuft eine Dauerstrich-Quelle, schwingt jeder Bildpunkt
     mit f0. Damit das auf dem Bildschirm hoechstens ~2 Hz sind, wird die Zahl
     der Zeitschritte je Frame begrenzt (bei ~60 fps: Periode >= 30 Frames).
     Zeitgemittelte Darstellungen (Intensitaet) flackern nicht und laufen frei. */
  function effSpf() {
    if (!S.view.antiflicker || !S.info || S.view.mode.startsWith("Intensit")) return S.spf;
    const sinus = S.scene.quellen.filter(q => q.aktiv && q.art === "sinus");
    if (!sinus.length) return S.spf;
    const fmax = Math.max(...sinus.map(q => q.f0));
    const perSteps = 1 / (fmax * S.info.dt);
    return Math.max(1, Math.min(S.spf, Math.floor(perSteps / 30)));
  }

  function pump() {
    if (!S.running || S.pending || !S.info) return;
    S.pending = true;
    post({ type: "run", n: effSpf(), intens: true, wantIntens: S.view.mode.startsWith("Intensit"), schirm: schirmI() });
  }
  /* Schirmspalte wie im Desktop-Studio: aus der Szene, sonst 20 Zellen vor
     dem rechten Rand. */
  function schirmI() {
    if (!S.info) return 0;
    const i = S.scene.schirm > 0 ? S.scene.schirm : S.info.nx - 20;
    return Math.min(Math.max(Math.round(i), 0), S.info.nx - 1);
  }
  function setRunning(on) {
    S.running = on && !!S.info;
    const b = $("#btn-play"); b.classList.toggle("on", S.running); $(".lbl", b).textContent = S.running ? "Pause" : "Start";
    if (S.running) { S.lastFrame = 0; pump(); } else requestSpectra();
  }
  function singleStep() {
    if (!S.info || S.pending) return;
    setRunning(false); S.pending = true;
    post({ type: "run", n: 1, intens: true, wantIntens: S.view.mode.startsWith("Intensit"), schirm: schirmI() });
  }
  function requestSpectra() {
    if (!S.info) return;
    S.scene.sonden.forEach((_, k) => post({ type: "spectrum", k }));
  }

  // ------------------------------------------------------------ Feld
  const cvF = $("#cv-field"), cvO = $("#cv-over"), box = $("#canvas-box");
  const ctxF = cvF.getContext("2d"), ctxO = cvO.getContext("2d");
  let img = null, lut = F.makeLut(S.view.cmap, 512), matLayer = null;

  function fit() {
    if (!S.info) return;
    const fld = $("#field"), cs = getComputedStyle(fld);
    const aw = fld.clientWidth - parseFloat(cs.paddingLeft) - parseFloat(cs.paddingRight);
    const ah = fld.clientHeight - parseFloat(cs.paddingTop) - parseFloat(cs.paddingBottom);
    const { nx, ny } = S.info, sc = Math.max(0.5, Math.min(aw / nx, ah / ny));
    const w = Math.floor(nx * sc), h = Math.floor(ny * sc);
    box.style.width = w + "px"; box.style.height = h + "px";
    if (cvF.width !== nx || cvF.height !== ny) { cvF.width = nx; cvF.height = ny; img = ctxF.createImageData(nx, ny); }
    cvO.width = Math.round(w * DPR()); cvO.height = Math.round(h * DPR());
    cvF.classList.toggle("pixel", !S.view.smooth);
    overlay();
  }
  new ResizeObserver(() => fit()).observe($("#field"));

  function fieldData() {
    const { mode } = S.view, ez = S.ez, N = ez.length;
    const out = new Float32Array(N);
    let signed = false, db = false;
    if (mode === "Ez") { out.set(ez); signed = true; }
    else if (mode === "|Ez|") for (let k = 0; k < N; k++) out[k] = Math.abs(ez[k]);
    else if (mode === "Ez in dB") { db = true; let p = 0; for (let k = 0; k < N; k++) p = Math.max(p, Math.abs(ez[k])); p = p || 1; for (let k = 0; k < N; k++) out[k] = Math.max(20 * Math.log10(Math.abs(ez[k]) / p + 1e-12), -60); }
    else if (mode.startsWith("Intensit")) {
      const I = S.intens || new Float32Array(N);
      if (mode.endsWith("dB")) { db = true; let p = 0; for (let k = 0; k < N; k++) p = Math.max(p, I[k]); p = p || 1; for (let k = 0; k < N; k++) out[k] = Math.max(10 * Math.log10(I[k] / p + 1e-12), -60); }
      else out.set(I);
    } else if (mode === "el. Energiedichte") for (let k = 0; k < N; k++) out[k] = 0.5 * EPS0 * S.eps[k] * ez[k] * ez[k];
    return { out, signed, db };
  }

  function render() {
    if (!S.ez || !img) return;
    const { out, signed, db } = fieldData();
    let vmin, vmax;
    if (db) { vmin = -60; vmax = 0; }
    else if (S.view.scale === "fest") { vmax = Math.abs(+S.view.vmax) || 1; vmin = signed ? -vmax : 0; }
    else {
      let p = 0; for (let k = 0; k < out.length; k++) { const a = Math.abs(out[k]); if (a > p) p = a; }
      // traege Skala: steigt zuegig, faellt langsam -> kein Helligkeitspumpen
      if (!S.peak || !isFinite(S.peak)) S.peak = p;
      else S.peak += (p - S.peak) * (p > S.peak ? 0.25 : 0.01);
      vmax = S.peak || 1; vmin = signed ? -vmax : 0;
    }
    F.renderField(img, out, S.info.nx, S.info.ny, lut, vmin, vmax);
    ctxF.putImageData(img, 0, 0);
    colorbar(vmin, vmax, db);
  }

  function colorbar(vmin, vmax, db) {
    const c = $("#cv-bar"), x = c.getContext("2d"), n = lut.length / 3;
    for (let y = 0; y < c.height; y++) {
      const l = Math.floor((1 - y / (c.height - 1)) * (n - 1)) * 3;
      x.fillStyle = `rgb(${lut[l]},${lut[l + 1]},${lut[l + 2]})`; x.fillRect(0, y, c.width, 1);
    }
    const f = v => db ? `${Math.round(v)} dB` : Math.abs(v) >= 1e3 || (Math.abs(v) < 1e-2 && v !== 0) ? v.toExponential(1) : v.toFixed(2);
    $("#bar-max").textContent = f(vmax); $("#bar-min").textContent = f(vmin);
  }

  function buildMaterialLayer() {
    const { nx, ny } = S.info;
    matLayer = document.createElement("canvas"); matLayer.width = nx; matLayer.height = ny;
    const mc = matLayer.getContext("2d"), id = mc.createImageData(nx, ny), px = id.data;
    let emax = 1; for (let k = 0; k < S.eps.length; k++) emax = Math.max(emax, S.eps[k]);
    const karten = F.karten(S.scene);          // fuer Verlustmaske (sigma)
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) {
      const k = i * ny + j, p = ((ny - 1 - j) * nx + i) * 4;
      if (S.pec[k]) { px[p] = 46; px[p + 1] = 52; px[p + 2] = 64; px[p + 3] = 255; }
      else if (karten.sig[k] > 0) { px[p] = 242; px[p + 1] = 166; px[p + 2] = 64; px[p + 3] = 80; }
      else if (S.eps[k] > 1 + 1e-9) { const s = (S.eps[k] - 1) / Math.max(emax - 1, 1e-9); px[p] = 51; px[p + 1] = 217; px[p + 2] = 178; px[p + 3] = Math.round(255 * (0.1 + 0.22 * s)); }
    }
    mc.putImageData(id, 0, 0);
  }

  // ------------------------------------------------------------ Overlay
  const cellToPx = (i, j) => [((i + 0.5) / S.info.nx) * cvO.width, (1 - (j + 0.5) / S.info.ny) * cvO.height];

  function overlay() {
    if (!S.info) return;
    const c = ctxO, W = cvO.width, H = cvO.height, sx = W / S.info.nx, sy = H / S.info.ny, d = DPR();
    c.clearRect(0, 0, W, H);
    if (S.view.mat && matLayer) { c.imageSmoothingEnabled = false; c.drawImage(matLayer, 0, 0, W, H); }
    const np = S.scene.npml;
    if (S.view.pml && np > 0) {
      c.fillStyle = "rgba(160,175,200,.07)";
      c.fillRect(0, 0, W, np * sy); c.fillRect(0, H - np * sy, W, np * sy);
      c.fillRect(0, np * sy, np * sx, H - 2 * np * sy); c.fillRect(W - np * sx, np * sy, np * sx, H - 2 * np * sy);
      c.setLineDash([5 * d, 5 * d]); c.strokeStyle = "rgba(255,255,255,.22)"; c.lineWidth = d;
      c.strokeRect(np * sx, np * sy, W - 2 * np * sx, H - 2 * np * sy); c.setLineDash([]);
      c.fillStyle = "rgba(255,255,255,.35)"; c.font = `500 ${10 * d}px JetBrains Mono`; c.fillText("PML", np * sx + 5 * d, np * sy + 14 * d);
    }
    if (S.view.schirm) {
      const x = (schirmI() + 0.5) * sx;
      c.strokeStyle = "rgba(255,209,102,.7)"; c.lineWidth = 1.2 * d; c.setLineDash([4 * d, 4 * d]);
      c.beginPath(); c.moveTo(x, 0); c.lineTo(x, H); c.stroke(); c.setLineDash([]);
      c.fillStyle = "rgba(255,209,102,.85)"; c.font = `500 ${10 * d}px JetBrains Mono`;
      c.textAlign = "right"; c.fillText("Schirm", x - 4 * d, H - 6 * d); c.textAlign = "left";
    }
    if (S.view.marks) {
      // Objektumrisse
      c.lineWidth = 1.2 * d;
      S.scene.objekte.forEach(o => {
        if (!o.aktiv) return;
        const [i0, j0, i1, j1] = F.objBox(o, S.info.nx, S.info.ny);
        c.strokeStyle = o.art === "pec" ? "rgba(200,210,230,.55)" : "rgba(51,217,178,.6)";
        shape(c, o.form, i0 * sx, H - j1 * sy, (i1 - i0) * sx, (j1 - j0) * sy); c.stroke();
      });
      // Quellen
      S.scene.quellen.forEach(q => {
        if (!q.aktiv) return;
        if (q.typ === "ebene") {
          const x = (q.i + 0.5) * sx;
          c.strokeStyle = "rgba(255,195,107,.75)"; c.lineWidth = 2 * d; c.setLineDash([7 * d, 5 * d]);
          c.beginPath(); c.moveTo(x, 0); c.lineTo(x, H); c.stroke(); c.setLineDash([]);
          const th = q.theta * Math.PI / 180, L = 34 * d, y = H / 2;
          c.beginPath(); c.moveTo(x, y); c.lineTo(x + L * Math.cos(th), y - L * Math.sin(th)); c.stroke();
          c.beginPath(); c.arc(x + L * Math.cos(th), y - L * Math.sin(th), 3.5 * d, 0, 7); c.fillStyle = "#ffc36b"; c.fill();
        } else {
          const [x, y] = cellToPx(q.i, q.j);
          c.beginPath(); c.arc(x, y, 5 * d, 0, 7); c.fillStyle = "#ffc36b"; c.fill();
          c.beginPath(); c.arc(x, y, 10 * d, 0, 7); c.strokeStyle = "rgba(255,195,107,.5)"; c.lineWidth = 1.5 * d; c.stroke();
        }
      });
      // Sonden
      S.scene.sonden.forEach((s, k) => {
        const [x, y] = cellToPx(s.i, s.j), col = PROBE_COLORS[k % PROBE_COLORS.length];
        c.beginPath(); c.arc(x, y, 7 * d, 0, 7); c.strokeStyle = col; c.lineWidth = 2 * d; c.stroke();
        c.beginPath(); c.arc(x, y, 1.8 * d, 0, 7); c.fillStyle = col; c.fill();
        c.font = `600 ${11 * d}px JetBrains Mono`; c.fillStyle = col; c.fillText(String(k + 1), x + 9 * d, y - 7 * d);
      });
    }
    // Zeichnen
    if (S.drag && (S.tool === "rechteck" || S.tool === "ellipse")) {
      const { a, b } = S.drag, i0 = Math.min(a.i, b.i), i1 = Math.max(a.i, b.i) + 1, j0 = Math.min(a.j, b.j), j1 = Math.max(a.j, b.j) + 1;
      c.fillStyle = S.material.art === "pec" ? "rgba(200,210,230,.25)" : "rgba(51,217,178,.2)";
      c.strokeStyle = "#fff"; c.lineWidth = 1.5 * d; c.setLineDash([5 * d, 4 * d]);
      shape(c, S.tool, i0 * sx, H - j1 * sy, (i1 - i0) * sx, (j1 - j0) * sy); c.fill(); c.stroke(); c.setLineDash([]);
      c.font = `500 ${11 * d}px JetBrains Mono`; c.fillStyle = "#fff";
      c.fillText(`${i1 - i0} × ${j1 - j0} Zellen`, i1 * sx + 6 * d, H - j1 * sy + 12 * d);
    }
    if (S.hover && S.tool !== "info") {
      const [x, y] = cellToPx(S.hover.i, S.hover.j);
      c.strokeStyle = "rgba(255,255,255,.5)"; c.lineWidth = d;
      c.beginPath(); c.moveTo(x - 8 * d, y); c.lineTo(x + 8 * d, y); c.moveTo(x, y - 8 * d); c.lineTo(x, y + 8 * d); c.stroke();
    }
  }
  function shape(c, form, x, y, w, h) {
    c.beginPath();
    if (form === "ellipse") c.ellipse(x + w / 2, y + h / 2, Math.max(w / 2, 1), Math.max(h / 2, 1), 0, 0, 7);
    else c.rect(x, y, w, h);
  }

  // ------------------------------------------------------------ Maus
  function evCell(e) {
    const r = cvO.getBoundingClientRect();
    const i = Math.floor((e.clientX - r.left) / r.width * S.info.nx);
    const j = Math.floor((1 - (e.clientY - r.top) / r.height) * S.info.ny);
    return { i: Math.min(Math.max(i, 0), S.info.nx - 1), j: Math.min(Math.max(j, 0), S.info.ny - 1) };
  }
  cvO.addEventListener("pointermove", e => {
    if (!S.info) return;
    const p = evCell(e); S.hover = p;
    if (S.drag) S.drag.b = p;
    const k = p.i * S.info.ny + p.j, d = S.info.d * 1e3;
    const ez = S.ez ? S.ez[k] : 0;
    $("#readout").textContent = `i ${p.i}, j ${p.j} · ${de(p.i * d, 1)} / ${de(p.j * d, 1)} mm · Ez ${ez.toExponential(2)} · εr ${de(S.eps[k], 2)}${S.pec[k] ? " · PEC" : ""}`;
    if (S.tool !== "info" || S.drag) overlay();
  });
  cvO.addEventListener("pointerleave", () => { S.hover = null; overlay(); $("#readout").textContent = "–"; });
  cvO.addEventListener("pointerdown", e => {
    if (!S.info) return;
    const p = evCell(e);
    if (S.tool === "rechteck" || S.tool === "ellipse") { S.drag = { a: p, b: p }; cvO.setPointerCapture(e.pointerId); overlay(); return; }
    if (S.tool === "sonde") { S.scene.sonden.push(F.Sonde({ i: p.i, j: p.j })); toast(`Sonde ${S.scene.sonden.length} bei (${p.i}, ${p.j})`); rebuild(); }
    else if (S.tool === "punkt") {
      const q0 = S.scene.quellen[0];
      S.scene.quellen.push(F.Quelle({ typ: "punkt", art: "sinus", i: p.i, j: p.j, f0: q0 ? q0.f0 : 15e9, zyklen: 3 }));
      toast(`Punktquelle bei (${p.i}, ${p.j})`); rebuild();
    } else if (S.tool === "loeschen") erase(p);
  });
  cvO.addEventListener("pointerup", () => {
    if (!S.drag) return;
    const { a, b } = S.drag; S.drag = null;
    const i0 = Math.min(a.i, b.i), i1 = Math.max(a.i, b.i) + 1, j0 = Math.min(a.j, b.j), j1 = Math.max(a.j, b.j) + 1;
    if (i1 - i0 < 2 && j1 - j0 < 2) { overlay(); return; }
    S.scene.objekte.push(F.Objekt({ form: S.tool, art: S.material.art, i0, j0, i1, j1, eps_r: S.material.eps_r, sigma: S.material.sigma }));
    S.open.objekte = new Set([S.scene.objekte.length - 1]);
    rebuild();
  });
  function erase(p) {
    const near = (i, j) => Math.hypot(i - p.i, j - p.j) < Math.max(4, S.info.nx / 60);
    let k = S.scene.sonden.findIndex(s => near(s.i, s.j));
    if (k >= 0) { S.scene.sonden.splice(k, 1); toast("Sonde entfernt"); return rebuild(); }
    k = S.scene.quellen.findIndex(q => q.typ === "punkt" ? near(q.i, q.j) : Math.abs(q.i - p.i) < 3);
    if (k >= 0) { S.scene.quellen.splice(k, 1); toast("Quelle entfernt"); return rebuild(); }
    for (k = S.scene.objekte.length - 1; k >= 0; k--)
      if (F.objEnthaelt(S.scene.objekte[k], p.i, p.j, S.info.nx, S.info.ny)) { S.scene.objekte.splice(k, 1); toast("Objekt entfernt"); return rebuild(); }
  }
  function setTool(t) {
    S.tool = t; $$(".st-tool").forEach(b => b.classList.toggle("on", b.dataset.tool === t));
    cvO.className = "t-" + t; overlay();
    if (t === "rechteck" || t === "ellipse") showTab("objekte");
  }
  $$(".st-tool").forEach(b => b.addEventListener("click", () => setTool(b.dataset.tool)));

  // ------------------------------------------------------------ Diagramme
  function prepChart(cv) {
    const r = cv.getBoundingClientRect(), d = DPR();
    if (cv.width !== Math.round(r.width * d) || cv.height !== Math.round(r.height * d)) { cv.width = Math.round(r.width * d); cv.height = Math.round(r.height * d); }
    const c = cv.getContext("2d"); c.setTransform(d, 0, 0, d, 0, 0); c.clearRect(0, 0, r.width, r.height);
    return { c, w: r.width, h: r.height };
  }
  function nice(lo, hi, n = 4) {
    const span = hi - lo || 1, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= step0);
    const t = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9 * span; v += step) t.push(v);
    return t;
  }
  const tickFmt = v => Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(0) : de(+v.toPrecision(3), 3);

  function chart(cv, series, o = {}) {
    const { c, w, h } = prepChart(cv);
    const L = 44, R = 10, T = 10, B = 22, pw = w - L - R, ph = h - T - B;
    let x0 = o.xmin, x1 = o.xmax, y0 = o.ymin, y1 = o.ymax;
    if (x0 === undefined) { x0 = Infinity; x1 = -Infinity; series.forEach(s => { if (s.x.length) { x0 = Math.min(x0, s.x[0]); x1 = Math.max(x1, s.x[s.x.length - 1]); } }); }
    if (y0 === undefined) { y0 = Infinity; y1 = -Infinity; series.forEach(s => { for (const v of s.y) if (isFinite(v)) { if (v < y0) y0 = v; if (v > y1) y1 = v; } }); const pad = (y1 - y0) * 0.08 || Math.abs(y1) * 0.1 || 1; y0 -= pad; y1 += pad; if (o.zero) y0 = Math.min(y0, 0); }
    if (!isFinite(x0) || x1 <= x0) { c.fillStyle = "#5d6679"; c.font = "500 11px JetBrains Mono"; c.textAlign = "center"; c.fillText(o.empty || "keine Daten", w / 2, h / 2); c.textAlign = "left"; return; }
    const X = v => L + (v - x0) / (x1 - x0) * pw, Y = v => T + (1 - (v - y0) / (y1 - y0)) * ph;
    c.font = "500 10px JetBrains Mono"; c.lineWidth = 1;
    nice(y0, y1).forEach(v => { c.strokeStyle = "rgba(255,255,255,.06)"; c.beginPath(); c.moveTo(L, Y(v)); c.lineTo(w - R, Y(v)); c.stroke(); c.fillStyle = "#5d6679"; c.textAlign = "right"; c.fillText(o.yfmt ? o.yfmt(v) : tickFmt(v), L - 6, Y(v) + 3); });
    nice(x0, x1, 5).forEach(v => { c.strokeStyle = "rgba(255,255,255,.04)"; c.beginPath(); c.moveTo(X(v), T); c.lineTo(X(v), T + ph); c.stroke(); c.fillStyle = "#5d6679"; c.textAlign = "center"; c.fillText(tickFmt(v), X(v), h - 7); });
    c.textAlign = "right"; c.fillStyle = "#8d97a9"; c.fillText(o.xlabel || "", w - R, h - 7 - 11); c.textAlign = "left";
    (o.marks || []).forEach((m, q) => {
      c.strokeStyle = "rgba(255,195,107,.5)"; c.setLineDash([3, 3]); c.beginPath(); c.moveTo(X(m), T); c.lineTo(X(m), T + ph); c.stroke(); c.setLineDash([]);
      if (o.markLabels) { c.fillStyle = "#8d97a9"; c.textAlign = "left"; c.fillText(o.markLabels[q], X(m) + 3, T + 9 + (q % 2) * 11); }
    });
    c.save(); c.beginPath(); c.rect(L, T, pw, ph); c.clip();
    series.forEach(s => {
      const n = s.x.length; if (!n) return;
      c.strokeStyle = s.color; c.lineWidth = 1.5; c.beginPath();
      // Min/Max je Pixelspalte, damit lange Reihen nicht zu Balken werden
      const cols = Math.max(2, Math.floor(pw)), per = n / cols;
      if (per <= 2) { for (let k = 0; k < n; k++) { const px = X(s.x[k]), py = Y(s.y[k]); k ? c.lineTo(px, py) : c.moveTo(px, py); } }
      else for (let q = 0; q < cols; q++) {
        const a = Math.floor(q * per), b = Math.min(n, Math.floor((q + 1) * per)); let mn = Infinity, mx = -Infinity;
        for (let k = a; k < b; k++) { const v = s.y[k]; if (v < mn) mn = v; if (v > mx) mx = v; }
        const px = X(s.x[a]); q ? c.lineTo(px, Y(mn)) : c.moveTo(px, Y(mn)); c.lineTo(px, Y(mx));
      }
      c.stroke();
    });
    c.restore();
  }

  function drawCharts() { drawOsc(); drawEnergy(); drawSchirm(); }

  /* Intensitaetsschirm wie im Desktop-Studio: <Ez^2> laengs der Spalte i,
     auf das Schirmmaximum normiert; im linearen Modus werden die Maxima
     (> 0,15, Randzellen ausgenommen) mit ihrem j markiert. */
  function drawSchirm() {
    if (!S.info) return;
    const i = schirmI(), db = $("#schirm-scale").value === "db", p = S.schirmProf;
    $("#schirm-info").textContent = `Spalte i = ${i}` + (S.nMittel ? ` · ${S.nMittel.toLocaleString("de-DE")} Schritte gemittelt` : "");
    if (!p || !p.length) { chart($("#ch-schirm"), [], { empty: "noch keine Werte" }); return; }
    let spitze = 0; for (const v of p) spitze = Math.max(spitze, v);
    const n = p.length, x = new Float64Array(n), y = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      const v = spitze > 0 ? p[j] / spitze : 0;
      x[j] = j; y[j] = db ? Math.max(10 * Math.log10(v + 1e-12), -50) : v;
    }
    const maxima = [];
    if (!db && spitze > 0)
      for (let j = 6; j < n - 6 && maxima.length < 9; j++)
        if (y[j] > y[j - 1] && y[j] >= y[j + 1] && y[j] > 0.15) maxima.push(j);
    chart($("#ch-schirm"), [{ x, y, color: "#ffd166" }], {
      xlabel: "j", xmin: 0, xmax: n - 1, ymin: db ? -52 : 0, ymax: db ? 3 : 1.12,
      marks: maxima, markLabels: maxima.map(String), yfmt: db ? v => `${Math.round(v)} dB` : undefined });
  }
  const schirmFeld = $("#schirm-i");
  schirmFeld.addEventListener("change", () => {
    if (!S.info) return;
    const v = parseInt(schirmFeld.value, 10);
    if (!isFinite(v)) { schirmFeld.value = schirmI(); return; }
    S.scene.schirm = Math.min(Math.max(v, 1), S.info.nx - 1); schirmFeld.value = S.scene.schirm;
    overlay(); post({ type: "schirm", i: S.scene.schirm });
  });
  $("#schirm-scale").addEventListener("change", drawSchirm);
  $("#schirm-reset").addEventListener("click", () => {
    if (!S.info) return;
    if (S.intens) S.intens = new Float32Array(S.intens.length);
    post({ type: "intensReset", i: schirmI() });
    if (!S.running) render();
    toast("Zeitmittelung neu gestartet");
  });
  function drawOsc() {
    if (!S.info) return;
    const dt = S.info.dt, f0 = S.scene.quellen.length ? S.scene.quellen[0].f0 : 15e9;
    const win = +$("#osc-win").value;
    const series = S.probes.map((y, k) => {
      let a = 0; if (win) a = Math.max(0, y.length - Math.round(win / (f0 * dt)));
      const n = y.length - a, x = new Float64Array(n), yy = new Float64Array(n);
      for (let q = 0; q < n; q++) { x[q] = (a + q + 1) * dt * 1e9; yy[q] = y[a + q]; }
      return { x, y: yy, color: PROBE_COLORS[k % PROBE_COLORS.length] };
    });
    chart($("#ch-osc"), series, { xlabel: "t / ns", empty: S.scene.sonden.length ? "noch keine Werte" : "keine Sonde gesetzt" });
  }
  function drawEnergy() {
    const log = $("#en-scale").value === "log";
    const x = S.energy.t.map(t => t * 1e9), y = log ? S.energy.w.map(w => Math.log10(Math.max(w, 1e-40))) : S.energy.w;
    chart($("#ch-energy"), [{ x, y, color: "#4cd38f" }], { xlabel: "t / ns", zero: !log, yfmt: log ? v => "1e" + Math.round(v) : undefined });
  }
  function drawSpectrum() {
    const db = $("#spec-scale").value === "db";
    const f0 = S.scene.quellen.length ? Math.max(...S.scene.quellen.map(q => q.f0)) : 15e9, fmax = 3 * f0 / 1e9;
    let gmax = 0; S.spectra.forEach(s => s && s.mag.forEach(v => { gmax = Math.max(gmax, v); }));
    const series = S.spectra.map((s, k) => {
      if (!s) return { x: [], y: [], color: "#fff" };
      const x = [], y = [];
      for (let q = 0; q < s.f.length; q++) { const g = s.f[q] / 1e9; if (g > fmax) break; x.push(g); y.push(db ? 20 * Math.log10(s.mag[q] / (gmax || 1) + 1e-9) : s.mag[q]); }
      return { x, y, color: PROBE_COLORS[k % PROBE_COLORS.length] };
    });
    chart($("#ch-spec"), series, { xlabel: "f / GHz", xmin: 0, xmax: fmax, marks: S.scene.quellen.filter(q => q.aktiv).map(q => q.f0 / 1e9),
      ymin: db ? -60 : 0, ymax: db ? 3 : undefined, zero: true, empty: "Spektrum ab 64 Werten" });
  }
  ["#osc-win", "#en-scale"].forEach(s => $(s).addEventListener("change", drawCharts));
  $("#spec-scale").addEventListener("change", drawSpectrum);
  addEventListener("resize", () => { drawCharts(); drawSpectrum(); });

  // ------------------------------------------------------------ Status
  const fmtT = t => t < 1e-9 ? de(t * 1e12, 1) + " ps" : de(t * 1e9, 3) + " ns";
  function updateStatus() {
    if (!S.info) return;
    $("#s-step").textContent = (S.step || 0).toLocaleString("de-DE");
    $("#s-time").textContent = fmtT(S.t || 0);
    $("#s-dt").textContent = de(S.info.dt * 1e12, 3) + " ps";
    const eff = effSpf();
    $("#s-rate").textContent = S.running ? `${Math.round(S.rate).toLocaleString("de-DE")} Schritte/s` + (eff < S.spf ? " · Flackerschutz" : "") : "pausiert";
    $("#s-ms").textContent = S.pyMs ? de(S.pyMs / Math.max(eff, 1), 2) + " ms/Schritt" : "–";
    $("#s-grid").textContent = `${S.info.nx} × ${S.info.ny} Zellen · Δ = ${de(S.info.d * 1e3, 2)} mm`;
  }

  // ------------------------------------------------------------ Formulare
  function field(parent, { label, unit, get, set, type = "number", options, step, full, change = rebuildSoon, min, max }) {
    const wrap = document.createElement("label"); wrap.className = "st-f" + (full ? " full" : "");
    if (type === "check") {
      wrap.className = "st-check" + (full ? " st-f full" : "");
      const inp = document.createElement("input"); inp.type = "checkbox"; inp.checked = !!get();
      inp.addEventListener("change", () => { set(inp.checked); change(); });
      wrap.append(inp, document.createTextNode(label)); parent.append(wrap); return inp;
    }
    const sp = document.createElement("span"); sp.innerHTML = label + (unit ? ` <i>${unit}</i>` : ""); wrap.append(sp);
    let inp;
    if (type === "select") {
      inp = document.createElement("select");
      options.forEach(([v, t]) => { const o = document.createElement("option"); o.value = v; o.textContent = t; inp.append(o); });
      inp.value = get();
      inp.addEventListener("change", () => { set(inp.value); change(); });
    } else {
      inp = document.createElement("input"); inp.type = "number"; inp.value = get(); inp.step = step || "any";
      if (min !== undefined) inp.min = min; if (max !== undefined) inp.max = max;
      inp.addEventListener("change", () => { const v = parseFloat(inp.value.replace(",", ".")); if (isFinite(v)) { set(v); change(); } else inp.value = get(); });
    }
    wrap.append(inp); parent.append(wrap); return inp;
  }

  const g = () => S.scene;
  function renderGitter() {
    const f = $("#form-gitter"); f.innerHTML = "";
    const nop = () => {};
    S.gitterDraft = { nx: g().nx, ny: g().ny, d: g().d, courant: g().courant, npml: g().npml, pml_m: g().pml_m, r0: g().r0 };
    const D = S.gitterDraft;
    field(f, { label: "nx", unit: "Zellen", get: () => D.nx, set: v => D.nx = Math.round(Math.min(Math.max(v, 20), 1200)), step: 1, change: nop });
    field(f, { label: "ny", unit: "Zellen", get: () => D.ny, set: v => D.ny = Math.round(Math.min(Math.max(v, 20), 1200)), step: 1, change: nop });
    field(f, { label: "Δ", unit: "mm", get: () => D.d * 1e3, set: v => D.d = Math.max(v, 1e-4) * 1e-3, change: nop });
    field(f, { label: "Courant S", get: () => D.courant, set: v => D.courant = v, step: 0.01, change: nop });
    field(f, { label: "PML", unit: "Zellen", get: () => D.npml, set: v => D.npml = Math.max(0, Math.round(v)), step: 1, change: nop });
    field(f, { label: "PML-Ordnung m", get: () => D.pml_m, set: v => D.pml_m = Math.max(1, Math.round(v)), step: 1, change: nop });
    field(f, { label: "R₀ (Soll-Reflexion)", get: () => D.r0, set: v => D.r0 = Math.min(Math.max(v, 1e-16), 0.5), full: true, change: nop });
  }
  $("#btn-gitter").addEventListener("click", () => {
    Object.assign(S.scene, S.gitterDraft);
    if (S.scene.courant > 1) toast("S > 1 verletzt die CFL-Bedingung: die Rechnung wird instabil.", "warn");
    if (S.scene.nx * S.scene.ny > 400000) toast("Großes Gitter: Python braucht dafür spürbar länger pro Schritt.", "warn");
    rebuild(); kpi();
  });

  function kpi() {
    const k = F.kennzahlen(S.scene), el = $("#kpi");
    const nl = k.n_lambda;
    el.innerHTML = `
      <dt>Gebiet</dt><dd>${de(k.breite * 1e3, 0)} × ${de(k.hoehe * 1e3, 0)} mm</dd>
      <dt>Zellen</dt><dd>${(S.scene.nx * S.scene.ny).toLocaleString("de-DE")}</dd>
      <dt>Zeitschritt Δt</dt><dd>${de(k.dt * 1e12, 3)} ps</dd>
      <dt>f₀ (Quelle 1)</dt><dd>${k.f0 ? de(k.f0 / 1e9, 2) + " GHz" : "–"}</dd>
      <dt>λ₀</dt><dd>${k.f0 ? de(C0 / k.f0 * 1e3, 1) + " mm" : "–"}</dd>
      <dt>Zellen je λ (N<sub>λ</sub>)</dt><dd class="${nl < 10 ? "warn" : ""}">${isFinite(nl) ? de(nl, 1) : "–"}</dd>
      <dt>Courant S</dt><dd class="${S.scene.courant > 1 ? "warn" : ""}">${de(S.scene.courant, 3)}</dd>`;
  }

  const ARTEN = [["sinus", "Sinus (Dauerstrich)"], ["ricker", "Ricker-Puls"], ["gauss", "Gauß-Puls"]];
  function item(list, key, idx, { color, title, sub, off, onDelete, body }) {
    const it = document.createElement("div"); it.className = "st-item" + (S.open[key].has(idx) ? " open" : "") + (off ? " off" : "");
    const h = document.createElement("div"); h.className = "st-item-h";
    h.innerHTML = `<span class="sw" style="background:${color}"></span><div style="min-width:0;display:flex;flex-direction:column"><b>${title}</b><small>${sub}</small></div>`;
    const x = document.createElement("button"); x.className = "x"; x.textContent = "✕"; x.title = "Entfernen";
    x.addEventListener("click", ev => { ev.stopPropagation(); onDelete(); });
    h.append(x);
    h.addEventListener("click", () => { S.open[key].has(idx) ? S.open[key].delete(idx) : S.open[key].add(idx); it.classList.toggle("open"); });
    const b = document.createElement("div"); b.className = "st-item-b"; body(b);
    it.append(h, b); list.append(it);
  }

  function renderLists() {
    // Quellen
    const lq = $("#list-quellen"); lq.innerHTML = "";
    if (!g().quellen.length) lq.innerHTML = `<div class="st-empty">Keine Quelle. Füge unten eine hinzu.</div>`;
    g().quellen.forEach((q, k) => item(lq, "quellen", k, {
      color: "#ffc36b", off: !q.aktiv,
      title: (q.typ === "ebene" ? "Ebene Welle" : "Punktquelle"),
      sub: `${ARTEN.find(a => a[0] === q.art)[1].split(" (")[0]} · ${de(q.f0 / 1e9, 2)} GHz · ` + (q.typ === "ebene" ? `i = ${q.i}, ${de(q.theta, 1)}°` : `(${q.i}, ${q.j})`),
      onDelete: () => { g().quellen.splice(k, 1); rebuild(); },
      body: b => {
        const re = () => { renderLists(); rebuildSoon(); };
        field(b, { label: "Typ", type: "select", options: [["ebene", "Ebene Welle"], ["punkt", "Punktquelle"]], get: () => q.typ, set: v => q.typ = v, change: re });
        field(b, { label: "Wellenform", type: "select", options: ARTEN, get: () => q.art, set: v => q.art = v, change: re });
        field(b, { label: "f₀", unit: "GHz", get: () => q.f0 / 1e9, set: v => q.f0 = Math.max(v, 0.01) * 1e9, change: () => { kpi(); re(); } });
        field(b, { label: "Amplitude", get: () => q.amp, set: v => q.amp = v, change: re });
        field(b, { label: "i", unit: "Spalte", get: () => q.i, set: v => q.i = Math.round(v), step: 1, change: re });
        if (q.typ === "punkt") field(b, { label: "j", unit: "Zeile", get: () => q.j, set: v => q.j = Math.round(v), step: 1, change: re });
        else field(b, { label: "Winkel θ", unit: "°", get: () => q.theta, set: v => q.theta = Math.min(Math.max(v, -85), 85), change: re });
        if (q.art === "sinus") field(b, { label: "Einschaltrampe", unit: "Perioden", get: () => q.zyklen, set: v => q.zyklen = Math.max(v, 0), change: re });
        field(b, { label: "Verzögerung", unit: "Perioden", get: () => q.verzoegerung, set: v => q.verzoegerung = Math.max(v, 0), change: re });
        field(b, { label: "Anregung", type: "select", options: [["strom", "Strom J_z"], ["weich", "weich (E_z)"]], get: () => q.modus, set: v => q.modus = v, change: re });
        if (q.typ === "ebene") {
          field(b, { label: "Arm entlang der Kante", type: "check", full: true, get: () => q.arm, set: v => q.arm = v, change: re });
          field(b, { label: "λ/4-Paar gegen Rückwärtswelle", type: "check", full: true, get: () => q.richtung, set: v => q.richtung = v, change: re });
        }
        field(b, { label: "aktiv", type: "check", full: true, get: () => q.aktiv, set: v => q.aktiv = v, change: re });
      }
    }));

    // Objekte
    const lo = $("#list-objekte"); lo.innerHTML = "";
    if (!g().objekte.length) lo.innerHTML = `<div class="st-empty">Keine Objekte. Zeichne eines ins Feld.</div>`;
    g().objekte.forEach((o, k) => item(lo, "objekte", k, {
      color: o.art === "pec" ? "#8a93a6" : "#33d9b2", off: !o.aktiv,
      title: `${o.form === "ellipse" ? "Ellipse" : "Rechteck"} · ${o.art === "pec" ? "PEC" : "εr = " + de(o.eps_r, 2)}`,
      sub: `i ${o.i0}…${o.i1}  j ${o.j0}…${o.j1}` + (o.art !== "pec" && o.sigma ? ` · σ = ${o.sigma}` : ""),
      onDelete: () => { g().objekte.splice(k, 1); S.open.objekte.clear(); rebuild(); },
      body: b => {
        const re = () => { renderLists(); rebuildSoon(); };
        field(b, { label: "Form", type: "select", options: [["rechteck", "Rechteck"], ["ellipse", "Ellipse"]], get: () => o.form, set: v => o.form = v, change: re });
        field(b, { label: "Material", type: "select", options: [["dielektrikum", "Dielektrikum"], ["pec", "PEC (Metall)"]], get: () => o.art, set: v => o.art = v, change: re });
        field(b, { label: "i₀", get: () => o.i0, set: v => o.i0 = Math.round(v), step: 1, change: re });
        field(b, { label: "i₁", get: () => o.i1, set: v => o.i1 = Math.round(v), step: 1, change: re });
        field(b, { label: "j₀", get: () => o.j0, set: v => o.j0 = Math.round(v), step: 1, change: re });
        field(b, { label: "j₁", get: () => o.j1, set: v => o.j1 = Math.round(v), step: 1, change: re });
        if (o.art !== "pec") {
          field(b, { label: "εr", get: () => o.eps_r, set: v => o.eps_r = Math.max(v, 1e-3), change: re });
          field(b, { label: "σ", unit: "S/m", get: () => o.sigma, set: v => o.sigma = Math.max(v, 0), change: re });
        }
        field(b, { label: "aktiv", type: "check", full: true, get: () => o.aktiv, set: v => o.aktiv = v, change: re });
      }
    }));

    // Sonden
    const ls = $("#list-sonden"); ls.innerHTML = "";
    if (!g().sonden.length) ls.innerHTML = `<div class="st-empty">Keine Sonden. Setze eine mit dem Sonden-Werkzeug.</div>`;
    g().sonden.forEach((s, k) => item(ls, "sonden", k, {
      color: PROBE_COLORS[k % PROBE_COLORS.length], title: `Sonde ${k + 1}${s.name ? " · " + s.name : ""}`, sub: `(${s.i}, ${s.j})`,
      onDelete: () => { g().sonden.splice(k, 1); rebuild(); },
      body: b => {
        field(b, { label: "i", get: () => s.i, set: v => s.i = Math.round(v), step: 1 });
        field(b, { label: "j", get: () => s.j, set: v => s.j = Math.round(v), step: 1 });
      }
    }));
  }

  $$("[data-add]").forEach(b => b.addEventListener("click", () => {
    const typ = b.dataset.add, sz = g();
    sz.quellen.push(F.Quelle(typ === "ebene" ? { typ, art: "sinus", i: Math.max(sz.npml + 10, 20), f0: 15e9 }
                                             : { typ, art: "sinus", i: sz.nx >> 1, j: sz.ny >> 1, f0: 15e9, zyklen: 3 }));
    S.open.quellen = new Set([sz.quellen.length - 1]); rebuild();
  }));
  $("#btn-obj-add").addEventListener("click", () => {
    const sz = g(), cx = sz.nx >> 1, cy = sz.ny >> 1, r = Math.round(Math.min(sz.nx, sz.ny) / 8);
    sz.objekte.push(F.Objekt({ form: "ellipse", art: S.material.art, i0: cx - r, i1: cx + r, j0: cy - r, j1: cy + r, eps_r: S.material.eps_r, sigma: S.material.sigma }));
    S.open.objekte = new Set([sz.objekte.length - 1]); rebuild();
  });
  $("#btn-obj-clear").addEventListener("click", () => { if (g().objekte.length) { g().objekte = []; rebuild(); } });

  function renderMaterial() {
    const f = $("#form-material"); f.innerHTML = "";
    const M = S.material, re = () => renderMaterial();
    field(f, { label: "Material", type: "select", options: [["dielektrikum", "Dielektrikum"], ["pec", "PEC (Metall)"]], get: () => M.art, set: v => M.art = v, change: re, full: M.art === "pec" });
    if (M.art !== "pec") {
      field(f, { label: "εr", get: () => M.eps_r, set: v => M.eps_r = Math.max(v, 1e-3), change: () => {} });
      field(f, { label: "σ", unit: "S/m", get: () => M.sigma, set: v => M.sigma = Math.max(v, 0), change: () => {} });
    }
  }

  const MODES = ["Ez", "|Ez|", "Ez in dB", "Intensität ⟨Ez²⟩", "Intensität in dB", "el. Energiedichte"];
  function renderAnzeige() {
    const f = $("#form-anzeige"); f.innerHTML = "";
    const V = S.view, upd = () => { saveView(); S.peak = 0; lut = F.makeLut(V.cmap, 512); fit(); render(); updateStatus(); };
    field(f, { label: "Darstellung", type: "select", full: true, options: MODES.map(m => [m.replace("⟨Ez²⟩", "").trim().replace("Intensität", "Intensität"), m]),
      get: () => V.mode, set: v => { V.mode = v; if (v.startsWith("Intensit") && !["Inferno", "Magma", "Viridis", "Turbo", "Grau"].includes(V.cmap)) V.cmap = "Inferno"; if (v === "Ez" && ["Inferno", "Magma", "Viridis", "Grau"].includes(V.cmap)) V.cmap = "Nacht"; renderAnzeige(); }, change: upd });
    field(f, { label: "Farbskala", type: "select", options: Object.keys(F.CMAP_STOPS).map(k => [k, k]), get: () => V.cmap, set: v => V.cmap = v, change: upd });
    field(f, { label: "Skalierung", type: "select", options: [["auto", "automatisch"], ["fest", "fest"]], get: () => V.scale, set: v => { V.scale = v; renderAnzeige(); }, change: upd });
    if (V.scale === "fest") field(f, { label: "Maximum", get: () => V.vmax, set: v => V.vmax = v, change: upd, full: true });
    field(f, { label: "Schritte je Frame", get: () => S.spf, set: v => { S.spf = Math.min(Math.max(Math.round(v), 1), 50); $("#spf").value = S.spf; $("#spf-out").textContent = S.spf; }, step: 1, min: 1, max: 50, full: true, change: updateStatus });
    field(f, { label: "Flackerschutz (empfohlen)", type: "check", full: true, get: () => V.antiflicker, set: v => { V.antiflicker = v; if (!v) toast("Flackerschutz aus: Dauerstrich-Quellen können schnell blinken.", "warn"); }, change: upd });
    field(f, { label: "Weich interpolieren", type: "check", full: true, get: () => V.smooth, set: v => V.smooth = v, change: upd });
    field(f, { label: "PML-Rand zeigen", type: "check", full: true, get: () => V.pml, set: v => V.pml = v, change: () => { saveView(); overlay(); } });
    field(f, { label: "Materialien einfärben", type: "check", full: true, get: () => V.mat, set: v => V.mat = v, change: () => { saveView(); overlay(); } });
    field(f, { label: "Quellen, Sonden, Umrisse", type: "check", full: true, get: () => V.marks, set: v => V.marks = v, change: () => { saveView(); overlay(); } });
    field(f, { label: "Schirmlinie zeigen", type: "check", full: true, get: () => V.schirm, set: v => V.schirm = v, change: () => { saveView(); overlay(); } });
    const p = document.createElement("p"); p.className = "st-hint st-f full";
    p.innerHTML = "Der <b>Flackerschutz</b> begrenzt bei Dauerstrich-Quellen die Rechenschritte pro Bild, sodass jeder Bildpunkt höchstens etwa zweimal pro Sekunde wechselt. Zeitgemittelte Darstellungen (<b>Intensität</b>) flackern nicht und laufen immer mit voller Geschwindigkeit.";
    f.append(p);
  }

  function showTab(t) {
    $$(".st-tabs button").forEach(b => b.classList.toggle("on", b.dataset.tab === t));
    $$(".st-pane").forEach(p => p.classList.toggle("on", p.dataset.pane === t));
  }
  $$(".st-tabs button").forEach(b => b.addEventListener("click", () => showTab(b.dataset.tab)));

  // ------------------------------------------------------------ Szenen
  const sel = $("#preset");
  Object.keys(F.PRESETS).forEach(k => { const o = document.createElement("option"); o.value = k; o.textContent = k; sel.append(o); });
  const custom = document.createElement("option"); custom.value = "__eigen"; custom.textContent = "Eigene Szene"; custom.hidden = true; sel.append(custom);
  function loadScene(sz, name) {
    S.scene = F.Szene(sz);
    S.open = { quellen: new Set([0]), objekte: new Set(), sonden: new Set() };
    sel.value = name && F.PRESETS[name] ? name : "__eigen";
    renderGitter(); renderLists(); kpi();
    // Intensitaetsbilder lohnen sich bei Beugungsszenen: Dauerstrich
    rebuild(true);
  }
  sel.addEventListener("change", () => { if (F.PRESETS[sel.value]) loadScene(F.PRESETS[sel.value](), sel.value); });

  // ------------------------------------------------------------ Export
  function download(name, blob) {
    const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name;
    document.body.append(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }
  $("#btn-save").addEventListener("click", () => download("szene.json", new Blob([JSON.stringify(S.scene, null, 2)], { type: "application/json" })));
  $("#file-load").addEventListener("change", e => {
    const f = e.target.files[0]; if (!f) return;
    f.text().then(t => { loadScene(JSON.parse(t)); toast("Szene geladen: " + f.name, "good"); })
      .catch(err => toast("Laden fehlgeschlagen: " + err.message, "bad"));
    e.target.value = "";
  });
  $("#btn-png").addEventListener("click", () => {
    if (!S.info) return;
    const W = cvO.width, H = cvO.height, c = document.createElement("canvas"); c.width = W; c.height = H;
    const x = c.getContext("2d"); x.imageSmoothingEnabled = S.view.smooth; x.drawImage(cvF, 0, 0, W, H); x.drawImage(cvO, 0, 0);
    c.toBlob(b => download(`feld_n${S.step || 0}.png`, b));
  });
  $("#btn-csv").addEventListener("click", () => {
    if (!S.probes.length || !S.probes[0].length) return toast("Noch keine Sondendaten.", "warn");
    const n = Math.max(...S.probes.map(p => p.length)), dt = S.info.dt;
    const rows = ["t_s," + S.scene.sonden.map((s, k) => `sonde${k + 1}_i${s.i}_j${s.j}`).join(",")];
    for (let q = 0; q < n; q++) rows.push([((q + 1) * dt).toExponential(9), ...S.probes.map(p => p[q] !== undefined ? p[q].toExponential(9) : "")].join(","));
    download("sonden.csv", new Blob([rows.join("\n")], { type: "text/csv" }));
  });

  // ------------------------------------------------------------ Konsole
  const conPre = $("#con-out"), conIn = $("#con-in"), hist = []; let hi = 0;
  function conOut(text, cls) {
    const s = document.createElement("span"); if (cls) s.className = cls; s.textContent = text + "\n";
    conPre.append(s); conPre.scrollTop = conPre.scrollHeight;
  }
  $("#con-form").addEventListener("submit", e => {
    e.preventDefault(); const code = conIn.value.trim(); if (!code || !S.ready) return;
    hist.push(code); hi = hist.length; conOut(">>> " + code, "in"); conIn.value = "";
    worker.postMessage({ type: "exec", code: "import numpy as np\n" + code });
  });
  conIn.addEventListener("keydown", e => {
    if (e.key === "ArrowUp" && hi > 0) { conIn.value = hist[--hi]; e.preventDefault(); }
    if (e.key === "ArrowDown") { hi = Math.min(hi + 1, hist.length); conIn.value = hist[hi] || ""; e.preventDefault(); }
  });
  $("#btn-src").addEventListener("click", e => { e.preventDefault(); e.stopPropagation(); if (S.coreSource) $("#dlg-src").showModal(); else toast("Python lädt noch …"); });
  $("#btn-help").addEventListener("click", () => $("#dlg-help").showModal());
  $$("[data-close]").forEach(b => b.addEventListener("click", () => b.closest("dialog").close()));
  $$("dialog").forEach(d => d.addEventListener("click", e => { if (e.target === d) d.close(); }));

  // ------------------------------------------------------------ Toasts
  function toast(text, kind = "") {
    const t = document.createElement("div"); t.className = "st-toast " + kind; t.textContent = text;
    $("#toasts").append(t); setTimeout(() => { t.style.transition = "opacity .4s"; t.style.opacity = 0; setTimeout(() => t.remove(), 400); }, 3200);
  }

  // ------------------------------------------------------------ Bedienung
  $("#btn-play").addEventListener("click", () => setRunning(!S.running));
  $("#btn-step").addEventListener("click", singleStep);
  $("#btn-reset").addEventListener("click", () => rebuild(false));
  const spf = $("#spf");
  spf.addEventListener("input", () => { S.spf = +spf.value; $("#spf-out").textContent = S.spf; spf.style.setProperty("--p", ((spf.value - spf.min) / (spf.max - spf.min) * 100) + "%"); updateStatus(); });
  spf.dispatchEvent(new Event("input"));

  addEventListener("keydown", e => {
    if (e.target.closest("input, select, textarea") || e.ctrlKey || e.metaKey || e.altKey) return;
    if (document.querySelector("dialog[open]")) return;
    const tools = ["info", "sonde", "rechteck", "ellipse", "punkt", "loeschen"];
    if (e.key === " ") { e.preventDefault(); setRunning(!S.running); }
    else if (e.key === ".") singleStep();
    else if (e.key === "r" || e.key === "R") rebuild(false);
    else if (e.key === "?") $("#dlg-help").showModal();
    else if (/^[1-6]$/.test(e.key)) setTool(tools[+e.key - 1]);
  });

  // ------------------------------------------------------------ Start
  renderMaterial(); renderAnzeige();
  const start = "Streuung am PEC-Zylinder";
  S.scene = F.PRESETS[start](); sel.value = start;
  renderGitter(); renderLists(); kpi();
  startWorker();
})();
