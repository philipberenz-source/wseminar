/* py-worker.js — fuehrt den Python-Kern (fdtd_core.py + bridge.py) mit
 * Pyodide in einem eigenen Thread aus. Die Oberflaeche schickt Befehle,
 * der Worker antwortet mit Feldern als Float32Array (transferiert, ohne Kopie).
 */
const PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/";
importScripts(PYODIDE + "pyodide.js");

let py = null, lauf = null, probeSeen = [], gen = 0;

const post = (msg, transfer) => { msg.gen = gen; self.postMessage(msg, transfer || []); };
const status = text => post({ type: "status", text });

function buf(proxy, Ctor) {
  // numpy -> JS: Speicher direkt ueber das Buffer-Protokoll lesen, dann kopieren.
  const b = proxy.getBuffer();
  const out = new Ctor(b.data.length); out.set(b.data);
  b.release(); proxy.destroy();
  return out;
}

async function init() {
  status("Lade Python-Laufzeit (Pyodide) …");
  py = await loadPyodide({ indexURL: PYODIDE });
  status("Lade NumPy …");
  await py.loadPackage("numpy");
  status("Lade fdtd_core.py …");
  for (const f of ["fdtd_core.py", "bridge.py"]) {
    const r = await fetch("../../py/" + f + "?v=" + Date.now());
    if (!r.ok) throw new Error(f + " nicht gefunden (" + r.status + ")");
    py.FS.writeFile(f, await r.text());
  }
  py.runPython("import sys; sys.path.insert(0, '.'); import bridge, numpy");
  const ver = py.runPython("import sys, numpy; f'Python {sys.version.split()[0]} · NumPy {numpy.__version__}'");
  const src = py.FS.readFile("fdtd_core.py", { encoding: "utf8" });
  post({ type: "ready", version: ver, coreSource: src });
}

function build(sceneJson) {
  if (lauf) { lauf.destroy(); lauf = null; }
  const Lauf = py.globals.get("bridge").Lauf;
  lauf = Lauf(sceneJson);
  const k = lauf.karten();
  const eps = buf(k.get(0), Float32Array), pec = buf(k.get(1), Uint8Array);
  k.destroy();
  const sim = lauf.sim;
  const probesProxy = sim._probes;
  const info = { nx: sim.nx, ny: sim.ny, dt: sim.dt, d: sim.d, nProbes: probesProxy.length };
  probesProxy.destroy(); sim.destroy();
  probeSeen = new Array(info.nProbes).fill(0);
  post({ type: "built", info, eps, pec, ez: buf(lauf.ez(), Float32Array) }, [eps.buffer, pec.buffer]);
}

function run({ n, intens, wantIntens }) {
  const t0 = performance.now();
  const ok = lauf.schritte(n, intens);
  const ms = performance.now() - t0;
  const sim = lauf.sim;
  const ez = buf(lauf.ez(), Float32Array);
  const transfer = [ez.buffer];
  let it = null;
  if (wantIntens) { it = buf(lauf.intensitaet(), Float32Array); transfer.push(it.buffer); }
  // Nur neue Sondenwerte schicken.
  const probes = [];
  for (let k = 0; k < probeSeen.length; k++) {
    const y = buf(lauf.sonde(k, probeSeen[k]), Float64Array);
    probeSeen[k] += y.length; probes.push(y); transfer.push(y.buffer);
  }
  const t = sim.t, step = sim.time_step; sim.destroy();
  post({ type: "frame", ok, ez, intens: it, probes, t, step, w: lauf.w, tw: lauf.t_w, n, ms }, transfer);
}

function spectrum({ k, fenster }) {
  const r = lauf.spektrum(k, fenster || 0);
  const f = buf(r.get(0), Float64Array), mag = buf(r.get(1), Float64Array);
  r.destroy();
  post({ type: "spectrum", k, f, mag }, [f.buffer, mag.buffer]);
}

self.onmessage = async e => {
  const m = e.data;
  if (m.gen !== undefined) gen = m.gen;
  try {
    if (m.type === "init") await init();
    else if (m.type === "build") build(m.scene);
    else if (m.type === "run") run(m);
    else if (m.type === "spectrum") spectrum(m);
    else if (m.type === "exec") {            // Python-Konsole
      let out = "";
      py.setStdout({ batched: s => { out += s + "\n"; } });
      py.setStderr({ batched: s => { out += s + "\n"; } });
      py.globals.set("lauf", lauf);
      let res;
      try { res = await py.runPythonAsync(m.code); }
      catch (err) { out += String(err.message || err).split("\n").slice(-4).join("\n"); }
      if (res !== undefined && res !== null) { out += String(res); if (res.destroy) res.destroy(); }
      post({ type: "exec", out: out.trimEnd() });
    }
  } catch (err) {
    post({ type: "error", text: String(err && err.message || err), during: m.type });
  }
};
