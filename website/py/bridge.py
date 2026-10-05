"""bridge.py - verbindet das Web-Studio mit dem Original-Kern fdtd_core.py.

Laeuft in Pyodide (Python/NumPy als WebAssembly im Browser). Die Szene kommt
als JSON im selben Format wie beim Desktop-Studio (Szene.als_dict). Die
Funktionen unten sind aus fdtd_studio.py uebernommen, ohne Tk/matplotlib.
"""
import json
import numpy as np
from fdtd_core import FDTD2D, C0, EPS0, MU0

ALTNAMEN = {"rampe": "sinus", "moduliert": "ricker", "nadel": "ricker"}


def mach_wellenform(art, f0, amp=1.0, zyklen=3.0, verzoegerung=0.0):
    art = ALTNAMEN.get(art, art)
    f0 = max(float(f0), 1.0)
    T = 1.0 / f0
    tv = float(verzoegerung) * T
    if art == "ricker":
        a, t0 = (np.pi * f0) ** 2, 1.2 * T + tv
        return lambda t: amp * (1.0 - 2 * a * (t - t0) ** 2) * np.exp(-a * (t - t0) ** 2)
    if art == "gauss":
        a, t0 = (np.pi * f0) ** 2, 1.2 * T + tv
        return lambda t: amp * np.exp(-a * (t - t0) ** 2)
    tr = max(float(zyklen), 0.0) * T

    def welle(t):
        if tr <= 0.0:
            return amp * np.sin(2 * np.pi * f0 * (t - tv)) * (t >= tv)
        x = np.clip((t - tv) / tr, 0.0, 1.0)
        return amp * 0.5 * (1.0 - np.cos(np.pi * x)) * np.sin(2 * np.pi * f0 * (t - tv)) * (t >= tv)
    return welle


def maske(o, nx, ny):
    i0, i1 = sorted((int(o["i0"]), int(o["i1"])))
    j0, j1 = sorted((int(o["j0"]), int(o["j1"])))
    i0, j0 = max(i0, 0), max(j0, 0)
    i1, j1 = min(max(i1, i0 + 1), nx), min(max(j1, j0 + 1), ny)
    ii, jj = np.ogrid[:nx, :ny]
    if o["form"] == "ellipse":
        ci, cj = 0.5 * (i0 + i1 - 1), 0.5 * (j0 + j1 - 1)
        ri, rj = max(0.5 * (i1 - i0), 0.5), max(0.5 * (j1 - j0), 0.5)
        return ((ii - ci) / ri) ** 2 + ((jj - cj) / rj) ** 2 <= 1.0
    return (ii >= i0) & (ii < i1) & (jj >= j0) & (jj < j1)


def setze_pec(sim, pec):
    for i in np.flatnonzero(pec.any(axis=1)):
        spalte = pec[i]
        j = 0
        while j < spalte.size:
            if spalte[j]:
                k = j
                while k < spalte.size and spalte[k]:
                    k += 1
                sim.set_pec_box(int(i), int(i) + 1, int(j), int(k))
                j = k
            else:
                j += 1


class LinienAnregung:
    def __init__(self, sim, idx_i, idx_j, gewicht, verzug, puls, modus):
        self.sim, self.idx_i, self.idx_j = sim, idx_i, idx_j
        self.gewicht, self.verzug, self.puls, self.modus = gewicht, verzug, puls, modus
        self.i, self.j = int(idx_i[0]), int(idx_j[0])

    def __call__(self, t):
        ziel = self.sim.Jz if self.modus == "strom" else self.sim.Ez
        ziel[self.idx_i, self.idx_j] += self.gewicht * self.puls(t - self.verzug)
        return 0.0


def linien_punkte(nx, ny, i_src, theta_grad, arm, richtung, f0, d, taper=0):
    th = np.deg2rad(float(theta_grad))
    kx, ky = np.cos(th), np.sin(th)
    versatz = [0]
    if richtung:
        lam = C0 / max(f0, 1.0)
        versatz.append(max(1, int(round(lam / (4.0 * max(abs(kx), 0.25) * d)))))
    alle_i, alle_j, alle_g = [], [], []
    for v in versatz:
        i0 = int(np.clip(i_src + v, 0, nx - 1))
        j = np.arange(ny)
        i = np.full(ny, i0)
        g = np.ones(ny)
        if taper > 0:
            r = np.clip(np.minimum(j, ny - 1 - j) / float(taper), 0.0, 1.0)
            g = np.sin(0.5 * np.pi * r) ** 2
        if arm and abs(ky) > 1e-9:
            j_kante = 0 if ky > 0 else ny - 1
            i_arm = np.arange(i0, nx)
            i = np.concatenate([i, i_arm])
            j = np.concatenate([j, np.full(i_arm.size, j_kante)])
            g = np.concatenate([g, np.ones(i_arm.size)])
        alle_i.append(i); alle_j.append(j); alle_g.append(g)
    i, j, g = np.concatenate(alle_i), np.concatenate(alle_j), np.concatenate(alle_g)
    _, einmalig = np.unique(i * ny + j, return_index=True)
    i, j, g = i[einmalig], j[einmalig], g[einmalig]
    verzug = (i * kx + j * ky) * d / C0
    return i, j, g, verzug - verzug.min()


def registriere_quelle(sim, q):
    puls = mach_wellenform(q["art"], q["f0"], q["amp"], q["zyklen"], q.get("verzoegerung", 0.0))
    if q["typ"] == "punkt":
        i = int(np.clip(q["i"], 0, sim.nx - 1))
        j = int(np.clip(q["j"], 0, sim.ny - 1))
        if q["modus"] == "weich":
            sim.add_soft_source(i, j, puls)
        else:
            sim.add_current_source(i, j, puls)
        return
    i, j, g, verzug = linien_punkte(sim.nx, sim.ny, q["i"], q["theta"], q["arm"],
                                    q["richtung"], q["f0"], sim.d, q.get("taper", 0))
    anreger = LinienAnregung(sim, i, j, g, verzug, puls, q["modus"])
    if q["modus"] == "weich":
        sim.add_soft_source(anreger.i, anreger.j, anreger)
    else:
        sim.add_current_source(anreger.i, anreger.j, anreger, amplitude=0.0)


class Lauf:
    """Eine laufende Simulation. Das Web-Studio ruft nur diese Methoden auf."""

    def __init__(self, szene_json):
        sz = json.loads(szene_json)
        nx, ny = sz["nx"], sz["ny"]
        eps = np.ones((nx, ny)); sig = np.zeros((nx, ny)); pec = np.zeros((nx, ny), dtype=bool)
        for o in sz["objekte"]:
            if not o.get("aktiv", True):
                continue
            m = maske(o, nx, ny)
            if o["art"] == "pec":
                pec |= m
            else:
                eps[m] = max(o["eps_r"], 1e-3)
                sig[m] = max(o["sigma"], 0.0)
        self.eps, self.sig, self.pec = eps, sig, pec
        self.sim = FDTD2D(nx, ny, sz["d"], courant=sz["courant"], eps_r=eps, sigma=sig,
                          npml=sz["npml"], pml_m=sz["pml_m"], r0=sz["r0"])
        setze_pec(self.sim, pec)
        for q in sz["quellen"]:
            if q.get("aktiv", True):
                registriere_quelle(self.sim, q)
        for s in sz["sonden"]:
            self.sim.add_probe(int(np.clip(s["i"], 0, nx - 1)), int(np.clip(s["j"], 0, ny - 1)))
        self.intens = np.zeros((nx, ny))
        self.intens_n = 0
        self.w = 0.0
        self.t_w = 0.0

    def karten(self):
        """eps_r und PEC-Maske flach (float32 / uint8) fuer die Anzeige."""
        return self.eps.astype(np.float32).ravel(), self.pec.astype(np.uint8).ravel()

    def schritte(self, n, intensitaet=True):
        """n Zeitschritte. Die Energie wird wie im Desktop-Studio zeitlich
        verschraenkt aus den Feldern vor dem ersten Schritt gebildet."""
        s = self.sim
        ez, hx, hy = s.Ez.copy(), s.Hx.copy(), s.Hy.copy()
        self.t_w = s.t
        for k in range(max(int(n), 1)):
            s.step()
            if k == 0:
                self.w = (0.5 * EPS0 * float(np.sum(self.eps * ez * ez))
                          + 0.5 * MU0 * (float(np.sum(hx * s.Hx)) + float(np.sum(hy * s.Hy)))) * s.d ** 2
            if intensitaet:
                self.intens += s.Ez * s.Ez
                self.intens_n += 1
        return bool(np.isfinite(s.Ez).all() and np.abs(s.Ez).max() < 1e20)

    def ez(self):
        return self.sim.Ez.astype(np.float32).ravel()

    def intensitaet(self):
        if self.intens_n == 0:
            return np.zeros(self.sim.nx * self.sim.ny, dtype=np.float32)
        return (self.intens / self.intens_n).astype(np.float32).ravel()

    def sonde(self, k, ab=0):
        p = self.sim._probes[k]
        return np.asarray(p["Ez"][int(ab):], dtype=np.float64)

    def spektrum(self, k, fenster=0):
        """Betragsspektrum einer Sonde: Hann-Fenster, bei langen Reihen nach
        Welch gemittelt (wie Studio._spektrum)."""
        y = np.asarray(self.sim._probes[k]["Ez"], dtype=float)
        if fenster and y.size > fenster:
            y = y[-int(fenster):]
        n = y.size
        dt = self.sim.dt
        if n < 64:
            return np.zeros(0), np.zeros(0)
        if n < 512:
            w = np.hanning(n)
            return np.fft.rfftfreq(n, dt), np.abs(np.fft.rfft(y * w)) * 2.0 / w.sum()
        laenge = max(128, n // 8)
        laenge -= laenge % 2
        w = np.hanning(laenge)
        summe = np.zeros(laenge // 2 + 1)
        anzahl = 0
        for a in range(0, n - laenge + 1, laenge // 2):
            summe += np.abs(np.fft.rfft(y[a:a + laenge] * w)) ** 2
            anzahl += 1
        return np.fft.rfftfreq(laenge, dt), np.sqrt(summe / max(anzahl, 1)) * 2.0 / w.sum()
