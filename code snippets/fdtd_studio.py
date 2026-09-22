"""fdtd_studio.py - Interaktive Oberflaeche fuer den FDTD-Kern.

Start:  python fdtd_studio.py

Der Kern (fdtd_core.FDTD2D) wird ausschliesslich importiert und ueber seine
oeffentliche Schnittstelle benutzt; diese Datei veraendert ihn nicht.

Aufbau
------
  Wellenformen : Sinus (Dauerstrich), Ricker-Puls, Gauss-Puls
  Szene        : Gitter + Objekte (Rechteck/Ellipse, Dielektrikum/PEC)
                 + Quellen (Punkt, ebene Welle) + Sonden  ->  baut ein FDTD2D
  Motor        : rechnet im Hintergrund-Thread, liefert Feldkopien
  Studio       : Tk-Oberflaeche mit Live-Animation, Analyse- und
                 Intensitaetsdiagrammen, Maus-Editor, Messwerkzeugen, Export

Energie
-------
Die angezeigte Feldenergie je Meter Tiefe wird zeitlich verschraenkt gebildet,

    W^n = Delta^2 * sum( eps0 eps_r (E_z^n)^2 / 2
                         + mu0 ( H^(n-1/2) . H^(n+1/2) ) / 2 ) ,

weil E und H im Yee-Schema einen halben Zeitschritt auseinanderliegen. So ist
W ohne Verluste exakt erhalten (Probe im PEC-Resonator: Schwankung 5e-16).
sim.total_energy() aus dem Kern nimmt dagegen H^(n+1/2) quadratisch zu E^(n+1)
und schwankt dadurch um einige Prozent mit der Schwingung.

Ebene Wellen
------------
Eine ebene Welle wird aus sehr vielen Punktquellen auf einer Gitterlinie
zusammengesetzt. Die Linie laeuft ueber die *gesamte* Kante, also bis in die
PML hinein - dadurch hat sie keine Enden, an denen gebeugt werden koennte.
Bei schraegem Einfall werden die Punkte zeitversetzt angeregt:

    t_d(i,j) = ( i*cos(theta) + j*sin(theta) ) * Delta / c ,

zusaetzlich kann ein Arm entlang der Kante angeregt werden, aus der die Welle
kommt (fuellt die Ecke), und ein zweites Blatt im Abstand lambda/(4 cos theta)
loescht die Rueckwaertswelle aus (Zweizeilen-Array).
"""

from __future__ import annotations

import copy
import json
import sys
import threading
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, Ellipse
from matplotlib.colors import Normalize

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fdtd_core import FDTD2D, C0, EPS0, MU0            # noqa: E402

# Als exe (PyInstaller) liegt __file__ im Auspackverzeichnis; fuer
# Dateidialoge ist der Ordner der exe gemeint.
if getattr(sys, "frozen", False):
    ORDNER = Path(sys.executable).resolve().parent
else:
    ORDNER = Path(__file__).resolve().parent
ABB = ORDNER.parent / "seminararbeit" / "abbildungen"

FARBE = {
    "bg": "#1b1e24", "panel": "#232830", "feld": "#12151a", "linie": "#39404d",
    "text": "#e3e7ee", "grau": "#98a2b3", "akzent": "#4aa3df",
    "gut": "#4cc38a", "warn": "#e5c07b", "fehler": "#e06c75",
}

CMAPS = ["RdBu_r", "seismic", "bwr", "coolwarm", "twilight_shifted",
         "inferno", "magma", "viridis", "turbo", "gray"]

DARSTELLUNGEN = ["Ez", "|Ez|", "Ez in dB", "Intensitaet <Ez^2>",
                 "Intensitaet in dB", "el. Energiedichte"]

VORZEICHENBEHAFTET = {"Ez"}

# Zeitausschnitt der Verlaufsdiagramme. Ohne Fenster wird ein langer Lauf zu
# einem ausgefuellten Balken, weil je Bildpunkt hunderte Schwingungen liegen.
ZEITFENSTER = ["letzte 10 Perioden", "letzte 25 Perioden", "letzte 100 Perioden",
               "letzte 500 Schritte", "letzte 5000 Schritte", "ganzer Lauf"]


# ============================================================== Wellenformen

WELLENFORMEN = {
    "sinus": "Sinus (Dauerstrich)",
    "ricker": "Ricker-Puls",
    "gauss": "Gauss-Puls",
}

# Frueher gespeicherte Szenen koennen noch andere Namen enthalten.
ALTNAMEN = {"rampe": "sinus", "moduliert": "ricker", "nadel": "ricker"}


def mach_wellenform(art, f0, amp=1.0, zyklen=3.0, dt=1e-12, verzoegerung=0.0):
    """Liefert w(t) fuer eine der drei Quellenarten; alle Formen gelten auch
    fuer Arrays von t.

    sinus   Dauerstrich sin(2 pi f0 t). Eingeschaltet wird ueber eine
            Kosinus-Rampe der Laenge `zyklen` Perioden. Ein hart geschalteter
            Sinus (zyklen = 0) ist ein Sprung und regt das ganze Band an -
            genau das macht Spektren unbrauchbar; die Rampe vermeidet das.
    ricker  zweite Ableitung der Gauss-Kurve, Spitzenfrequenz f0,
            mittelwertfrei - der uebliche FDTD-Breitbandpuls.
    gauss   reine Gauss-Kurve, Bandbreite ueber f0 gesetzt. Enthaelt einen
            Gleichanteil, der als statisches Restfeld stehen bleibt.
    """
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

    if art == "sinus":
        tr = max(float(zyklen), 0.0) * T

        def welle(t):
            if tr <= 0.0:
                return amp * np.sin(2 * np.pi * f0 * (t - tv)) * (t >= tv)
            x = np.clip((t - tv) / tr, 0.0, 1.0)
            rampe = 0.5 * (1.0 - np.cos(np.pi * x))
            return amp * rampe * np.sin(2 * np.pi * f0 * (t - tv)) * (t >= tv)
        return welle

    raise ValueError(f"unbekannte Wellenform: {art}")


# ==================================================================== Szene

@dataclass
class Objekt:
    """Rechteck oder einbeschriebene Ellipse, aus Material oder PEC."""
    form: str = "rechteck"            # rechteck | ellipse
    art: str = "dielektrikum"         # dielektrikum | pec
    i0: int = 0
    j0: int = 0
    i1: int = 10
    j1: int = 10
    eps_r: float = 2.25
    sigma: float = 0.0
    aktiv: bool = True

    def box(self, nx, ny):
        i0, i1 = sorted((int(self.i0), int(self.i1)))
        j0, j1 = sorted((int(self.j0), int(self.j1)))
        i0, j0 = max(i0, 0), max(j0, 0)
        return i0, j0, min(max(i1, i0 + 1), nx), min(max(j1, j0 + 1), ny)

    def maske(self, nx, ny):
        i0, j0, i1, j1 = self.box(nx, ny)
        ii, jj = np.ogrid[:nx, :ny]
        if self.form == "ellipse":
            ci, cj = 0.5 * (i0 + i1 - 1), 0.5 * (j0 + j1 - 1)
            ri, rj = max(0.5 * (i1 - i0), 0.5), max(0.5 * (j1 - j0), 0.5)
            return ((ii - ci) / ri) ** 2 + ((jj - cj) / rj) ** 2 <= 1.0
        return (ii >= i0) & (ii < i1) & (jj >= j0) & (jj < j1)

    def enthaelt(self, i, j, nx, ny):
        return bool(self.maske(nx, ny)[int(i), int(j)])

    def zeile(self):
        form = "Ellipse" if self.form == "ellipse" else "Rechteck"
        stoff = "PEC" if self.art == "pec" else f"eps_r = {self.eps_r:g}"
        if self.art != "pec" and self.sigma:
            stoff += f" | sigma = {self.sigma:g}"
        return (("" if self.aktiv else "(aus) ") + form,
                f"i {self.i0}..{self.i1}  j {self.j0}..{self.j1}", stoff)


@dataclass
class Quelle:
    typ: str = "ebene"                # punkt | ebene
    art: str = "sinus"                # Wellenform
    i: int = 30
    j: int = 100
    f0: float = 15e9
    amp: float = 1.0
    zyklen: float = 4.0
    theta: float = 0.0                # Grad; nur bei ebener Welle
    modus: str = "strom"              # strom (J_z) | weich (E_z)
    taper: int = 0                    # Kantenglaettung in Zellen (0 = aus)
    arm: bool = True                  # Arm entlang der Einfallskante
    richtung: bool = False            # lambda/4-Paar gegen Rueckwaertswelle
    verzoegerung: float = 0.0         # Startverzoegerung in Perioden
    aktiv: bool = True

    def __post_init__(self):
        self.art = ALTNAMEN.get(self.art, self.art)

    def zeile(self):
        typ = "Ebene Welle" if self.typ == "ebene" else "Punktquelle"
        ort = f"Linie i = {self.i}" if self.typ == "ebene" else f"({self.i}, {self.j})"
        name = WELLENFORMEN.get(self.art, self.art).split(" (")[0]
        zus = f"{name}, {self.f0 / 1e9:g} GHz"
        if self.typ == "ebene":
            zus += f", {self.theta:g} Grad"
        return ("" if self.aktiv else "(aus) ") + typ, ort, zus


@dataclass
class Sonde:
    i: int = 100
    j: int = 100
    name: str = ""


@dataclass
class Szene:
    nx: int = 300
    ny: int = 240
    d: float = 1e-3
    courant: float = 0.99
    npml: int = 12
    pml_m: int = 3
    r0: float = 1e-6
    schirm: int = 0                   # Spalte fuer das Schirmprofil (0 = automatisch)
    objekte: list = field(default_factory=list)
    quellen: list = field(default_factory=list)
    sonden: list = field(default_factory=list)

    def karten(self):
        eps = np.ones((self.nx, self.ny))
        sig = np.zeros((self.nx, self.ny))
        pec = np.zeros((self.nx, self.ny), dtype=bool)
        for o in self.objekte:
            if not o.aktiv:
                continue
            m = o.maske(self.nx, self.ny)
            if o.art == "pec":
                pec |= m
            else:
                eps[m] = max(o.eps_r, 1e-3)
                sig[m] = max(o.sigma, 0.0)
        return eps, sig, pec

    def baue(self):
        eps, sig, pec = self.karten()
        sim = FDTD2D(self.nx, self.ny, self.d, courant=self.courant, eps_r=eps,
                     sigma=sig, npml=self.npml, pml_m=self.pml_m, r0=self.r0)
        setze_pec(sim, pec)
        for q in self.quellen:
            if q.aktiv:
                registriere_quelle(sim, q)
        for s in self.sonden:
            sim.add_probe(int(np.clip(s.i, 0, self.nx - 1)),
                          int(np.clip(s.j, 0, self.ny - 1)))
        return sim, eps, sig, pec

    # -- Kennzahlen -------------------------------------------------------
    def kennzahlen(self):
        dt = self.courant * self.d / (C0 * np.sqrt(2.0))
        f0 = self.quellen[0].f0 if self.quellen else 0.0
        n_lambda = C0 / (f0 * self.d) if f0 > 0 else float("inf")
        return {
            "dt": dt,
            "schritte_quer": int(round(max(self.nx, self.ny) * self.d / (C0 * dt))),
            "breite": self.nx * self.d,
            "hoehe": self.ny * self.d,
            "n_lambda": n_lambda,
            "f0": f0,
        }

    def als_dict(self):
        d = asdict(self)
        return d

    @staticmethod
    def aus_dict(d):
        s = Szene(**{k: v for k, v in d.items()
                     if k not in ("objekte", "quellen", "sonden")})
        s.objekte = [Objekt(**o) for o in d.get("objekte", [])]
        s.quellen = [Quelle(**q) for q in d.get("quellen", [])]
        s.sonden = [Sonde(**p) for p in d.get("sonden", [])]
        return s


def setze_pec(sim, pec):
    """Traegt eine beliebige PEC-Maske ueber die oeffentliche Methode
    set_pec_box ein - je Spalte werden die zusammenhaengenden Laeufe gesetzt."""
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


def sondendaten(sim, k):
    """Liest die vom Kern aufgezeichneten Sondenwerte (nur lesend)."""
    p = sim._probes[k]
    return np.asarray(p["t"], dtype=float), np.asarray(p["Ez"], dtype=float)


# ===================================================== Quellen registrieren

class LinienAnregung:
    """Traegt eine ganze Quellenlinie in einem Rutsch ein.

    Der Kern ruft je Zeitschritt jede registrierte Wellenform genau einmal
    auf. Diese "Wellenform" schreibt dabei alle Punkte der Linie vektorisiert
    in J_z bzw. E_z und gibt 0 zurueck, damit der Einzelpunkt-Beitrag des
    Kerns wegfaellt. So kosten auch 600 Punktquellen nur einen Array-Zugriff
    pro Zeitschritt.
    """

    def __init__(self, sim, idx_i, idx_j, gewicht, verzug, puls, modus):
        self.sim = sim
        self.idx_i = idx_i
        self.idx_j = idx_j
        self.gewicht = gewicht
        self.verzug = verzug
        self.puls = puls
        self.modus = modus
        self.i, self.j = int(idx_i[0]), int(idx_j[0])

    def __call__(self, t):
        ziel = self.sim.Jz if self.modus == "strom" else self.sim.Ez
        ziel[self.idx_i, self.idx_j] += self.gewicht * self.puls(t - self.verzug)
        return 0.0


def linien_punkte(nx, ny, i_src, theta_grad, arm, richtung, f0, d, taper=0):
    """Gitterpunkte, Gewichte und Zeitversatz einer ebenen Welle.

    Ausbreitungsrichtung k = (cos theta, sin theta); positives theta laeuft
    also nach schraeg oben. Die Startlinie belegt alle j (auch in der PML),
    der optionale Arm belegt die Kante, aus der die Welle kommt.
    """
    th = np.deg2rad(float(theta_grad))
    kx, ky = np.cos(th), np.sin(th)

    versatz = [0]
    if richtung:                                   # zweites Blatt: lambda/4
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
        alle_i.append(i)
        alle_j.append(j)
        alle_g.append(g)

    i = np.concatenate(alle_i)
    j = np.concatenate(alle_j)
    g = np.concatenate(alle_g)

    # Doppelte Gitterpunkte entfernen: Linie und Arm teilen sich eine Zelle,
    # und bei += mit doppelten Indizes zaehlt numpy nur einmal.
    _, einmalig = np.unique(i * ny + j, return_index=True)
    i, j, g = i[einmalig], j[einmalig], g[einmalig]

    verzug = (i * kx + j * ky) * d / C0
    return i, j, g, verzug - verzug.min()


def registriere_quelle(sim, q):
    puls = mach_wellenform(q.art, q.f0, q.amp, q.zyklen, sim.dt, q.verzoegerung)

    if q.typ == "punkt":
        i = int(np.clip(q.i, 0, sim.nx - 1))
        j = int(np.clip(q.j, 0, sim.ny - 1))
        if q.modus == "weich":
            sim.add_soft_source(i, j, puls)
        else:
            sim.add_current_source(i, j, puls)
        return 1

    i, j, g, verzug = linien_punkte(sim.nx, sim.ny, q.i, q.theta, q.arm,
                                    q.richtung, q.f0, sim.d, q.taper)
    anreger = LinienAnregung(sim, i, j, g, verzug, puls, q.modus)
    if q.modus == "weich":
        sim.add_soft_source(anreger.i, anreger.j, anreger)
    else:
        sim.add_current_source(anreger.i, anreger.j, anreger, amplitude=0.0)
    return i.size


# ================================================================== Presets

def _wand(i0, i1, ny, schlitze):
    """PEC-Wand von i0 bis i1 mit Schlitzen [(j_von, j_bis), ...]."""
    kanten = [0]
    for a, b in schlitze:
        kanten += [a, b]
    kanten.append(ny)
    bloecke = []
    for k in range(0, len(kanten), 2):
        j0, j1 = kanten[k], kanten[k + 1]
        if j1 > j0:
            bloecke.append(Objekt("rechteck", "pec", i0, j0, i1, j1))
    return bloecke


def preset_freiraum():
    sz = Szene(nx=320, ny=260)
    sz.quellen = [Quelle(typ="ebene", art="ricker", i=30, f0=15e9, theta=0.0)]
    sz.sonden = [Sonde(120, 130), Sonde(220, 130)]
    return sz


def preset_schraege_welle():
    sz = Szene(nx=320, ny=320)
    sz.quellen = [Quelle(typ="ebene", art="sinus", i=25, f0=15e9, theta=30.0,
                         arm=True, richtung=True, zyklen=4)]
    sz.sonden = [Sonde(200, 200)]
    return sz


def preset_doppelspalt():
    """Spaltabstand 3 lambda, Spaltbreite lambda/2: Maxima 0. und 1. Ordnung
    liegen bei j = 150 und j = 150 +- 52, die 2. Ordnung bei j = 150 +- 131."""
    sz = Szene(nx=340, ny=300, schirm=280)
    sz.quellen = [Quelle(typ="ebene", art="sinus", i=30, f0=15e9, zyklen=4)]
    sz.objekte = _wand(130, 134, sz.ny, [(115, 125), (175, 185)])
    sz.sonden = [Sonde(280, 150), Sonde(280, 202)]
    return sz


def preset_einzelspalt():
    """Spaltbreite 2 lambda: erstes Minimum bei sin theta = 1/2."""
    sz = Szene(nx=340, ny=300, schirm=280)
    sz.quellen = [Quelle(typ="ebene", art="sinus", i=30, f0=15e9, zyklen=4)]
    sz.objekte = _wand(130, 134, sz.ny, [(130, 170)])
    sz.sonden = [Sonde(280, 150)]
    return sz


def preset_gitter():
    """Fuenf Spalte im Abstand 3 lambda: scharfe Ordnungen bei 19.5 und 41.8 Grad."""
    sz = Szene(nx=360, ny=320, schirm=300)
    sz.quellen = [Quelle(typ="ebene", art="sinus", i=30, f0=15e9, zyklen=4)]
    schlitze = [(j - 5, j + 5) for j in range(40, 281, 60)]
    sz.objekte = _wand(140, 144, sz.ny, schlitze)
    sz.sonden = [Sonde(300, 160)]
    return sz


def preset_brechung():
    sz = Szene(nx=360, ny=300)
    sz.quellen = [Quelle(typ="ebene", art="sinus", i=25, f0=15e9, theta=30.0,
                         arm=True, richtung=True, zyklen=4)]
    sz.objekte = [Objekt("rechteck", "dielektrikum", 190, 0, 360, 300, eps_r=2.25)]
    sz.sonden = [Sonde(120, 150), Sonde(280, 150)]
    return sz


def preset_totalreflexion():
    sz = Szene(nx=360, ny=300)
    sz.quellen = [Quelle(typ="ebene", art="sinus", i=25, f0=15e9, theta=50.0,
                         arm=True, richtung=True, zyklen=4)]
    sz.objekte = [Objekt("rechteck", "dielektrikum", 0, 0, 190, 300, eps_r=2.25)]
    sz.sonden = [Sonde(210, 150)]
    return sz


def preset_linse():
    sz = Szene(nx=360, ny=280)
    sz.quellen = [Quelle(typ="ebene", art="sinus", i=25, f0=15e9, zyklen=4)]
    sz.objekte = [Objekt("ellipse", "dielektrikum", 140, 70, 200, 210, eps_r=2.25)]
    sz.sonden = [Sonde(250, 140), Sonde(300, 140)]
    return sz


def preset_zylinder():
    sz = Szene(nx=340, ny=300)
    sz.quellen = [Quelle(typ="ebene", art="ricker", i=30, f0=15e9)]
    sz.objekte = [Objekt("ellipse", "pec", 150, 120, 200, 180)]
    sz.sonden = [Sonde(280, 150), Sonde(80, 150)]
    return sz


def preset_hohlleiter():
    sz = Szene(nx=360, ny=260)
    sz.quellen = [Quelle(typ="punkt", art="sinus", i=60, j=130, f0=15e9, zyklen=3)]
    sz.objekte = [Objekt("rechteck", "pec", 40, 104, 360, 110),
                  Objekt("rechteck", "pec", 40, 150, 360, 156)]
    sz.sonden = [Sonde(200, 130), Sonde(320, 130)]
    return sz


def preset_bragg():
    sz = Szene(nx=380, ny=200)
    sz.quellen = [Quelle(typ="ebene", art="ricker", i=30, f0=7.5e9)]
    i = 140
    for _ in range(8):                       # lambda/4-Schichtstapel bei 7.5 GHz
        sz.objekte.append(Objekt("rechteck", "dielektrikum", i, 0, i + 5, 200, eps_r=4.0))
        i += 15
    sz.sonden = [Sonde(90, 100), Sonde(330, 100)]
    return sz


def preset_punktquelle():
    sz = Szene(nx=280, ny=280)
    sz.quellen = [Quelle(typ="punkt", art="ricker", i=140, j=140, f0=15e9)]
    sz.sonden = [Sonde(200, 140), Sonde(140, 200)]
    return sz


def preset_resonator():
    sz = Szene(nx=141, ny=141, npml=0)
    sz.quellen = [Quelle(typ="punkt", art="ricker", i=70, j=70, f0=15e9)]
    sz.sonden = [Sonde(40, 95)]
    return sz


PRESETS = {
    "Freiraum: ebene Welle": preset_freiraum,
    "Schraege ebene Welle (30 Grad)": preset_schraege_welle,
    "Doppelspalt (Beugung)": preset_doppelspalt,
    "Einzelspalt": preset_einzelspalt,
    "Beugungsgitter (5 Spalte)": preset_gitter,
    "Brechung an Grenzflaeche": preset_brechung,
    "Totalreflexion (50 Grad)": preset_totalreflexion,
    "Dielektrische Linse": preset_linse,
    "Streuung am PEC-Zylinder": preset_zylinder,
    "Parallelplatten-Hohlleiter": preset_hohlleiter,
    "Bragg-Spiegel (Schichtstapel)": preset_bragg,
    "Punktquelle im Freiraum": preset_punktquelle,
    "PEC-Resonator (ohne PML)": preset_resonator,
}


# ==================================================================== Motor

class Motor:
    """Rechnet die Simulation in einem Hintergrund-Thread."""

    FILM_MAX = 600

    def __init__(self, protokoll=None):
        self.lock = threading.RLock()
        self.protokoll = protokoll or (lambda *a, **k: None)

        self.szene = None
        self.sim = None
        self.eps = self.sig = self.pec = None

        self.feld = None
        self.intens = None
        self.intens_n = 0
        self.zeit_h = []
        self.energie_h = []

        self.chunk = 4
        self.max_rate = 0.0          # Schritte/s, 0 = unbegrenzt
        self.grenze = 0              # Stoppschritt, 0 = kein Limit
        self.intens_an = True
        self.rate = 0.0
        self.instabil = False

        self.aufnahme = False
        self.film_jede = 4
        self.filme = []

        self._laeuft = threading.Event()
        self._ende = threading.Event()
        self._thread = threading.Thread(target=self._schleife, daemon=True)
        self._thread.start()

    # -- Steuerung ---------------------------------------------------------
    @property
    def laeuft(self):
        return self._laeuft.is_set()

    def start(self):
        if self.sim is not None and not self.instabil:
            self._laeuft.set()

    def pause(self):
        self._laeuft.clear()

    def beenden(self):
        self._laeuft.clear()
        self._ende.set()

    def setze_szene(self, szene):
        with self.lock:
            self._laeuft.clear()
            self.szene = szene
            self.sim, self.eps, self.sig, self.pec = szene.baue()
            self.feld = self.sim.Ez.copy()
            self.intens = np.zeros((szene.nx, szene.ny))
            self.intens_n = 0
            self.zeit_h, self.energie_h = [], []
            self.instabil = False
            self.filme = []
        return self.sim

    def zuruecksetzen(self):
        if self.szene is not None:
            self.setze_szene(self.szene)

    def einzelschritt(self, n=1):
        with self.lock:
            if self.sim is None:
                return
            for _ in range(n):
                zeit, ez, hx, hy = self._vorher()
                self.sim.step()
                if self.intens_an:
                    self.intens += self.sim.Ez * self.sim.Ez
                    self.intens_n += 1
                self._notiere(zeit, self.energie(ez, hx, hy))
            self.feld = self.sim.Ez.copy()

    # -- interne Schleife --------------------------------------------------
    def _vorher(self):
        """Zeit und Felder vor dem naechsten Schritt (fuer die Energie)."""
        return (self.sim.t, self.sim.Ez.copy(),
                self.sim.Hx.copy(), self.sim.Hy.copy())

    def energie(self, ez_alt, hx_alt, hy_alt):
        """Feldenergie je Meter Tiefe, zeitlich richtig verschraenkt:

            W^n = Delta^2 * sum( eps0 eps_r (E_z^n)^2 / 2
                                 + mu0 ( H^(n-1/2) . H^(n+1/2) ) / 2 )

        E und H liegen im Yee-Schema einen halben Zeitschritt auseinander.
        Der magnetische Anteil wird deshalb aus dem Produkt der beiden
        benachbarten H-Halbschritte gebildet - so ist W ohne Verluste exakt
        erhalten. sim.total_energy() nimmt dagegen H^(n+1/2) quadratisch zu
        E^(n+1) und schwankt dadurch sichtbar mit der Schwingung.

        ez_alt, hx_alt, hy_alt sind die Felder vor dem letzten Schritt.
        """
        flaeche = self.sim.d ** 2
        w_e = 0.5 * EPS0 * float(np.sum(self.eps * ez_alt * ez_alt))
        w_h = 0.5 * MU0 * (float(np.sum(hx_alt * self.sim.Hx))
                           + float(np.sum(hy_alt * self.sim.Hy)))
        return (w_e + w_h) * flaeche

    def _notiere(self, zeit=None, wert=None):
        self.zeit_h.append(self.sim.t if zeit is None else zeit)
        self.energie_h.append(self.sim.total_energy() if wert is None else wert)
        if self.aufnahme and len(self.filme) < self.FILM_MAX:
            if self.sim.time_step % max(self.film_jede, 1) == 0:
                s = max(1, int(np.ceil(max(self.sim.nx, self.sim.ny) / 320)))
                self.filme.append(self.sim.Ez[::s, ::s].copy())

    def _schleife(self):
        while not self._ende.is_set():
            if not self._laeuft.wait(0.05):
                continue
            t0 = time.perf_counter()
            with self.lock:
                if self.sim is None:
                    self._laeuft.clear()
                    continue
                n = 0
                # Nur vor dem ersten Schritt des Blocks kopieren: daraus wird
                # die verschraenkte Energie gebildet (einmal je Block statt
                # je Zeitschritt).
                zeit, ez, hx, hy = self._vorher()
                for _ in range(max(self.chunk, 1)):
                    self.sim.step()
                    n += 1
                    if n == 1:
                        w = self.energie(ez, hx, hy)
                    if self.intens_an:
                        self.intens += self.sim.Ez * self.sim.Ez
                        self.intens_n += 1
                    if self.grenze and self.sim.time_step >= self.grenze:
                        self._laeuft.clear()
                        break
                self.feld = self.sim.Ez.copy()
                self._notiere(zeit, w)
                spitze = float(np.abs(self.feld).max())
            verbraucht = time.perf_counter() - t0
            self.rate = n / verbraucht if verbraucht > 0 else 0.0

            if not np.isfinite(spitze) or spitze > 1e20:
                self._laeuft.clear()
                self.instabil = True
                self.protokoll("Rechnung divergiert - Simulation gestoppt. "
                               "Courant-Zahl S <= 1 waehlen.", "fehler")

            if self.max_rate > 0:
                soll = n / self.max_rate
                if soll > verbraucht:
                    time.sleep(soll - verbraucht)
            else:
                time.sleep(0.0005)          # der Oberflaeche Luft lassen


# =================================================================== Studio

def zahl(var, standard=0.0, ganz=False):
    try:
        w = float(str(var.get()).replace(",", "."))
    except (ValueError, tk.TclError):
        return standard
    if not np.isfinite(w):
        return standard
    return int(round(w)) if ganz else w


class Studio(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("FDTD-Studio  -  2D-TM-Simulation (fdtd_core)")
        self.geometry("1620x980")
        self.minsize(1180, 740)
        self.configure(bg=FARBE["bg"])

        self._stil()
        self._mpl_stil()

        self.motor = Motor(protokoll=self._log_thread)
        self.szene = preset_freiraum()

        self._bg = None
        self._spitze = 1e-12
        self._letzte_cb = 0.0
        self._letzte_analyse = 0.0
        self._letzte_intens = 0.0
        self._zieh_start = None
        self._meldungen = []

        self._variablen()
        self._menue()
        self._werkzeugleiste()

        rumpf = ttk.Frame(self, style="Fenster.TFrame")
        rumpf.pack(fill="both", expand=True, padx=6, pady=(0, 4))
        self._linkes_panel(rumpf)
        self._rechtes_panel(rumpf)
        self._statusleiste()

        self._tastatur()
        self.protocol("WM_DELETE_WINDOW", self._beenden)

        self._szene_anwenden(self.szene, "Preset 'Freiraum: ebene Welle' geladen.")
        self.after(30, self._takt)

    # ------------------------------------------------------------ Aussehen
    def _stil(self):
        st = ttk.Style(self)
        st.theme_use("clam")
        p, t, g, f = FARBE["panel"], FARBE["text"], FARBE["grau"], FARBE["feld"]
        st.configure(".", background=p, foreground=t, fieldbackground=f,
                     bordercolor=FARBE["linie"], lightcolor=p, darkcolor=p,
                     focuscolor=FARBE["akzent"])
        st.configure("TFrame", background=p)
        st.configure("Fenster.TFrame", background=FARBE["bg"])
        st.configure("Leiste.TFrame", background=FARBE["bg"])
        st.configure("TLabel", background=p, foreground=t, font=("Segoe UI", 9))
        st.configure("Leiste.TLabel", background=FARBE["bg"], foreground=g)
        st.configure("Kopf.TLabel", foreground=FARBE["akzent"],
                     font=("Segoe UI Semibold", 10))
        st.configure("Grau.TLabel", foreground=g)
        st.configure("Wert.TLabel", foreground=FARBE["gut"],
                     font=("Consolas", 9))
        st.configure("TButton", background="#2e3540", foreground=t,
                     borderwidth=0, padding=(9, 5), font=("Segoe UI", 9))
        st.map("TButton", background=[("active", "#3b4350"), ("pressed", "#1d222b")])
        st.configure("Akzent.TButton", background=FARBE["akzent"],
                     foreground="#10131a", font=("Segoe UI Semibold", 9))
        st.map("Akzent.TButton", background=[("active", "#66b6ea")])
        st.configure("Klein.TButton", padding=(5, 3), font=("Segoe UI", 8))
        st.configure("TEntry", fieldbackground=f, foreground=t,
                     insertcolor=t, borderwidth=1, padding=2)
        st.configure("TSpinbox", fieldbackground=f, foreground=t,
                     arrowcolor=t, borderwidth=1)
        st.configure("TCombobox", fieldbackground=f, background="#2e3540",
                     foreground=t, arrowcolor=t, borderwidth=1, padding=2)
        st.map("TCombobox", fieldbackground=[("readonly", f)],
               foreground=[("readonly", t)], background=[("readonly", "#2e3540")])
        st.configure("TCheckbutton", background=p, foreground=t,
                     indicatorcolor=f, font=("Segoe UI", 9))
        st.map("TCheckbutton", background=[("active", p)],
               indicatorcolor=[("selected", FARBE["akzent"])])
        st.configure("TRadiobutton", background=p, foreground=t, indicatorcolor=f)
        st.map("TRadiobutton", background=[("active", p)],
               indicatorcolor=[("selected", FARBE["akzent"])])
        st.configure("TNotebook", background=FARBE["bg"], borderwidth=0, tabmargins=0)
        st.configure("TNotebook.Tab", background="#272d36", foreground=g,
                     padding=(9, 6), borderwidth=0, font=("Segoe UI", 9))
        st.map("TNotebook.Tab", background=[("selected", p)],
               foreground=[("selected", t)])
        st.configure("Treeview", background=f, fieldbackground=f, foreground=t,
                     borderwidth=0, rowheight=21, font=("Segoe UI", 9))
        st.configure("Treeview.Heading", background="#272d36", foreground=g,
                     borderwidth=0, font=("Segoe UI", 8))
        st.map("Treeview", background=[("selected", FARBE["akzent"])],
               foreground=[("selected", "#10131a")])
        st.configure("TLabelframe", background=p, bordercolor=FARBE["linie"])
        st.configure("TLabelframe.Label", background=p, foreground=FARBE["akzent"],
                     font=("Segoe UI Semibold", 9))
        st.configure("TScale", background=p, troughcolor=f)
        st.configure("TSeparator", background=FARBE["linie"])

    @staticmethod
    def _mpl_stil():
        matplotlib.rcParams.update({
            "figure.facecolor": FARBE["panel"], "axes.facecolor": FARBE["feld"],
            "axes.edgecolor": FARBE["linie"], "axes.labelcolor": FARBE["text"],
            "axes.titlecolor": FARBE["text"], "text.color": FARBE["text"],
            "xtick.color": FARBE["grau"], "ytick.color": FARBE["grau"],
            "grid.color": FARBE["linie"], "grid.alpha": 0.5,
            "font.size": 8.5, "axes.titlesize": 9.5, "axes.labelsize": 8.5,
            "legend.fontsize": 7.5, "figure.dpi": 100,
            "savefig.facecolor": FARBE["panel"], "savefig.dpi": 200,
        })

    # ----------------------------------------------------------- Variablen
    def _variablen(self):
        sz = self.szene
        self.v_nx = tk.StringVar(value=str(sz.nx))
        self.v_ny = tk.StringVar(value=str(sz.ny))
        self.v_d = tk.StringVar(value=f"{sz.d * 1e3:g}")
        self.v_S = tk.StringVar(value=f"{sz.courant:g}")
        self.v_npml = tk.StringVar(value=str(sz.npml))
        self.v_pml_m = tk.StringVar(value=str(sz.pml_m))
        self.v_r0 = tk.StringVar(value=f"{sz.r0:g}")

        self.v_q_typ = tk.StringVar(value="ebene")
        self.v_q_art = tk.StringVar(value="rampe")
        self.v_q_i = tk.StringVar(value="30")
        self.v_q_j = tk.StringVar(value="120")
        self.v_q_f0 = tk.StringVar(value="15")
        self.v_q_amp = tk.StringVar(value="1")
        self.v_q_zyk = tk.StringVar(value="4")
        self.v_q_theta = tk.StringVar(value="0")
        self.v_q_modus = tk.StringVar(value="strom")
        self.v_q_taper = tk.StringVar(value="0")
        self.v_q_verz = tk.StringVar(value="0")
        self.v_q_arm = tk.BooleanVar(value=True)
        self.v_q_richt = tk.BooleanVar(value=False)
        self.v_q_aktiv = tk.BooleanVar(value=True)

        self.v_o_form = tk.StringVar(value="rechteck")
        self.v_o_art = tk.StringVar(value="dielektrikum")
        self.v_o_i0 = tk.StringVar(value="120")
        self.v_o_i1 = tk.StringVar(value="160")
        self.v_o_j0 = tk.StringVar(value="100")
        self.v_o_j1 = tk.StringVar(value="160")
        self.v_o_eps = tk.StringVar(value="2.25")
        self.v_o_sig = tk.StringVar(value="0")

        self.v_s_i = tk.StringVar(value="200")
        self.v_s_j = tk.StringVar(value="130")

        self.v_darstellung = tk.StringVar(value="Ez")
        self.v_cmap = tk.StringVar(value="RdBu_r")
        self.v_skala = tk.StringVar(value="auto")
        self.v_vmax = tk.StringVar(value="1.0")
        self.v_ov_mat = tk.BooleanVar(value=True)
        self.v_ov_pml = tk.BooleanVar(value=True)
        self.v_ov_qu = tk.BooleanVar(value=True)
        self.v_ov_so = tk.BooleanVar(value=True)
        self.v_ov_schirm = tk.BooleanVar(value=True)
        self.v_chunk = tk.IntVar(value=4)
        self.v_rate = tk.IntVar(value=0)
        self.v_grenze = tk.StringVar(value="0")
        self.v_intens = tk.BooleanVar(value=True)
        self.v_analyse_ms = tk.IntVar(value=400)
        self.v_schirm = tk.StringVar(value="0")
        self.v_zeitfenster = tk.StringVar(value="letzte 25 Perioden")
        self.v_intens_norm = tk.StringVar(value="Gesamtmaximum")

        self.v_werkzeug = tk.StringVar(value="zeiger")
        self.v_preset = tk.StringVar(value="Freiraum: ebene Welle")
        self.v_aufnahme = tk.BooleanVar(value=False)

        self.st_schritt = tk.StringVar(value="-")
        self.st_zeit = tk.StringVar(value="-")
        self.st_rate = tk.StringVar(value="-")
        self.st_spitze = tk.StringVar(value="-")
        self.st_energie = tk.StringVar(value="-")
        self.st_cursor = tk.StringVar(value="")
        self.st_hinweis = tk.StringVar(value="")

    # ------------------------------------------------------------ Menueband
    def _menue(self):
        leiste = tk.Menu(self, tearoff=0, bg=FARBE["panel"], fg=FARBE["text"],
                         activebackground=FARBE["akzent"], activeforeground="#10131a",
                         borderwidth=0)

        datei = tk.Menu(leiste, tearoff=0, bg=FARBE["panel"], fg=FARBE["text"],
                        activebackground=FARBE["akzent"], activeforeground="#10131a")
        datei.add_command(label="Szene laden ...", command=self._szene_laden)
        datei.add_command(label="Szene speichern ...", command=self._szene_speichern)
        datei.add_separator()
        datei.add_command(label="Aktuelle Ansicht als Bild ...", command=self._bild_speichern)
        datei.add_command(label="Sondendaten als CSV ...", command=self._csv_speichern)
        datei.add_command(label="Feld als NPY ...", command=lambda: self._npy_speichern("feld"))
        datei.add_command(label="Intensitaet als NPY ...", command=lambda: self._npy_speichern("intens"))
        datei.add_command(label="Aufnahme als GIF ...", command=self._gif_speichern)
        datei.add_separator()
        datei.add_command(label="Beenden", command=self._beenden)
        leiste.add_cascade(label="Datei", menu=datei)

        sim = tk.Menu(leiste, tearoff=0, bg=FARBE["panel"], fg=FARBE["text"],
                      activebackground=FARBE["akzent"], activeforeground="#10131a")
        sim.add_command(label="Start / Pause        [Leertaste]", command=self._umschalten)
        sim.add_command(label="Einzelschritt        [S]", command=lambda: self._schritt(1))
        sim.add_command(label="10 Schritte", command=lambda: self._schritt(10))
        sim.add_command(label="Zuruecksetzen        [R]", command=self._reset)
        sim.add_command(label="Neu aufbauen         [N]", command=self._neu_aufbauen)
        leiste.add_cascade(label="Simulation", menu=sim)

        mess = tk.Menu(leiste, tearoff=0, bg=FARBE["panel"], fg=FARBE["text"],
                       activebackground=FARBE["akzent"], activeforeground="#10131a")
        mess.add_command(label="Phasengeschwindigkeit (Sonde 1 -> 2)",
                         command=self._mess_phasengeschwindigkeit)
        mess.add_command(label="Uebertragung Sonde 2 / Sonde 1",
                         command=self._mess_uebertragung)
        mess.add_command(label="Sondenkennzahlen", command=self._mess_kennzahlen)
        mess.add_command(label="Energiebilanz", command=self._mess_energie)
        leiste.add_cascade(label="Messen", menu=mess)

        hilfe = tk.Menu(leiste, tearoff=0, bg=FARBE["panel"], fg=FARBE["text"],
                        activebackground=FARBE["akzent"], activeforeground="#10131a")
        hilfe.add_command(label="Kurzanleitung", command=self._hilfe)
        leiste.add_cascade(label="Hilfe", menu=hilfe)
        self.config(menu=leiste)

    def _werkzeugleiste(self):
        leiste = ttk.Frame(self, style="Leiste.TFrame")
        leiste.pack(fill="x", padx=8, pady=6)

        self.btn_start = ttk.Button(leiste, text="▶  Start", width=10,
                                    style="Akzent.TButton", command=self._umschalten)
        self.btn_start.pack(side="left")
        ttk.Button(leiste, text="Schritt", width=8,
                   command=lambda: self._schritt(1)).pack(side="left", padx=3)
        ttk.Button(leiste, text="+10", width=5,
                   command=lambda: self._schritt(10)).pack(side="left")
        ttk.Button(leiste, text="Reset", width=7,
                   command=self._reset).pack(side="left", padx=3)
        ttk.Button(leiste, text="Neu aufbauen", width=13,
                   command=self._neu_aufbauen).pack(side="left")

        ttk.Separator(leiste, orient="vertical").pack(side="left", fill="y", padx=10)
        ttk.Label(leiste, text="Preset", style="Leiste.TLabel").pack(side="left", padx=(0, 4))
        box = ttk.Combobox(leiste, textvariable=self.v_preset, values=list(PRESETS),
                           state="readonly", width=30)
        box.pack(side="left")
        box.bind("<<ComboboxSelected>>", lambda e: self._preset_laden())

        ttk.Separator(leiste, orient="vertical").pack(side="left", fill="y", padx=10)
        ttk.Checkbutton(leiste, text="Aufnahme (GIF)", variable=self.v_aufnahme,
                        command=self._aufnahme_schalten).pack(side="left")
        self.lb_film = ttk.Label(leiste, text="0 Bilder", style="Leiste.TLabel")
        self.lb_film.pack(side="left", padx=6)

        ttk.Label(leiste, textvariable=self.st_hinweis,
                  style="Leiste.TLabel").pack(side="right")

    # --------------------------------------------------------- linkes Panel
    def _linkes_panel(self, eltern):
        links = ttk.Frame(eltern, style="Fenster.TFrame", width=360)
        links.pack(side="left", fill="y", padx=(0, 6))
        links.pack_propagate(False)

        nb = ttk.Notebook(links)
        nb.pack(fill="both", expand=True)
        self._tab_gitter(nb)
        self._tab_quellen(nb)
        self._tab_objekte(nb)
        self._tab_sonden(nb)
        self._tab_anzeige(nb)

        rahmen = ttk.Frame(links)
        rahmen.pack(fill="both", pady=(6, 0))
        ttk.Label(rahmen, text="Protokoll", style="Kopf.TLabel").pack(anchor="w", padx=6, pady=(4, 2))
        self.txt_log = tk.Text(rahmen, height=9, bg=FARBE["feld"], fg=FARBE["grau"],
                               bd=0, highlightthickness=0, wrap="word",
                               font=("Consolas", 8), padx=6, pady=4)
        self.txt_log.pack(fill="both", expand=True)
        self.txt_log.tag_config("info", foreground=FARBE["grau"])
        self.txt_log.tag_config("gut", foreground=FARBE["gut"])
        self.txt_log.tag_config("warn", foreground=FARBE["warn"])
        self.txt_log.tag_config("fehler", foreground=FARBE["fehler"])
        self.txt_log.configure(state="disabled")

    def _zeile(self, eltern, text, var, breite=10, einheit="", width_label=17):
        z = ttk.Frame(eltern)
        z.pack(fill="x", pady=2)
        ttk.Label(z, text=text, width=width_label).pack(side="left")
        e = ttk.Entry(z, textvariable=var, width=breite)
        e.pack(side="left")
        if einheit:
            ttk.Label(z, text=einheit, style="Grau.TLabel").pack(side="left", padx=4)
        return e

    def _auswahl(self, eltern, text, var, werte, breite=18, width_label=17, cb=None):
        z = ttk.Frame(eltern)
        z.pack(fill="x", pady=2)
        ttk.Label(z, text=text, width=width_label).pack(side="left")
        b = ttk.Combobox(z, textvariable=var, values=werte, state="readonly", width=breite)
        b.pack(side="left", fill="x", expand=True)
        if cb:
            b.bind("<<ComboboxSelected>>", lambda e: cb())
        return b

    # -- Tab Gitter --------------------------------------------------------
    def _tab_gitter(self, nb):
        tab = ttk.Frame(nb, padding=8)
        nb.add(tab, text="Gitter")

        g = ttk.Labelframe(tab, text="Rechengebiet", padding=8)
        g.pack(fill="x")
        self._zeile(g, "Zellen in x (nx)", self.v_nx, 8, "Zellen")
        self._zeile(g, "Zellen in y (ny)", self.v_ny, 8, "Zellen")
        self._zeile(g, "Zellbreite Delta", self.v_d, 8, "mm")
        self._zeile(g, "Courant-Zahl S", self.v_S, 8, "(S <= 1 stabil)")

        p = ttk.Labelframe(tab, text="Absorbierender Rand (CPML)", padding=8)
        p.pack(fill="x", pady=8)
        self._zeile(p, "Dicke npml", self.v_npml, 8, "Zellen")
        self._zeile(p, "Ordnung m", self.v_pml_m, 8)
        self._zeile(p, "Restreflexion R0", self.v_r0, 8)

        ttk.Button(tab, text="Uebernehmen und neu aufbauen", style="Akzent.TButton",
                   command=self._gitter_uebernehmen).pack(fill="x", pady=(2, 8))

        info = ttk.Labelframe(tab, text="Kennzahlen", padding=8)
        info.pack(fill="both", expand=True)
        self.lb_info = tk.Text(info, height=11, bg=FARBE["feld"], fg=FARBE["text"],
                               bd=0, highlightthickness=0, font=("Consolas", 8),
                               padx=6, pady=4)
        self.lb_info.pack(fill="both", expand=True)
        self.lb_info.tag_config("warn", foreground=FARBE["warn"])
        self.lb_info.tag_config("gut", foreground=FARBE["gut"])
        self.lb_info.configure(state="disabled")

    # -- Tab Quellen -------------------------------------------------------
    def _tab_quellen(self, nb):
        tab = ttk.Frame(nb, padding=8)
        nb.add(tab, text="Quellen")

        self.tv_q = ttk.Treeview(tab, columns=("a", "b", "c"), show="headings", height=5)
        for spalte, text, br in (("a", "Typ", 95), ("b", "Ort", 80), ("c", "Signal", 150)):
            self.tv_q.heading(spalte, text=text)
            self.tv_q.column(spalte, width=br, anchor="w")
        self.tv_q.pack(fill="x")
        self.tv_q.bind("<<TreeviewSelect>>", lambda e: self._quelle_gewaehlt())

        kn = ttk.Frame(tab)
        kn.pack(fill="x", pady=4)
        ttk.Button(kn, text="Neu", style="Klein.TButton", width=7,
                   command=self._quelle_neu).pack(side="left")
        ttk.Button(kn, text="Kopie", style="Klein.TButton", width=7,
                   command=self._quelle_kopie).pack(side="left", padx=3)
        ttk.Button(kn, text="Loeschen", style="Klein.TButton", width=9,
                   command=self._quelle_loeschen).pack(side="left")
        ttk.Button(kn, text="An/Aus", style="Klein.TButton", width=8,
                   command=self._quelle_schalten).pack(side="left", padx=3)

        e = ttk.Labelframe(tab, text="Eigenschaften", padding=8)
        e.pack(fill="both", expand=True)
        self._auswahl(e, "Typ", self.v_q_typ, ["ebene", "punkt"], 16,
                      cb=self._quelle_typ_wechsel)
        self._auswahl(e, "Wellenform", self.v_q_art, list(WELLENFORMEN), 16)
        self._zeile(e, "Frequenz f0", self.v_q_f0, 8, "GHz")
        self._zeile(e, "Amplitude", self.v_q_amp, 8)
        self._zeile(e, "Einschaltrampe", self.v_q_zyk, 8, "Perioden (Sinus)")
        self._zeile(e, "Start bei i", self.v_q_i, 8, "Zelle")
        self.zl_j = self._zeile(e, "Punkt bei j", self.v_q_j, 8, "Zelle")
        self._zeile(e, "Winkel theta", self.v_q_theta, 8, "Grad")
        self._zeile(e, "Kantenglaettung", self.v_q_taper, 8, "Zellen")
        self._zeile(e, "Verzoegerung", self.v_q_verz, 8, "Perioden")
        self._auswahl(e, "Einspeisung", self.v_q_modus, ["strom", "weich"], 16)
        ttk.Checkbutton(e, text="Arm entlang der Einfallskante (gegen Eckluecke)",
                        variable=self.v_q_arm).pack(anchor="w", pady=(4, 0))
        ttk.Checkbutton(e, text="Rueckwaertswelle unterdruecken (lambda/4-Paar)",
                        variable=self.v_q_richt).pack(anchor="w")
        ttk.Checkbutton(e, text="Quelle aktiv", variable=self.v_q_aktiv).pack(anchor="w")
        ttk.Label(e, text="Ebene Welle = Punktquellen ueber die ganze Kante,\n"
                          "bis in die PML hinein; schraeg = zeitversetzt.\n"
                          "Rampe 0 schaltet hart ein und regt das ganze Band an.",
                  style="Grau.TLabel", justify="left").pack(anchor="w", pady=(6, 4))
        ttk.Button(e, text="Uebernehmen und neu aufbauen", style="Akzent.TButton",
                   command=self._quelle_uebernehmen).pack(fill="x")
        self._quelle_typ_wechsel()

    # -- Tab Objekte -------------------------------------------------------
    def _tab_objekte(self, nb):
        tab = ttk.Frame(nb, padding=8)
        nb.add(tab, text="Objekte")

        self.tv_o = ttk.Treeview(tab, columns=("a", "b", "c"), show="headings", height=7)
        for spalte, text, br in (("a", "Form", 80), ("b", "Bereich", 130), ("c", "Material", 110)):
            self.tv_o.heading(spalte, text=text)
            self.tv_o.column(spalte, width=br, anchor="w")
        self.tv_o.pack(fill="x")
        self.tv_o.bind("<<TreeviewSelect>>", lambda e: self._objekt_gewaehlt())

        kn = ttk.Frame(tab)
        kn.pack(fill="x", pady=4)
        ttk.Button(kn, text="Neu", style="Klein.TButton", width=7,
                   command=self._objekt_neu).pack(side="left")
        ttk.Button(kn, text="Kopie", style="Klein.TButton", width=7,
                   command=self._objekt_kopie).pack(side="left", padx=3)
        ttk.Button(kn, text="Loeschen", style="Klein.TButton", width=9,
                   command=self._objekt_loeschen).pack(side="left")
        ttk.Button(kn, text="An/Aus", style="Klein.TButton", width=8,
                   command=self._objekt_schalten).pack(side="left", padx=3)
        ttk.Button(kn, text="Alle weg", style="Klein.TButton", width=9,
                   command=self._objekte_leeren).pack(side="left")

        e = ttk.Labelframe(tab, text="Eigenschaften", padding=8)
        e.pack(fill="x", pady=(4, 0))
        self._auswahl(e, "Form", self.v_o_form, ["rechteck", "ellipse"], 14)
        self._auswahl(e, "Material", self.v_o_art, ["dielektrikum", "pec"], 14)
        z = ttk.Frame(e)
        z.pack(fill="x", pady=2)
        ttk.Label(z, text="Bereich i", width=17).pack(side="left")
        ttk.Entry(z, textvariable=self.v_o_i0, width=6).pack(side="left")
        ttk.Label(z, text="..", style="Grau.TLabel").pack(side="left", padx=3)
        ttk.Entry(z, textvariable=self.v_o_i1, width=6).pack(side="left")
        z = ttk.Frame(e)
        z.pack(fill="x", pady=2)
        ttk.Label(z, text="Bereich j", width=17).pack(side="left")
        ttk.Entry(z, textvariable=self.v_o_j0, width=6).pack(side="left")
        ttk.Label(z, text="..", style="Grau.TLabel").pack(side="left", padx=3)
        ttk.Entry(z, textvariable=self.v_o_j1, width=6).pack(side="left")
        self._zeile(e, "eps_r", self.v_o_eps, 8)
        self._zeile(e, "sigma", self.v_o_sig, 8, "S/m")
        ttk.Button(e, text="Uebernehmen und neu aufbauen", style="Akzent.TButton",
                   command=self._objekt_uebernehmen).pack(fill="x", pady=(6, 0))

        w = ttk.Labelframe(tab, text="Bausteine", padding=8)
        w.pack(fill="x", pady=8)
        ttk.Button(w, text="Wand mit Doppelspalt einfuegen",
                   command=self._baustein_doppelspalt).pack(fill="x", pady=1)
        ttk.Button(w, text="Halbraum (Grenzflaeche) einfuegen",
                   command=self._baustein_halbraum).pack(fill="x", pady=1)
        ttk.Button(w, text="Schichtstapel einfuegen",
                   command=self._baustein_stapel).pack(fill="x", pady=1)

    # -- Tab Sonden --------------------------------------------------------
    def _tab_sonden(self, nb):
        tab = ttk.Frame(nb, padding=8)
        nb.add(tab, text="Sonden")

        self.tv_s = ttk.Treeview(tab, columns=("a", "b"), show="headings", height=7)
        self.tv_s.heading("a", text="Nr.")
        self.tv_s.column("a", width=45, anchor="w")
        self.tv_s.heading("b", text="Position (i, j)")
        self.tv_s.column("b", width=180, anchor="w")
        self.tv_s.pack(fill="x")

        z = ttk.Frame(tab)
        z.pack(fill="x", pady=4)
        ttk.Label(z, text="i", width=2).pack(side="left")
        ttk.Entry(z, textvariable=self.v_s_i, width=6).pack(side="left", padx=(0, 6))
        ttk.Label(z, text="j", width=2).pack(side="left")
        ttk.Entry(z, textvariable=self.v_s_j, width=6).pack(side="left", padx=(0, 6))
        ttk.Button(z, text="Setzen", style="Klein.TButton",
                   command=lambda: self._sonde_hinzu(zahl(self.v_s_i, 0, True),
                                                     zahl(self.v_s_j, 0, True))).pack(side="left")
        ttk.Button(z, text="Loeschen", style="Klein.TButton",
                   command=self._sonde_loeschen).pack(side="left", padx=3)
        ttk.Label(tab, text="Tipp: Werkzeug 'Sonde' waehlen und ins Feld klicken.",
                  style="Grau.TLabel").pack(anchor="w")

        m = ttk.Labelframe(tab, text="Messwerkzeuge", padding=8)
        m.pack(fill="x", pady=8)
        ttk.Button(m, text="Phasengeschwindigkeit (Sonde 1 -> 2)",
                   command=self._mess_phasengeschwindigkeit).pack(fill="x", pady=1)
        ttk.Button(m, text="Uebertragung Sonde 2 / Sonde 1",
                   command=self._mess_uebertragung).pack(fill="x", pady=1)
        ttk.Button(m, text="Sondenkennzahlen", command=self._mess_kennzahlen).pack(fill="x", pady=1)
        ttk.Button(m, text="Energiebilanz", command=self._mess_energie).pack(fill="x", pady=1)

        x = ttk.Labelframe(tab, text="Export", padding=8)
        x.pack(fill="x")
        ttk.Button(x, text="Sondendaten als CSV", command=self._csv_speichern).pack(fill="x", pady=1)
        ttk.Button(x, text="Ansicht als Bild (PNG/PDF)", command=self._bild_speichern).pack(fill="x", pady=1)
        ttk.Button(x, text="Aufnahme als GIF", command=self._gif_speichern).pack(fill="x", pady=1)
        ttk.Button(x, text="Szene speichern", command=self._szene_speichern).pack(fill="x", pady=1)
        ttk.Button(x, text="Szene laden", command=self._szene_laden).pack(fill="x", pady=1)

    # -- Tab Anzeige -------------------------------------------------------
    def _tab_anzeige(self, nb):
        tab = ttk.Frame(nb, padding=8)
        nb.add(tab, text="Anzeige")

        d = ttk.Labelframe(tab, text="Feldbild", padding=8)
        d.pack(fill="x")
        self._auswahl(d, "Groesse", self.v_darstellung, DARSTELLUNGEN, 16,
                      cb=self._darstellung_wechsel)
        self._auswahl(d, "Farbskala", self.v_cmap, CMAPS, 16, cb=self._cmap_wechsel)
        self._auswahl(d, "Skalierung", self.v_skala, ["auto", "fest"], 16)
        self._zeile(d, "fester Maximalwert", self.v_vmax, 8)
        for text, var in (("Material / PEC zeigen", self.v_ov_mat),
                          ("PML-Rand zeigen", self.v_ov_pml),
                          ("Quellen zeigen", self.v_ov_qu),
                          ("Sonden zeigen", self.v_ov_so),
                          ("Schirmlinie zeigen", self.v_ov_schirm)):
            ttk.Checkbutton(d, text=text, variable=var,
                            command=self._overlay_neu).pack(anchor="w")

        r = ttk.Labelframe(tab, text="Rechenlauf", padding=8)
        r.pack(fill="x", pady=8)
        z = ttk.Frame(r)
        z.pack(fill="x", pady=2)
        ttk.Label(z, text="Schritte je Bild", width=17).pack(side="left")
        ttk.Scale(z, from_=1, to=40, variable=self.v_chunk, orient="horizontal",
                  command=lambda e: self._lauf_uebernehmen()).pack(side="left", fill="x", expand=True)
        self.lb_chunk = ttk.Label(z, text="4", style="Wert.TLabel", width=4)
        self.lb_chunk.pack(side="left")
        z = ttk.Frame(r)
        z.pack(fill="x", pady=2)
        ttk.Label(z, text="Tempolimit", width=17).pack(side="left")
        ttk.Scale(z, from_=0, to=4000, variable=self.v_rate, orient="horizontal",
                  command=lambda e: self._lauf_uebernehmen()).pack(side="left", fill="x", expand=True)
        self.lb_rate = ttk.Label(z, text="aus", style="Wert.TLabel", width=6)
        self.lb_rate.pack(side="left")
        self._zeile(r, "Stopp bei Schritt", self.v_grenze, 8, "(0 = nie)")
        ttk.Checkbutton(r, text="Intensitaet mitteln (<Ez^2>)", variable=self.v_intens,
                        command=self._lauf_uebernehmen).pack(anchor="w")
        ttk.Button(r, text="Lauf-Einstellungen uebernehmen",
                   command=self._lauf_uebernehmen).pack(fill="x", pady=(4, 0))

        a = ttk.Labelframe(tab, text="Diagramme", padding=8)
        a.pack(fill="x")
        self._auswahl(a, "Zeitfenster", self.v_zeitfenster, ZEITFENSTER, 14,
                      cb=self._analyse_zeichnen)
        self._zeile(a, "Schirm bei i", self.v_schirm, 8, "Zelle")
        self._auswahl(a, "Intensitaet bezogen auf", self.v_intens_norm,
                      ["Gesamtmaximum", "Schirmmaximum", "logarithmisch (dB)"],
                      14, cb=self._intens_zeichnen)
        z = ttk.Frame(a)
        z.pack(fill="x", pady=2)
        ttk.Label(z, text="Auffrischung", width=17).pack(side="left")
        ttk.Scale(z, from_=100, to=2000, variable=self.v_analyse_ms,
                  orient="horizontal").pack(side="left", fill="x", expand=True)
        ttk.Label(z, text="ms", style="Grau.TLabel").pack(side="left")
        ttk.Button(a, text="Diagramme jetzt auffrischen",
                   command=self._diagramme_auffrischen).pack(fill="x", pady=(4, 0))

    # -------------------------------------------------------- rechtes Panel
    def _rechtes_panel(self, eltern):
        rechts = ttk.Frame(eltern, style="Fenster.TFrame")
        rechts.pack(side="left", fill="both", expand=True)

        leiste = ttk.Frame(rechts, style="Leiste.TFrame")
        leiste.pack(fill="x", pady=(0, 4))
        ttk.Label(leiste, text="Maus:", style="Leiste.TLabel").pack(side="left", padx=(2, 6))
        werkzeuge = [("Zeiger", "zeiger"), ("Sonde", "sonde"), ("Punktquelle", "quelle"),
                     ("Rechteck", "rechteck"), ("Ellipse", "ellipse"),
                     ("PEC-Rechteck", "pec_rechteck"), ("PEC-Ellipse", "pec_ellipse"),
                     ("Radierer", "radierer")]
        for text, wert in werkzeuge:
            ttk.Radiobutton(leiste, text=text, value=wert, variable=self.v_werkzeug,
                            style="TRadiobutton").pack(side="left", padx=3)
        ttk.Label(leiste, textvariable=self.st_cursor, style="Leiste.TLabel",
                  font=("Consolas", 9)).pack(side="right", padx=6)

        self.nb_ansicht = ttk.Notebook(rechts)
        self.nb_ansicht.pack(fill="both", expand=True)

        # Feld ------------------------------------------------------------
        tab_feld = ttk.Frame(self.nb_ansicht)
        self.nb_ansicht.add(tab_feld, text="  Feld  ")
        self.fig_feld = Figure(figsize=(8, 6))
        self.kanvas_feld = FigureCanvasTkAgg(self.fig_feld, master=tab_feld)
        self.kanvas_feld.get_tk_widget().pack(fill="both", expand=True)
        nav = NavigationToolbar2Tk(self.kanvas_feld, tab_feld, pack_toolbar=False)
        nav.update()
        nav.pack(fill="x")
        self._toolbar_faerben(nav)
        self.kanvas_feld.mpl_connect("draw_event", self._auf_draw)
        self.kanvas_feld.mpl_connect("button_press_event", self._maus_druck)
        self.kanvas_feld.mpl_connect("motion_notify_event", self._maus_bewegung)
        self.kanvas_feld.mpl_connect("button_release_event", self._maus_los)

        # Analyse ---------------------------------------------------------
        tab_an = ttk.Frame(self.nb_ansicht)
        self.nb_ansicht.add(tab_an, text="  Analyse  ")
        self.fig_an = Figure(figsize=(8, 6), layout="constrained")
        raster = self.fig_an.add_gridspec(2, 2, height_ratios=[1.15, 1.0])
        self.ax_zeit = self.fig_an.add_subplot(raster[0, :])
        self.ax_spek = self.fig_an.add_subplot(raster[1, 0])
        self.ax_ener = self.fig_an.add_subplot(raster[1, 1])
        self.kanvas_an = FigureCanvasTkAgg(self.fig_an, master=tab_an)
        self.kanvas_an.get_tk_widget().pack(fill="both", expand=True)
        nav2 = NavigationToolbar2Tk(self.kanvas_an, tab_an, pack_toolbar=False)
        nav2.update()
        nav2.pack(fill="x")
        self._toolbar_faerben(nav2)

        # Intensitaet -----------------------------------------------------
        tab_in = ttk.Frame(self.nb_ansicht)
        self.nb_ansicht.add(tab_in, text="  Intensitaet  ")
        self.fig_in = Figure(figsize=(8, 6), layout="constrained")
        self.ax_in_bild = self.fig_in.add_subplot(121)
        self.ax_in_prof = self.fig_in.add_subplot(122)
        self.kanvas_in = FigureCanvasTkAgg(self.fig_in, master=tab_in)
        self.kanvas_in.get_tk_widget().pack(fill="both", expand=True)
        nav3 = NavigationToolbar2Tk(self.kanvas_in, tab_in, pack_toolbar=False)
        nav3.update()
        nav3.pack(fill="x")
        self._toolbar_faerben(nav3)

    @staticmethod
    def _toolbar_faerben(nav):
        try:
            nav.config(background=FARBE["panel"])
            for kind in nav.winfo_children():
                kind.config(background=FARBE["panel"])
                if isinstance(kind, tk.Label):
                    kind.config(foreground=FARBE["grau"])
        except tk.TclError:
            pass

    def _statusleiste(self):
        s = ttk.Frame(self, style="Leiste.TFrame")
        s.pack(fill="x", padx=8, pady=(0, 6))
        for text, var in (("Schritt", self.st_schritt), ("t", self.st_zeit),
                          ("Tempo", self.st_rate), ("max|Ez|", self.st_spitze),
                          ("Energie", self.st_energie)):
            ttk.Label(s, text=text + ":", style="Leiste.TLabel").pack(side="left", padx=(10, 3))
            ttk.Label(s, textvariable=var, style="Leiste.TLabel",
                      font=("Consolas", 9)).pack(side="left")

    def _tastatur(self):
        self.bind("<space>", lambda e: self._umschalten())
        self.bind("<KeyPress-s>", lambda e: self._schritt(1))
        self.bind("<KeyPress-S>", lambda e: self._schritt(10))
        self.bind("<KeyPress-r>", lambda e: self._reset())
        self.bind("<KeyPress-n>", lambda e: self._neu_aufbauen())
        self.bind("<Control-s>", lambda e: self._szene_speichern())
        self.bind("<Control-o>", lambda e: self._szene_laden())

    # ================================================== Protokoll / Status
    def _log(self, text, art="info"):
        self.txt_log.configure(state="normal")
        self.txt_log.insert("end", time.strftime("%H:%M:%S  ") + text + "\n", art)
        self.txt_log.see("end")
        self.txt_log.configure(state="disabled")

    def _log_thread(self, text, art="info"):
        self._meldungen.append((text, art))

    # ====================================================== Szene anwenden
    def _szene_anwenden(self, szene, meldung=None):
        self.szene = szene
        try:
            self.motor.setze_szene(szene)
        except Exception as fehler:                     # noqa: BLE001
            messagebox.showerror("Aufbau fehlgeschlagen", str(fehler))
            self._log(f"Aufbau fehlgeschlagen: {fehler}", "fehler")
            return
        self._spitze = 1e-12
        self._bg = None
        if szene.schirm:
            self.v_schirm.set(str(int(np.clip(szene.schirm, 0, szene.nx - 1))))
        elif not self.v_schirm.get().strip() or zahl(self.v_schirm, 0, True) >= szene.nx:
            self.v_schirm.set(str(max(szene.nx - 20, 1)))
        self._listen_fuellen()
        self._info_fuellen()
        self._feld_aufbauen()
        self._diagramme_auffrischen()
        self.btn_start.config(text="▶  Start")
        if meldung:
            self._log(meldung, "gut")

    def _neu_aufbauen(self):
        self._szene_anwenden(self.szene, "Szene neu aufgebaut, Felder zurueckgesetzt.")

    def _preset_laden(self):
        name = self.v_preset.get()
        if name in PRESETS:
            sz = PRESETS[name]()
            self._gitter_vars_setzen(sz)
            self._szene_anwenden(sz, f"Preset '{name}' geladen.")

    def _gitter_vars_setzen(self, sz):
        self.v_nx.set(str(sz.nx))
        self.v_ny.set(str(sz.ny))
        self.v_d.set(f"{sz.d * 1e3:g}")
        self.v_S.set(f"{sz.courant:g}")
        self.v_npml.set(str(sz.npml))
        self.v_pml_m.set(str(sz.pml_m))
        self.v_r0.set(f"{sz.r0:g}")

    def _gitter_uebernehmen(self):
        sz = self.szene
        sz.nx = max(zahl(self.v_nx, sz.nx, True), 20)
        sz.ny = max(zahl(self.v_ny, sz.ny, True), 20)
        sz.d = max(zahl(self.v_d, sz.d * 1e3) * 1e-3, 1e-9)
        sz.courant = max(zahl(self.v_S, sz.courant), 0.01)
        sz.npml = max(zahl(self.v_npml, sz.npml, True), 0)
        sz.pml_m = max(zahl(self.v_pml_m, sz.pml_m, True), 1)
        sz.r0 = max(zahl(self.v_r0, sz.r0), 1e-14)
        if sz.courant > 1.0:
            self._log(f"S = {sz.courant:g} > 1: Die Rechnung wird instabil "
                      f"(das ist als Demonstration gewollt).", "warn")
        self._szene_anwenden(sz, "Gitter uebernommen.")

    def _info_fuellen(self):
        sz, k = self.szene, self.szene.kennzahlen()
        eps, sig, pec = self.motor.eps, self.motor.sig, self.motor.pec
        eps_max = float(eps.max())
        n_lam = k["n_lambda"] / np.sqrt(eps_max)
        zeilen = [
            ("Gitter", f"{sz.nx} x {sz.ny} Zellen = {k['breite']*1e3:.0f} x {k['hoehe']*1e3:.0f} mm"),
            ("Zellbreite", f"{sz.d*1e3:.3f} mm"),
            ("Zeitschritt dt", f"{k['dt']*1e12:.4f} ps"),
            ("Courant S", f"{sz.courant:.3f}" + ("   INSTABIL" if sz.courant > 1 else "   stabil")),
            ("Durchlauf", f"{k['schritte_quer']} Schritte quer durchs Gebiet"),
            ("f0 (Quelle 1)", f"{k['f0']/1e9:.3f} GHz" if k["f0"] else "-"),
            ("lambda0", f"{C0/k['f0']*1e3:.2f} mm" if k["f0"] else "-"),
            ("Zellen je lambda", f"{n_lam:.1f}" + ("   (grob!)" if n_lam < 10 else "")),
            ("PML", f"{sz.npml} Zellen, m = {sz.pml_m}, R0 = {sz.r0:g}"),
            ("Material", f"eps_r 1 .. {eps_max:g}, PEC-Zellen: {int(pec.sum())}"),
            ("Objekte/Quellen", f"{len(sz.objekte)} / {len(sz.quellen)}, Sonden: {len(sz.sonden)}"),
        ]
        self.lb_info.configure(state="normal")
        self.lb_info.delete("1.0", "end")
        for name, wert in zeilen:
            tag = "warn" if ("INSTABIL" in wert or "grob" in wert) else "gut"
            self.lb_info.insert("end", f"{name:<17}", "")
            self.lb_info.insert("end", f"{wert}\n", tag if tag == "warn" else "")
        self.lb_info.configure(state="disabled")

    # ====================================================== Listen / Editor
    def _listen_fuellen(self):
        for zeile in self.tv_q.get_children():
            self.tv_q.delete(zeile)
        for k, q in enumerate(self.szene.quellen):
            self.tv_q.insert("", "end", iid=str(k), values=q.zeile())
        for zeile in self.tv_o.get_children():
            self.tv_o.delete(zeile)
        for k, o in enumerate(self.szene.objekte):
            self.tv_o.insert("", "end", iid=str(k), values=o.zeile())
        for zeile in self.tv_s.get_children():
            self.tv_s.delete(zeile)
        for k, s in enumerate(self.szene.sonden):
            self.tv_s.insert("", "end", iid=str(k), values=(k + 1, f"({s.i}, {s.j})"))

    @staticmethod
    def _index(tv):
        wahl = tv.selection()
        return int(wahl[0]) if wahl else None

    # -- Quellen -----------------------------------------------------------
    def _quelle_typ_wechsel(self):
        """Bei der ebenen Welle belegt die Quelle die ganze Spalte i -
        ein einzelnes j gibt es dann nicht."""
        punkt = self.v_q_typ.get() == "punkt"
        self.zl_j.configure(state="normal" if punkt else "disabled")

    def _quelle_gewaehlt(self):
        k = self._index(self.tv_q)
        if k is None:
            return
        q = self.szene.quellen[k]
        self.v_q_typ.set(q.typ)
        self.v_q_art.set(q.art)
        self.v_q_i.set(str(q.i))
        self.v_q_j.set(str(q.j))
        self.v_q_f0.set(f"{q.f0 / 1e9:g}")
        self.v_q_amp.set(f"{q.amp:g}")
        self.v_q_zyk.set(f"{q.zyklen:g}")
        self.v_q_theta.set(f"{q.theta:g}")
        self.v_q_modus.set(q.modus)
        self.v_q_taper.set(str(q.taper))
        self.v_q_verz.set(f"{q.verzoegerung:g}")
        self.v_q_arm.set(q.arm)
        self.v_q_richt.set(q.richtung)
        self.v_q_aktiv.set(q.aktiv)
        self._quelle_typ_wechsel()

    def _quelle_aus_vars(self):
        return Quelle(
            typ=self.v_q_typ.get(), art=self.v_q_art.get(),
            i=zahl(self.v_q_i, 30, True), j=zahl(self.v_q_j, 100, True),
            f0=max(zahl(self.v_q_f0, 15) * 1e9, 1.0),
            amp=zahl(self.v_q_amp, 1.0), zyklen=max(zahl(self.v_q_zyk, 4.0), 0.25),
            theta=zahl(self.v_q_theta, 0.0), modus=self.v_q_modus.get(),
            taper=max(zahl(self.v_q_taper, 0, True), 0),
            arm=bool(self.v_q_arm.get()), richtung=bool(self.v_q_richt.get()),
            verzoegerung=zahl(self.v_q_verz, 0.0), aktiv=bool(self.v_q_aktiv.get()))

    def _quelle_uebernehmen(self):
        k = self._index(self.tv_q)
        neu = self._quelle_aus_vars()
        if k is None:
            self.szene.quellen.append(neu)
        else:
            self.szene.quellen[k] = neu
        self._szene_anwenden(self.szene, "Quelle uebernommen.")

    def _quelle_neu(self):
        self.szene.quellen.append(self._quelle_aus_vars())
        self._szene_anwenden(self.szene, "Quelle hinzugefuegt.")

    def _quelle_kopie(self):
        k = self._index(self.tv_q)
        if k is None:
            return
        self.szene.quellen.append(copy.deepcopy(self.szene.quellen[k]))
        self._szene_anwenden(self.szene, "Quelle kopiert.")

    def _quelle_loeschen(self):
        k = self._index(self.tv_q)
        if k is None:
            return
        self.szene.quellen.pop(k)
        self._szene_anwenden(self.szene, "Quelle geloescht.")

    def _quelle_schalten(self):
        k = self._index(self.tv_q)
        if k is None:
            return
        self.szene.quellen[k].aktiv = not self.szene.quellen[k].aktiv
        self._szene_anwenden(self.szene, "Quelle umgeschaltet.")

    # -- Objekte -----------------------------------------------------------
    def _objekt_gewaehlt(self):
        k = self._index(self.tv_o)
        if k is None:
            return
        o = self.szene.objekte[k]
        self.v_o_form.set(o.form)
        self.v_o_art.set(o.art)
        self.v_o_i0.set(str(o.i0))
        self.v_o_i1.set(str(o.i1))
        self.v_o_j0.set(str(o.j0))
        self.v_o_j1.set(str(o.j1))
        self.v_o_eps.set(f"{o.eps_r:g}")
        self.v_o_sig.set(f"{o.sigma:g}")

    def _objekt_aus_vars(self):
        return Objekt(form=self.v_o_form.get(), art=self.v_o_art.get(),
                      i0=zahl(self.v_o_i0, 0, True), i1=zahl(self.v_o_i1, 10, True),
                      j0=zahl(self.v_o_j0, 0, True), j1=zahl(self.v_o_j1, 10, True),
                      eps_r=max(zahl(self.v_o_eps, 2.25), 1e-3),
                      sigma=max(zahl(self.v_o_sig, 0.0), 0.0))

    def _objekt_uebernehmen(self):
        k = self._index(self.tv_o)
        neu = self._objekt_aus_vars()
        if k is None:
            self.szene.objekte.append(neu)
        else:
            neu.aktiv = self.szene.objekte[k].aktiv
            self.szene.objekte[k] = neu
        self._szene_anwenden(self.szene, "Objekt uebernommen.")

    def _objekt_neu(self):
        self.szene.objekte.append(self._objekt_aus_vars())
        self._szene_anwenden(self.szene, "Objekt hinzugefuegt.")

    def _objekt_kopie(self):
        k = self._index(self.tv_o)
        if k is None:
            return
        self.szene.objekte.append(copy.deepcopy(self.szene.objekte[k]))
        self._szene_anwenden(self.szene, "Objekt kopiert.")

    def _objekt_loeschen(self):
        k = self._index(self.tv_o)
        if k is None:
            return
        self.szene.objekte.pop(k)
        self._szene_anwenden(self.szene, "Objekt geloescht.")

    def _objekt_schalten(self):
        k = self._index(self.tv_o)
        if k is None:
            return
        self.szene.objekte[k].aktiv = not self.szene.objekte[k].aktiv
        self._szene_anwenden(self.szene, "Objekt umgeschaltet.")

    def _objekte_leeren(self):
        if self.szene.objekte and messagebox.askyesno("Loeschen", "Alle Objekte entfernen?"):
            self.szene.objekte.clear()
            self._szene_anwenden(self.szene, "Alle Objekte entfernt.")

    def _baustein_doppelspalt(self):
        sz = self.szene
        i0 = sz.nx // 3
        mitte, spalt, abstand = sz.ny // 2, 10, 30
        schlitze = [(mitte - abstand // 2 - spalt, mitte - abstand // 2),
                    (mitte + abstand // 2, mitte + abstand // 2 + spalt)]
        sz.objekte.extend(_wand(i0, i0 + 4, sz.ny, schlitze))
        self._szene_anwenden(sz, "Doppelspalt eingefuegt.")

    def _baustein_halbraum(self):
        sz = self.szene
        sz.objekte.append(Objekt("rechteck", "dielektrikum", sz.nx // 2, 0,
                                 sz.nx, sz.ny, eps_r=max(zahl(self.v_o_eps, 2.25), 1e-3),
                                 sigma=max(zahl(self.v_o_sig, 0.0), 0.0)))
        self._szene_anwenden(sz, "Halbraum eingefuegt.")

    def _baustein_stapel(self):
        sz = self.szene
        eps = max(zahl(self.v_o_eps, 4.0), 1e-3)
        f0 = sz.quellen[0].f0 if sz.quellen else 15e9
        lam = C0 / f0
        dicke_hoch = max(int(round(lam / (4 * np.sqrt(eps) * sz.d))), 1)
        dicke_tief = max(int(round(lam / (4 * sz.d))), 1)
        i = sz.nx // 3
        for _ in range(8):
            if i + dicke_hoch >= sz.nx:
                break
            sz.objekte.append(Objekt("rechteck", "dielektrikum", i, 0,
                                     i + dicke_hoch, sz.ny, eps_r=eps))
            i += dicke_hoch + dicke_tief
        self._szene_anwenden(sz, f"Schichtstapel eingefuegt "
                                 f"({dicke_hoch}/{dicke_tief} Zellen je Schicht).")

    # -- Sonden ------------------------------------------------------------
    def _sonde_hinzu(self, i, j):
        """Sonden duerfen im laufenden Betrieb dazukommen - der Kern nimmt
        add_probe jederzeit an, also bleibt das Feld erhalten."""
        i = int(np.clip(i, 0, self.szene.nx - 1))
        j = int(np.clip(j, 0, self.szene.ny - 1))
        self.szene.sonden.append(Sonde(i, j))
        with self.motor.lock:
            if self.motor.sim is not None:
                self.motor.sim.add_probe(i, j)
        self._listen_fuellen()
        self._info_fuellen()
        if hasattr(self, "ln_sonden"):
            self.ln_sonden.set_data([s.i for s in self.szene.sonden],
                                    [s.j for s in self.szene.sonden])
            self.kanvas_feld.draw()
        self._log(f"Sonde {len(self.szene.sonden)} bei ({i}, {j}) gesetzt.", "gut")

    def _sonde_loeschen(self):
        k = self._index(self.tv_s)
        if k is None:
            return
        self.szene.sonden.pop(k)
        self._szene_anwenden(self.szene, "Sonde geloescht.")

    # ========================================================== Lauf/Takt
    def _lauf_uebernehmen(self):
        self.motor.chunk = max(int(self.v_chunk.get()), 1)
        self.motor.max_rate = float(self.v_rate.get())
        self.motor.grenze = max(zahl(self.v_grenze, 0, True), 0)
        self.motor.intens_an = bool(self.v_intens.get())
        self.lb_chunk.config(text=str(self.motor.chunk))
        self.lb_rate.config(text="aus" if self.motor.max_rate <= 0
                            else f"{int(self.motor.max_rate)}/s")

    def _umschalten(self):
        if self.motor.laeuft:
            self.motor.pause()
            self.btn_start.config(text="▶  Start")
        else:
            if self.motor.instabil:
                self._log("Erst zuruecksetzen - die Rechnung ist divergiert.", "warn")
                return
            self._lauf_uebernehmen()
            self.motor.start()
            self.btn_start.config(text="⏸  Pause")

    def _schritt(self, n=1):
        self.motor.pause()
        self.btn_start.config(text="▶  Start")
        self.motor.einzelschritt(n)

    def _reset(self):
        self.motor.zuruecksetzen()
        self._spitze = 1e-12
        self.btn_start.config(text="▶  Start")
        self._log("Felder zurueckgesetzt.", "info")
        self._diagramme_auffrischen()

    def _aufnahme_schalten(self):
        self.motor.aufnahme = bool(self.v_aufnahme.get())
        if self.motor.aufnahme:
            self._log(f"Aufnahme laeuft (jedes {self.motor.film_jede}. Bild, "
                      f"max. {Motor.FILM_MAX}).", "info")
        else:
            self._log(f"Aufnahme gestoppt: {len(self.motor.filme)} Bilder.", "info")

    def _takt(self):
        try:
            self._ansicht_auffrischen()
        except Exception as fehler:                      # noqa: BLE001
            self._log(f"Anzeige: {fehler}", "fehler")
        self.after(30, self._takt)

    def _ansicht_auffrischen(self):
        while self._meldungen:
            text, art = self._meldungen.pop(0)
            self._log(text, art)

        m = self.motor
        if m.sim is None:
            return
        with m.lock:
            schritt, t = m.sim.time_step, m.sim.t
            energie = m.energie_h[-1] if m.energie_h else 0.0
        spitze = float(np.abs(m.feld).max()) if m.feld is not None else 0.0
        self.st_schritt.set(f"{schritt:>7d}")
        self.st_zeit.set(f"{t * 1e9:8.3f} ns")
        self.st_rate.set(f"{m.rate:7.0f}/s")
        self.st_spitze.set(f"{spitze:9.3e}")
        self.st_energie.set(f"{energie:9.3e} J/m")
        self.lb_film.config(text=f"{len(m.filme)} Bilder")
        if m.instabil:
            self.st_hinweis.set("DIVERGIERT - S > 1?")
            self.btn_start.config(text="▶  Start")
        elif m.grenze and schritt >= m.grenze:
            self.st_hinweis.set(f"Stopp bei Schritt {m.grenze}")
        else:
            self.st_hinweis.set("")

        aktiv = self.nb_ansicht.index(self.nb_ansicht.select())
        jetzt = time.perf_counter()
        if aktiv == 0:
            self._feld_auffrischen()
        elif aktiv == 1 and jetzt - self._letzte_analyse > self.v_analyse_ms.get() / 1000:
            self._letzte_analyse = jetzt
            self._analyse_zeichnen()
        elif aktiv == 2 and jetzt - self._letzte_intens > self.v_analyse_ms.get() / 1000:
            self._letzte_intens = jetzt
            self._intens_zeichnen()

    def _diagramme_auffrischen(self):
        self._analyse_zeichnen()
        self._intens_zeichnen()

    # =============================================================== Feld
    def feld_daten(self):
        """Die gerade gewaehlte Feldgroesse als Array (nx, ny)."""
        m = self.motor
        art = self.v_darstellung.get()
        with m.lock:
            if m.feld is None:
                return None
            ez = m.feld
            if art.startswith("Intensitaet"):
                d = (m.intens / m.intens_n).copy() if m.intens_n else np.zeros_like(ez)
            elif art == "el. Energiedichte":
                d = 0.5 * EPS0 * m.eps * ez * ez
            elif art == "|Ez|":
                d = np.abs(ez)
            else:
                d = ez
        if art.endswith("dB"):
            spitze = float(np.abs(d).max())
            spitze = spitze if spitze > 0 else 1.0
            d = 20.0 * np.log10(np.abs(d) / spitze + 1e-12)
            d = np.clip(d, -60.0, 0.0)
        return d

    def _skala(self, daten):
        art = self.v_darstellung.get()
        if art.endswith("dB"):
            return -60.0, 0.0
        if self.v_skala.get() == "fest":
            v = abs(zahl(self.v_vmax, 1.0)) or 1.0
        else:
            v = float(np.nanmax(np.abs(daten)))
            self._spitze = max(v, self._spitze * 0.995)
            v = self._spitze if self._spitze > 0 else 1.0
        return (-v, v) if art in VORZEICHENBEHAFTET else (0.0, v)

    def _overlay(self):
        sz = self.szene
        nx, ny = sz.nx, sz.ny
        rgba = np.zeros((nx, ny, 4))
        if self.v_ov_pml.get() and sz.npml > 0:
            rand = np.zeros((nx, ny), dtype=bool)
            rand[:sz.npml, :] = rand[-sz.npml:, :] = True
            rand[:, :sz.npml] = rand[:, -sz.npml:] = True
            rgba[rand] = (0.62, 0.68, 0.78, 0.13)
        if self.v_ov_mat.get():
            eps, sig, pec = self.motor.eps, self.motor.sig, self.motor.pec
            if eps is not None:
                hoch = eps > 1.0 + 1e-9
                if hoch.any():
                    staerke = (eps - 1.0) / max(float(eps.max()) - 1.0, 1e-9)
                    rgba[hoch, 0] = 0.20
                    rgba[hoch, 1] = 0.85
                    rgba[hoch, 2] = 0.70
                    rgba[hoch, 3] = 0.10 + 0.25 * staerke[hoch]
                verlust = sig > 0
                if verlust.any():
                    rgba[verlust] = (0.95, 0.65, 0.25, 0.30)
                if pec is not None and pec.any():
                    rgba[pec] = (0.24, 0.27, 0.33, 1.0)
        return np.transpose(rgba, (1, 0, 2))

    def _feld_aufbauen(self):
        sz = self.szene
        self.fig_feld.clf()
        self.ax_feld = self.fig_feld.add_subplot(111)
        ax = self.ax_feld
        daten = self.feld_daten()
        if daten is None:
            return
        vmin, vmax = self._skala(daten)
        self.im_feld = ax.imshow(daten.T, origin="lower", extent=(0, sz.nx, 0, sz.ny),
                                 cmap=self.v_cmap.get(), vmin=vmin, vmax=vmax,
                                 interpolation="nearest", animated=True, zorder=1)
        self.im_over = ax.imshow(self._overlay(), origin="lower",
                                 extent=(0, sz.nx, 0, sz.ny),
                                 interpolation="nearest", animated=True, zorder=3)
        self.cb_feld = self.fig_feld.colorbar(self.im_feld, ax=ax, fraction=0.045, pad=0.02)
        self.cb_feld.ax.tick_params(labelsize=7)

        qi, qj = [], []
        for q in self.szene.quellen:
            if not q.aktiv:
                continue
            if q.typ == "punkt":
                qi.append(q.i)
                qj.append(q.j)
            else:
                qi.extend([q.i] * 9)
                qj.extend(np.linspace(0, sz.ny - 1, 9))
        self.ln_quellen, = ax.plot(qi, qj, "x", color="#67e8f9", ms=5, mew=1.3,
                                   animated=True, zorder=5,
                                   visible=bool(self.v_ov_qu.get()))
        si = [s.i for s in self.szene.sonden]
        sj = [s.j for s in self.szene.sonden]
        self.ln_sonden, = ax.plot(si, sj, "o", mfc="none", mec="#ffd166", ms=8, mew=1.5,
                                  animated=True, zorder=5,
                                  visible=bool(self.v_ov_so.get()))
        self.ln_schirm = ax.axvline(zahl(self.v_schirm, sz.nx - 20, True),
                                    color=FARBE["gut"], ls="--", lw=1.1, animated=True,
                                    zorder=5, visible=bool(self.v_ov_schirm.get()))
        self.band_r = Rectangle((0, 0), 0, 0, fill=False, ec=FARBE["akzent"],
                                lw=1.2, ls="--", animated=True, zorder=6, visible=False)
        self.band_e = Ellipse((0, 0), 0, 0, fill=False, ec=FARBE["akzent"],
                              lw=1.2, ls="--", animated=True, zorder=6, visible=False)
        ax.add_patch(self.band_r)
        ax.add_patch(self.band_e)

        ax.set_xlabel(f"Gitterindex i   (Delta = {sz.d*1e3:g} mm)")
        ax.set_ylabel("Gitterindex j")
        ax.set_title(self.v_darstellung.get())
        self.fig_feld.tight_layout()
        self.kanvas_feld.draw()

    def _auf_draw(self, event):
        try:
            self._bg = self.kanvas_feld.copy_from_bbox(self.ax_feld.bbox)
        except Exception:                                # noqa: BLE001
            self._bg = None

    def _feld_auffrischen(self):
        daten = self.feld_daten()
        if daten is None or not hasattr(self, "im_feld"):
            return
        if daten.shape != self.im_feld.get_array().shape[::-1]:
            self._feld_aufbauen()
            return
        vmin, vmax = self._skala(daten)
        a, b = self.im_feld.get_clim()
        voll = False
        spanne = max(abs(b - a), 1e-30)
        if abs(vmin - a) / spanne > 0.04 or abs(vmax - b) / spanne > 0.04:
            self.im_feld.set_clim(vmin, vmax)
            jetzt = time.perf_counter()
            if jetzt - self._letzte_cb > 0.4:
                self._letzte_cb = jetzt
                voll = True
        self.im_feld.set_data(daten.T)
        self.ln_schirm.set_xdata([zahl(self.v_schirm, 0, True)] * 2)

        if voll or self._bg is None:
            self.kanvas_feld.draw()
        if self._bg is None:
            return
        self.kanvas_feld.restore_region(self._bg)
        for kuenstler in (self.im_feld, self.im_over, self.ln_quellen,
                          self.ln_sonden, self.ln_schirm, self.band_r, self.band_e):
            self.ax_feld.draw_artist(kuenstler)
        self.kanvas_feld.blit(self.ax_feld.bbox)

    def _darstellung_wechsel(self):
        art = self.v_darstellung.get()
        if art in VORZEICHENBEHAFTET:
            self.v_cmap.set("RdBu_r")
        elif self.v_cmap.get() in ("RdBu_r", "seismic", "bwr", "coolwarm"):
            self.v_cmap.set("inferno")
        self._spitze = 1e-12
        self._feld_aufbauen()

    def _cmap_wechsel(self):
        if hasattr(self, "im_feld"):
            self.im_feld.set_cmap(self.v_cmap.get())
            self.kanvas_feld.draw()

    def _overlay_neu(self):
        if not hasattr(self, "im_over"):
            return
        self.im_over.set_data(self._overlay())
        self.ln_quellen.set_visible(bool(self.v_ov_qu.get()))
        self.ln_sonden.set_visible(bool(self.v_ov_so.get()))
        self.ln_schirm.set_visible(bool(self.v_ov_schirm.get()))
        self.kanvas_feld.draw()

    # ============================================================= Analyse
    @staticmethod
    def _spektrum(y, dt, mitteln=True, segmente=8):
        """Betragsspektrum eines Sondensignals.

        Ohne Fensterung leckt ein abgeschnittener Dauerstrich ueber das ganze
        Band (das ist der Zappel in einem rohen FFT-Bild), deshalb immer ein
        Hann-Fenster. Bei langen Aufzeichnungen wird zusaetzlich nach Welch
        gemittelt: das Signal wird in ueberlappende Abschnitte zerlegt, jedes
        Betragsquadrat gemittelt. Das kostet Frequenzaufloesung, macht die
        Kurve aber glatt und die Pegel ablesbar."""
        n = y.size
        if n < 64:
            return None, None
        if not mitteln or n < 512:
            fenster = np.hanning(n)
            betrag = np.abs(np.fft.rfft(y * fenster)) * 2.0 / fenster.sum()
            return np.fft.rfftfreq(n, dt), betrag
        laenge = max(128, n // max(int(segmente), 1))
        laenge -= laenge % 2
        fenster = np.hanning(laenge)
        summe = np.zeros(laenge // 2 + 1)
        anzahl = 0
        for a in range(0, n - laenge + 1, laenge // 2):
            summe += np.abs(np.fft.rfft(y[a:a + laenge] * fenster)) ** 2
            anzahl += 1
        betrag = np.sqrt(summe / max(anzahl, 1)) * 2.0 / fenster.sum()
        return np.fft.rfftfreq(laenge, dt), betrag

    def _fenster_start(self, t):
        """Anfangszeit des gewaehlten Zeitfensters (None = ganzer Lauf)."""
        wahl = self.v_zeitfenster.get()
        if t.size == 0 or wahl.startswith("ganzer"):
            return None
        f0 = self.szene.quellen[0].f0 if self.szene.quellen else 0.0
        teile = wahl.split()
        try:
            anzahl = float(teile[1])
        except (IndexError, ValueError):
            return None
        if wahl.endswith("Perioden"):
            if f0 <= 0:
                return None
            spanne = anzahl / f0
        else:
            spanne = anzahl * self.motor.sim.dt
        return t[-1] - spanne if t[-1] - spanne > t[0] else None

    @staticmethod
    def _ausschnitt(t, y, t0):
        if t0 is None or t.size == 0:
            return t, y
        k = int(np.searchsorted(t, t0))
        return t[k:], y[k:]

    @staticmethod
    def _zeitkurve(ax, t, y, farbe, label=None, punkte=900, log=False,
                   block_min=1, mittel=True, rand=False):
        """Zeichnet y(t). Liegen mehr Abtastwerte vor als Bildpunkte, wird
        statt einer (aliasenden) Ausduennung die Huellkurve aus Minimum und
        Maximum je Block gezeichnet - so bleibt auch ein langer Lauf lesbar.

        block_min erzwingt eine Mindestblocklaenge (bei Schwingungen eine
        volle Periode, sonst wellt der Rand der Huellkurve).
        Rueckgabe: True, wenn die Huellkurve benutzt wurde."""
        if t.size == 0:
            return False
        zeichner = ax.semilogy if log else ax.plot
        k = max(t.size // punkte, int(block_min), 1)
        if t.size <= 2 * punkte or t.size <= 4 * k:
            zeichner(t * 1e9, y, lw=1.0, color=farbe, label=label)
            return False
        anzahl = t.size // k
        bloecke = y[:anzahl * k].reshape(anzahl, k)
        unten, oben = bloecke.min(axis=1), bloecke.max(axis=1)
        mitte = bloecke.mean(axis=1)
        schwanz = y[(anzahl - 1) * k:]            # letzter Block nimmt den Rest
        unten[-1], oben[-1] = schwanz.min(), schwanz.max()
        mitte[-1] = schwanz.mean()
        # Stuetzstelle am Ende anhaengen, damit die Kurve bis t[-1] laeuft
        tt = np.concatenate([t[:anzahl * k:k], t[-1:]]) * 1e9
        unten = np.append(unten, unten[-1])
        oben = np.append(oben, oben[-1])
        mitte = np.append(mitte, mitte[-1])
        ax.fill_between(tt, unten, oben, color=farbe,
                        alpha=0.20 if rand else 0.40, lw=0, label=label)
        if rand:                                  # Raender zeichnen, damit sich
            zeichner(tt, oben, lw=0.9, color=farbe)    # mehrere Sonden nicht
            zeichner(tt, unten, lw=0.9, color=farbe)   # zu Matsch ueberlagern
        if mittel:
            zeichner(tt, mitte, lw=0.9, color=farbe)
        return True

    def _analyse_zeichnen(self):
        """Drei Diagramme: Oszilloskop, Balkenspektrum, Energie."""
        m = self.motor
        if m.sim is None:
            return
        farben = ["#4aa3df", "#e5c07b", "#4cc38a", "#e06c75", "#c678dd", "#56b6c2"]
        with m.lock:
            dt = m.sim.dt
            daten = [sondendaten(m.sim, k) for k in range(len(self.szene.sonden))]
        f0 = self.szene.quellen[0].f0 if self.szene.quellen else 0.0

        self._zeichne_oszilloskop(self.ax_zeit, daten, dt, f0, farben)
        self._zeichne_spektrum(self.ax_spek, daten, dt, f0, farben)
        self._zeichne_energie(self.ax_ener)
        self.kanvas_an.draw_idle()

    # -- 1. Oszilloskop ----------------------------------------------------
    def _zeichne_oszilloskop(self, ax, daten, dt, f0, farben):
        ax.clear()
        t0 = self._fenster_start(daten[0][0]) if daten else None
        pro_periode = int(round(1.0 / (f0 * dt))) if f0 > 0 else 1
        gehuellt = False
        for k, (t, ez) in enumerate(daten):
            tt, ee = self._ausschnitt(t, ez, t0)
            gehuellt |= self._zeitkurve(ax, tt, ee, farben[k % len(farben)],
                                        f"Sonde {k+1}", block_min=pro_periode,
                                        mittel=False, rand=True)
        ax.axhline(0.0, color=FARBE["linie"], lw=0.8)
        # Ablenkung am Lauf messen, nicht nur am Fenster: sonst wird nach dem
        # Durchlauf das numerische Restrauschen formatfuellend aufgezogen.
        # Der Boden liegt bei 2 % der groessten je gemessenen Feldstaerke.
        gesamt = max((float(np.abs(ez).max()) for _, ez in daten if ez.size),
                     default=0.0)
        sichtbar = 0.0
        for t, ez in daten:
            _, ee = self._ausschnitt(t, ez, t0)
            if ee.size:
                sichtbar = max(sichtbar, float(np.abs(ee).max()))
        ablenkung = max(sichtbar, 0.02 * gesamt)
        if ablenkung > 0:
            ax.set_ylim(-1.15 * ablenkung, 1.15 * ablenkung)
        ax.set_xlabel("Zeit t in ns")
        ax.set_ylabel("elektrisches Feld $E_z$")
        ax.set_title("Oszilloskop: Feldstaerke an den Sonden  -  "
                     + self.v_zeitfenster.get()
                     + (" (Huellkurve)" if gehuellt else ""))
        ax.grid(True, alpha=0.3)
        if daten and daten[0][0].size:
            zeit = self._ausschnitt(daten[0][0], daten[0][1], t0)[0]
            if zeit.size > 1:
                ax.set_xlim(zeit[0] * 1e9, zeit[-1] * 1e9)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="upper right", ncol=max(1, len(daten)))

    # -- 2. Balkenspektrum -------------------------------------------------
    def _zeichne_spektrum(self, ax, daten, dt, f0, farben, balken=28):
        """Klassisches Balkenspektrum: das Band wird in gleich breite
        Streifen geteilt, je Streifen wird die mittlere Amplitude als Saeule
        gezeichnet. Das ist gegen Rauschen unempfindlich und direkt ablesbar."""
        ax.clear()
        # Beim Dauerstrich sitzt alles bei f0, beim Puls reicht das Band bis
        # etwa 3 f0 - der Ausschnitt richtet sich danach.
        dauerstrich = bool(self.szene.quellen) and self.szene.quellen[0].art == "sinus"
        weite = 2.0 if dauerstrich else 3.0
        f_max = min(weite * f0, 0.5 / dt) if f0 > 0 else 0.5 / dt
        kanten = np.linspace(0.0, f_max, balken + 1)
        mitte = 0.5 * (kanten[:-1] + kanten[1:]) / 1e9
        breite = (kanten[1] - kanten[0]) / 1e9

        saeulen = []
        f_spitze, staerkste = 0.0, 0.0
        for k, (t, ez) in enumerate(daten[:4]):        # mehr wird unleserlich
            f, betrag = self._spektrum(ez, dt, mitteln=True)
            if f is None or betrag.max() <= 0:
                continue
            # Lage des Maximums aus dem feinen Spektrum der staerksten
            # Sonde ablesen, nicht aus der Balkenmitte.
            if float(betrag.max()) > staerkste:
                staerkste = float(betrag.max())
                f_spitze = float(f[int(np.argmax(betrag))])
            fach = np.digitize(f, kanten) - 1
            gilt = (fach >= 0) & (fach < balken)
            summe = np.bincount(fach[gilt], weights=betrag[gilt] ** 2,
                                minlength=balken)
            zahl = np.bincount(fach[gilt], minlength=balken)
            saeulen.append((k, np.sqrt(summe / np.maximum(zahl, 1))))

        if not saeulen:
            ax.set_title("Frequenzspektrum")
            ax.set_xlabel("Frequenz f in GHz")
            ax.set_ylabel("Amplitude (normiert)")
            return

        bezug = max(float(h.max()) for _, h in saeulen) or 1.0
        n = len(saeulen)
        for j, (k, hoehe) in enumerate(saeulen):
            versatz = (j - (n - 1) / 2.0) * breite * 0.82 / n
            ax.bar(mitte + versatz, hoehe / bezug, width=breite * 0.82 / n,
                   color=farben[k % len(farben)], edgecolor="none",
                   label=f"Sonde {k+1}")

        if f0:
            ax.axvline(f0 / 1e9, color=FARBE["grau"], ls=":", lw=1.2)
            ax.annotate("Quelle", (f0 / 1e9, 1.1), ha="center", fontsize=7.5,
                        color=FARBE["grau"])
        ax.set_ylim(0, 1.2)
        ax.set_xlim(0, f_max / 1e9)
        ax.set_xlabel("Frequenz f in GHz")
        ax.set_ylabel("Amplitude (normiert)")
        ax.set_title(f"Frequenzspektrum  -  Maximum bei {f_spitze/1e9:.1f} GHz")
        ax.grid(True, axis="y", alpha=0.3)
        if n > 1:
            ax.legend(loc="upper right", fontsize=7.5)

    # -- 3. Gesamtenergie --------------------------------------------------
    def _zeichne_energie(self, ax):
        ax.clear()
        m = self.motor
        with m.lock:
            t = np.asarray(m.zeit_h)
            w = np.asarray(m.energie_h)
        if t.size > 1:
            tt, ww = self._ausschnitt(t, np.maximum(w, 1e-30),
                                      self._fenster_start(t))
            self._zeitkurve(ax, tt, ww, "#4cc38a", log=True)
            gross, klein = float(np.max(ww)), float(np.min(ww))
            if gross > 0:
                ax.set_ylim(max(klein * 0.5, gross * 1e-8), gross * 2.0)
            if tt.size > 1:
                ax.set_xlim(tt[0] * 1e9, tt[-1] * 1e9)
        ax.set_xlabel("Zeit t in ns")
        ax.set_ylabel("Energie W in J/m")
        ax.set_title("Gesamtenergie im Rechengebiet")
        ax.grid(True, alpha=0.3, which="both")

    def _intens_zeichnen(self):
        m = self.motor
        if m.sim is None:
            return
        sz = self.szene
        with m.lock:
            bild = (m.intens / m.intens_n).copy() if m.intens_n else np.zeros((sz.nx, sz.ny))
            n_mittel = m.intens_n
        i_schirm = int(np.clip(zahl(self.v_schirm, sz.nx - 20, True), 0, sz.nx - 1))
        bezug = self.v_intens_norm.get()

        if bezug.startswith("Schirm"):
            # Auf die Umgebung des Schirms normieren: sonst erschlaegt die
            # einlaufende Welle das viel schwaechere Beugungsbild dahinter.
            fenster = bild[max(i_schirm - 15, 0):i_schirm + 15, :]
            spitze = float(fenster.max())
        else:
            spitze = float(bild.max())
        if spitze <= 0:
            spitze = 1.0
        norm = bild / spitze
        log = bezug.startswith("log")
        if log:
            norm = np.clip(10 * np.log10(norm + 1e-12), -50.0, 0.0)
            vmin, vmax = -50.0, 0.0
        else:
            vmin, vmax = 0.0, 1.0

        ax = self.ax_in_bild
        ax.clear()
        bildchen = ax.imshow(norm.T, origin="lower", extent=(0, sz.nx, 0, sz.ny),
                             cmap="inferno", vmin=vmin, vmax=vmax,
                             interpolation="nearest")
        ax.imshow(self._overlay(), origin="lower", extent=(0, sz.nx, 0, sz.ny),
                  interpolation="nearest")
        ax.axvline(i_schirm, color=FARBE["gut"], ls="--", lw=1.1)
        ax.set_xlabel("Gitterindex i")
        ax.set_ylabel("Gitterindex j")
        ax.set_title("zeitgemittelte Intensitaet <Ez^2>\n"
                     f"{n_mittel} Zeitschritte, bezogen auf {bezug}", fontsize=9)
        if getattr(self, "_cb_in", None) is None:
            self._cb_in = self.fig_in.colorbar(bildchen, ax=ax, fraction=0.045, pad=0.02)
            self._cb_in.ax.tick_params(labelsize=7)
        else:
            self._cb_in.update_normal(bildchen)

        ax = self.ax_in_prof
        ax.clear()
        schnitt = bild[i_schirm, :]
        spitze_schirm = float(schnitt.max()) or 1.0
        profil = schnitt / spitze_schirm
        if log:
            profil = np.clip(10 * np.log10(profil + 1e-12), -50.0, 0.0)
        ax.plot(profil, np.arange(sz.ny), lw=1.3, color="#ffd166")
        ax.set_ylim(0, sz.ny)
        if log:
            ax.set_xlim(-52.0, 3.0)
        else:
            ax.set_xlim(-0.02, 1.18)
        ax.set_xlabel("Intensitaet in dB" if log else "Intensitaet (auf Schirmmaximum normiert)")
        ax.set_ylabel("Gitterindex j")
        ax.set_title(f"Schirmprofil bei i = {i_schirm}")
        ax.grid(True, alpha=0.3)

        # Lage der Maxima ausweisen - das ist beim Doppelspalt die eigentliche
        # Aussage des Diagramms.
        if not log:
            innen = profil[5:-5]
            spitzen = [k + 5 for k in range(1, innen.size - 1)
                       if innen[k] > innen[k - 1] and innen[k] >= innen[k + 1]
                       and innen[k] > 0.15]
            for j in spitzen[:9]:
                ax.annotate(f"j = {j}", (profil[j], j), fontsize=7,
                            color=FARBE["grau"], xytext=(-4, 4),
                            textcoords="offset points", ha="right")

        self.kanvas_in.draw_idle()

    # =============================================================== Maus
    def _maus_indizes(self, event):
        if event.inaxes is not self.ax_feld or event.xdata is None:
            return None
        i = int(np.clip(event.xdata, 0, self.szene.nx - 1))
        j = int(np.clip(event.ydata, 0, self.szene.ny - 1))
        return i, j

    def _maus_druck(self, event):
        pos = self._maus_indizes(event)
        if pos is None or event.button != 1:
            return
        if self.kanvas_feld.toolbar is not None and self.kanvas_feld.toolbar.mode:
            return                                   # Zoom/Pan der mpl-Leiste
        i, j = pos
        wz = self.v_werkzeug.get()
        if wz == "sonde":
            self._sonde_hinzu(i, j)
        elif wz == "quelle":
            q = self._quelle_aus_vars()
            q.typ, q.i, q.j = "punkt", i, j
            self.szene.quellen.append(q)
            self._szene_anwenden(self.szene, f"Punktquelle bei ({i}, {j}).")
        elif wz == "radierer":
            for k in range(len(self.szene.objekte) - 1, -1, -1):
                if self.szene.objekte[k].enthaelt(i, j, self.szene.nx, self.szene.ny):
                    self.szene.objekte.pop(k)
                    self._szene_anwenden(self.szene, "Objekt entfernt.")
                    return
            self._log("Hier liegt kein Objekt.", "info")
        elif wz in ("rechteck", "ellipse", "pec_rechteck", "pec_ellipse"):
            self._zieh_start = (i, j)

    def _maus_bewegung(self, event):
        pos = self._maus_indizes(event)
        if pos is None:
            self.st_cursor.set("")
            return
        i, j = pos
        m = self.motor
        wert = float(m.feld[i, j]) if m.feld is not None else 0.0
        eps = float(m.eps[i, j]) if m.eps is not None else 1.0
        self.st_cursor.set(f"i={i:4d} j={j:4d} | x={i*self.szene.d*1e3:7.2f} mm "
                           f"y={j*self.szene.d*1e3:7.2f} mm | Ez={wert:+.3e} | eps_r={eps:g}")
        if self._zieh_start is None:
            return
        i0, j0 = self._zieh_start
        wz = self.v_werkzeug.get()
        if wz.endswith("rechteck"):
            self.band_r.set_bounds(min(i0, i), min(j0, j), abs(i - i0), abs(j - j0))
            self.band_r.set_visible(True)
        else:
            self.band_e.set_center((0.5 * (i0 + i), 0.5 * (j0 + j)))
            self.band_e.set_width(abs(i - i0))
            self.band_e.set_height(abs(j - j0))
            self.band_e.set_visible(True)

    def _maus_los(self, event):
        if self._zieh_start is None:
            return
        pos = self._maus_indizes(event)
        i0, j0 = self._zieh_start
        self._zieh_start = None
        self.band_r.set_visible(False)
        self.band_e.set_visible(False)
        if pos is None:
            return
        i, j = pos
        if abs(i - i0) < 2 or abs(j - j0) < 2:
            self._log("Zu kleiner Bereich - Objekt nicht angelegt.", "info")
            return
        wz = self.v_werkzeug.get()
        o = Objekt(form="ellipse" if wz.endswith("ellipse") else "rechteck",
                   art="pec" if wz.startswith("pec") else "dielektrikum",
                   i0=min(i0, i), i1=max(i0, i), j0=min(j0, j), j1=max(j0, j),
                   eps_r=max(zahl(self.v_o_eps, 2.25), 1e-3),
                   sigma=max(zahl(self.v_o_sig, 0.0), 0.0))
        self.szene.objekte.append(o)
        self._szene_anwenden(self.szene, f"{o.zeile()[0]} angelegt: {o.zeile()[1]}.")

    # ============================================================= Messen
    def _sonden_pruefen(self, anzahl=2):
        if len(self.szene.sonden) < anzahl:
            messagebox.showinfo("Sonden fehlen",
                                f"Dafuer werden mindestens {anzahl} Sonden gebraucht.")
            return False
        with self.motor.lock:
            t, _ = sondendaten(self.motor.sim, 0)
        if t.size < 32:
            messagebox.showinfo("Zu wenig Daten", "Erst ein Stueck weit rechnen lassen.")
            return False
        return True

    def _mess_phasengeschwindigkeit(self):
        if not self._sonden_pruefen(2):
            return
        s1, s2 = self.szene.sonden[0], self.szene.sonden[1]
        abstand = self.szene.d * np.hypot(s2.i - s1.i, s2.j - s1.j)
        if abstand <= 0:
            self._log("Die beiden Sonden liegen aufeinander.", "warn")
            return
        with self.motor.lock:
            f, Y1 = self.motor.sim.get_spectrum(0, pad=4)
            _, Y2 = self.motor.sim.get_spectrum(1, pad=4)
        gut = np.where(np.abs(Y1) > 0.05 * np.abs(Y1).max())[0]
        if gut.size < 8:
            self._log("Spektrum zu schwach fuer eine Auswertung.", "warn")
            return
        aus = slice(gut[0], gut[-1] + 1)
        fb = f[aus]
        phi = -np.unwrap(np.angle(Y2[aus] / Y1[aus]))
        phi += 2 * np.pi * np.round((2 * np.pi * fb[0] / C0 * abstand - phi[0]) / (2 * np.pi))
        with np.errstate(divide="ignore", invalid="ignore"):
            vp = 2 * np.pi * fb * abstand / phi
        n_lambda = C0 / (fb * self.szene.d)
        f0 = self.szene.quellen[0].f0 if self.szene.quellen else 0.0
        self._log(f"Phasengeschwindigkeit ueber {abstand*1e3:.1f} mm:", "gut")
        if f0 and fb[0] <= f0 <= fb[-1]:
            self._log(f"   bei f0 = {f0/1e9:g} GHz:  v_p/c = "
                      f"{np.interp(f0, fb, vp)/C0:.5f}")
        for ziel in (10.0, 20.0):
            if n_lambda.min() <= ziel <= n_lambda.max():
                v = np.interp(ziel, n_lambda[::-1], vp[::-1])
                self._log(f"   bei N_lambda = {ziel:g}:   v_p/c = {v/C0:.5f}")

    def _mess_uebertragung(self):
        if not self._sonden_pruefen(2):
            return
        with self.motor.lock:
            f, Y1 = self.motor.sim.get_spectrum(0, pad=4)
            _, Y2 = self.motor.sim.get_spectrum(1, pad=4)
        gut = np.abs(Y1) > 0.05 * np.abs(Y1).max()
        if gut.sum() < 8:
            self._log("Spektrum zu schwach fuer eine Auswertung.", "warn")
            return
        betrag = np.abs(Y2[gut] / Y1[gut])
        f0 = self.szene.quellen[0].f0 if self.szene.quellen else 0.0
        self._log("Uebertragung |E2/E1|:", "gut")
        if f0:
            self._log(f"   bei f0 = {f0/1e9:g} GHz: "
                      f"{np.interp(f0, f[gut], betrag):.4f} "
                      f"({20*np.log10(np.interp(f0, f[gut], betrag) + 1e-12):+.2f} dB)")
        self._log(f"   Minimum {betrag.min():.4f}, Maximum {betrag.max():.4f} "
                  f"im Band {f[gut][0]/1e9:.1f} .. {f[gut][-1]/1e9:.1f} GHz")

    def _mess_kennzahlen(self):
        if not self.szene.sonden:
            messagebox.showinfo("Sonden fehlen", "Es ist keine Sonde gesetzt.")
            return
        self._log("Sondenkennzahlen:", "gut")
        with self.motor.lock:
            for k, s in enumerate(self.szene.sonden):
                t, ez = sondendaten(self.motor.sim, k)
                if t.size == 0:
                    continue
                effektiv = float(np.sqrt(np.mean(ez ** 2)))
                self._log(f"   Sonde {k+1} ({s.i:>4d},{s.j:>4d}): "
                          f"Spitze {np.abs(ez).max():.4e}, Effektivwert {effektiv:.4e}, "
                          f"Spitze bei t = {t[int(np.argmax(np.abs(ez)))]*1e9:.3f} ns")

    def _mess_energie(self):
        m = self.motor
        with m.lock:
            w = np.asarray(m.energie_h)
            t = np.asarray(m.zeit_h)
        if w.size < 3:
            messagebox.showinfo("Zu wenig Daten", "Erst ein Stueck weit rechnen lassen.")
            return
        self._log("Energiebilanz  W = Delta^2 sum( eps0 eps_r (E^n)^2 / 2 "
                  "+ mu0 H^(n-1/2).H^(n+1/2) / 2 ):", "gut")
        self._log(f"   Maximum {w.max():.4e} J/m bei t = {t[int(np.argmax(w))]*1e9:.3f} ns")
        self._log(f"   jetzt   {w[-1]:.4e} J/m  = {w[-1]/max(w.max(),1e-30)*100:.3f} % "
                  f"des Maximums")
        spanne = (w.max() - w.min()) / max(w.mean(), 1e-30)
        if self.szene.npml > 0:
            self._log("   Der Abfall ist die von der PML aufgenommene Leistung.")
        else:
            self._log(f"   Schwankung {spanne:.2e} (ohne PML und nach dem Abschalten "
                      f"der Quelle muss W erhalten bleiben).")
        self._log("   H wird zwischen den beiden Halbschritten verschraenkt; "
                  "sim.total_energy() nimmt H^(n+1/2) quadratisch zu E^(n+1) "
                  "und schwankt daher mit der Schwingung.")

    # ============================================================= Export
    def _bild_speichern(self):
        aktiv = self.nb_ansicht.index(self.nb_ansicht.select())
        pfad = filedialog.asksaveasfilename(
            title="Ansicht speichern", initialdir=str(ABB if ABB.exists() else ORDNER),
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if not pfad:
            return
        if aktiv == 0:
            fig = self._feld_standbild()
            fig.savefig(pfad, bbox_inches="tight")
        else:
            fig = self.fig_an if aktiv == 1 else self.fig_in
            fig.savefig(pfad, bbox_inches="tight")
        self._log(f"Bild gespeichert: {pfad}", "gut")

    def _feld_standbild(self):
        """Eigene Figur fuer den Export (die Live-Figur nutzt animated-Artists)."""
        sz = self.szene
        daten = self.feld_daten()
        vmin, vmax = self._skala(daten)
        fig = Figure(figsize=(7.5, 6.0))
        ax = fig.add_subplot(111)
        bild = ax.imshow(daten.T, origin="lower", extent=(0, sz.nx, 0, sz.ny),
                         cmap=self.v_cmap.get(), vmin=vmin, vmax=vmax,
                         interpolation="nearest")
        ax.imshow(self._overlay(), origin="lower", extent=(0, sz.nx, 0, sz.ny),
                  interpolation="nearest")
        if self.v_ov_so.get():
            for k, s in enumerate(sz.sonden):
                ax.plot(s.i, s.j, "o", mfc="none", mec="#ffd166", ms=8, mew=1.5)
                ax.annotate(str(k + 1), (s.i, s.j), color="#ffd166",
                            fontsize=7, xytext=(5, 4), textcoords="offset points")
        fig.colorbar(bild, ax=ax, fraction=0.045, pad=0.02)
        ax.set_xlabel(f"Gitterindex i   (Delta = {sz.d*1e3:g} mm)")
        ax.set_ylabel("Gitterindex j")
        ax.set_title(f"{self.v_darstellung.get()}   -   Schritt "
                     f"{self.motor.sim.time_step}, t = {self.motor.sim.t*1e9:.3f} ns")
        fig.tight_layout()
        return fig

    def _csv_speichern(self):
        if not self.szene.sonden:
            messagebox.showinfo("Sonden fehlen", "Es ist keine Sonde gesetzt.")
            return
        pfad = filedialog.asksaveasfilename(
            title="Sondendaten speichern", defaultextension=".csv",
            initialdir=str(ORDNER), filetypes=[("CSV", "*.csv")])
        if not pfad:
            return
        with self.motor.lock:
            reihen = [sondendaten(self.motor.sim, k) for k in range(len(self.szene.sonden))]
        n = max(r[0].size for r in reihen)
        if n == 0:
            messagebox.showinfo("Keine Daten", "Die Sonden haben noch nichts aufgezeichnet.")
            return
        # Spaeter gesetzte Sonden haben kuerzere Reihen; sie enden aber zum
        # gleichen Zeitschritt, also hinten buendig ausrichten und vorn mit
        # NaN auffuellen.
        zeit = next(r[0] for r in reihen if r[0].size == n)
        spalten = [zeit]
        for _, ez in reihen:
            spalte = np.full(n, np.nan)
            spalte[n - ez.size:] = ez
            spalten.append(spalte)
        kopf = "t_s," + ",".join(f"Ez_S{k+1}_i{s.i}_j{s.j}"
                                 for k, s in enumerate(self.szene.sonden))
        np.savetxt(pfad, np.column_stack(spalten), delimiter=",", header=kopf,
                   comments="", fmt="%.8e")
        self._log(f"{n} Zeitpunkte gespeichert: {pfad}", "gut")

    def _npy_speichern(self, was):
        pfad = filedialog.asksaveasfilename(
            title="Array speichern", defaultextension=".npy",
            initialdir=str(ORDNER), filetypes=[("NumPy", "*.npy")])
        if not pfad:
            return
        with self.motor.lock:
            if was == "feld":
                daten = self.motor.feld.copy()
            else:
                daten = ((self.motor.intens / self.motor.intens_n).copy()
                         if self.motor.intens_n else np.zeros_like(self.motor.feld))
        np.save(pfad, daten)
        self._log(f"Array {daten.shape} gespeichert: {pfad}", "gut")

    def _gif_speichern(self):
        filme = list(self.motor.filme)
        if len(filme) < 2:
            messagebox.showinfo("Keine Aufnahme",
                                "Erst 'Aufnahme (GIF)' einschalten und rechnen lassen.")
            return
        pfad = filedialog.asksaveasfilename(
            title="Animation speichern", defaultextension=".gif",
            initialdir=str(ORDNER), filetypes=[("GIF", "*.gif")])
        if not pfad:
            return
        try:
            from PIL import Image
        except ImportError:
            messagebox.showerror("Pillow fehlt", "Fuer GIF-Export wird Pillow gebraucht.")
            return
        spitze = max(float(np.abs(f).max()) for f in filme) or 1.0
        karte = matplotlib.colormaps[self.v_cmap.get()]
        norm = Normalize(-spitze, spitze)
        bilder = []
        for f in filme:
            rgb = (karte(norm(f.T[::-1]))[..., :3] * 255).astype(np.uint8)
            bilder.append(Image.fromarray(rgb))
        bilder[0].save(pfad, save_all=True, append_images=bilder[1:],
                       duration=40, loop=0)
        self._log(f"{len(bilder)} Bilder als GIF gespeichert: {pfad}", "gut")

    def _szene_speichern(self):
        pfad = filedialog.asksaveasfilename(
            title="Szene speichern", defaultextension=".json",
            initialdir=str(ORDNER), filetypes=[("JSON", "*.json")])
        if not pfad:
            return
        Path(pfad).write_text(json.dumps(self.szene.als_dict(), indent=2), encoding="utf-8")
        self._log(f"Szene gespeichert: {pfad}", "gut")

    def _szene_laden(self):
        pfad = filedialog.askopenfilename(title="Szene laden", initialdir=str(ORDNER),
                                          filetypes=[("JSON", "*.json")])
        if not pfad:
            return
        try:
            sz = Szene.aus_dict(json.loads(Path(pfad).read_text(encoding="utf-8")))
        except Exception as fehler:                      # noqa: BLE001
            messagebox.showerror("Laden fehlgeschlagen", str(fehler))
            return
        self._gitter_vars_setzen(sz)
        self._szene_anwenden(sz, f"Szene geladen: {pfad}")

    # =============================================================== Hilfe
    def _hilfe(self):
        text = (
            "FDTD-Studio - Kurzanleitung\n"
            "\n"
            "Tastatur\n"
            "   Leertaste   Start / Pause\n"
            "   S / Shift-S Einzelschritt / 10 Schritte\n"
            "   R           Felder zuruecksetzen\n"
            "   N           Szene neu aufbauen\n"
            "\n"
            "Maus im Feldbild\n"
            "   Werkzeug oben waehlen, dann klicken oder ziehen:\n"
            "   Sonde / Punktquelle setzen, Rechteck oder Ellipse aus\n"
            "   Dielektrikum oder PEC aufziehen, Radierer loescht ein Objekt.\n"
            "   eps_r und sigma des Pinsels stehen im Tab 'Objekte'.\n"
            "\n"
            "Ebene Wellen\n"
            "   Eine ebene Welle besteht aus je einer Punktquelle pro Zelle\n"
            "   der Startlinie - ueber die ganze Kante, bis in die PML hinein,\n"
            "   damit an den Enden nichts gebeugt wird. Bei theta != 0 werden\n"
            "   die Punkte gemaess t_d = (i cos t + j sin t) * Delta / c\n"
            "   zeitversetzt angeregt; der Arm entlang der Einfallskante\n"
            "   fuellt die Ecke. Mit 'lambda/4-Paar' wird ein zweites Blatt\n"
            "   gesetzt, dessen Rueckwaertswelle die erste ausloescht -\n"
            "   das wirkt sauber nur fuer schmalbandige Signale.\n"
            "\n"
            "Quellenarten\n"
            "   sinus    Dauerstrich, eingeschaltet ueber eine Kosinus-Rampe\n"
            "            von 'Einschaltrampe' Perioden. Rampe 0 heisst harter\n"
            "            Sprung - der regt das ganze Band an und macht das\n"
            "            Spektrum unbrauchbar.\n"
            "   ricker   zweite Ableitung der Gauss-Kurve, mittelwertfrei,\n"
            "            Spitzenfrequenz f0 - der Breitbandpuls fuer Spektren.\n"
            "   gauss    reine Gauss-Kurve; enthaelt einen Gleichanteil, der\n"
            "            als statisches Restfeld stehen bleibt.\n"
            "\n"
            "Energie\n"
            "   W = Delta^2 sum( eps0 eps_r (E^n)^2 / 2\n"
            "                    + mu0 H^(n-1/2).H^(n+1/2) / 2 )\n"
            "   E und H liegen einen halben Zeitschritt auseinander, deshalb\n"
            "   der Produktterm aus beiden H-Halbschritten. Ohne Verluste ist\n"
            "   W damit exakt erhalten. sim.total_energy() quadriert dagegen\n"
            "   H^(n+1/2) zu E^(n+1) und schwankt um einige Prozent.\n"
            "\n"
            "Stabilitaet\n"
            "   S <= 1 ist die Grenze. S > 1 laeuft absichtlich weiter und\n"
            "   divergiert sichtbar - die Anzeige meldet das und stoppt.\n"
            "\n"
            "Diagramme\n"
            "   Feld        Live-Bild, Farbskala und Groesse im Tab 'Anzeige'\n"
            "   Analyse     drei Diagramme:\n"
            "               Oszilloskop     E_z(t) an den Sonden\n"
            "               Spektrum        Balken je Frequenzband, normiert\n"
            "               Gesamtenergie   W(t) im ganzen Gebiet\n"
            "               Das Zeitfenster im Tab 'Anzeige' legt den gezeigten\n"
            "               Ausschnitt fest. Bei 'ganzer Lauf' gibt es mehr\n"
            "               Abtastwerte als Bildpunkte; dann wird die Huellkurve\n"
            "               (Minimum/Maximum je Block, Block = eine Periode)\n"
            "               gezeichnet, statt das Signal aliasend auszuduennen.\n"
            "   Intensitaet zeitgemittelte Intensitaet plus Schirmprofil\n"
            "               (Interferenzmuster hinter dem Doppelspalt)\n"
        )
        fenster = tk.Toplevel(self)
        fenster.title("Kurzanleitung")
        fenster.configure(bg=FARBE["panel"])
        fenster.geometry("620x620")
        feld = tk.Text(fenster, bg=FARBE["feld"], fg=FARBE["text"], bd=0,
                       font=("Consolas", 9), padx=12, pady=10, wrap="word")
        feld.pack(fill="both", expand=True, padx=10, pady=10)
        feld.insert("1.0", text)
        feld.configure(state="disabled")

    def _beenden(self):
        self.motor.beenden()
        self.destroy()


def selbsttest():
    """Baut das Fenster, rechnet ein Stueck und zeichnet alle drei Reiter.
    Dient dazu, die gepackte exe ohne Konsole pruefen zu koennen - das
    Ergebnis landet in selbsttest.log neben der exe."""
    zeilen = []
    try:
        app = Studio()
        app.update()
        app._schritt(60)
        app._sonde_hinzu(app.szene.nx // 2, app.szene.ny // 2)
        for reiter in range(3):
            app.nb_ansicht.select(reiter)
            app._ansicht_auffrischen()
            app.update()
        app._analyse_zeichnen()
        app._intens_zeichnen()
        zeilen.append(f"OK  -  {app.szene.nx}x{app.szene.ny} Zellen, "
                      f"{app.motor.sim.time_step} Schritte gerechnet, "
                      f"max|Ez| = {np.abs(app.motor.feld).max():.4e}")
        app._beenden()
    except Exception as fehler:                      # noqa: BLE001
        import traceback
        zeilen.append("FEHLER: " + traceback.format_exc())
    text = "\n".join(zeilen)
    print(text)
    try:
        (ORDNER / "selbsttest.log").write_text(text, encoding="utf-8")
    except OSError:
        pass
    return 0 if text.startswith("OK") else 1


def main():
    if "--selbsttest" in sys.argv:
        raise SystemExit(selbsttest())
    Studio().mainloop()


if __name__ == "__main__":
    main()
