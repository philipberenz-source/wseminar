"""
FDTD-Wellen-Baukasten
=====================
Interaktives 2D-FDTD-Spielfeld (TMz), das auf fdtd_core.FDTD2D aufbaut
(dieselbe Physik wie in der Seminararbeit). Quellen, dielektrische Bloecke,
PEC-Waende und ganze Quellen-Arrays (ebene Welle nach dem Huygensschen
Prinzip) werden per Drag-and-Drop platziert; das Feld laeuft bei "Play"
live mit -- als Matplotlib-Heatmap mit Achsen in cm, Colorbar und Zoom/Pan.
Rechts: Matplotlib-Graphen fuer Zeitsignal/Spektrum der Sonden.

Start:  python "fdtd_playground.py"   (im Ordner "code snippets")
"""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
import time
from collections import deque
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np
import pygame

import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from fdtd_core import FDTD2D, EPS0, MU0, C0

# ----------------------------------------------------------------------------
# Solver: dieselbe Physik wie fdtd_core, nur um harte Quellen erweitert
# (fdtd_core kennt nur weiche Quellen und reale Stromquellen)
# ----------------------------------------------------------------------------

class PlaygroundFDTD(FDTD2D):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._hard_sources = []

    def add_hard_source(self, i, j, waveform):
        self._hard_sources.append((i, j, waveform))

    def step(self):
        super().step()
        for i, j, wf in self._hard_sources:
            self.Ez[i, j] = wf(self.t)


def energy_density(sim: PlaygroundFDTD) -> np.ndarray:
    """Naeherung der Energiedichte auf dem Ez-Gitter (H auf benachbarte
    Ez-Punkte gemittelt). Nur fuer die Anzeige gedacht, nicht fuer total_energy()."""
    hx, hy = sim.Hx, sim.Hy
    hx_full = np.empty((sim.nx, sim.ny))
    hx_full[:, 1:-1] = 0.5 * (hx[:, :-1] + hx[:, 1:])
    hx_full[:, 0] = hx[:, 0]
    hx_full[:, -1] = hx[:, -1]
    hy_full = np.empty((sim.nx, sim.ny))
    hy_full[1:-1, :] = 0.5 * (hy[:-1, :] + hy[1:, :])
    hy_full[0, :] = hy[0, :]
    hy_full[-1, :] = hy[-1, :]
    u_e = 0.5 * EPS0 * sim.eps_r * sim.Ez ** 2
    u_h = 0.5 * MU0 * (hx_full ** 2 + hy_full ** 2)
    return u_e + u_h


# ----------------------------------------------------------------------------
# Szene: Quellen / Bloecke / Arrays / Sonden als eigene Datenhaltung. sim wird
# bei jeder strukturellen Aenderung frisch aus der Szene aufgebaut (billig
# genug bei diesen Gittergroessen) statt die (nur anfuegbaren) internen
# Listen von FDTD2D zu verwalten.
# ----------------------------------------------------------------------------

@dataclass
class SourceItem:
    kind: str        # 'soft' | 'hard' | 'current'
    i: int
    j: int
    waveform: str = "ricker"   # 'gaussian' | 'ricker' | 'cw'
    freq: float = 5.0e9        # Hz
    amp: float = 1.0


@dataclass
class BlockItem:
    kind: str         # 'dielectric' | 'pec'
    i0: int; i1: int; j0: int; j1: int
    eps_r: float = 4.0
    sigma: float = 0.0


@dataclass
class ArrayItem:
    """Reihe eng benachbarter weicher Quellen entlang einer Linie. Nach dem
    Huygensschen Prinzip ueberlagern sich die vielen kugelfoermigen
    Elementarwellen zu einer naeherungsweise ebenen Welle, die senkrecht
    zur Linie in beide Richtungen ausstrahlt (reine Broadside-Anregung, alle
    Elemente exakt gleichzeitig)."""
    i0: float; j0: float; i1: float; j1: float
    waveform: str = "cw"
    freq: float = 5.0e9
    amp: float = 1.0
    extend_to_edge: bool = True


@dataclass
class ProbeItem:
    i: int
    j: int
    t_buf: deque = field(default_factory=lambda: deque(maxlen=4096))
    e_buf: deque = field(default_factory=lambda: deque(maxlen=4096))


ARRAY_SPACING_CELLS = 1.0
MAX_ARRAY_ELEMENTS = 220


def make_waveform(src) -> callable:
    """Baut die Zeitfunktion einer Quelle. src braucht nur .amp/.freq/.waveform
    (funktioniert fuer SourceItem und ArrayItem gleichermassen)."""
    amp, freq = src.amp, max(src.freq, 1e6)
    if src.waveform == "gaussian":
        tau = 0.32 / freq
        t0 = 6.0 * tau
        return lambda t: amp * FDTD2D.gaussian_pulse(t, t0, tau)
    if src.waveform == "ricker":
        t0 = 3.0 / freq
        return lambda t: amp * FDTD2D.ricker_wavelet(t, t0, freq)
    # cw: kontinuierliche Sinuswelle mit weichem Einschalt-Ramp (raised cosine)
    period = 1.0 / freq
    ramp_periods = 4.0

    def cw(t):
        ramp = min(1.0, max(0.0, t / (ramp_periods * period)))
        env = 0.5 * (1.0 - math.cos(math.pi * ramp))
        return amp * env * math.sin(2.0 * math.pi * freq * t)

    return cw


def extend_segment_to_bounds(mid_i, mid_j, tx, ty, nx, ny, margin=1):
    """Schnittparameter t_min/t_max, an denen die Linie durch (mid_i,mid_j)
    mit Richtung (tx,ty) den Innenbereich [margin, dim-1-margin] verlaesst
    (Ray-Box-Schnitt). Wird genutzt, um ein Array bis in die PML zu ziehen."""
    lo_i, hi_i = margin, nx - 1 - margin
    lo_j, hi_j = margin, ny - 1 - margin

    def t_range(m, t, lo, hi):
        if abs(t) < 1e-9:
            return -math.inf, math.inf
        t0, t1 = (lo - m) / t, (hi - m) / t
        return min(t0, t1), max(t0, t1)

    ti0, ti1 = t_range(mid_i, tx, lo_i, hi_i)
    tj0, tj1 = t_range(mid_j, ty, lo_j, hi_j)
    return max(ti0, tj0), min(ti1, tj1)


def effective_array_endpoints(arr: ArrayItem, nx: int, ny: int):
    """Liefert die tatsaechlich zur Quellenerzeugung genutzten Endpunkte. Ein
    endlich langes Quellen-Array beugt an seinen Enden (sichtbare Kreiswellen
    von den Kanten) -- deshalb wird die Linie standardmaessig bis in die PML
    hineingezogen, dort absorbiert die PML die Enden weg und es bleibt eine
    saubere ebene Welle im sichtbaren Innenbereich uebrig."""
    if not arr.extend_to_edge:
        return arr.i0, arr.j0, arr.i1, arr.j1
    di, dj = arr.i1 - arr.i0, arr.j1 - arr.j0
    length = math.hypot(di, dj)
    if length < 1e-9:
        return arr.i0, arr.j0, arr.i1, arr.j1
    tx, ty = di / length, dj / length
    mid_i, mid_j = (arr.i0 + arr.i1) / 2.0, (arr.j0 + arr.j1) / 2.0
    t_min, t_max = extend_segment_to_bounds(mid_i, mid_j, tx, ty, nx, ny, margin=1)
    if not (math.isfinite(t_min) and math.isfinite(t_max)) or t_max <= t_min:
        return arr.i0, arr.j0, arr.i1, arr.j1
    return (mid_i + t_min * tx, mid_j + t_min * ty, mid_i + t_max * tx, mid_j + t_max * ty)


def array_elements(arr: ArrayItem, nx: int, ny: int):
    """Liefert die Gitterzellen [(i, j), ...] der Elemente eines Quellen-
    Arrays (alle exakt gleichzeitig angeregt, reine Broadside-Abstrahlung)."""
    ai0, aj0, ai1, aj1 = effective_array_endpoints(arr, nx, ny)
    di, dj = ai1 - ai0, aj1 - aj0
    length_cells = math.hypot(di, dj)
    if length_cells < 1e-9:
        return [(int(round(ai0)), int(round(aj0)))]
    n = int(min(max(round(length_cells / ARRAY_SPACING_CELLS) + 1, 2), MAX_ARRAY_ELEMENTS))
    elems, seen = [], set()
    for k in range(n):
        t = k / (n - 1)
        i, j = ai0 + di * t, aj0 + dj * t
        ii, jj = int(round(i)), int(round(j))
        if (ii, jj) in seen:
            continue
        seen.add((ii, jj))
        elems.append((ii, jj))
    return elems


class Scene:
    def __init__(self, nx, ny, d, npml=10, courant=0.99):
        self.nx, self.ny, self.d, self.npml, self.courant = nx, ny, d, npml, courant
        self.blocks: list[BlockItem] = []
        self.sources: list[SourceItem] = []
        self.arrays: list[ArrayItem] = []
        self.probes: list[ProbeItem] = []
        self.sim: PlaygroundFDTD = None
        self.rebuild()

    def rebuild(self):
        sim = PlaygroundFDTD(self.nx, self.ny, self.d, courant=self.courant, npml=self.npml)
        for b in self.blocks:
            if b.kind == "dielectric":
                sim.set_material_box(b.i0, b.i1, b.j0, b.j1, eps_r=b.eps_r, sigma=b.sigma)
            else:
                sim.set_pec_box(b.i0, b.i1, b.j0, b.j1)
        for s in self.sources:
            wf = make_waveform(s)
            if s.kind == "soft":
                sim.add_soft_source(s.i, s.j, wf)
            elif s.kind == "hard":
                sim.add_hard_source(s.i, s.j, wf)
            else:
                sim.add_current_source(s.i, s.j, wf, amplitude=1.0)
        for arr in self.arrays:
            base_wf = make_waveform(arr)
            for i, j in array_elements(arr, self.nx, self.ny):
                sim.add_soft_source(i, j, base_wf)
        self.sim = sim
        for p in self.probes:
            p.t_buf.clear()
            p.e_buf.clear()

    def clear_all(self):
        self.blocks.clear()
        self.sources.clear()
        self.arrays.clear()
        self.probes.clear()
        self.rebuild()

    def sample_probes(self):
        sim = self.sim
        for p in self.probes:
            p.t_buf.append(sim.t)
            p.e_buf.append(sim.Ez[p.i, p.j])

    def to_dict(self):
        return {
            "blocks": [asdict(b) for b in self.blocks],
            "sources": [asdict(s) for s in self.sources],
            "arrays": [asdict(a) for a in self.arrays],
            "probes": [{"i": p.i, "j": p.j} for p in self.probes],
        }

    def load_dict(self, data):
        self.blocks = [BlockItem(**b) for b in data.get("blocks", [])]
        self.sources = [SourceItem(**s) for s in data.get("sources", [])]
        self.arrays = [ArrayItem(**a) for a in data.get("arrays", [])]
        self.probes = [ProbeItem(i=p["i"], j=p["j"]) for p in data.get("probes", [])]
        self.rebuild()


# ----------------------------------------------------------------------------
# Theme
# ----------------------------------------------------------------------------

BG = (18, 20, 26)
PANEL_BG = (26, 29, 37)
PANEL_BG_LIGHT = (33, 37, 47)
BORDER = (45, 49, 60)
TEXT = (225, 228, 235)
TEXT_DIM = (150, 155, 168)
ACCENT = (86, 160, 255)
ACCENT_LIGHT = (140, 195, 255)
ACCENT_SOFT = (80, 220, 200)
ACCENT_HARD = (240, 90, 160)
ACCENT_CUR = (250, 200, 70)
ACCENT_PROBE = (255, 255, 255)
ACCENT_DIEL = (90, 210, 140)
ACCENT_PEC = (210, 90, 90)
ACCENT_ERASE = (200, 90, 90)
ACCENT_PLANE = (175, 130, 255)
DANGER = (220, 80, 80)
TRACK_BG = (52, 57, 70)

PROBE_COLORS = [
    (255, 255, 255), (255, 200, 80), (120, 220, 255), (255, 140, 200),
    (140, 255, 170), (255, 120, 120), (180, 160, 255), (255, 220, 120),
]


def probe_color(idx):
    return PROBE_COLORS[idx % len(PROBE_COLORS)]


def rgb01(c):
    return (c[0] / 255.0, c[1] / 255.0, c[2] / 255.0)


def lerp(a, b, t):
    return a + (b - a) * t


DEFAULT_WIN_W, DEFAULT_WIN_H = 1560, 1140
MIN_WIN_W, MIN_WIN_H = 1220, 820
HEADER_H = 148
FOOTER_H = 96
LEFT_W = 250
RIGHT_W = 380

GRID_PRESETS = [
    ("Klein", 120, 90),
    ("Mittel", 170, 110),
    ("Gross", 220, 160),
    ("XL", 240, 240),
]

SCENARIOS = [
    ("leer", "Leer"),
    ("doppelspalt", "Doppelspalt"),
    ("resonator", "Resonator"),
]
CELL_D = 2.0e-3
NPML = 10
MIN_ZOOM_SPAN = 8.0     # kleinster sichtbarer Ausschnitt beim Reinzoomen (Zellen)
FIELD_FLOOR = 4e-3      # Kontrast-Skala faellt nie unter dieses Niveau. Verhindert zwei
                        # Effekte: (1) unipolare Quellen (Gauss) hinterlassen ein
                        # Restfeld, das sonst die ganze Anzeige einfaerben wuerde,
                        # (2) minimale numerische PML-Restreflexion (typ. -100 dB und
                        # tiefer, physikalisch bedeutungslos) wuerde sonst, sobald das
                        # eigentliche Signal abgeklungen ist, zum neuen "Maximum"
                        # werden und durch die schnelle Skalierung wie eine echte,
                        # aus der PML zurücklaufende Welle aussehen
FIELD_SCALE_DECAY = 0.75    # zieht die Farbskala schnell auf den aktuellen
                            # Maximalwert zusammen (statt laenger auf einem vergangenen
                            # Peak zu bleiben), damit auch schwaechere lokale Maxima
                            # (z. B. hinter einem Doppelspalt) kontrastreich bleiben
FIELD_REFRESH_DT = 0.035
PLOT_REFRESH_DT = 0.09


def rrect(surf, color, rect, radius=10, width=0):
    pygame.draw.rect(surf, color, rect, width=width, border_radius=radius)


def point_segment_dist(px, py, ax_, ay_, bx_, by_):
    dx, dy = bx_ - ax_, by_ - ay_
    l2 = dx * dx + dy * dy
    if l2 < 1e-12:
        return math.hypot(px - ax_, py - ay_)
    t = max(0.0, min(1.0, ((px - ax_) * dx + (py - ay_) * dy) / l2))
    cx, cy = ax_ + t * dx, ay_ + t * dy
    return math.hypot(px - cx, py - cy)


class Slider:
    def __init__(self, x, y, w, lo, hi, value, label, fmt="{:.2f}", unit="", disp_scale=1.0):
        self.rect = pygame.Rect(x, y, w, 40)
        self.track = pygame.Rect(x, y + 24, w, 6)
        self.lo, self.hi = lo, hi
        self.value = value
        self.label = label
        self.fmt = fmt
        self.unit = unit
        self.disp_scale = disp_scale   # value ist intern in Basiseinheiten (z. B. Hz), Anzeige skaliert (z. B. GHz)
        self.dragging = False

    def handle_x(self):
        t = (self.value - self.lo) / (self.hi - self.lo)
        return self.track.x + t * self.track.w

    def handle_event(self, ev) -> bool:
        if ev.type == pygame.MOUSEBUTTONDOWN and ev.button == 1:
            if self.track.inflate(0, 16).collidepoint(ev.pos):
                self.dragging = True
                self._set_from_x(ev.pos[0])
                return True
        elif ev.type == pygame.MOUSEBUTTONUP:
            if self.dragging:
                self.dragging = False
                return True
        elif ev.type == pygame.MOUSEMOTION and self.dragging:
            self._set_from_x(ev.pos[0])
            return True
        return False

    def _set_from_x(self, x):
        t = min(1.0, max(0.0, (x - self.track.x) / self.track.w))
        self.value = self.lo + t * (self.hi - self.lo)

    def draw(self, surf, font):
        txt = font.render(f"{self.label}: {self.fmt.format(self.value / self.disp_scale)}{self.unit}", True, TEXT_DIM)
        surf.blit(txt, (self.rect.x, self.rect.y))
        rrect(surf, TRACK_BG, self.track, radius=3)
        t = (self.value - self.lo) / (self.hi - self.lo)
        fill = pygame.Rect(self.track.x, self.track.y, max(6, int(self.track.w * t)), self.track.h)
        rrect(surf, ACCENT, fill, radius=3)
        hx = self.handle_x()
        pygame.draw.circle(surf, ACCENT_LIGHT, (int(hx), self.track.centery), 8)
        pygame.draw.circle(surf, PANEL_BG, (int(hx), self.track.centery), 8, 2)


class TextInput:
    """Kleines numerisches Eingabefeld (Klick = Fokus, Ziffern tippen, Enter/
    Klick daneben = bestaetigen). Solange fokussiert, muessen Tastendruecke
    hier ankommen statt als Werkzeug-Hotkeys interpretiert zu werden."""
    def __init__(self, x, y, w, h, text="", max_len=4):
        self.rect = pygame.Rect(x, y, w, h)
        self.text = text
        self.max_len = max_len
        self.focused = False

    def handle_event(self, ev) -> bool:
        if ev.type == pygame.MOUSEBUTTONDOWN and ev.button == 1:
            self.focused = self.rect.collidepoint(ev.pos)
            return self.focused
        if self.focused and ev.type == pygame.KEYDOWN:
            if ev.key == pygame.K_BACKSPACE:
                self.text = self.text[:-1]
            elif ev.key in (pygame.K_RETURN, pygame.K_KP_ENTER, pygame.K_TAB, pygame.K_ESCAPE):
                self.focused = False
            elif ev.unicode.isdigit() and len(self.text) < self.max_len:
                self.text += ev.unicode
            return True
        return False

    def draw(self, surf, font):
        rrect(surf, PANEL_BG_LIGHT, self.rect, radius=6)
        rrect(surf, ACCENT if self.focused else BORDER, self.rect, radius=6, width=2 if self.focused else 1)
        t = font.render(self.text, True, TEXT)
        surf.blit(t, (self.rect.x + 8, self.rect.y + (self.rect.h - t.get_height()) // 2))
        if self.focused and int(time.time() * 2) % 2 == 0:
            cx = self.rect.x + 8 + t.get_width() + 2
            pygame.draw.line(surf, TEXT, (cx, self.rect.y + 4), (cx, self.rect.bottom - 4), 1)


TOOLS = [
    ("select", "Auswahl"),
    ("soft", "Weiche Quelle"),
    ("hard", "Harte Quelle"),
    ("current", "Stromquelle J"),
    ("probe", "Sonde"),
    ("dielectric", "Dielektrikum"),
    ("pec", "PEC-Wand"),
    ("planewave", "Ebene Welle"),
    ("erase", "Radierer"),
]
TOOL_COLORS = {
    "soft": ACCENT_SOFT, "hard": ACCENT_HARD, "current": ACCENT_CUR,
    "probe": ACCENT_PROBE, "dielectric": ACCENT_DIEL, "pec": ACCENT_PEC,
    "erase": ACCENT_ERASE, "select": ACCENT, "planewave": ACCENT_PLANE,
}
POINT_TOOLS = {"soft", "hard", "current", "probe"}
RECT_TOOLS = {"dielectric", "pec"}
LINE_TOOLS = {"planewave"}


def draw_tool_icon(surf, kind, center, r, color):
    cx, cy = center
    if kind == "select":
        pts = [(cx - r, cy - r), (cx - r, cy + r * 0.9), (cx - r * 0.25, cy + r * 0.35),
               (cx + r * 0.15, cy + r * 0.95), (cx + r * 0.45, cy + r * 0.75),
               (cx + r * 0.05, cy + r * 0.05), (cx + r * 0.75, cy - r * 0.05)]
        pygame.draw.polygon(surf, color, pts)
    elif kind == "soft":
        pygame.draw.circle(surf, color, (cx, cy), r, 2)
        pygame.draw.circle(surf, color, (cx, cy), int(r * 0.55), 2)
    elif kind == "hard":
        pygame.draw.circle(surf, color, (cx, cy), r)
    elif kind == "current":
        pts = [(cx - r * 0.2, cy - r), (cx - r * 0.9, cy + r * 0.15), (cx - r * 0.1, cy + r * 0.15),
               (cx - r * 0.5, cy + r), (cx + r * 0.9, cy - r * 0.25), (cx + r * 0.05, cy - r * 0.25)]
        pygame.draw.polygon(surf, color, pts)
    elif kind == "probe":
        pygame.draw.circle(surf, color, (cx, cy - r * 0.3), r * 0.5)
        pygame.draw.polygon(surf, color, [(cx, cy - r * 0.05), (cx - r * 0.35, cy + r), (cx + r * 0.35, cy + r)])
    elif kind == "dielectric":
        rect = pygame.Rect(0, 0, int(r * 1.8), int(r * 1.8))
        rect.center = (cx, cy)
        rrect(surf, color, rect, radius=4)
        for k in range(1, 4):
            y = rect.y + k * rect.h / 4
            pygame.draw.line(surf, PANEL_BG, (rect.x, y), (rect.right, y), 1)
        rrect(surf, color, rect, radius=4, width=2)
    elif kind == "pec":
        rect = pygame.Rect(0, 0, int(r * 1.8), int(r * 1.8))
        rect.center = (cx, cy)
        rrect(surf, color, rect, radius=3)
        for k in range(4):
            x0 = rect.x + k * rect.w / 4
            pygame.draw.line(surf, PANEL_BG, (x0, rect.bottom), (x0 + rect.w / 4, rect.top), 1)
    elif kind == "erase":
        rect = pygame.Rect(0, 0, int(r * 1.7), int(r * 1.2))
        rect.center = (cx, cy)
        rrect(surf, color, rect, radius=4)
        pygame.draw.line(surf, PANEL_BG, rect.topleft, rect.bottomright, 2)
    elif kind == "planewave":
        for dx in (-r * 0.9, 0, r * 0.9):
            pygame.draw.circle(surf, color, (int(cx + dx), int(cy + r * 0.5)), max(2, int(r * 0.22)))
        rc = pygame.Rect(0, 0, int(r * 3.2), int(r * 3.2))
        rc.center = (cx, cy + r * 2.0)
        pygame.draw.arc(surf, color, rc, math.radians(200), math.radians(340), 2)


WAVEFORMS = ["ricker", "gaussian", "cw"]
WAVE_LABELS = {"ricker": "Ricker", "gaussian": "Gauss", "cw": "CW"}


# ----------------------------------------------------------------------------
# Einfache Schwarz-Blau-Farbskala fuer die Feldanzeige (nur Betrag, kein
# Vorzeichen). Eine Log/dB-Normalisierung sorgt dafuer, dass auch lokale
# Maxima weit unterhalb des globalen Maximums (z. B. hinter einem Doppel-
# spalt) noch als eigene Helligkeitsstufe sichtbar bleiben, statt zusammen
# mit dem Hintergrund zu verschwinden -- "decades" (per Regler einstellbar)
# legt fest, wie viele Zehnerpotenzen Dynamikbereich dargestellt werden,
# also wie stark niedrige Werte ausgeblendet werden.
# ----------------------------------------------------------------------------

BLUE_CMAP = LinearSegmentedColormap.from_list("fdtd_blue", [
    rgb01((7, 8, 11)), rgb01((70, 170, 255)),
])


class LogNormPositive(Normalize):
    """Log/dB-Normalisierung auf [0, 1]."""
    def __init__(self, decades=2.75):
        super().__init__(vmin=0.0, vmax=1.0, clip=True)
        self.decades = decades

    def __call__(self, value, clip=None):
        v = np.clip(np.asarray(value, dtype=float), 0.0, 1.0)
        floor = 10.0 ** (-self.decades - 1)
        log_v = np.log10(np.where(v > floor, v, floor))
        return np.ma.masked_array(np.clip(1.0 + log_v / self.decades, 0.0, 1.0))


# ----------------------------------------------------------------------------
# Matplotlib-Einbettung: Figures werden einmal erzeugt und per ax.cla() +
# neu plotten aktualisiert (billiger als Figures neu zu bauen), dann als
# RGBA-Puffer in eine pygame-Surface kopiert.
# ----------------------------------------------------------------------------

def style_axes(ax):
    ax.set_facecolor(rgb01(PANEL_BG_LIGHT))
    for spine in ax.spines.values():
        spine.set_color(rgb01(BORDER))
    ax.tick_params(colors=rgb01(TEXT_DIM), labelsize=8)
    ax.xaxis.label.set_color(rgb01(TEXT_DIM))
    ax.yaxis.label.set_color(rgb01(TEXT_DIM))
    ax.grid(True, color=rgb01(BORDER), linewidth=0.6, alpha=0.8)


def make_figure(w_px, h_px, dpi=100):
    fig = Figure(figsize=(w_px / dpi, h_px / dpi), dpi=dpi)
    fig.patch.set_facecolor(rgb01(PANEL_BG_LIGHT))
    ax = fig.add_subplot(111)
    style_axes(ax)
    canvas = FigureCanvasAgg(fig)
    return fig, ax, canvas


def figure_to_surface(canvas):
    canvas.draw()
    w, h = canvas.get_width_height()
    buf = bytes(canvas.buffer_rgba())
    return pygame.image.frombuffer(buf, (w, h), "RGBA")


TIME_PLOT_SAMPLES = 450     # Zeitsignal zeigt immer nur die letzten N Rohschritte, nicht
                            # den ganzen (bis zu 4096 Schritte tiefen) Sonden-Puffer. Sonst
                            # waechst die Zahl der auf die x-Achse gequetschten Perioden mit
                            # der Laufzeit immer weiter -- bei laenger schwingenden Quellen
                            # wird die Linie dann dichter und "verzerrt" aussehend. Mit einer
                            # festen Fensterbreite bleibt die Anzeige immer gleich klar, egal
                            # wie lange die Simulation schon laeuft (Spektrum/CSV-Export
                            # nutzen weiterhin den vollen Puffer fuer gute Aufloesung).


class App:
    def __init__(self):
        pygame.init()
        pygame.display.set_caption("FDTD Wellen-Baukasten")
        self.win_w, self.win_h = DEFAULT_WIN_W, DEFAULT_WIN_H
        self.screen = pygame.display.set_mode((self.win_w, self.win_h), pygame.RESIZABLE)
        self._pending_resize = None
        self.clock = pygame.time.Clock()
        self.font = pygame.font.SysFont("Segoe UI", 15)
        self.font_small = pygame.font.SysFont("Segoe UI", 12)
        self.font_bold = pygame.font.SysFont("Segoe UI", 15, bold=True)
        self.font_title = pygame.font.SysFont("Segoe UI", 21, bold=True)

        self.grid_preset = 1  # "Mittel"
        self.scenario = "leer"
        self.nx, self.ny = GRID_PRESETS[self.grid_preset][1], GRID_PRESETS[self.grid_preset][2]
        self.npml = NPML
        self.courant = 0.99
        self.scene = Scene(self.nx, self.ny, CELL_D, self.npml, self.courant)
        self.canvas_rect = pygame.Rect(LEFT_W, HEADER_H,
                                        self.win_w - LEFT_W - RIGHT_W, self.win_h - HEADER_H - FOOTER_H)
        self.vmax_running = 1e-6

        self.tool = "soft"
        self.playing = False
        self.speed = 4
        self.view_mode = "ez"
        # Feste, nicht ueber die UI einstellbare Werte fuer die Farbskala (siehe
        # BLUE_CMAP/LogNormPositive): color_decades steuert, wie stark niedrige
        # Werte ausgeblendet werden (kleiner = schaerfer), gain ist ein reiner
        # Skalierungsfaktor auf die Autoskala. Beide wurden experimentell auf
        # gut lesbare Werte eingestellt.
        self.color_decades = 2.0
        self.gain = 1.0
        self.show_grid = False
        self.layer_field = True
        self.layer_material = True
        self.layer_sources = True
        self.layer_probes = True
        self.layer_arrays = True

        self.new_freq = 5.0e9
        self.new_amp = 1.0
        self.new_waveform = "ricker"
        self.new_eps = 4.0
        self.new_sigma = 0.0

        self.dragging_new = None
        self.rect_start = None
        self.rect_cur = None
        self.line_start = None
        self.line_cur = None
        self.selected = None       # (kind, index), kind in {'source','probe','block','array'}
        self.moving = False

        self._right_down = False
        self._pan_dragged = False
        self._pan_start_pos = None
        self._pan_start_view = None

        self.clear_confirm_until = 0.0
        self.toasts = []
        self.t_last = time.time()
        self.fps_smooth = 60.0
        self._wave_button_rects = {}

        self._build_sliders()
        self._build_figures()
        self._last_plot_refresh = 0.0
        self._force_plot_refresh = True
        self._last_field_refresh = 0.0
        self._force_field_refresh = True
        self._time_surf = self._spec_surf = self._field_surf = None

    # -------------------------------------------------------------- setup
    def _rebuild_domain(self, nx, ny, npml, courant):
        self.nx, self.ny, self.npml, self.courant = nx, ny, npml, courant
        self.scene = Scene(nx, ny, CELL_D, npml, courant)
        self.selected = None
        self.vmax_running = 1e-6
        self._build_field_figure(reset_view=True)
        self._force_field_refresh = True
        self._force_plot_refresh = True

    def set_grid_size(self, preset_idx):
        self.grid_preset = preset_idx
        _, nx, ny = GRID_PRESETS[preset_idx]
        self._rebuild_domain(nx, ny, self.npml, self.courant)
        self.nx_input.text, self.ny_input.text = str(nx), str(ny)
        self.toast(f"Neues Gitter {nx}×{ny} – Szene geleert")

    def apply_domain_settings(self):
        try:
            nx = int(self.nx_input.text)
            ny = int(self.ny_input.text)
        except ValueError:
            self.toast("Ungültige Gittergröße")
            return
        nx = max(30, min(800, nx))
        ny = max(30, min(800, ny))
        npml = int(round(self.s_pml.value))
        npml = max(2, min(npml, min(nx, ny) // 2 - 5))
        courant = round(self.s_courant.value, 4)
        self.grid_preset = -1
        self._rebuild_domain(nx, ny, npml, courant)
        self.nx_input.text, self.ny_input.text = str(nx), str(ny)
        self.toast(f"Neue Domäne {nx}×{ny}, Courant={courant:.3f}, PML={npml} Zellen")

    def load_scenario(self, key):
        """Baut ein vorgefertigtes Setup, skaliert auf die aktuelle Gittergroesse."""
        self.scenario = key
        self.scene.clear_all()
        nx, ny, d = self.nx, self.ny, self.scene.d
        margin = self.npml + 6

        if key == "doppelspalt":
            cx = nx // 2
            wall_j = int(ny * 0.42)
            thick = max(2, int(ny * 0.02))
            gap_w = max(3, int(nx * 0.018))
            gap_off = max(gap_w, int(nx * 0.09))
            s1c, s2c = cx - gap_off, cx + gap_off
            # Die aeusseren Wandsegmente reichen bis an den Gitterrand (in die
            # PML hinein) statt nur bis "margin" -- sonst kann die Welle seitlich
            # an der Wand vorbei laufen statt nur durch die beiden Spalte.
            for i0, i1 in ((1, s1c - gap_w // 2),
                           (s1c + gap_w // 2, s2c - gap_w // 2),
                           (s2c + gap_w // 2, nx - 1)):
                if i1 > i0:
                    self.scene.blocks.append(BlockItem("pec", i0, i1, wall_j, wall_j + thick))
            src_j = int(ny * 0.18)
            self.scene.arrays.append(ArrayItem(margin, src_j, nx - margin, src_j,
                                                 waveform="cw", freq=8.0e9, amp=1.0))
            self.scene.probes.append(ProbeItem(cx, int(ny * 0.8)))
            self.toast("Szenario: Doppelspalt geladen")
        elif key == "resonator":
            cx, cy = nx // 2, ny // 2
            side = max(20, int(min(nx, ny) * 0.5))
            thick = max(2, int(min(nx, ny) * 0.02))
            half = side // 2
            i0, i1 = cx - half, cx + half
            j0, j1 = cy - half, cy + half
            # Vier PEC-Waende schliessen einen quadratischen Hohlraum-Resonator
            # (statt nur zweier paralleler Spiegel) -- echte 2D-Kavitaetsmoden.
            self.scene.blocks.append(BlockItem("pec", i0 - thick, i0, j0 - thick, j1 + thick))  # links
            self.scene.blocks.append(BlockItem("pec", i1, i1 + thick, j0 - thick, j1 + thick))  # rechts
            self.scene.blocks.append(BlockItem("pec", i0, i1, j0 - thick, j0))                  # oben
            self.scene.blocks.append(BlockItem("pec", i0, i1, j1, j1 + thick))                  # unten
            L = side * d
            m, n = 3, 2
            freq = (C0 / 2.0) * math.sqrt((m / L) ** 2 + (n / L) ** 2)
            # Quelle auf den ersten Schwingungsbauch des (m,n)-Modus setzen,
            # sonst wuerde eine Quelle auf einer Knotenlinie den Modus kaum anregen.
            src_i = i0 + max(1, int(round(side / (2 * m))))
            src_j = j0 + max(1, int(round(side / (2 * n))))
            self.scene.sources.append(SourceItem("soft", src_i, src_j, "cw", freq, 1.0))
            self.scene.probes.append(ProbeItem(min(i1 - 1, src_i + 3), min(j1 - 1, src_j + 2)))
            self.toast(f"Szenario: Resonator geladen ({freq / 1e9:.2f} GHz, Modus m={m}, n={n})")
        else:
            self.toast("Szene geleert")

        self.scene.rebuild()
        self.selected = None
        self.vmax_running = 1e-6
        self._force_field_refresh = True
        self._force_plot_refresh = True

    def _apply_resize(self, new_w, new_h):
        new_w = max(int(new_w), MIN_WIN_W)
        new_h = max(int(new_h), MIN_WIN_H)
        if new_w == self.win_w and new_h == self.win_h:
            return
        self.win_w, self.win_h = new_w, new_h
        self.screen = pygame.display.set_mode((self.win_w, self.win_h), pygame.RESIZABLE)
        self.canvas_rect = pygame.Rect(LEFT_W, HEADER_H,
                                        self.win_w - LEFT_W - RIGHT_W, self.win_h - HEADER_H - FOOTER_H)
        self._build_field_figure(reset_view=False)
        self._force_field_refresh = True

    def _build_sliders(self):
        x = 20
        w = LEFT_W - 40
        self.s_freq = Slider(x, 0, w, 0.5e9, 20e9, self.new_freq, "Frequenz", "{:.2f}", " GHz", disp_scale=1e9)
        self.s_amp = Slider(x, 0, w, 0.1, 3.0, self.new_amp, "Amplitude", "{:.2f}", " a.u.")
        self.s_eps = Slider(x, 0, w, 1.0, 12.0, self.new_eps, "eps_r", "{:.1f}")
        self.s_sigma = Slider(x, 0, w, 0.0, 2.0, self.new_sigma, "sigma", "{:.2f}", " S/m")
        self.s_speed = Slider(0, 0, 170, 1, 32, self.speed, "Tempo", "{:.0f}", "x")
        row4_y = 96
        self.s_courant = Slider(270, row4_y, 130, 0.5, 0.999, self.courant, "Courant", "{:.3f}")
        self.s_pml = Slider(430, row4_y, 110, 2, 40, self.npml, "PML", "{:.0f}", " Zellen")
        self.nx_input = TextInput(70, row4_y + 2, 55, 24, str(self.nx))
        self.ny_input = TextInput(142, row4_y + 2, 55, 24, str(self.ny))

    def _build_figures(self):
        pad = 18
        plot_w = RIGHT_W - 2 * pad
        self.fig_time, self.ax_time, self.canvas_time = make_figure(plot_w, 236)
        self.fig_spec, self.ax_spec, self.canvas_spec = make_figure(plot_w, 236)
        self._norm_field = LogNormPositive(decades=self.color_decades)
        self._build_field_figure(reset_view=True)

    def _build_field_figure(self, reset_view):
        """(Neu-)Erstellt die Feld-Figure fuer die aktuelle Canvas-Groesse.
        Wird beim Start, bei Gittergroessen-Wechsel und beim Fenster-Resize
        aufgerufen; reset_view=False erhaelt Zoom/Pan (z. B. beim Resize)."""
        fw, fh = self.canvas_rect.w, self.canvas_rect.h
        self.fig_field, self.ax_field, self.canvas_field = make_figure(fw, fh)
        self.ax_field.set_aspect("equal")
        self.ax_field.grid(False)
        dummy = np.zeros((self.ny, self.nx))
        self.im_field = self.ax_field.imshow(dummy, cmap=BLUE_CMAP, norm=self._norm_field,
                                              extent=(0, self.nx, self.ny, 0), origin="upper",
                                              interpolation="nearest", zorder=1)
        d_cm = self.scene.d * 100
        fmt = FuncFormatter(lambda val, pos: f"{val * d_cm:.1f}")
        self.ax_field.xaxis.set_major_formatter(fmt)
        self.ax_field.yaxis.set_major_formatter(fmt)
        self.ax_field.set_xlabel("x (cm)", fontsize=8)
        self.ax_field.set_ylabel("y (cm)", fontsize=8)
        self.cbar_field = self.fig_field.colorbar(self.im_field, ax=self.ax_field, fraction=0.046, pad=0.03)
        self.cbar_field.ax.tick_params(colors=rgb01(TEXT_DIM), labelsize=8)
        self.cbar_field.outline.set_edgecolor(rgb01(BORDER))

        if reset_view or not hasattr(self, "view_i0"):
            self.view_i0, self.view_i1 = 0.0, float(self.nx)
            self.view_j0, self.view_j1 = 0.0, float(self.ny)
        self._sync_field_axes_limits()
        self.fig_field.tight_layout(pad=0.6)
        self._field_px_size = self.canvas_field.get_width_height()
        self._field_artists = []

    def toast(self, msg):
        self.toasts.append([msg, 2.4])

    # -------------------------------------------------------------- zoom/pan
    def _sync_field_axes_limits(self):
        self.ax_field.set_xlim(self.view_i0, self.view_i1)
        self.ax_field.set_ylim(self.view_j1, self.view_j0)

    def _clamp_span(self, a, b, full):
        w = b - a
        if a < 0:
            a, b = 0.0, w
        if b > full:
            a, b = full - w, full
        return a, b

    def zoom_at(self, data_i, data_j, factor):
        cur_w = self.view_i1 - self.view_i0
        cur_h = self.view_j1 - self.view_j0
        full_w, full_h = float(self.nx), float(self.ny)
        new_w = cur_w * factor
        if new_w > full_w:
            factor = full_w / cur_w
        elif new_w < MIN_ZOOM_SPAN:
            factor = MIN_ZOOM_SPAN / cur_w
        new_w = cur_w * factor
        new_h = cur_h * factor
        ni0 = data_i - (data_i - self.view_i0) * factor
        nj0 = data_j - (data_j - self.view_j0) * factor
        ni1, nj1 = ni0 + new_w, nj0 + new_h
        if ni0 < 0:
            ni1 -= ni0; ni0 = 0.0
        if ni1 > full_w:
            ni0 -= (ni1 - full_w); ni1 = full_w
        if nj0 < 0:
            nj1 -= nj0; nj0 = 0.0
        if nj1 > full_h:
            nj0 -= (nj1 - full_h); nj1 = full_h
        self.view_i0, self.view_i1 = ni0, ni1
        self.view_j0, self.view_j1 = nj0, nj1
        self._sync_field_axes_limits()

    def reset_zoom(self):
        self.view_i0, self.view_i1 = 0.0, float(self.nx)
        self.view_j0, self.view_j1 = 0.0, float(self.ny)
        self._sync_field_axes_limits()

    def pixels_per_cell(self):
        span = self.view_i1 - self.view_i0
        return self.canvas_rect.w / span if span > 1e-9 else 1.0

    # -------------------------------------------------------------- coords
    def screen_to_grid(self, pos):
        fig_w, fig_h = self._field_px_size
        fx = pos[0] - self.canvas_rect.x
        fy_mpl = fig_h - (pos[1] - self.canvas_rect.y)
        data_x, data_y = self.ax_field.transData.inverted().transform((fx, fy_mpl))
        return float(data_x), float(data_y)

    def grid_to_screen(self, i, j):
        fig_w, fig_h = self._field_px_size
        fx, fy_mpl = self.ax_field.transData.transform((i, j))
        fy_top = fig_h - fy_mpl
        return self.canvas_rect.x + fx, self.canvas_rect.y + fy_top

    def clamp_cell(self, i, j):
        return int(min(max(round(i), 1), self.nx - 2)), int(min(max(round(j), 1), self.ny - 2))

    # -------------------------------------------------------------- hit test
    def find_item_at(self, pos):
        if not self.canvas_rect.collidepoint(pos):
            return None
        gi, gj = self.screen_to_grid(pos)
        best, best_d = None, 10.0 / self.pixels_per_cell()
        for k, s in enumerate(self.scene.sources):
            d = math.hypot(gi - s.i, gj - s.j)
            if d < best_d:
                best, best_d = ("source", k), d
        for k, p in enumerate(self.scene.probes):
            d = math.hypot(gi - p.i, gj - p.j)
            if d < best_d:
                best, best_d = ("probe", k), d
        if best is not None:
            return best
        for k, a in enumerate(self.scene.arrays):
            ai0, aj0, ai1, aj1 = effective_array_endpoints(a, self.nx, self.ny)
            d = point_segment_dist(gi, gj, ai0, aj0, ai1, aj1)
            if d < best_d:
                return ("array", k)
        for k, b in enumerate(self.scene.blocks):
            if b.i0 <= gi <= b.i1 and b.j0 <= gj <= b.j1:
                return ("block", k)
        return None

    def delete_item(self, item):
        kind, idx = item
        if kind == "source":
            del self.scene.sources[idx]
        elif kind == "probe":
            del self.scene.probes[idx]
        elif kind == "array":
            del self.scene.arrays[idx]
        else:
            del self.scene.blocks[idx]
        self.scene.rebuild()
        self.selected = None

    # -------------------------------------------------------------- events
    def run(self):
        while True:
            now = time.time()
            dt_wall = now - self.t_last
            self.t_last = now
            for ev in pygame.event.get():
                if ev.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit(0)
                if ev.type == pygame.VIDEORESIZE:
                    self._pending_resize = (ev.w, ev.h, time.time())
                    continue
                self.handle_event(ev)

            if self._pending_resize is not None:
                w, h, ts = self._pending_resize
                if time.time() - ts > 0.15:
                    self._apply_resize(w, h)
                    self._pending_resize = None

            if self.playing:
                for _ in range(int(round(self.speed))):
                    self.scene.sim.step()
                    self.scene.sample_probes()

            for t in self.toasts:
                t[1] -= dt_wall
            self.toasts = [t for t in self.toasts if t[1] > 0]

            self.fps_smooth = lerp(self.fps_smooth, self.clock.get_fps() or 60.0, 0.1)
            self.draw()
            self.clock.tick(60)

    def active_sliders(self):
        """Sliders je nach Tool/Auswahl im linken Panel."""
        if self.selected and self.selected[0] == "source":
            s = self.scene.sources[self.selected[1]]
            self.s_freq.value, self.s_amp.value = s.freq, s.amp
            return [self.s_freq, self.s_amp]
        if self.selected and self.selected[0] == "array":
            a = self.scene.arrays[self.selected[1]]
            self.s_freq.value, self.s_amp.value = a.freq, a.amp
            return [self.s_freq, self.s_amp]
        if self.selected and self.selected[0] == "block" and self.scene.blocks[self.selected[1]].kind == "dielectric":
            b = self.scene.blocks[self.selected[1]]
            self.s_eps.value, self.s_sigma.value = b.eps_r, b.sigma
            return [self.s_eps, self.s_sigma]
        if self.selected:
            return []
        if self.tool in {"soft", "hard", "current", "planewave"}:
            return [self.s_freq, self.s_amp]
        if self.tool == "dielectric":
            return [self.s_eps, self.s_sigma]
        return []

    def handle_event(self, ev):
        if ev.type == pygame.KEYDOWN:
            if self.nx_input.focused:
                self.nx_input.handle_event(ev)
                return
            if self.ny_input.focused:
                self.ny_input.handle_event(ev)
                return
            self.handle_key(ev)
            return

        if ev.type == pygame.MOUSEBUTTONDOWN and ev.button == 1:
            nx_focus = self.nx_input.handle_event(ev)
            ny_focus = self.ny_input.handle_event(ev)
            if nx_focus or ny_focus:
                return
            if self.domain_apply_rect().collidepoint(ev.pos):
                self.apply_domain_settings()
                return
            if self.td3d_button_rect().collidepoint(ev.pos):
                self.export_3d_snapshot()
                return

        for sl in self.active_sliders():
            if sl.handle_event(ev):
                self._apply_slider_change(sl)
                return
        if self.s_speed.handle_event(ev):
            self.speed = self.s_speed.value
            return
        if self.s_courant.handle_event(ev):
            return
        if self.s_pml.handle_event(ev):
            return

        if ev.type == pygame.MOUSEWHEEL:
            mpos = pygame.mouse.get_pos()
            if self.canvas_rect.collidepoint(mpos):
                data_i, data_j = self.screen_to_grid(mpos)
                factor = 0.85 if ev.y > 0 else (1.0 / 0.85)
                self.zoom_at(data_i, data_j, factor)
        elif ev.type == pygame.MOUSEBUTTONDOWN and ev.button == 1:
            self.on_mouse_down(ev.pos)
        elif ev.type == pygame.MOUSEBUTTONDOWN and ev.button == 3:
            self._right_down = True
            self._pan_dragged = False
            self._pan_start_pos = ev.pos
            self._pan_start_view = (self.view_i0, self.view_i1, self.view_j0, self.view_j1)
        elif ev.type == pygame.MOUSEBUTTONUP and ev.button == 1:
            self.on_mouse_up(ev.pos)
        elif ev.type == pygame.MOUSEBUTTONUP and ev.button == 3:
            if self._right_down and not self._pan_dragged:
                self.on_right_click(ev.pos)
            self._right_down = False
        elif ev.type == pygame.MOUSEMOTION:
            if self._right_down:
                self._drag_pan(ev.pos)
            self.on_mouse_move(ev.pos)

    def _drag_pan(self, pos):
        dx = pos[0] - self._pan_start_pos[0]
        dy = pos[1] - self._pan_start_pos[1]
        if not self._pan_dragged and (dx * dx + dy * dy) > 16:
            self._pan_dragged = True
        if not self._pan_dragged:
            return
        i0, i1, j0, j1 = self._pan_start_view
        ppc_x = self.canvas_rect.w / (i1 - i0)
        ppc_y = self.canvas_rect.h / (j1 - j0)
        di, dj = dx / ppc_x, dy / ppc_y
        ni0, ni1 = self._clamp_span(i0 - di, i1 - di, self.nx)
        nj0, nj1 = self._clamp_span(j0 - dj, j1 - dj, self.ny)
        self.view_i0, self.view_i1 = ni0, ni1
        self.view_j0, self.view_j1 = nj0, nj1
        self._sync_field_axes_limits()

    def _apply_slider_change(self, sl):
        if self.selected and self.selected[0] == "source" and sl in (self.s_freq, self.s_amp):
            s = self.scene.sources[self.selected[1]]
            s.freq, s.amp = self.s_freq.value, self.s_amp.value
            self.scene.rebuild()
        elif self.selected and self.selected[0] == "array" and sl in (self.s_freq, self.s_amp):
            a = self.scene.arrays[self.selected[1]]
            a.freq, a.amp = self.s_freq.value, self.s_amp.value
            self.scene.rebuild()
        elif self.selected and self.selected[0] == "block" and sl in (self.s_eps, self.s_sigma):
            b = self.scene.blocks[self.selected[1]]
            b.eps_r, b.sigma = self.s_eps.value, self.s_sigma.value
            self.scene.rebuild()
        else:
            self.new_freq, self.new_amp = self.s_freq.value, self.s_amp.value
            self.new_eps, self.new_sigma = self.s_eps.value, self.s_sigma.value

    def handle_key(self, ev):
        mods = pygame.key.get_mods()
        if ev.key == pygame.K_SPACE:
            self.playing = not self.playing
        elif ev.key == pygame.K_s and mods & pygame.KMOD_CTRL:
            self.save_scene()
        elif ev.key == pygame.K_o and mods & pygame.KMOD_CTRL:
            self.load_scene()
        elif ev.key == pygame.K_s:
            if not self.playing:
                self.scene.sim.step()
                self.scene.sample_probes()
                self._force_plot_refresh = True
                self._force_field_refresh = True
        elif ev.key == pygame.K_r:
            self.scene.rebuild()
            self.vmax_running = 1e-6
            self._force_plot_refresh = True
            self._force_field_refresh = True
            self.toast("Felder zurückgesetzt")
        elif ev.key == pygame.K_g:
            self.show_grid = not self.show_grid
        elif ev.key == pygame.K_v:
            self.view_mode = "energy" if self.view_mode == "ez" else "ez"
        elif ev.key == pygame.K_HOME:
            self.reset_zoom()
        elif ev.key == pygame.K_p:
            self.export_png()
        elif ev.key == pygame.K_e:
            self.export_probe_csv()
        elif ev.key == pygame.K_d:
            self.export_3d_snapshot()
        elif ev.key == pygame.K_ESCAPE:
            self.selected = None
            self.dragging_new = None
            self.rect_start = None
            self.line_start = None
        elif ev.key in (pygame.K_DELETE, pygame.K_BACKSPACE):
            if self.selected:
                self.delete_item(self.selected)
        elif pygame.K_1 <= ev.key <= pygame.K_9:
            idx = ev.key - pygame.K_1
            if idx < len(TOOLS):
                self.tool = TOOLS[idx][0]
                self.selected = None

    # -------------------------------------------------------------- mouse
    def tools_top(self):
        return HEADER_H + 90

    def palette_rect_for(self, idx):
        return pygame.Rect(14, self.tools_top() + idx * 46, LEFT_W - 28, 40)

    def grid_preset_rects(self):
        rects = {}
        x = 20
        for idx, (label, nx, ny) in enumerate(GRID_PRESETS):
            w = 90
            rects[idx] = pygame.Rect(x, 44, w, 24)
            x += w + 8
        return rects

    def zoom_reset_rect(self):
        rects = self.grid_preset_rects()
        last = rects[len(GRID_PRESETS) - 1]
        return pygame.Rect(last.right + 16, 44, 150, 24)

    def scenario_rects(self):
        rects = {}
        x = 20
        for key, label in SCENARIOS:
            w = 14 + 8 * len(label)
            rects[key] = pygame.Rect(x, 72, w, 24)
            x += w + 8
        return rects

    def domain_apply_rect(self):
        return pygame.Rect(560, 96, 100, 28)

    def layer_rects(self):
        labels = [("lay_field", "Feld"), ("lay_mat", "Material"), ("lay_src", "Quellen"),
                  ("lay_probe", "Sonden"), ("lay_arr", "Arrays")]
        rects = {}
        x = 700
        for key, label in labels:
            w = 22 + 8 * len(label)
            rects[key] = pygame.Rect(x, 96, w, 28)
            x += w + 8
        return rects

    def td3d_button_rect(self):
        return pygame.Rect(1110, 96, 130, 28)

    def on_mouse_down(self, pos):
        if self.handle_wave_button_click(pos):
            return
        for idx, rect in self.grid_preset_rects().items():
            if rect.collidepoint(pos) and idx != self.grid_preset:
                self.set_grid_size(idx)
                return
        if self.zoom_reset_rect().collidepoint(pos):
            self.reset_zoom()
            return
        for key, rect in self.scenario_rects().items():
            if rect.collidepoint(pos):
                self.load_scenario(key)
                return
        layer_attr = {"lay_field": "layer_field", "lay_mat": "layer_material", "lay_src": "layer_sources",
                      "lay_probe": "layer_probes", "lay_arr": "layer_arrays"}
        for key, rect in self.layer_rects().items():
            if rect.collidepoint(pos):
                attr = layer_attr[key]
                setattr(self, attr, not getattr(self, attr))
                self._force_field_refresh = True
                return

        for idx, (kind, _label) in enumerate(TOOLS):
            if self.palette_rect_for(idx).collidepoint(pos):
                self.tool = kind
                self.selected = None
                if kind in POINT_TOOLS:
                    self.dragging_new = kind
                return

        if self.footer_button_rects()["clear"].collidepoint(pos):
            self.on_clear_clicked()
            return
        for name, rect in self.footer_button_rects().items():
            if name == "clear":
                continue
            if rect.collidepoint(pos):
                self.on_footer_button(name)
                return

        if not self.canvas_rect.collidepoint(pos):
            return

        if self.tool == "select":
            item = self.find_item_at(pos)
            self.selected = item
            if item is not None:
                self.moving = True
        elif self.tool in POINT_TOOLS:
            self.place_point_tool(pos, self.tool)
        elif self.tool in RECT_TOOLS:
            gi, gj = self.screen_to_grid(pos)
            self.rect_start = self.clamp_cell(gi, gj)
            self.rect_cur = self.rect_start
        elif self.tool in LINE_TOOLS:
            gi, gj = self.screen_to_grid(pos)
            self.line_start = self.clamp_cell(gi, gj)
            self.line_cur = self.line_start
        elif self.tool == "erase":
            item = self.find_item_at(pos)
            if item is not None:
                self.delete_item(item)

    def on_right_click(self, pos):
        """Rechtsklick waehlt ein Element unabhaengig vom aktiven Werkzeug aus,
        damit man Quellen/Arrays/Bloecke jederzeit nachtraeglich anklicken und
        ihre Parameter im linken Panel anpassen kann. (Rechtsklick-Ziehen
        verschiebt stattdessen den Bildausschnitt, siehe _drag_pan.)"""
        item = self.find_item_at(pos)
        self.selected = item
        self.moving = False

    def on_mouse_move(self, pos):
        if self.rect_start is not None and self.tool in RECT_TOOLS:
            gi, gj = self.screen_to_grid(pos)
            gi = min(max(gi, 0), self.nx - 1)
            gj = min(max(gj, 0), self.ny - 1)
            self.rect_cur = self.clamp_cell(gi, gj)
        elif self.line_start is not None and self.tool in LINE_TOOLS:
            gi, gj = self.screen_to_grid(pos)
            gi = min(max(gi, 0), self.nx - 1)
            gj = min(max(gj, 0), self.ny - 1)
            self.line_cur = self.clamp_cell(gi, gj)
        elif self.moving and self.selected is not None:
            gi, gj = self.screen_to_grid(pos)
            kind, idx = self.selected
            if kind == "source":
                i, j = self.clamp_cell(gi, gj)
                self.scene.sources[idx].i, self.scene.sources[idx].j = i, j
            elif kind == "probe":
                i, j = self.clamp_cell(gi, gj)
                self.scene.probes[idx].i, self.scene.probes[idx].j = i, j
            elif kind == "block":
                b = self.scene.blocks[idx]
                w, h = b.i1 - b.i0, b.j1 - b.j0
                i0 = int(min(max(round(gi - w / 2), 1), self.nx - 2 - w))
                j0 = int(min(max(round(gj - h / 2), 1), self.ny - 2 - h))
                b.i0, b.i1, b.j0, b.j1 = i0, i0 + w, j0, j0 + h
            elif kind == "array":
                a = self.scene.arrays[idx]
                hw, hh = (a.i1 - a.i0) / 2, (a.j1 - a.j0) / 2
                a.i0, a.i1 = gi - hw, gi + hw
                a.j0, a.j1 = gj - hh, gj + hh

    def on_mouse_up(self, pos):
        if self.dragging_new is not None:
            kind = self.dragging_new
            self.dragging_new = None
            if self.canvas_rect.collidepoint(pos):
                self.place_point_tool(pos, kind)
        if self.rect_start is not None and self.tool in RECT_TOOLS:
            i0, j0 = self.rect_start
            i1, j1 = self.rect_cur if self.rect_cur else self.rect_start
            i0, i1 = sorted((i0, i1))
            j0, j1 = sorted((j0, j1))
            i1, j1 = max(i1, i0 + 1), max(j1, j0 + 1)
            if self.tool == "dielectric":
                self.scene.blocks.append(BlockItem("dielectric", i0, i1, j0, j1,
                                                     eps_r=self.new_eps, sigma=self.new_sigma))
                self.toast(f"Dielektrikum eps_r={self.new_eps:.1f} platziert")
            else:
                self.scene.blocks.append(BlockItem("pec", i0, i1, j0, j1))
                self.toast("PEC-Wand platziert")
            self.scene.rebuild()
            self.rect_start = None
            self.rect_cur = None
        if self.line_start is not None and self.tool in LINE_TOOLS:
            i0, j0 = self.line_start
            i1, j1 = self.line_cur if self.line_cur else self.line_start
            if math.hypot(i1 - i0, j1 - j0) < 2:
                self.toast("Linie zu kurz für ein Array")
            else:
                self.scene.arrays.append(ArrayItem(i0, j0, i1, j1, self.new_waveform,
                                                     self.new_freq, self.new_amp))
                self.scene.rebuild()
                self.toast("Ebene Welle (Quellen-Array) platziert")
            self.line_start = None
            self.line_cur = None
        self.moving = False

    def place_point_tool(self, pos, kind):
        gi, gj = self.screen_to_grid(pos)
        i, j = self.clamp_cell(gi, gj)
        if kind == "probe":
            self.scene.probes.append(ProbeItem(i, j))
            self.toast("Sonde platziert")
        else:
            self.scene.sources.append(SourceItem(kind, i, j, self.new_waveform, self.new_freq, self.new_amp))
            label = {"soft": "Weiche Quelle", "hard": "Harte Quelle", "current": "Stromquelle"}[kind]
            self.toast(f"{label} platziert")
        self.scene.rebuild()

    def on_clear_clicked(self):
        if time.time() < self.clear_confirm_until:
            self.scene.clear_all()
            self.selected = None
            self.clear_confirm_until = 0.0
            self._force_plot_refresh = True
            self._force_field_refresh = True
            self.toast("Szene geleert")
        else:
            self.clear_confirm_until = time.time() + 2.0
            self.toast("Nochmal klicken zum Bestätigen")

    def on_footer_button(self, name):
        if name == "play":
            self.playing = not self.playing
        elif name == "step":
            if not self.playing:
                self.scene.sim.step()
                self.scene.sample_probes()
                self._force_plot_refresh = True
                self._force_field_refresh = True
        elif name == "reset":
            self.scene.rebuild()
            self.vmax_running = 1e-6
            self._force_plot_refresh = True
            self._force_field_refresh = True
            self.toast("Felder zurückgesetzt")
        elif name == "view":
            self.view_mode = "energy" if self.view_mode == "ez" else "ez"
        elif name == "grid":
            self.show_grid = not self.show_grid
        elif name == "png":
            self.export_png()
        elif name == "csv":
            self.export_probe_csv()

    def save_scene(self):
        path = Path(__file__).with_name("playground_scene.json")
        data = self.scene.to_dict()
        data["grid_preset"] = self.grid_preset
        data["nx"], data["ny"] = self.nx, self.ny
        data["npml"], data["courant"] = self.npml, self.courant
        path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        self.toast(f"Szene gespeichert ({path.name})")

    def load_scene(self):
        path = Path(__file__).with_name("playground_scene.json")
        if not path.exists():
            self.toast("Keine gespeicherte Szene gefunden")
            return
        data = json.loads(path.read_text(encoding="utf-8"))
        nx, ny = data.get("nx"), data.get("ny")
        npml = data.get("npml", self.npml)
        courant = data.get("courant", self.courant)
        if nx and ny and (nx != self.nx or ny != self.ny or npml != self.npml or courant != self.courant):
            self.grid_preset = data.get("grid_preset", -1)
            self._rebuild_domain(nx, ny, npml, courant)
            self.nx_input.text, self.ny_input.text = str(nx), str(ny)
            self.s_courant.value, self.s_pml.value = courant, npml
        self.scene.load_dict(data)
        self.selected = None
        self._force_plot_refresh = True
        self._force_field_refresh = True
        self.toast(f"Szene geladen ({path.name})")

    def export_png(self):
        self.refresh_field_figure()
        path = Path(__file__).with_name(f"fdtd_snapshot_{int(time.time())}.png")
        pygame.image.save(self._field_surf, str(path))
        self.toast(f"Bild gespeichert ({path.name})")

    def export_probe_csv(self):
        idx = self.active_probe_index()
        if idx < 0:
            self.toast("Keine Sonde vorhanden")
            return
        p = self.scene.probes[idx]
        path = Path(__file__).with_name(f"probe_{idx + 1}_{int(time.time())}.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["t_s", "Ez"])
            for t, e in zip(p.t_buf, p.e_buf):
                w.writerow([t, e])
        self.toast(f"CSV exportiert ({path.name})")

    def export_3d_snapshot(self):
        """Oeffnet den aktuellen Feld-Schnappschuss als interaktiven 3D-
        Oberflaechenplot (x, y, z = Feldwert). Laeuft im selben Skript
        (fdtd_playground.py --3d-view <datei>), aber als eigener Prozess mit
        eigenem interaktiven Matplotlib-Backend (TkAgg) -- das Agg-Backend
        der Haupt-App erlaubt kein interaktives Fenster im selben Prozess."""
        sim = self.scene.sim
        if self.view_mode == "ez":
            z = sim.Ez.copy()
            label = "Ez (a.u.)"
        else:
            z = energy_density(sim)
            label = "Energiedichte (a.u.)"
        path = Path(__file__).with_name(f"_3d_snapshot_{int(time.time() * 1000)}.npz")
        np.savez(path, Z=z, d_cm=self.scene.d * 100, t_ns=sim.t * 1e9, label=label)
        try:
            subprocess.Popen([sys.executable, __file__, "--3d-view", str(path)])
            self.toast("3D-Ansicht wird geöffnet …")
        except OSError:
            self.toast("3D-Ansicht konnte nicht gestartet werden")

    # -------------------------------------------------------------- drawing
    def footer_button_rects(self):
        y = self.win_h - FOOTER_H + 14
        rx = self.win_w - RIGHT_W - 20
        gain_x = rx - 170
        grid_x = gain_x - 14 - 70
        view_x = grid_x - 10 - 90
        csv_x = view_x - 14 - 50
        png_x = csv_x - 10 - 50
        return {
            "play": pygame.Rect(LEFT_W + 20, y, 46, 46),
            "step": pygame.Rect(LEFT_W + 76, y, 46, 46),
            "reset": pygame.Rect(LEFT_W + 132, y, 46, 46),
            "png": pygame.Rect(png_x, y, 50, 46),
            "csv": pygame.Rect(csv_x, y, 50, 46),
            "view": pygame.Rect(view_x, y, 90, 46),
            "grid": pygame.Rect(grid_x, y, 70, 46),
            "clear": pygame.Rect(20, 12, LEFT_W - 40, 30),
        }

    def active_probe_index(self):
        if self.selected and self.selected[0] == "probe":
            return self.selected[1]
        return len(self.scene.probes) - 1 if self.scene.probes else -1

    def draw(self):
        screen = self.screen
        screen.fill(BG)
        self.draw_header()
        self.draw_palette()
        self.draw_param_panel()
        self.draw_canvas()
        self.draw_plots()
        self.draw_footer()
        self.draw_toasts()
        pygame.display.flip()

    def draw_header(self):
        screen = self.screen
        rrect(screen, PANEL_BG, pygame.Rect(0, 0, self.win_w, HEADER_H), radius=0)
        pygame.draw.line(screen, BORDER, (0, HEADER_H), (self.win_w, HEADER_H))
        screen.blit(self.font_title.render("FDTD Wellen-Baukasten", True, TEXT), (20, 8))
        sim = self.scene.sim
        domain_cm = (self.nx * self.scene.d * 100, self.ny * self.scene.d * 100)
        stats = (f"Schritt {sim.time_step}   t = {sim.t * 1e9:8.2f} ns   "
                 f"Energie {sim.total_energy():.2e} J/m   {domain_cm[0]:.0f}×{domain_cm[1]:.0f} cm   "
                 f"{self.fps_smooth:4.0f} FPS")
        txt = self.font.render(stats, True, TEXT_DIM)
        screen.blit(txt, (self.win_w - txt.get_width() - 20, 12))

        screen.blit(self.font_small.render("Gitter:", True, TEXT_DIM), (20, 24))
        for idx, rect in self.grid_preset_rects().items():
            label, nx, ny = GRID_PRESETS[idx]
            active = idx == self.grid_preset
            rrect(screen, ACCENT if active else PANEL_BG_LIGHT, rect, radius=6)
            t = self.font_small.render(f"{label} {nx}×{ny}", True, (10, 12, 16) if active else TEXT_DIM)
            screen.blit(t, t.get_rect(center=rect.center))

        zr = self.zoom_reset_rect()
        rrect(screen, PANEL_BG_LIGHT, zr, radius=6)
        t = self.font_small.render("Zoom zurücksetzen (Pos1)", True, TEXT_DIM)
        screen.blit(t, t.get_rect(center=zr.center))

        screen.blit(self.font_small.render("Szenario:", True, TEXT_DIM), (20, 78))
        for key, rect in self.scenario_rects().items():
            label = dict(SCENARIOS)[key]
            active = key == self.scenario
            rrect(screen, ACCENT if active else PANEL_BG_LIGHT, rect, radius=6)
            t = self.font_small.render(label, True, (10, 12, 16) if active else TEXT_DIM)
            screen.blit(t, t.get_rect(center=rect.center))
        hint = "Fenster kann per Rand/Maximieren vergrößert werden – mehr Platz fürs Feld"
        t = self.font_small.render(hint, True, TEXT_DIM)
        screen.blit(t, (self.win_w - t.get_width() - 20, 82))

        screen.blit(self.font_small.render("Größe:", True, TEXT_DIM), (20, 104))
        self.nx_input.draw(screen, self.font_small)
        screen.blit(self.font_small.render("×", True, TEXT_DIM), (128, 104))
        self.ny_input.draw(screen, self.font_small)
        screen.blit(self.font_small.render("Zellen", True, TEXT_DIM), (200, 104))
        self.s_courant.draw(screen, self.font_small)
        self.s_pml.draw(screen, self.font_small)
        ar = self.domain_apply_rect()
        rrect(screen, ACCENT, ar, radius=6)
        t = self.font_small.render("Anwenden", True, (10, 12, 16))
        screen.blit(t, t.get_rect(center=ar.center))

        flags = {"lay_field": "layer_field", "lay_mat": "layer_material", "lay_src": "layer_sources",
                 "lay_probe": "layer_probes", "lay_arr": "layer_arrays"}
        for key, rect in self.layer_rects().items():
            active = getattr(self, flags[key])
            rrect(screen, ACCENT if active else PANEL_BG_LIGHT, rect, radius=6)
            label = {"lay_field": "Feld", "lay_mat": "Material", "lay_src": "Quellen",
                     "lay_probe": "Sonden", "lay_arr": "Arrays"}[key]
            t = self.font_small.render(label, True, (10, 12, 16) if active else TEXT_DIM)
            screen.blit(t, t.get_rect(center=rect.center))

        td3d_rect = self.td3d_button_rect()
        rrect(screen, PANEL_BG_LIGHT, td3d_rect, radius=6)
        t = self.font_small.render("3D-Ansicht", True, TEXT_DIM)
        screen.blit(t, t.get_rect(center=td3d_rect.center))

    def draw_palette(self):
        screen = self.screen
        rrect(screen, PANEL_BG, pygame.Rect(0, HEADER_H, LEFT_W, self.win_h - HEADER_H), radius=0)
        pygame.draw.line(screen, BORDER, (LEFT_W, HEADER_H), (LEFT_W, self.win_h))

        btn = self.footer_button_rects()["clear"]
        btn = pygame.Rect(btn.x, HEADER_H + 14, btn.w, btn.h)
        confirming = time.time() < self.clear_confirm_until
        rrect(screen, DANGER if confirming else PANEL_BG_LIGHT, btn, radius=8)
        rrect(screen, BORDER, btn, radius=8, width=1)
        label = "Sicher? Nochmal klicken" if confirming else "Szene leeren"
        t = self.font_small.render(label, True, TEXT if not confirming else (255, 255, 255))
        screen.blit(t, t.get_rect(center=btn.center))

        screen.blit(self.font_bold.render("Werkzeuge", True, TEXT_DIM), (18, HEADER_H + 60))
        for idx, (kind, label) in enumerate(TOOLS):
            rect = self.palette_rect_for(idx)
            active = self.tool == kind
            hover = rect.collidepoint(pygame.mouse.get_pos())
            bg = PANEL_BG_LIGHT if (active or hover) else PANEL_BG
            rrect(screen, bg, rect, radius=8)
            if active:
                rrect(screen, TOOL_COLORS[kind], rect, radius=8, width=2)
            draw_tool_icon(screen, kind, (rect.x + 24, rect.centery), 10, TOOL_COLORS[kind])
            txt = self.font.render(label, True, TEXT if active else TEXT_DIM)
            screen.blit(txt, (rect.x + 48, rect.centery - txt.get_height() // 2))
            screen.blit(self.font_small.render(str(idx + 1), True, TEXT_DIM), (rect.right - 16, rect.y + 4))

    def draw_param_panel(self):
        screen = self.screen
        y0 = self.tools_top() + len(TOOLS) * 46 + 20
        pygame.draw.line(screen, BORDER, (16, y0), (LEFT_W - 16, y0))
        y0 += 14

        heading = "Ausgewähltes Element" if self.selected else "Neue Elemente"
        screen.blit(self.font_bold.render(heading, True, TEXT_DIM), (20, y0))
        y0 += 26

        sliders = self.active_sliders()
        want_wave = (self.selected and self.selected[0] in ("source", "array")) or \
            self.tool in {"soft", "hard", "current", "planewave"}
        if want_wave:
            if self.selected and self.selected[0] == "source":
                wf = self.scene.sources[self.selected[1]].waveform
            elif self.selected and self.selected[0] == "array":
                wf = self.scene.arrays[self.selected[1]].waveform
            else:
                wf = self.new_waveform
            screen.blit(self.font.render("Signalform", True, TEXT_DIM), (20, y0))
            y0 += 20
            bw = (LEFT_W - 40 - 12) / 3
            self._wave_button_rects = {}
            for k, wave in enumerate(WAVEFORMS):
                rect = pygame.Rect(20 + k * (bw + 6), y0, bw, 28)
                active = wf == wave
                rrect(screen, ACCENT if active else PANEL_BG_LIGHT, rect, radius=6)
                t = self.font_small.render(WAVE_LABELS[wave], True, (10, 12, 16) if active else TEXT_DIM)
                screen.blit(t, t.get_rect(center=rect.center))
                self._wave_button_rects[wave] = rect
            y0 += 40
            if wf == "gaussian":
                note = self.font_small.render("Gauss hat DC-Anteil: Restfeld", True, TEXT_DIM)
                note2 = self.font_small.render("kann lange nachleuchten.", True, TEXT_DIM)
                screen.blit(note, (20, y0)); screen.blit(note2, (20, y0 + 15))
                y0 += 34

        is_array_context = (self.selected and self.selected[0] == "array") or self.tool == "planewave"
        if is_array_context:
            for line in ("Strahlt senkrecht zur Linie ab (beide",
                         "Richtungen, s. Pfeile). Reicht bis in",
                         "die PML, damit die Enden nicht beugen.",
                         "Nah am Array: einzelne Kreiswellen –",
                         "erst mit Abstand ebene Welle (Huygens)."):
                screen.blit(self.font_small.render(line, True, TEXT_DIM), (20, y0))
                y0 += 15
            y0 += 5

        for sl in sliders:
            sl.rect.y = y0
            sl.track.y = y0 + 24
            sl.draw(screen, self.font)
            y0 += 46

        if self.selected and self.selected[0] == "array":
            n = len(array_elements(self.scene.arrays[self.selected[1]], self.nx, self.ny))
            info = self.font_small.render(f"{n} Einzelquellen im Array", True, TEXT_DIM)
            screen.blit(info, (20, y0))
            y0 += 20

        if self.selected:
            hint = self.font_small.render("Entf löscht  ·  Esc hebt Auswahl auf", True, TEXT_DIM)
            screen.blit(hint, (20, y0 + 4))
        elif not is_array_context:
            hint_lines = [
                "Punkte: klicken/ziehen zum Platzieren.",
                "Flächen/Ebene Welle: aufziehen.",
                "Rechtsklick wählt aus, Ziehen verschiebt",
                "den Ausschnitt. Mausrad zoomt.",
            ]
            for k, line in enumerate(hint_lines):
                screen.blit(self.font_small.render(line, True, TEXT_DIM), (20, y0 + 4 + k * 16))

    def handle_wave_button_click(self, pos):
        for wave, rect in self._wave_button_rects.items():
            if rect.collidepoint(pos):
                if self.selected and self.selected[0] == "source":
                    self.scene.sources[self.selected[1]].waveform = wave
                    self.scene.rebuild()
                elif self.selected and self.selected[0] == "array":
                    self.scene.arrays[self.selected[1]].waveform = wave
                    self.scene.rebuild()
                else:
                    self.new_waveform = wave
                return True
        return False

    # -------------------------------------------------------------- field view
    def refresh_field_figure(self):
        self._sync_field_axes_limits()
        sim = self.scene.sim
        self._norm_field.decades = self.color_decades
        if self.view_mode == "ez":
            mag = np.abs(sim.Ez)
            cur_max = float(np.max(mag)) + 1e-30
            self.vmax_running = max(cur_max, self.vmax_running * FIELD_SCALE_DECAY, FIELD_FLOOR)
            data = np.clip(mag / (self.vmax_running / self.gain), 0, 1)
        else:
            e = energy_density(sim)
            cur_max = float(np.max(e)) + 1e-30
            self.vmax_running = max(cur_max, self.vmax_running * FIELD_SCALE_DECAY, FIELD_FLOOR)
            data = np.clip(e / (self.vmax_running / self.gain), 0, 1)
        self.im_field.set_data(data.T)
        self.im_field.set_visible(self.layer_field)
        self._update_field_colorbar()
        self._rebuild_field_overlays()
        self._field_surf = figure_to_surface(self.canvas_field)
        self._field_px_size = self.canvas_field.get_width_height()

    def _update_field_colorbar(self):
        cbar = self.cbar_field
        cbar.update_normal(self.im_field)
        cbar.set_ticks([0.0, 0.001, 0.01, 0.1, 1.0])
        cbar.set_ticklabels(["0", "0.001", "0.01", "0.1", "1"])
        label = "|Ez| / Skalenmaximum (log, a.u.)" if self.view_mode == "ez" else "Energiedichte / Max. (log, a.u.)"
        cbar.set_label(label, color=rgb01(TEXT_DIM), fontsize=8)
        cbar.ax.tick_params(colors=rgb01(TEXT_DIM), labelsize=8)

    def _rebuild_field_overlays(self):
        ax = self.ax_field
        for artist in self._field_artists:
            artist.remove()
        self._field_artists = []

        def add(artist):
            self._field_artists.append(artist)
            return artist

        nx, ny, npml = self.nx, self.ny, self.npml
        alpha = 0.55
        add(ax.add_patch(Rectangle((0, 0), nx, npml, facecolor="black", edgecolor="none", alpha=alpha, zorder=3)))
        add(ax.add_patch(Rectangle((0, ny - npml), nx, npml, facecolor="black", edgecolor="none", alpha=alpha, zorder=3)))
        add(ax.add_patch(Rectangle((0, 0), npml, ny, facecolor="black", edgecolor="none", alpha=alpha, zorder=3)))
        add(ax.add_patch(Rectangle((nx - npml, 0), npml, ny, facecolor="black", edgecolor="none", alpha=alpha, zorder=3)))

        if self.show_grid:
            step = 10
            xs = np.arange(0, nx + 1, step)
            ys = np.arange(0, ny + 1, step)
            add(ax.vlines(xs, 0, ny, colors=[(1, 1, 1, 0.12)], linewidth=0.6, zorder=4))
            add(ax.hlines(ys, 0, nx, colors=[(1, 1, 1, 0.12)], linewidth=0.6, zorder=4))

        if self.layer_material:
            for k, b in enumerate(self.scene.blocks):
                w, h = b.i1 - b.i0, b.j1 - b.j0
                if b.kind == "dielectric":
                    shade = min(1.0, (b.eps_r - 1) / 11)
                    col = (lerp(0.24, 0.08, shade), lerp(0.72, 0.36, shade), lerp(0.55, 0.28, shade))
                    rect = Rectangle((b.i0, b.j0), w, h, facecolor=col, edgecolor=rgb01(ACCENT_DIEL),
                                      linewidth=1.6, alpha=0.6, zorder=5)
                else:
                    rect = Rectangle((b.i0, b.j0), w, h, facecolor=(0.24, 0.12, 0.12), edgecolor=rgb01(ACCENT_PEC),
                                      linewidth=1.6, alpha=0.9, zorder=5)
                add(ax.add_patch(rect))
                if self.selected == ("block", k):
                    hi = Rectangle((b.i0 - 0.6, b.j0 - 0.6), w + 1.2, h + 1.2, facecolor="none",
                                    edgecolor=rgb01(TEXT), linewidth=1.6, linestyle="--", zorder=6)
                    add(ax.add_patch(hi))

        if self.layer_arrays:
            for k, a in enumerate(self.scene.arrays):
                ai0, aj0, ai1, aj1 = effective_array_endpoints(a, self.nx, self.ny)
                lw = 4 if self.selected == ("array", k) else 2.6
                add(ax.plot([ai0, ai1], [aj0, aj1], color=rgb01(ACCENT_PLANE),
                            linewidth=lw, solid_capstyle="round", zorder=6)[0])
                elems = array_elements(a, self.nx, self.ny)
                add(ax.scatter([e[0] for e in elems], [e[1] for e in elems],
                               s=16, color=rgb01(ACCENT_PLANE), zorder=7))
                di, dj = a.i1 - a.i0, a.j1 - a.j0
                length = math.hypot(di, dj) or 1.0
                nx_, ny_ = -dj / length, di / length    # Breitseiten-Normale
                mx, my = (a.i0 + a.i1) / 2, (a.j0 + a.j1) / 2
                arrow_len = max(4.0, length * 0.18)
                for side in (1, -1):
                    add(ax.annotate("", xy=(mx + nx_ * arrow_len * side, my + ny_ * arrow_len * side),
                                    xytext=(mx, my),
                                    arrowprops=dict(arrowstyle="-|>", color=rgb01(ACCENT_LIGHT), linewidth=1.6),
                                    zorder=7))
                if self.selected == ("array", k):
                    add(ax.plot([ai0, ai1], [aj0, aj1], "o", markersize=9, markerfacecolor="none",
                                markeredgecolor=rgb01(TEXT), markeredgewidth=1.8, zorder=8)[0])

        if self.layer_sources:
            for k, s in enumerate(self.scene.sources):
                marker = "^" if s.kind == "current" else "o"
                filled = s.kind != "soft"
                add(ax.plot(s.i, s.j, marker=marker, markersize=9,
                            markerfacecolor=rgb01(TOOL_COLORS[s.kind]) if filled else "none",
                            markeredgecolor=rgb01(TOOL_COLORS[s.kind]), markeredgewidth=1.8, zorder=8)[0])
                if self.selected == ("source", k):
                    add(ax.plot(s.i, s.j, "o", markersize=16, markerfacecolor="none",
                                markeredgecolor=rgb01(TEXT), markeredgewidth=1.8, zorder=9)[0])

        if self.layer_probes:
            active_probe = self.active_probe_index()
            pulse = 3 + 2 * math.sin(time.time() * 4)
            for k, p in enumerate(self.scene.probes):
                add(ax.plot(p.i, p.j, marker="s", markersize=8, markerfacecolor=rgb01(probe_color(k)),
                            markeredgecolor=(0, 0, 0), markeredgewidth=0.8, zorder=8)[0])
                if k == active_probe:
                    add(ax.plot(p.i, p.j, "o", markersize=14 + pulse, markerfacecolor="none",
                                markeredgecolor=rgb01(ACCENT_LIGHT), markeredgewidth=1.4, zorder=8)[0])
                if self.selected == ("probe", k):
                    add(ax.plot(p.i, p.j, "o", markersize=17, markerfacecolor="none",
                                markeredgecolor=rgb01(TEXT), markeredgewidth=1.8, zorder=9)[0])

        hover_item = self.find_item_at(pygame.mouse.get_pos())
        if hover_item and hover_item != self.selected and hover_item[0] in ("source", "probe") and \
                ((hover_item[0] == "source" and self.layer_sources) or (hover_item[0] == "probe" and self.layer_probes)):
            items = self.scene.sources if hover_item[0] == "source" else self.scene.probes
            it = items[hover_item[1]]
            add(ax.plot(it.i, it.j, "o", markersize=13, markerfacecolor="none",
                        markeredgecolor=rgb01(TEXT_DIM), markeredgewidth=1.2, zorder=8)[0])

    def draw_canvas(self):
        screen = self.screen
        rrect(screen, (10, 11, 15), self.canvas_rect, radius=0)

        due = (time.time() - self._last_field_refresh > FIELD_REFRESH_DT) or self._force_field_refresh
        if due or self._field_surf is None:
            self.refresh_field_figure()
            self._last_field_refresh = time.time()
            self._force_field_refresh = False
        if self._field_surf is not None:
            screen.blit(self._field_surf, self.canvas_rect.topleft)

        if self.rect_start is not None and self.rect_cur is not None:
            i0, j0 = self.rect_start
            i1, j1 = self.rect_cur
            i0, i1 = sorted((i0, i1))
            j0, j1 = sorted((j0, j1))
            x0, y0 = self.grid_to_screen(i0, j0)
            x1, y1 = self.grid_to_screen(i1, j1)
            col = ACCENT_DIEL if self.tool == "dielectric" else ACCENT_PEC
            rect = pygame.Rect(min(x0, x1), min(y0, y1), abs(x1 - x0), abs(y1 - y0))
            pygame.draw.rect(screen, col, rect, width=2)

        if self.line_start is not None and self.line_cur is not None:
            x0, y0 = self.grid_to_screen(*self.line_start)
            x1, y1 = self.grid_to_screen(*self.line_cur)
            pygame.draw.line(screen, ACCENT_PLANE, (x0, y0), (x1, y1), 3)

        if self.dragging_new is not None:
            mx, my = pygame.mouse.get_pos()
            if self.canvas_rect.collidepoint((mx, my)):
                draw_tool_icon(screen, self.dragging_new, (mx, my), 11, TOOL_COLORS[self.dragging_new])

        mx, my = pygame.mouse.get_pos()
        readout = ""
        if self.canvas_rect.collidepoint((mx, my)) and not self._right_down:
            gi, gj = self.screen_to_grid((mx, my))
            if 0 <= gi <= self.nx and 0 <= gj <= self.ny:
                ii = min(max(int(round(gi)), 0), self.nx - 1)
                jj = min(max(int(round(gj)), 0), self.ny - 1)
                ez_val = self.scene.sim.Ez[ii, jj]
                x_cm, y_cm = gi * self.scene.d * 100, gj * self.scene.d * 100
                readout = f"x={x_cm:5.1f} cm  y={y_cm:5.1f} cm  Ez={ez_val:+.3f} a.u."
        legend = (f"{'|Ez| (a.u.)' if self.view_mode == 'ez' else 'Energiedichte (a.u.)'}  ·  "
                  f"PML schattiert  ·  Mausrad = Zoom  ·  Rechtsklick-Ziehen = Verschieben")
        screen.blit(self.font_small.render(legend, True, TEXT_DIM), (self.canvas_rect.x + 8, self.canvas_rect.bottom - 36))
        if readout:
            screen.blit(self.font_small.render(readout, True, TEXT_DIM), (self.canvas_rect.x + 8, self.canvas_rect.bottom - 18))

    # -------------------------------------------------------------- side plots
    def draw_plots(self):
        screen = self.screen
        rx = self.win_w - RIGHT_W
        rrect(screen, PANEL_BG, pygame.Rect(rx, HEADER_H, RIGHT_W, self.win_h - HEADER_H), radius=0)
        pygame.draw.line(screen, BORDER, (rx, HEADER_H), (rx, self.win_h))

        pad = 18
        screen.blit(self.font_bold.render("Sonden", True, TEXT), (rx + pad, HEADER_H + 14))

        due = (time.time() - self._last_plot_refresh > PLOT_REFRESH_DT) or self._force_plot_refresh
        if due:
            self.refresh_time_figure()
            self.refresh_spec_figure()
            self._last_plot_refresh = time.time()
            self._force_plot_refresh = False

        content_y = HEADER_H + 14 + 26
        if self._time_surf is not None:
            screen.blit(self._time_surf, (rx + pad, content_y))
        if self._spec_surf is not None:
            screen.blit(self._spec_surf, (rx + pad, content_y + 236 + 14))
        idx = self.active_probe_index()
        if idx >= 0:
            info = f"Aktive Sonde: Sonde {idx + 1} bei Zelle ({self.scene.probes[idx].i}, {self.scene.probes[idx].j})"
            screen.blit(self.font_small.render(info, True, TEXT_DIM), (rx + pad, content_y + 236 + 14 + 236 + 10))

    def refresh_time_figure(self):
        ax = self.ax_time
        ax.cla()
        style_axes(ax)
        active = self.active_probe_index()
        any_data = False
        for idx, p in enumerate(self.scene.probes):
            if len(p.e_buf) < 2:
                continue
            any_data = True
            t = np.array(p.t_buf)[-TIME_PLOT_SAMPLES:] * 1e9
            e = np.array(p.e_buf)[-TIME_PLOT_SAMPLES:]
            lw = 2.2 if idx == active else 1.1
            alpha = 1.0 if idx == active else 0.55
            ax.plot(t, e, color=rgb01(probe_color(idx)), linewidth=lw, alpha=alpha, label=f"Sonde {idx + 1}")
        ax.set_xlabel("t (ns)", fontsize=8)
        ax.set_ylabel("Ez (a.u.)", fontsize=8)
        ax.set_title(f"Zeitsignal Ez(t)  (letzte {TIME_PLOT_SAMPLES} Schritte)", color=rgb01(TEXT), fontsize=10)
        if any_data:
            leg = ax.legend(fontsize=7, facecolor=rgb01(PANEL_BG_LIGHT), edgecolor=rgb01(BORDER), loc="upper right")
            for txt in leg.get_texts():
                txt.set_color(rgb01(TEXT))
        else:
            ax.text(0.5, 0.5, "Sonden-Werkzeug wählen\nund ins Feld klicken",
                    ha="center", va="center", color=rgb01(TEXT_DIM), fontsize=9, transform=ax.transAxes)
        self.fig_time.tight_layout(pad=0.5)
        self._time_surf = figure_to_surface(self.canvas_time)

    def refresh_spec_figure(self):
        ax = self.ax_spec
        ax.cla()
        style_axes(ax)
        active = self.active_probe_index()
        any_data = False
        fmax_ghz = 25.0
        for idx, p in enumerate(self.scene.probes):
            if len(p.e_buf) < 16:
                continue
            any_data = True
            e = np.array(p.e_buf)
            n = len(e)
            spec = np.abs(np.fft.rfft(e * np.hanning(n)))
            freqs = np.fft.rfftfreq(n, self.scene.sim.dt) / 1e9
            mask = freqs <= fmax_ghz
            spec_m = spec[mask]
            if spec_m.max() > 0:
                spec_m = spec_m / spec_m.max()
            lw = 2.2 if idx == active else 1.1
            alpha = 1.0 if idx == active else 0.55
            ax.plot(freqs[mask], spec_m, color=rgb01(probe_color(idx)), linewidth=lw, alpha=alpha, label=f"Sonde {idx + 1}")
        ax.set_xlabel("f (GHz)", fontsize=8)
        ax.set_ylabel("|E(f)| (norm.)", fontsize=8)
        ax.set_title("Spektrum", color=rgb01(TEXT), fontsize=10)
        ax.set_xlim(0, fmax_ghz)
        if any_data:
            leg = ax.legend(fontsize=7, facecolor=rgb01(PANEL_BG_LIGHT), edgecolor=rgb01(BORDER), loc="upper right")
            for txt in leg.get_texts():
                txt.set_color(rgb01(TEXT))
        else:
            ax.text(0.5, 0.5, "Keine Sondendaten", ha="center", va="center",
                    color=rgb01(TEXT_DIM), fontsize=9, transform=ax.transAxes)
        self.fig_spec.tight_layout(pad=0.5)
        self._spec_surf = figure_to_surface(self.canvas_spec)

    def draw_footer(self):
        screen = self.screen
        y = self.win_h - FOOTER_H
        rrect(screen, PANEL_BG, pygame.Rect(0, y, self.win_w, FOOTER_H), radius=0)
        pygame.draw.line(screen, BORDER, (0, y), (self.win_w, y))

        rects = self.footer_button_rects()
        play_rect = rects["play"]
        rrect(screen, ACCENT if self.playing else PANEL_BG_LIGHT, play_rect, radius=10)
        cx, cy = play_rect.center
        if self.playing:
            pygame.draw.rect(screen, (10, 12, 16), (cx - 8, cy - 8, 6, 16))
            pygame.draw.rect(screen, (10, 12, 16), (cx + 2, cy - 8, 6, 16))
        else:
            pygame.draw.polygon(screen, TEXT, [(cx - 7, cy - 9), (cx - 7, cy + 9), (cx + 9, cy)])

        step_rect = rects["step"]
        rrect(screen, PANEL_BG_LIGHT, step_rect, radius=10)
        cx, cy = step_rect.center
        pygame.draw.polygon(screen, TEXT, [(cx - 8, cy - 8), (cx - 8, cy + 8), (cx + 3, cy)])
        pygame.draw.rect(screen, TEXT, (cx + 4, cy - 8, 3, 16))

        reset_rect = rects["reset"]
        rrect(screen, PANEL_BG_LIGHT, reset_rect, radius=10)
        pygame.draw.circle(screen, TEXT, reset_rect.center, 9, 2)
        pygame.draw.polygon(screen, TEXT, [(reset_rect.centerx + 6, reset_rect.centery - 9),
                                            (reset_rect.centerx + 12, reset_rect.centery - 9),
                                            (reset_rect.centerx + 9, reset_rect.centery - 3)])

        png_rect = rects["png"]
        rrect(screen, PANEL_BG_LIGHT, png_rect, radius=8)
        pygame.draw.rect(screen, TEXT, pygame.Rect(png_rect.centerx - 9, png_rect.centery - 6, 18, 13), 2, border_radius=2)
        pygame.draw.circle(screen, TEXT, png_rect.center, 4, 1)
        t = self.font_small.render("P", True, TEXT_DIM)
        screen.blit(t, (png_rect.right - 12, png_rect.bottom - 14))

        csv_rect = rects["csv"]
        rrect(screen, PANEL_BG_LIGHT, csv_rect, radius=8)
        for row in range(3):
            yy = csv_rect.centery - 7 + row * 6
            pygame.draw.line(screen, TEXT, (csv_rect.centerx - 8, yy), (csv_rect.centerx + 8, yy), 1)
        t = self.font_small.render("E", True, TEXT_DIM)
        screen.blit(t, (csv_rect.right - 12, csv_rect.bottom - 14))

        self.s_speed.rect.x, self.s_speed.track.x = LEFT_W + 200, LEFT_W + 200
        self.s_speed.rect.y, self.s_speed.track.y = y + 8, y + 30
        self.s_speed.draw(screen, self.font)

        view_rect = rects["view"]
        rrect(screen, PANEL_BG_LIGHT, view_rect, radius=8)
        label = "Ez" if self.view_mode == "ez" else "Energie"
        t = self.font_small.render(f"Ansicht: {label}", True, TEXT)
        screen.blit(t, t.get_rect(center=view_rect.center))

        grid_rect = rects["grid"]
        rrect(screen, ACCENT if self.show_grid else PANEL_BG_LIGHT, grid_rect, radius=8)
        t = self.font_small.render("Gitter", True, (10, 12, 16) if self.show_grid else TEXT)
        screen.blit(t, t.get_rect(center=grid_rect.center))

        hint = "Leertaste Play/Pause · S Schritt · R Reset · G Gitter · V Ansicht · Mausrad Zoom · Rechtsklick-Ziehen Pan · P PNG · E CSV · D 3D-Ansicht"
        t = self.font_small.render(hint, True, TEXT_DIM)
        screen.blit(t, (self.win_w / 2 - t.get_width() / 2, y + 68))

    def draw_toasts(self):
        for k, (msg, life) in enumerate(reversed(self.toasts[-4:])):
            alpha = min(1.0, life / 0.4) if life < 0.4 else 1.0
            t = self.font.render(msg, True, TEXT)
            w = t.get_width() + 28
            rect = pygame.Rect(0, 0, w, 34)
            rect.centerx = LEFT_W + (self.win_w - LEFT_W - RIGHT_W) / 2
            rect.bottom = self.canvas_rect.bottom - 16 - k * 42
            surf = pygame.Surface(rect.size, pygame.SRCALPHA)
            a = int(220 * alpha)
            pygame.draw.rect(surf, (*PANEL_BG_LIGHT, a), surf.get_rect(), border_radius=8)
            pygame.draw.rect(surf, (*ACCENT, a), surf.get_rect(), width=1, border_radius=8)
            self.screen.blit(surf, rect.topleft)
            self.screen.blit(t, (rect.x + 14, rect.y + 8))


def show_3d_snapshot(path):
    """Zeigt einen per export_3d_snapshot() gespeicherten Feld-Schnappschuss
    als interaktiven 3D-Oberflaechenplot (x, y, z = Feldwert). Laeuft als
    eigener Prozess (siehe export_3d_snapshot) mit eigenem interaktivem
    Backend -- muss VOR jedem "import matplotlib.pyplot" in diesem Prozess
    gesetzt werden, deshalb hier lokal und nicht am Modulanfang."""
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

    with np.load(path, allow_pickle=True) as data:
        z = data["Z"]
        d_cm = float(data["d_cm"])
        t_ns = float(data["t_ns"])
        label = str(data["label"])
    try:
        os.remove(path)
    except OSError:
        pass

    nx, ny = z.shape
    x = np.arange(nx) * d_cm
    y = np.arange(ny) * d_cm
    xx, yy = np.meshgrid(x, y, indexing="ij")

    fig = plt.figure(figsize=(9, 7), facecolor="#12141a")
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("#12141a")
    surf = ax.plot_surface(xx, yy, z, cmap="viridis", linewidth=0, antialiased=True,
                            rcount=min(nx, 200), ccount=min(ny, 200))
    ax.set_xlabel("x (cm)")
    ax.set_ylabel("y (cm)")
    ax.set_zlabel(label)
    ax.set_title(f"{label}  bei  t = {t_ns:.2f} ns")
    fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.08, label=label)
    fig.tight_layout()
    plt.show()


def main():
    App().run()


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--3d-view":
        show_3d_snapshot(sys.argv[2])
    else:
        main()
