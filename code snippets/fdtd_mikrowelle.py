"""
fdtd_mikrowelle.py -- Anwendung: Erwaermung im Mikrowellenherd (Kap. 6)
================================================================================
Ein Mikrowellenherd ist ein Hohlraumresonator mit ideal leitenden Waenden. Das
Magnetron speist bei 2,45 GHz ein, die Welle wird an den Metallwaenden
vollstaendig reflektiert und ueberlagert sich mit sich selbst zu einer stehenden
Welle mit ortsfesten Baeuchen und Knoten.

Randbedingung: hier wird ausdruecklich KEINE PML verwendet. Der absorbierende
Rand waere gerade falsch -- die stehende Welle entsteht ja erst dadurch, dass
die Waende vollstaendig reflektieren. Der Solver laeuft deshalb mit
boundary="pec" (E_tan = 0 an der Wand); Split-Felder werden gar nicht angelegt.

Physik der Erwaermung
--------------------------------------------------------------------------------
Im verlustbehafteten Dielektrikum fliesst nach dem Ohmschen Gesetz die
Stromdichte J = sigma*E. Die in Waerme umgesetzte Leistungsdichte ist

    p(x,y) = <J*E> = sigma * <E_z^2>            [W/m^3]

Aus ihr folgt mit der volumetrischen Waermekapazitaet rho*c_p unmittelbar die
anfaengliche Aufheizrate

    dT/dt = p / (rho*c_p)                        [K/s]

Welche Stelle am heissesten wird, haengt von der Objektgroesse ab. Im Medium ist
die Wellenlaenge um sqrt(eps_r) verkuerzt (in Wasser von 12,2 cm auf rund
1,4 cm), und die gekruemmte Oberflaeche wirkt wie eine dielektrische Linse:

  * kleine Objekte (bis etwa 5 Wellenlaengen im Medium) buendeln die Welle in
    ihrer Mitte -- das Maximum liegt im Kern
  * grosse Objekte werden von der Eindringtiefe bestimmt -- aussen heiss, innen
    kalt

Warum das Feldmuster unregelmaessig aussieht
--------------------------------------------------------------------------------
Der leere Garraum besitzt diskrete Eigenmoden f_mn = (c/2)*sqrt((m/Lx)^2+(n/Ly)^2).
Bei 30 x 20 cm liegen um 2,45 GHz gleich zwei davon dicht beieinander (TM_23 bei
2,461 GHz und TM_42 bei 2,498 GHz). Eine Anregung bei 2,45 GHz trifft keine
davon exakt, regt aber beide teilweise an -- das Ergebnis ist eine Ueberlagerung
und kein sauberes Schachbrettmuster. Zusaetzlich verzerrt das Objekt mit
eps_r = 78 die Moden erheblich. Beides ist physikalisch korrekt und genau der
Grund, weshalb reale Geraete ungleichmaessig garen und einen Drehteller haben.
Die Anzeige nennt die naechstgelegene Eigenmode samt Verstimmung; stellt man die
Frequenz genau darauf, wird das Muster deutlich regelmaessiger.

Einheiten -- durchgaengig SI
--------------------------------------------------------------------------------
    E_z [V/m]   sigma [S/m]   J_z [A/m^2]   p [W/m^3]   dT/dt [K/s]
Die Regler arbeiten in cm und GHz und rechnen beim Eintritt sofort in m und Hz um.
Die Quelle ist eine eingepraegte Stromdichte J_z mit fest vorgegebenem
Scheitelwert; die Feldstaerken sind damit echte V/m und nicht normiert. Wie viel
Leistung tatsaechlich aufgenommen wird, ist dann ein Ergebnis der Rechnung und
keine Vorgabe -- nur so bleibt der Vergleich verschiedener Speisen ehrlich
(Eis nimmt um Groessenordnungen weniger auf als Wasser). Die 2D-Rechnung liefert
Leistung je Meter Tiefe [W/m]; erst die Annahme einer Geraetetiefe TIEFE macht
daraus absolute Watt.

Grenzen des Modells (Kap. 3.4)
--------------------------------------------------------------------------------
  * 2D-Schnitt statt des dreidimensionalen Garraums
  * linear, isotrop, nichtdispersiv -- bei Wasser ist Letzteres falsch
    (Debye-Relaxation); die Materialwerte gelten nur nahe 2,45 GHz
  * keine Rueckkopplung: real steigt sigma mit der Temperatur, warme Stellen
    heizen sich dann noch schneller auf
  * keine Waermeleitung, kein Phasenuebergang, kein Drehteller. Gezeigt wird die
    momentane Aufheizrate, nicht die Temperatur nach Minuten

Aufruf:  python fdtd_mikrowelle.py
"""

from __future__ import annotations

import math
import tkinter as tk
from tkinter import ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.figure import Figure

from fdtd_core import C0, EPS0, MU0, FDTD2D, Source, grid_step, sigma_aus_epsr2

# --- Garraum und Betrieb ------------------------------------------------------
BREITE, HOEHE = 0.30, 0.20   # Innenmasse des Garraums in m
TIEFE = 0.20                 # angenommene Tiefe fuer die Umrechnung nach Watt
J_QUELLE = 43.0e3            # Scheitelwert der Quellstromdichte J_0 in A/m^2.
                             # Fest vorgegeben, NICHT auf eine Sollleistung
                             # nachgeregelt: nur so bleibt der Vergleich der
                             # Materialien aussagekraeftig. Der Wert ist so
                             # gewaehlt, dass die Standardgeometrie mit Wasser
                             # rund 800 W aufnimmt -- die Groessenordnung eines
                             # Haushaltsgeraets. Wie viel ein anderes Material
                             # aufnimmt, ergibt die Rechnung dann von selbst.
N_LAMBDA = 10                # Zellen je Wellenlaenge im dichtesten Medium
EPS_MAX = 80.0               # Auslegung des Gitters auf Wasser
S_COURANT = 0.5              # Courant-Zahl S = c_0*dt/dx (CFL in 2D: <= 0,707)
EINSCHWINGEN = 25            # Perioden bis zum eingeschwungenen Zustand
MESSEN = 20                  # Perioden, ueber die gemittelt wird
SCHRITTE_PRO_BILD = 40
# Kontrastanhebung: dargestellt wird (Wert/Maximum)**GAMMA. Werte > 1 druecken
# die mittleren Helligkeiten nach unten, sodass der helle Bereich um die
# Wellenbaeuche schmaler wird und der Untergrund schnell in Schwarz uebergeht.
GAMMA_FELD = 2.0
GAMMA_WAERME = 1.8

# Materialdaten bei 2,45 GHz: (eps_r', eps_r'', rho*c_p in J/(m^3 K)).
# eps_r = eps' - j*eps'' ist die uebliche Angabeform; der Imaginaerteil wird
# ueber sigma = omega*eps_0*eps'' in die Konvention des Kerns umgerechnet.
# Groessenordnungen aus der Literatur zur dielektrischen Erwaermung; sie
# streuen je nach Temperatur, Wasser- und Salzgehalt erheblich.
MATERIAL = {
    "Wasser (25 C)":  (78.0, 9.5, 4.18e6),
    "Kartoffel":      (57.0, 16.0, 3.5e6),
    "Fleisch, mager": (52.0, 17.0, 3.6e6),
    "Butter / Fett":  (4.5, 0.6, 1.9e6),
    "Eis (-10 C)":    (3.2, 0.003, 1.9e6),
}

# --- Farbschema: gruen auf schwarz -------------------------------------------
# Feld: Knoten schwarz, Baeuche hell. Die Skala ist symmetrisch, das Vorzeichen
# also nicht ablesbar -- fuer eine stehende Welle zaehlt allein der Betrag.
CMAP_FELD = LinearSegmentedColormap.from_list(
    "feld", ["#eaffef", "#5cf09a", "#12b866", "#053d21", "#000000",
             "#053d21", "#12b866", "#5cf09a", "#eaffef"])
CMAP_WAERME = LinearSegmentedColormap.from_list(
    "waerme", ["#000000", "#022a16", "#0b7a3c", "#3ee584", "#eaffef"])
BG, FG = "#0b0f0c", "#b6ffd4"


class Mikrowelle:
    """Hohlraumresonator mit Speise, Objekt und Auswertung der Verlustleistung."""

    def __init__(self, material="Wasser (25 C)", f_ghz=2.45, radius=0.04,
                 x_obj=0.15, y_obj=0.10):
        self.f0 = f_ghz * 1e9
        self.lam = C0 / self.f0
        self.material, self.radius = material, radius
        eps_r1, eps_r2, self.rho_cp = MATERIAL[material]

        self.dx = grid_step(2.45e9, N_LAMBDA, EPS_MAX)   # festes Gitter
        self.Nx = int(round(BREITE / self.dx)) + 1
        self.Ny = int(round(HOEHE / self.dx)) + 1

        # Objekt: Kreisquerschnitt (in 2D ein unendlich langer Zylinder)
        ii, jj = np.meshgrid(np.arange(self.Nx), np.arange(self.Ny), indexing="ij")
        self.maske = ((ii - x_obj / self.dx) ** 2 + (jj - y_obj / self.dx) ** 2
                      <= (radius / self.dx) ** 2)

        eps_r = np.ones((self.Nx, self.Ny))
        sigma = np.zeros((self.Nx, self.Ny))
        eps_r[self.maske] = eps_r1
        self.sigma_obj = sigma_aus_epsr2(eps_r2, self.f0)
        sigma[self.maske] = self.sigma_obj

        # Ideal leitende Waende, keine PML -- siehe Modulkopf
        self.s = FDTD2D(eps_r, sigma, self.dx, S=S_COURANT, boundary="pec")

        # Magnetron: kurze Apertur eine Viertelwellenlaenge vor der linken Wand.
        # Direkt an der Wand erzwingt E_tan = 0 einen Knoten, dort wuerde kaum
        # Leistung einkoppeln; bei lambda/4 liegt ein Bauch. Aussermittig in y,
        # damit nicht nur symmetrische Moden angeregt werden.
        self.src = Source(kind="sinus", f0=self.f0, amplitude=J_QUELLE,
                          ramp_cycles=5.0)
        self.i_src = int(round(0.25 * self.lam / self.dx))
        j0 = int(round(0.70 * self.Ny))
        self.ap = slice(j0 - 6, j0 + 7)

        self.n_per = max(int(round(1.0 / (self.f0 * self.s.dt))), 1)
        self.n_settle = EINSCHWINGEN * self.n_per
        self.n_mess = MESSEN * self.n_per
        self.acc = np.zeros((self.Nx, self.Ny))    # Summe von E_z^2
        self.n_acc = 0

    # -- Simulation -----------------------------------------------------------
    def step(self) -> None:
        self.s.update_h()
        self.s.update_e()
        self.s._deposit((self.i_src, self.ap), self.src(self.s.t))   # J_z in A/m^2
        self.s.t += self.s.dt
        self.s.n += 1
        if self.s.n > self.n_settle and self.n_acc < self.n_mess:
            self.acc += self.s.Ez ** 2
            self.n_acc += 1

    @property
    def phase(self) -> str:
        if self.s.n <= self.n_settle:
            return "einschwingen"
        return "messen" if self.n_acc < self.n_mess else "fertig"

    @property
    def fortschritt(self) -> float:
        return min((self.s.n) / (self.n_settle + self.n_mess), 1.0)

    # -- Auswertung -----------------------------------------------------------
    def leistungsdichte(self) -> np.ndarray:
        """Absorbierte Leistungsdichte p = sigma*<E_z^2> in W/m^3."""
        if self.n_acc == 0:
            return np.zeros_like(self.acc)
        return self.s.sigma * (self.acc / self.n_acc)

    def leistung(self) -> float:
        """Insgesamt aufgenommene Leistung in W.

        Die 2D-Rechnung liefert zunaechst Leistung je Meter Tiefe [W/m]; erst
        die Annahme einer Geraetetiefe macht daraus absolute Watt.
        """
        return float(self.leistungsdichte().sum() * self.dx ** 2 * TIEFE)

    def aufheizrate(self) -> np.ndarray:
        """Anfaengliche Temperaturaenderung dT/dt = p/(rho*c_p) in K/s."""
        return self.leistungsdichte() / self.rho_cp

    def eindringtiefe(self) -> float:
        """Leistungs-Eindringtiefe D_p = 1/(2*alpha) in m.

        Die Amplitude faellt wie exp(-alpha*x), die Leistung wie exp(-2*alpha*x).
        Da die Heatmap eine Leistungsdichte zeigt, ist D_p die passende Groesse;
        sie ist die in der Mikrowellentechnik uebliche Angabe.
        """
        w = 2.0 * math.pi * self.f0
        eps = EPS0 * MATERIAL[self.material][0]
        p = self.sigma_obj / (w * eps)
        alpha = w * math.sqrt(MU0 * eps / 2.0) * math.sqrt(math.sqrt(1 + p * p) - 1)
        return float("inf") if alpha <= 0 else 0.5 / alpha

    def naechste_mode(self):
        """Naechstgelegene Eigenmode des leeren Garraums und Verstimmung."""
        best = min(((0.5 * C0 * math.hypot(m / BREITE, n / HOEHE), m, n)
                    for m in range(1, 9) for n in range(1, 9)),
                   key=lambda t: abs(t[0] - self.f0))
        return best[0], best[1], best[2], 100.0 * (self.f0 - best[0]) / best[0]

    def kennzahlen(self) -> dict:
        """Gleichmaessigkeit der Erwaermung innerhalb des Objekts."""
        r = self.aufheizrate()[self.maske]
        if r.size == 0 or r.max() <= 0:
            return {"kalt": 0.0, "gleich": 0.0, "rmax": 0.0, "rmit": 0.0}
        return {"kalt": 100.0 * float((r < 0.2 * r.max()).mean()),
                "gleich": 100.0 * float(r.mean() / r.max()),
                "rmax": float(r.max()), "rmit": float(r.mean())}


# ------------------------------------------------------------------------------
class App:
    """Tkinter-Oberflaeche: erst laeuft die Simulation, am Ende die Heatmap."""

    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("FDTD - Erwaermung im Mikrowellenherd")
        root.configure(bg=BG)
        self.mw = Mikrowelle()
        self.laeuft = True

        stil = ttk.Style()
        stil.theme_use("clam")
        stil.configure("D.TLabel", background=BG, foreground=FG)
        stil.configure("D.TFrame", background=BG)
        stil.configure("D.TButton", background="#123", foreground=FG)

        panel = ttk.Frame(root, style="D.TFrame", padding=10)
        panel.grid(row=0, column=0, sticky="ns")

        def titel(txt, r):
            ttk.Label(panel, text=txt, style="D.TLabel").grid(
                row=r, column=0, sticky="w", pady=(10, 2))

        titel("Speise", 0)
        self.v_mat = tk.StringVar(value="Wasser (25 C)")
        box = ttk.Combobox(panel, textvariable=self.v_mat, state="readonly",
                           values=list(MATERIAL), width=18)
        box.grid(row=1, column=0, sticky="we")
        box.bind("<<ComboboxSelected>>", lambda e: self.neu())

        self.regler = {}
        for r, (name, txt, lo, hi, init, res) in enumerate([
                ("f", "Frequenz / GHz", 2.0, 3.0, 2.45, 0.005),
                ("radius", "Radius / cm", 1.0, 7.0, 4.0, 0.05),
                ("x", "Position x / cm", 6.0, 24.0, 15.0, 0.05),
                ("y", "Position y / cm", 4.0, 16.0, 10.0, 0.05)], start=1):
            titel(txt, 2 * r)
            v = tk.DoubleVar(value=init)
            sc = tk.Scale(panel, from_=lo, to=hi, resolution=res, variable=v,
                          orient="horizontal", length=190, bg=BG, fg=FG,
                          troughcolor="#123", highlightthickness=0, bd=0)
            sc.grid(row=2 * r + 1, column=0)
            # erst beim Loslassen neu rechnen, sonst startet die Simulation
            # waehrend des Ziehens bei jedem Pixel neu
            sc.bind("<ButtonRelease-1>", lambda e: self.neu())
            self.regler[name] = v

        knoepfe = ttk.Frame(panel, style="D.TFrame")
        knoepfe.grid(row=12, column=0, pady=12)
        ttk.Button(knoepfe, text="Start / Pause", style="D.TButton",
                   command=self.pause).grid(row=0, column=0, padx=2)
        ttk.Button(knoepfe, text="Neu starten", style="D.TButton",
                   command=self.neu).grid(row=0, column=1, padx=2)

        self.info = ttk.Label(panel, style="D.TLabel", justify="left")
        self.info.grid(row=13, column=0, sticky="w")
        ttk.Label(panel, style="D.TLabel", justify="left", wraplength=205,
                  text="Modell: 2D-Schnitt, lineares, isotropes und nicht "
                       "dispersives Medium. Ohne Waermeleitung, ohne "
                       "Temperatur-Rueckkopplung auf sigma, ohne Drehteller. "
                       "Gezeigt wird die anfaengliche Aufheizrate, nicht die "
                       "Temperatur nach Minuten."
                  ).grid(row=14, column=0, sticky="w", pady=(14, 0))

        # -- Grafik ------------------------------------------------------------
        self.fig = Figure(figsize=(10.4, 4.6), facecolor=BG)
        self.fig.subplots_adjust(left=0.06, right=0.88, bottom=0.14, top=0.90,
                                 wspace=0.22)
        self.ax1 = self.fig.add_subplot(1, 2, 1)
        self.ax2 = self.fig.add_subplot(1, 2, 2)
        ext = (0, BREITE * 100, 0, HOEHE * 100)
        for ax in (self.ax1, self.ax2):
            ax.set_facecolor("black")
            ax.set_xlabel("x / cm", color=FG); ax.set_ylabel("y / cm", color=FG)
            ax.tick_params(colors=FG, labelsize=8)
            for sp in ax.spines.values():
                sp.set_color(FG)
        self.im1 = self.ax1.imshow(np.zeros((5, 5)), origin="lower", extent=ext,
                                   cmap=CMAP_FELD, vmin=-1, vmax=1,
                                   interpolation="bilinear")
        self.im2 = self.ax2.imshow(np.zeros((5, 5)), origin="lower", extent=ext,
                                   cmap=CMAP_WAERME, interpolation="bilinear",
                                   norm=PowerNorm(1.0 / GAMMA_WAERME, 0, 1))
        cax = self.fig.add_axes((0.895, 0.14, 0.016, 0.76))
        self.cbar = self.fig.colorbar(self.im2, cax=cax, label="dT/dt in K/s")
        self.cbar.ax.tick_params(colors=FG, labelsize=7)
        self.cbar.ax.yaxis.label.set_color(FG)
        self.cbar.outline.set_edgecolor(FG)
        self.t1 = self.ax1.set_title("", color=FG, fontsize=10)
        self.t2 = self.ax2.set_title("", color=FG, fontsize=10)
        self.kontur = []
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=0, column=1, sticky="nsew")
        root.columnconfigure(1, weight=1); root.rowconfigure(0, weight=1)
        self.tick()

    # -- Steuerung ------------------------------------------------------------
    def pause(self):
        self.laeuft = not self.laeuft

    def neu(self):
        self.mw = Mikrowelle(self.v_mat.get(), self.regler["f"].get(),
                             self.regler["radius"].get() / 100,
                             self.regler["x"].get() / 100,
                             self.regler["y"].get() / 100)
        self.im2.set_data(np.zeros((5, 5)))
        self.laeuft = True

    def tick(self):
        mw = self.mw
        if self.laeuft and mw.phase != "fertig":
            for _ in range(SCHRITTE_PRO_BILD):
                mw.step()

        # -- links: laufende Simulation des Feldes ----------------------------
        E = mw.s.Ez
        v = max(float(np.abs(E).max()), 1e-12)
        # Kontrastanhebung, Vorzeichen bleibt erhalten
        self.im1.set_data((np.sign(E) * (np.abs(E) / v) ** GAMMA_FELD).T)
        f_mode, m, n, verst = mw.naechste_mode()
        self.t1.set_text(f"stehende Welle $E_z$   (max {v:.0f} V/m)")

        # -- rechts: Heatmap erst, wenn die Mittelung abgeschlossen ist --------
        if mw.phase == "fertig":
            if self.laeuft:                       # einmalig beim Fertigwerden
                self.laeuft = False
            r = mw.aufheizrate()
            rm = max(float(r.max()), 1e-12)
            self.im2.set_data(r.T); self.im2.set_norm(
                PowerNorm(1.0 / GAMMA_WAERME, 0, rm))
            self.t2.set_text(f"Aufheizrate - {mw.leistung():.0f} W aufgenommen")
        else:
            self.t2.set_text(f"Messung laeuft - {100 * mw.fortschritt:.0f} % "
                             f"({mw.phase})")

        for c in self.kontur:                     # Objektumriss in beide Karten
            c.remove()
        self.kontur = [ax.contour(np.linspace(0, BREITE * 100, mw.Nx),
                                  np.linspace(0, HOEHE * 100, mw.Ny),
                                  mw.maske.T.astype(float), levels=[0.5],
                                  colors=[FG], linewidths=0.9)
                       for ax in (self.ax1, self.ax2)]

        k = mw.kennzahlen()
        self.info.configure(
            text=(f"sigma = {mw.sigma_obj:.3f} S/m\n"
                  f"Eindringtiefe = {mw.eindringtiefe() * 1e3:.1f} mm\n"
                  f"lambda: {mw.lam * 100:.1f} cm Luft, "
                  f"{mw.lam * 100 / math.sqrt(MATERIAL[mw.material][0]):.1f} cm im Medium\n"
                  f"Gitter {mw.Nx} x {mw.Ny}, dx = {mw.dx * 1e3:.2f} mm\n\n"
                  f"naechste Mode: TM_{m}{n}\n"
                  f"  bei {f_mode / 1e9:.3f} GHz ({verst:+.1f} %)\n\n"
                  f"aufgenommen: {mw.leistung():.0f} W\n"
                  f"  max {k['rmax']:.2f} K/s, im Mittel {k['rmit']:.2f} K/s\n"
                  f"  kalte Zonen: {k['kalt']:.0f} % der Flaeche\n"
                  f"  Gleichmaessigkeit: {k['gleich']:.0f} %"))
        self.canvas.draw_idle()
        self.root.after(30, self.tick)


if __name__ == "__main__":
    wurzel = tk.Tk()
    App(wurzel)
    wurzel.mainloop()
