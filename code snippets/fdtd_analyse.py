"""
fdtd_analyse.py -- Untersuchung der Fehlerquellen des Verfahrens (Kap. 5.6)
================================================================================
Drei voneinander unabhaengige Funktionen, die jeweils eine der in der Arbeit
diskutierten Grenzen des FDTD-Verfahrens quantifizieren:

    stabilitaet()      Kap. 5.6.2  -- CFL-Bedingung
    dispersion()       Kap. 5.6.1  -- numerische Dispersion
    pml_reflexion()    Kap. 5.6.4  -- Restreflexion des absorbierenden Randes

Aufruf:  python fdtd_analyse.py          # alle drei, Tabellen + analyse.png
"""

from __future__ import annotations

import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import C0, FDTD2D, Source, grid_step


F0 = 2.4e9


# ------------------------------------------------------------------------------
# 1) Stabilitaet: die CFL-Bedingung
# ------------------------------------------------------------------------------
def stabilitaet(S_werte=(0.35, 0.60, 0.70, 0.72, 0.80), n_steps=600):
    """Verfolgt max|E_z| ueber die Zeit fuer verschiedene Courant-Zahlen
    S = c_0*dt/dx.

    Erwartung: unterhalb der CFL-Grenze S_max = 1/sqrt(2) ~ 0,7071 bleibt die
    Amplitude beschraenkt, oberhalb waechst sie exponentiell an. Die
    Instabilitaet ist kein Rundungsproblem, sondern folgt aus der
    von-Neumann-Analyse: fuer S > S_max wird der Verstaerkungsfaktor einzelner
    Fourier-Moden dem Betrage nach groesser als eins. Die Werte 0,70 und 0,72
    liegen bewusst knapp beiderseits der Grenze.
    """
    dx = grid_step(F0, 20)
    N = 90
    verlauf = {}
    for S in S_werte:
        s = FDTD2D(np.ones((N, N)), np.zeros((N, N)), dx, S=S)
        src = Source(kind="ricker", f0=F0)
        v = np.full(n_steps, np.nan)
        # Im instabilen Fall laufen die Feldwerte ueber den Zahlenbereich hinaus;
        # der Overflow ist hier das erwartete Ergebnis und keine Fehlfunktion.
        with np.errstate(over="ignore", invalid="ignore"):
            for k in range(n_steps):
                s.step(src, N // 2, N // 2)
                v[k] = np.max(np.abs(s.Ez))
                if v[k] > 1e12:              # divergiert eindeutig, abbrechen
                    break
        verlauf[S] = v
        letzt = v[np.isfinite(v)][-1]
        print(f"  S = {S:4.2f}:  max|E_z| = {letzt:10.3e} nach {k + 1:4d} Schritten"
              f"{'   -> divergiert' if letzt > 1e3 else '   -> stabil'}")
    return verlauf


# ------------------------------------------------------------------------------
# 2) Numerische Dispersion
# ------------------------------------------------------------------------------
def dispersion_theorie(n_lambda: float, S: float) -> float:
    """Analytische Abweichung der numerischen Phasengeschwindigkeit in Prozent.

    Aus der Dispersionsrelation des 2D-Yee-Schemas folgt fuer Ausbreitung
    entlang einer Gitterachse (k_y = 0)

        sin(k_num*dx/2) = (1/S) * sin(pi*S/N_lambda),   S = c_0*dt/dx,

    also v_p = omega/k_num < c_0. Der Fehler waechst mit grober Aufloesung.
    """
    k_num_dx = 2.0 * math.asin(min(math.sin(math.pi * S / n_lambda) / S, 1.0))
    return 100.0 * ((2.0 * math.pi / n_lambda) / k_num_dx - 1.0)


def _messe_k(n_lambda: float, S: float, Nx=400, npml=12) -> float:
    """Misst die Wellenzahl einer ebenen Welle im Vakuum per Lock-in.

    Durch Nullsetzen von H_x wird d/dy = 0 erzwungen; die Rechnung ist damit
    effektiv eindimensional und beschreibt eine unendlich ausgedehnte ebene
    Welle ohne Beugungsartefakte einer endlichen Linienquelle.
    """
    Ny, dx = 6, grid_step(F0, n_lambda)
    s = FDTD2D(np.ones((Nx, Ny)), np.zeros((Nx, Ny)), dx, S=S, boundary="pml", npml=npml)
    src = Source(kind="sinus", f0=F0, ramp_cycles=4.0)
    i_a, i_b = 120, 340
    n_per = int(round(1.0 / (F0 * s.dt)))

    def schritt():
        s.update_h()
        s.Hx[:] = 0.0                       # y-Invarianz erzwingen
        s.update_e()
        s.inject(src, 40)
        s.t += s.dt

    for _ in range(40 * n_per):             # einschwingen
        schritt()
    co = np.zeros(i_b - i_a); si = np.zeros(i_b - i_a)
    for _ in range(20 * n_per):             # Phase entlang x bestimmen
        schritt()
        w = 2.0 * math.pi * F0 * s.t
        co += s.Ez[i_a:i_b, Ny // 2] * math.cos(w)
        si += s.Ez[i_a:i_b, Ny // 2] * math.sin(w)
    phase = np.unwrap(np.arctan2(si, co))
    return abs(np.polyfit(np.arange(i_b - i_a) * dx, phase, 1)[0])


def dispersion(n_lambdas=(8, 10, 15, 20, 30, 40), S=0.5):
    """Vergleicht gemessene und theoretische Phasengeschwindigkeitsfehler."""
    sim, theo = [], []
    for nl in n_lambdas:
        v = 2.0 * math.pi * F0 / _messe_k(nl, S)
        sim.append(100.0 * (v - C0) / C0)
        theo.append(dispersion_theorie(nl, S))
        print(f"  N_lambda = {nl:3.0f}:  simuliert {sim[-1]:+7.4f} %   "
              f"Theorie {theo[-1]:+7.4f} %   Differenz {abs(sim[-1]-theo[-1]):.4f} %-Pkt.")
    return np.array(n_lambdas, float), np.array(sim), np.array(theo)


# ------------------------------------------------------------------------------
# 3) Restreflexion der PML
# ------------------------------------------------------------------------------
def pml_reflexion(npml_werte=(4, 8, 12, 16)):
    """Misst die vom Rand zurueckgeworfene Feldamplitude.

    Eine 2D-Punktquelle erzeugt hinter der eigentlichen Wellenfront einen
    physikalischen Nachlauf (die 2D-Greensfunktion klingt nur wie
    1/sqrt(t^2 - r^2/c^2) ab). Ein blosses Zeitfenster wuerde diesen Nachlauf
    faelschlich als Reflexion zaehlen. Verglichen wird deshalb gegen eine
    Rechnung in einem so grossen Gebiet, dass dessen Raender im
    Beobachtungsfenster noch nicht wirken koennen.
    """
    dx = grid_step(F0, 20)
    d_probe = 30

    def spur(N, npml, n_steps):
        s = FDTD2D(np.ones((N, N)), np.zeros((N, N)), dx, boundary="pml", npml=npml)
        src = Source(kind="ricker", f0=F0)
        c = N // 2
        return np.array([(s.step(src, c, c), s.Ez[c + d_probe, c])[1]
                         for _ in range(n_steps)])

    N_test, N_ref = 140, 380
    # Zeitschritte, die die Front fuer eine Zelle braucht: dx/(c_0*dt) = 1/S
    n_steps = int((N_ref // 2 - 25) / 0.5)
    ref = spur(N_ref, 16, n_steps)
    werte = []
    for npml in npml_werte:
        tst = spur(N_test, npml, n_steps)
        werte.append(100.0 * np.max(np.abs(tst - ref)) / np.max(np.abs(ref)))
        print(f"  {npml:2d} Zellen PML:  Restreflexion {werte[-1]:.4f} %")
    return np.array(npml_werte, float), np.array(werte)


# ------------------------------------------------------------------------------
def main(out="analyse.png"):
    print("1) Stabilitaet (CFL-Bedingung)")
    verlauf = stabilitaet()
    print("\n2) Numerische Dispersion")
    nl, sim, theo = dispersion()
    print("\n3) Restreflexion der PML")
    npml, refl = pml_reflexion()

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for S, v in verlauf.items():
        ax[0].semilogy(np.maximum(v, 1e-12), label=f"S = {S:g}")
    ax[0].set_ylim(1e-4, 1e13)      # feste Grenzen: die abgebrochenen Kurven
    ax[0].set_xlim(0, 300)          # wuerden die Autoskalierung verzerren
    ax[0].set_xlabel("Zeitschritt"); ax[0].set_ylabel(r"max$|E_z|$ / (V/m)")
    ax[0].set_title("Stabilitaet: CFL-Bedingung"); ax[0].legend(fontsize=8)

    ax[1].plot(nl, sim, "o-", label="simuliert")
    ax[1].plot(nl, theo, "s--", label="Theorie (Yee)")
    ax[1].set_xlabel(r"$N_\lambda = \lambda/\Delta x$")
    ax[1].set_ylabel(r"$(v_{sim}-c_0)/c_0$ in %")
    ax[1].set_title("Numerische Dispersion"); ax[1].legend(fontsize=8)

    ax[2].semilogy(npml, refl, "o-")
    ax[2].set_xlabel("PML-Dicke in Zellen"); ax[2].set_ylabel("Restreflexion in %")
    ax[2].set_title("Absorbierender Rand")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"\ngeschrieben: {out}")


if __name__ == "__main__":
    main()
