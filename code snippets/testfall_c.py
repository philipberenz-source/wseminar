import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import FDTD2D, C0
from testfall_a import ABB, D, FARBEN, ricker, schritte

# Punktquelle J_z (Ricker), Sonde vor der rechten PML. Lauf A: Gebiet, dessen rechter Rand
# waehrend der Messdauer nicht erreicht wird (Referenz), Lauf B: Gebiet mit PML.
# Gemessen wird |E_B - E_A| / |E_A| im Frequenzbereich (wie Testfall 3), ausgewertet bei N_lambda = 10.
# Die reflektierte Zylinderwelle legt einen laengeren Weg zurueck als die einfallende; ihre
# Amplitude faellt wie r^(-1/2), deshalb wird mit sqrt(r_R / r_E) auf gleiche Weglaenge korrigiert.
L_KLEIN, L_GROSS = 200, 600          # Innengebiet in Zellen (ohne PML), senkrechter Einfall
ABSTAND = 95                         # Quelle -> Sonde (senkrechter Einfall)
LAUF = 210                           # Laufweg in Zellen; Reflexionen der anderen Raender kommen spaeter
FC = C0 / (15 * D)                   # Ricker-Spitze bei N_lambda = 15
S = 0.99
M = 3
N_BEZUG = 10                         # Auswertung bei N_lambda = 10 wie in den anderen Testfaellen

DICKEN = ((5, FARBEN[0]), (10, FARBEN[1]), (20, FARBEN[2]))
R0_WERTE = np.logspace(-1, -20, 39)

# schraeger Einfall
DX_W = 60                            # Quelle -> Sondenspalte in x
WINKEL = np.array([0, 10, 20, 30, 40, 50, 55, 60, 65, 70, 75])
FAELLE_W = ((20, 1e-2, "tab:purple", r"20 Zellen, $|R(0)| = 10^{-2}$"),
            (14, 1e-6, "tab:red", r"14 Zellen, $|R(0)| = 10^{-6}$ (Standard)"))


def lauf(innen, npml, r0, innen_y=None, quelle=None, sonden=((ABSTAND, 0),), laufweg=LAUF):
    nx, ny = innen + 2 * npml, (innen if innen_y is None else innen_y) + 2 * npml
    if quelle is None:
        xs, ys = nx // 2, ny // 2
    else:
        xs, ys = quelle[0] + npml, quelle[1] + npml
    sim = FDTD2D(nx, ny, D, courant=S, npml=npml, pml_m=M, r0=r0)
    sim.add_current_source(xs, ys, ricker(FC))
    for dx, dy in sonden:
        sim.add_probe(xs + dx, ys + dy)
    sim.run(schritte(sim, laufweg))
    return sim


def restreflexion(ref, sim, k=0, pad=8):
    """|R_Sim| ueber N_lambda aus Referenzlauf A und Lauf B (Sonde k)."""
    f, A = ref.get_spectrum(k, pad)
    B = sim.get_spectrum(k, pad)[1]
    m = (f > 0) & (np.abs(A) > 0.02 * np.abs(A).max())
    return C0 / (f[m] * D), np.abs(B - A)[m] / np.abs(A[m])


def bei(nl, r, n):
    return np.interp(n, nl[::-1], r[::-1])


def r0_text(r0):
    return rf"$|R(0)| = 10^{{{int(round(np.log10(r0)))}}}$"


def senkrecht(ref, npml, r0):
    korr = np.sqrt((ABSTAND + 2 * (5 + npml)) / ABSTAND)
    return bei(*restreflexion(ref, lauf(L_KLEIN, npml, r0)), N_BEZUG) * korr


def schraeg(npml, r0):
    """|R_Sim| ueber dem Einfallswinkel; Spiegelpunkt auf dem PEC-Rand hinter der Schicht."""
    a = DX_W + 5 + npml                       # Quelle -> PEC-Rand
    weg = 2 * a - DX_W                        # x-Weg Quelle -> Spiegelpunkt -> Sonde
    dys = [int(round(weg * np.tan(np.radians(w)))) for w in WINKEL]
    winkel = np.degrees(np.arctan(np.array(dys) / weg))
    korr = np.sqrt(np.hypot(weg, dys) / np.hypot(DX_W, dys))
    sonden = [(DX_W, dy) for dy in dys]
    quelle, nx, ny = (300, 120), 300 + DX_W + 5, 720
    sim = lauf(nx, npml, r0, ny, quelle, sonden, 500)
    ref = lauf(nx + 500, npml, r0, ny, quelle, sonden, 500)
    r = np.array([bei(*restreflexion(ref, sim, k), N_BEZUG) for k in range(len(sonden))])
    return winkel, r * korr


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    ref = lauf(L_GROSS, 10, 1e-6)

    # Abbildung: links Vorgabe |R(0)| fuer drei Dicken (N_lambda = 10), rechts Dicke fuer vier Aufloesungen
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for npml, farbe in DICKEN:
        werte = np.array([senkrecht(ref, npml, r0) for r0 in R0_WERTE])
        a1.loglog(R0_WERTE, werte, color=farbe, lw=1.6, label=f"{npml} PML-Zellen")
        k = np.argmin(werte)
        print(f"  {npml} Zellen: Minimum {werte[k]:.1e} bei r0 = {R0_WERTE[k]:.0e}, "
              f"bei 1e-6: {werte[np.argmin(abs(np.log(R0_WERTE / 1e-6)))]:.1e}, bei 1e-20: {werte[-1]:.1e}")
    a1.loglog(R0_WERTE, R0_WERTE, ":", color="k", lw=1.1, label=r"$|R_\mathrm{Sim}| = |R(0)|$")
    a1.axvline(1e-6, color="0.5", lw=0.8, ls="--")
    a1.text(1.3e-6, 2.2e-1, "Standardwert", fontsize=8, color="0.35", ha="right")
    a1.set_xlim(1e-1, 1e-20)
    a1.set_ylim(1e-8, 1)
    a1.set_xlabel(r"vorgegebene Restreflexion $|R(0)|$")
    a1.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    a1.grid(True, which="major", alpha=0.3)
    a1.legend(fontsize=7, loc="lower right", ncol=2)

    # rechts: kleinste Reflexion je Dicke (beste Vorgabe) fuer vier Aufloesungen, Daten aus testfall_c_karten.py
    from testfall_c_karten import laden, N_PML, NL
    best = laden().min(axis=1)
    # eigene, einfarbige Abstufung: N_lambda ist eine Groesse, keine Kategorie wie die Dicken links
    for nl, farbe in zip((7, 10, 20, 30), plt.cm.Purples(np.linspace(0.45, 1.0, 4))):
        werte = np.array([np.interp(nl, NL, best[i]) for i in range(len(N_PML))])
        a2.semilogy(N_PML, werte, "o-", ms=3.5, lw=1.5, color=farbe, label=rf"$N_\lambda = {nl}$")
        print(f"  N = {nl}: beste Reflexion bei 5, 10, 30 Zellen: "
              + ", ".join(f"{werte[list(N_PML).index(n)]:.1e}" for n in (5, 10, 30)))
    a2.set_xlim(0, 31)
    a2.set_ylim(1e-7, 1)
    a2.grid(True, which="major", alpha=0.3)
    a2.set_xlabel("Dicke der PML in Zellen")
    a2.set_ylabel(r"kleinste gemessene Reflexion $|R_\mathrm{Sim}|$")
    a2.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_c_restreflexion.pdf")
    plt.close(fig)

    # Schraeger Einfall (nur Zahlen fuer den Text)
    for npml, r0, _, _ in FAELLE_W:
        w, r = schraeg(npml, r0)
        print(f"  schraeg {npml} Zellen, r0 = {r0:.0e}: "
              + ", ".join(f"{a:.0f}:{b:.1e} ({r0 ** np.cos(np.radians(a)):.1e})" for a, b in zip(w, r)))


if __name__ == "__main__":
    main()
    print(f"Abbildungen in {ABB}")
