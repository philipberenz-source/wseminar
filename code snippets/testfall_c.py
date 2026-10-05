
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import FDTD2D, C0
from testfall_a import ABB, D, FARBEN, ricker, schritte, log_achse

# Punktquelle in der Mitte, Sonde auf der Achse 5 Zellen vor der rechten PML (senkrechter Einfall).
# Lauf A: grosses Gebiet (Referenz, kein Rand erreicht die Sonde), Lauf B: kleines Gebiet.
# Gemessen wird |E_B - E_A| / |E_A| im Frequenzbereich (wie Testfall 4).
L_KLEIN, L_GROSS = 200, 600          # Innengebiet in Zellen (ohne PML)
ABSTAND = 95                         # Quelle -> Sonde
LAUF = 210                           # Laufweg in Zellen; Reflexionen der anderen Raender kommen spaeter
FC = C0 / (15 * D)                   # Ricker-Spitze bei N_lambda = 15
S = 0.99
N_BEZUG = 10                         # Auswertung bei N_lambda = 10 wie in den anderen Testfaellen
R0_DREI = (1e-3, 1e-6, 1e-16)        # 1e-6 ist der Standardwert des Cores
FARBEN_R0 = ("tab:purple", "tab:red", "tab:brown")
DICKEN = (4, 6, 8, 10, 12, 15, 20, 25, 30)
R0_WERTE = np.logspace(-1, -20, 39)


def lauf(innen, npml, r0):
    n = innen + 2 * npml
    c = n // 2
    sim = FDTD2D(n, n, D, courant=S, npml=npml, r0=r0)
    sim.add_current_source(c, c, ricker(FC))
    sim.add_probe(c + ABSTAND, c)
    sim.run(schritte(sim, LAUF))
    return sim


def restreflexion(ref, sim, pad=8):
    """|R_Sim| ueber N_lambda aus Referenzlauf A und Lauf B."""
    f, A = ref.get_spectrum(0, pad)
    B = sim.get_spectrum(0, pad)[1]
    m = (f > 0) & (np.abs(A) > 0.02 * np.abs(A).max())
    return C0 / (f[m] * D), np.abs(B - A)[m] / np.abs(A[m])


def bei(nl, r, n):
    return np.interp(n, nl[::-1], r[::-1])


def r0_text(r0):
    return rf"$|R(0)| = 10^{{{int(round(np.log10(r0)))}}}$"


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    ref = lauf(L_GROSS, 10, 1e-6)

    # Abbildung 1: links Schichtdicke, rechts Vorgabe r0 (jeweils N_lambda = 10)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for r0, farbe in zip(R0_DREI, FARBEN_R0):
        werte = [bei(*restreflexion(ref, lauf(L_KLEIN, n, r0)), N_BEZUG) for n in DICKEN]
        a1.semilogy(DICKEN, werte, "o-", color=farbe, lw=1.5, ms=3.5, label=r0_text(r0))
        a1.axhline(r0, color=farbe, lw=0.9, ls=":")
        print(f"  r0 = {r0:.0e}: " + ", ".join(f"{n}:{w:.1e}" for n, w in zip(DICKEN, werte)))
    a1.set_xlim(3, 31)
    a1.set_ylim(1e-7, 1e-1)
    a1.set_xlabel("Dicke der PML in Zellen")
    a1.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    a1.grid(True, which="major", alpha=0.3)
    a1.legend(fontsize=7.5, loc="upper right")

    for npml, farbe in ((10, FARBEN[0]), (20, FARBEN[1])):
        werte = np.array([bei(*restreflexion(ref, lauf(L_KLEIN, npml, r0)), N_BEZUG) for r0 in R0_WERTE])
        a2.loglog(R0_WERTE, werte, color=farbe, lw=1.6, label=f"{npml} PML-Zellen")
        k = np.argmin(werte)
        print(f"  npml = {npml}: Minimum {werte[k]:.1e} bei r0 = {R0_WERTE[k]:.0e}, bei 1e-20: {werte[-1]:.1e}")
    a2.loglog(R0_WERTE, R0_WERTE, ":", color="k", lw=1.1, label=r"$|R_\mathrm{Sim}| = |R(0)|$")
    a2.axvline(1e-6, color="0.5", lw=0.8, ls="--")
    a2.text(1.3e-6, 2.2e-1, "Standardwert", fontsize=8, color="0.35", ha="right")
    a2.set_xlim(1e-1, 1e-20)
    a2.set_ylim(1e-7, 5e-1)
    a2.set_xlabel(r"vorgegebene Restreflexion $|R(0)|$")
    a2.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    a2.grid(True, which="major", alpha=0.3)
    a2.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_c_restreflexion.pdf")
    plt.close(fig)

    # Abbildung 2: Konvergenz ueber N_lambda (10 PML-Zellen)
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for r0, farbe in zip(R0_DREI, FARBEN_R0):
        nl, r = restreflexion(ref, lauf(L_KLEIN, 10, r0))
        ax.loglog(nl, r, color=farbe, lw=1.6, label=r0_text(r0))
        print(f"  r0 = {r0:.0e}, 10 Zellen: |R| bei N = 10, 20, 30: "
              + ", ".join(f"{bei(nl, r, n):.1e}" for n in (10, 20, 30)))
    log_achse(ax)
    ax.set_xlim(7, 30)
    ax.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    ax.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_c_konvergenz.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
    print(f"Abbildungen in {ABB}")
