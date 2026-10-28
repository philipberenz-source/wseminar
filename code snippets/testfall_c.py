import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import FDTD2D, C0
from testfall_a import ABB, D, FARBEN, linienquelle, schritte

BREITE_KLEIN = 200
BREITE_GROSS = 600
HOEHE = 500
ABSTAND = 95
LAUFWEG = 210
FC = C0 / (15 * D)
S = 0.99
M = 3
R0_WERTE = np.logspace(-1, -20, 20)
DICKEN = (5, 10, 20)


def lauf(breite, npml, r0):
    nx = breite + 2 * npml
    ny = HOEHE + 2 * npml
    sim = FDTD2D(nx, ny, D, courant=S, npml=npml, pml_m=M, r0=r0)
    linienquelle(sim, nx // 2, FC)
    sim.add_probe(nx // 2 + ABSTAND, ny // 2)
    sim.run(schritte(sim, LAUFWEG))
    return sim


def referenz():
    return lauf(BREITE_GROSS, 10, 1e-6)


def reflexion(ref, npml, r0):
    sim = lauf(BREITE_KLEIN, npml, r0)
    f, A = ref.get_spectrum(0, 8)
    B = sim.get_spectrum(0, 8)[1]
    m = (f > 0) & (np.abs(A) > 0.02 * np.abs(A).max())
    return C0 / (f[m] * D), np.abs(B - A)[m] / np.abs(A[m])


def bei(nl, r, wert):
    return np.interp(wert, nl[::-1], r[::-1])


def main():
    from testfall_c_karten import laden, N_PML, NL

    ABB.mkdir(parents=True, exist_ok=True)
    ref = referenz()
    fig, (links, rechts) = plt.subplots(1, 2, figsize=(9.5, 3.8))

    for npml, farbe in zip(DICKEN, FARBEN):
        werte = np.array([bei(*reflexion(ref, npml, r0), 10) for r0 in R0_WERTE])
        links.loglog(R0_WERTE, werte, color=farbe, lw=1.6, label=f"{npml} PML-Zellen")
        k = np.argmin(werte)
        print(f"{npml} Zellen: Minimum {werte[k]:.1e} bei r0 = {R0_WERTE[k]:.0e}, "
              f"bei 1e-6: {werte[5]:.1e}, bei 1e-20: {werte[-1]:.1e}")
    links.loglog(R0_WERTE, R0_WERTE, ":", color="k", lw=1.1, label=r"$|R_\mathrm{Sim}| = |R(0)|$")
    links.axvline(1e-6, color="0.5", lw=0.8, ls="--")
    links.text(1.3e-6, 2.2e-1, "Standardwert", fontsize=8, color="0.35", ha="right")
    links.set_xlim(1e-1, 1e-20)
    links.set_ylim(1e-8, 1)
    links.set_xlabel(r"vorgegebene Restreflexion $|R(0)|$")
    links.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    links.grid(True, which="major", alpha=0.3)
    links.legend(fontsize=7, loc="lower right", ncol=2)

    beste = laden().min(axis=1)
    farben = plt.cm.Purples(np.linspace(0.45, 1.0, 4))
    for nl, farbe in zip((7, 10, 20, 30), farben):
        werte = np.array([np.interp(nl, NL, beste[i]) for i in range(len(N_PML))])
        rechts.semilogy(N_PML, werte, "o-", ms=3.5, lw=1.5, color=farbe, label=rf"$N_\lambda = {nl}$")
        print(f"N = {nl}: beste Reflexion bei 5, 10, 30 Zellen: "
              f"{werte[3]:.1e}, {werte[6]:.1e}, {werte[12]:.1e}")
    rechts.set_xlim(0, 31)
    rechts.set_ylim(1e-7, 1)
    rechts.grid(True, which="major", alpha=0.3)
    rechts.set_xlabel("Dicke der PML in Zellen")
    rechts.set_ylabel(r"kleinste gemessene Reflexion $|R_\mathrm{Sim}|$")
    rechts.legend(fontsize=7.5, loc="upper right")

    fig.tight_layout()
    fig.savefig(ABB / "testfall_c_restreflexion.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
