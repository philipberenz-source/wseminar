import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from fdtd_core import FDTD2D, C0
from testfall_a import ABB, D, FARBEN, linienquelle
from testfall_c import R0_WERTE, DICKEN, M, referenz, reflexion, bei

S = 0.99
FC = C0 / (15 * D)
NPML = 14
WINKEL = (0.0, 15.0, 30.0, 45.0, 60.0)
VORGABEN = (1e-3, 1e-6)
ABSTAND = 95
J_SONDE = 300
RAND = 60
DATEI = Path(__file__).resolve().parent / "testfall_c_winkel.npz"


def fenster(theta):
    # Laufweg bis die PML-Reflexion vollstaendig an der Sonde vorbei ist
    return (2 * (ABSTAND + 5) - ABSTAND + 2 * NPML) / np.cos(theta) + 50


def lauf(breite, theta, r0, faktor=1.2):
    w = fenster(theta)
    # Hoehe oberhalb der Sonde so, dass die Welle vom oberen Ende der
    # Quelllinie die Sonde erst nach der Messdauer erreicht
    oben = (w - ABSTAND * np.sin(theta) * np.tan(theta)) / (1 - np.sin(theta))
    oben = max(300, int(faktor * oben) + RAND)
    nx = breite + 2 * NPML
    ny = J_SONDE + oben + 2 * NPML
    sim = FDTD2D(nx, ny, D, courant=S, npml=NPML, pml_m=M, r0=r0)
    i_quelle = NPML + 100
    linienquelle(sim, i_quelle, FC, theta, rand=RAND)
    j = NPML + J_SONDE
    sim.add_probe(i_quelle + ABSTAND, j)
    verzug = (ny - 1 - j - ABSTAND * np.tan(theta)) * np.sin(theta)
    sim.run(int(round((verzug + w) * D / (C0 * sim.dt))))
    return sim


def reflexion_winkel(theta, r0, faktor=1.2):
    f, A = lauf(600, theta, r0, faktor).get_spectrum(0, 8)
    m = (f > 0) & (np.abs(A) > 0.02 * np.abs(A).max())
    nl = C0 / (f[m] * D)
    B = lauf(200, theta, r0, faktor).get_spectrum(0, 8)[1]
    return bei(nl, np.abs(B - A)[m] / np.abs(A[m]), 10)


def messen():
    R = np.empty((len(VORGABEN), len(WINKEL)))
    for k, grad in enumerate(WINKEL):
        for i, r0 in enumerate(VORGABEN):
            R[i, k] = reflexion_winkel(np.deg2rad(grad), r0)
        print(f"{grad:.0f} Grad: " + ", ".join(f"{r:.1e}" for r in R[:, k]))
    np.savez(DATEI, R=R)
    return R


def kontrolle():
    # gleiche Messung mit groesserer Gebietshoehe; Werte duerfen sich nicht aendern
    R = laden()
    for grad in (45.0, 60.0):
        k = WINKEL.index(grad)
        for i, r0 in enumerate(VORGABEN):
            r = reflexion_winkel(np.deg2rad(grad), r0, faktor=1.6)
            print(f"{grad:.0f} Grad, r0 = {r0:.0e}: {R[i, k]:.2e} (Hoehe x1,2), {r:.2e} (x1,6)")


def laden():
    if DATEI.exists():
        return np.load(DATEI)["R"]
    return messen()


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    winkel = laden()
    ref = referenz()
    fig, (links, rechts) = plt.subplots(1, 2, figsize=(9.5, 3.8))

    for npml, farbe in zip(DICKEN, FARBEN):
        werte = np.array([bei(*reflexion(ref, npml, r0), 10) for r0 in R0_WERTE])
        links.loglog(R0_WERTE, werte, color=farbe, lw=1.6, label=f"{npml} PML-Zellen")
    links.loglog(R0_WERTE, R0_WERTE, ":", color="k", lw=1.1, label=r"$|R_\mathrm{Sim}| = |R(0)|$")
    links.axvline(1e-6, color="0.5", lw=0.8, ls="--")
    links.text(1.3e-6, 2.2e-1, "Standardwert", fontsize=8, color="0.35", ha="right")
    links.set_xlim(1e-1, 1e-20)
    links.set_ylim(1e-8, 1)
    links.set_xlabel(r"vorgegebene Restreflexion $|R(0)|$")
    links.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    links.grid(True, which="major", alpha=0.3)
    links.legend(fontsize=7, loc="lower right", ncol=2)

    th = np.linspace(0, 62, 200)
    for r0, werte, farbe in zip(VORGABEN, winkel, ("tab:red", "tab:purple")):
        rechts.semilogy(th, r0 ** np.cos(np.deg2rad(th)), "--", color=farbe, lw=1.1)
        rechts.semilogy(WINKEL, werte, "o-", color=farbe, ms=4, lw=1.4,
                        label=rf"$|R(0)| = 10^{{{int(np.log10(r0))}}}$")
    rechts.plot([], [], "--", color="0.4", lw=1.1, label=r"$|R(0)|^{\cos\theta}$")
    rechts.set_xlim(0, 62)
    rechts.set_ylim(1e-7, 1)
    rechts.set_xticks(WINKEL)
    rechts.set_xlabel(r"Einfallswinkel $\theta$ in Grad")
    rechts.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    rechts.grid(True, which="major", alpha=0.3)
    rechts.legend(fontsize=7.5, loc="upper left")

    fig.tight_layout()
    fig.savefig(ABB / "testfall_c_winkel.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
