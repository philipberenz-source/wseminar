"""Testfall 3: Reflexion der PML fuer eine ebene Welle.

Gemessen wird |R_Sim| nach dem Zwei-Lauf-Verfahren (kleines Gebiet mit PML
gegen breite Referenz, Taflove Abb. 7.2) in Abhaengigkeit von der Vorgabe
|R(0)|, der Schichtdicke, der Aufloesung N_lambda und dem Einfallswinkel.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from fdtd_core import FDTD2D, C0
from testfall_1_vakuum import ABB, D, FARBEN, linienquelle, schritte

S = 0.99
M = 3
FC = C0 / (15 * D)            # Ricker-Puls, Spektrum etwa N_lambda = 6 ... 100
NL_STANDARD = 20              # Aufloesung fuer (a) und (c)
DICKEN = (5, 10, 20)          # gemessen
DICKEN_ABB = (5, 20)          # in der Abbildung gezeigt
R0_WERTE = np.logspace(-1, -16, 16)
R0_STANDARD = 1e-6            # Standardwert von fdtd_core.py
WINKEL = (0.0, 15.0, 30.0, 45.0, 60.0)
WINKEL_DICKE = 20
WINKEL_VORGABEN = (1e-3, 1e-6)
HOEHE, ABSTAND, LAUFWEG = 500, 95, 210
DATEI = Path(__file__).resolve().parent / "testfall_3_pml.npz"


def lauf_senkrecht(breite, npml, r0, courant=S):
    nx, ny = breite + 2 * npml, HOEHE + 2 * npml
    sim = FDTD2D(nx, ny, D, courant=courant, npml=npml, pml_m=M, r0=r0)
    linienquelle(sim, nx // 2, FC)
    sim.add_probe(nx // 2 + ABSTAND, ny // 2)   # 5 Zellen vor der Schicht
    sim.run(schritte(sim, LAUFWEG))
    return sim


def lauf_schraeg(breite, npml, theta, r0, rand=60, j_sonde=300):
    w = (2 * (ABSTAND + 5) - ABSTAND + 2 * npml) / np.cos(theta) + 50
    oben = (w - ABSTAND * np.sin(theta) * np.tan(theta)) / (1 - np.sin(theta))
    oben = max(300, int(1.2 * oben) + rand)
    nx, ny = breite + 2 * npml, j_sonde + oben + 2 * npml
    sim = FDTD2D(nx, ny, D, courant=S, npml=npml, pml_m=M, r0=r0)
    i_quelle = npml + 100
    linienquelle(sim, i_quelle, FC, theta, rand=rand)
    j = npml + j_sonde
    sim.add_probe(i_quelle + ABSTAND, j)
    verzug = (ny - 1 - j - ABSTAND * np.tan(theta)) * np.sin(theta)
    sim.run(int(round((verzug + w) * D / (C0 * sim.dt))))
    return sim


def reflexion(ref, sim):
    """|R_Sim| ueber N_lambda aus Referenz (A) und Lauf mit PML (B)."""
    f, A = ref.get_spectrum(0, 8)
    B = sim.get_spectrum(0, 8)[1]
    m = (f > 0) & (np.abs(A) > 0.02 * np.abs(A).max())
    nl = C0 / (f[m] * D)
    return nl[::-1], (np.abs(B - A)[m] / np.abs(A[m]))[::-1]


def bei(nl, r, wert):
    return np.interp(wert, nl, r)


def messen():
    ref = lauf_senkrecht(600, 10, 1e-6)
    vorgabe = np.empty((len(DICKEN), len(R0_WERTE)))
    for a, npml in enumerate(DICKEN):
        for b, r0 in enumerate(R0_WERTE):
            vorgabe[a, b] = bei(*reflexion(ref, lauf_senkrecht(200, npml, r0)), NL_STANDARD)
    nl_achse = np.linspace(6, 30, 97)
    aufl = np.array([bei(*reflexion(ref, lauf_senkrecht(200, npml, R0_STANDARD)), nl_achse)
                     for npml in DICKEN])
    winkel = np.empty((len(WINKEL_VORGABEN), len(WINKEL)))
    for k, grad in enumerate(WINKEL):
        th = np.deg2rad(grad)
        for i, r0 in enumerate(WINKEL_VORGABEN):
            A = lauf_schraeg(600, WINKEL_DICKE, th, r0)
            B = lauf_schraeg(200, WINKEL_DICKE, th, r0)
            winkel[i, k] = bei(*reflexion(A, B), NL_STANDARD)
    # Kontrolle: kleinere Courant-Zahl -> kleinere Daempfung je Zeitschritt
    ref05 = lauf_senkrecht(600, 10, 1e-6, courant=0.5)
    kontrolle = {S_: bei(*reflexion(r_, lauf_senkrecht(200, 10, 1e-3, courant=S_)), NL_STANDARD)
                 for S_, r_ in ((S, ref), (0.5, ref05))}
    np.savez(DATEI, vorgabe=vorgabe, nl_achse=nl_achse, aufl=aufl, winkel=winkel,
             kontrolle=np.array([kontrolle[S], kontrolle[0.5]]))


def laden():
    if not DATEI.exists():
        messen()
    return np.load(DATEI)


def main():
    d = laden()
    print("Kontrolle d=10, |R(0)|=1e-3: S=0.99 -> {:.2e}, S=0.5 -> {:.2e}".format(*d["kontrolle"]))
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(11.5, 3.6))

    zeigen = [DICKEN.index(n) for n in DICKEN_ABB]
    for werte, npml, farbe in zip(d["vorgabe"][zeigen], DICKEN_ABB, (FARBEN[0], FARBEN[2])):
        a.loglog(R0_WERTE, werte, "o-", ms=2.5, color=farbe, lw=1.5, label=f"{npml} Zellen")
    a.loglog(R0_WERTE, R0_WERTE, ":", color="k", lw=1.1, label=r"$|R_\mathrm{Sim}| = |R(0)|$")
    a.set_xlim(1e-1, 1e-16)
    a.set_ylim(1e-7, 1)
    a.set_xlabel(r"vorgegebene Restreflexion $|R(0)|$")
    a.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    a.set_title(rf"(a) $\theta = 0$, $N_\lambda = {NL_STANDARD}$", fontsize=9)

    nl = d["nl_achse"]
    for werte, npml, farbe in zip(d["aufl"][zeigen], DICKEN_ABB, (FARBEN[0], FARBEN[2])):
        b.loglog(nl, werte, color=farbe, lw=1.5, label=f"{npml} Zellen")
    r10 = np.interp(10, nl, d["aufl"][DICKEN.index(DICKEN_ABB[0])])
    # tatsaechliche Steigung der 5-Zellen-Kurve im doppelt-logarithmischen Diagramm
    duenn = d["aufl"][DICKEN.index(DICKEN_ABB[0])]
    bereich = (nl >= 6) & (nl <= 20)
    steigung = np.polyfit(np.log(nl[bereich]), np.log(duenn[bereich]), 1)[0]
    print(f"Steigung {DICKEN_ABB[0]} Zellen, N_lambda 6...20: {steigung:.2f}")
    b.loglog(nl, 1.6 * r10 * (nl / 10) ** -2, "--", color="0.45", lw=1.0, label="Steigung $-2$")
    b.axhline(R0_STANDARD, color="k", ls=":", lw=1.1, label=r"$|R(0)| = 10^{-6}$")
    b.set_xlim(6, 30)
    b.set_xticks([6, 8, 10, 15, 20, 30])
    b.set_xticklabels(["6", "8", "10", "15", "20", "30"])
    b.minorticks_off()
    b.set_ylim(1e-7, 1)
    b.set_xlabel(r"Zellen pro Wellenlänge $N_\lambda$")
    b.set_title(r"(b) $\theta = 0$, $|R(0)| = 10^{-6}$", fontsize=9)

    for r0, werte, farbe in zip(WINKEL_VORGABEN, d["winkel"], ("tab:red", "tab:purple")):
        c.semilogy(WINKEL, werte, "o-", color=farbe, ms=4, lw=1.4,
                   label=rf"$|R(0)| = 10^{{{int(np.log10(r0))}}}$")
    c.set_xlim(0, 62)
    c.set_ylim(1e-7, 1)
    c.set_xticks(WINKEL)
    c.set_xlabel(r"Einfallswinkel $\theta$ in Grad")
    c.set_title(rf"(c) {WINKEL_DICKE} Zellen, $N_\lambda = {NL_STANDARD}$", fontsize=9)

    for ax in (a, b, c):
        ax.grid(True, which="major", alpha=0.3)
        ax.legend(fontsize=7, loc="upper right" if ax is not c else "upper left")
    fig.tight_layout()
    ABB.mkdir(parents=True, exist_ok=True)
    fig.savefig(ABB / "testfall_pml.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
