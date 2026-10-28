import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from pathlib import Path

from fdtd_core import FDTD2D, C0

ABB = Path(__file__).resolve().parent.parent / "seminararbeit" / "abbildungen"
D = 1e-3
NPML = 14
COURANT = (0.1, 0.5, 0.99)
FARBEN = ("tab:blue", "tab:orange", "tab:green")


def ricker(fc):
    a = (np.pi * fc) ** 2
    t0 = 1.2 / fc

    def puls(t):
        return (1 - 2 * a * (t - t0) ** 2) * np.exp(-a * (t - t0) ** 2)
    return puls


def verzoegert(puls, td):
    return lambda t: puls(t - td)


def linienquelle(sim, i, fc, theta=0.0, rand=0):
    puls = ricker(fc)
    for j in range(sim.ny):
        if rand == 0:
            gewicht = 1.0
        else:
            abstand = min(j, sim.ny - 1 - j)
            gewicht = np.sin(0.5 * np.pi * np.clip(abstand / rand, 0, 1)) ** 2
        if gewicht > 0:
            td = (sim.ny - 1 - j) * sim.d * np.sin(theta) / C0
            sim.add_current_source(i, j, verzoegert(puls, td), amplitude=gewicht)


def schritte(sim, zellen):
    return int(round(zellen * sim.d / (C0 * sim.dt)))


def phasengeschwindigkeit(sim, sonde1, sonde2, abstand):
    f, y1 = sim.get_spectrum(sonde1, pad=4)
    _, y2 = sim.get_spectrum(sonde2, pad=4)
    idx = np.where(np.abs(y1) > 0.05 * np.abs(y1).max())[0]
    f = f[idx[0]:idx[-1] + 1]
    phase = -np.unwrap(np.angle(y2[idx[0]:idx[-1] + 1] / y1[idx[0]:idx[-1] + 1]))
    phase += 2 * np.pi * np.round((2 * np.pi * f[0] / C0 * abstand - phase[0]) / (2 * np.pi))
    return f, 2 * np.pi * f * abstand / phase


def log_achse(ax):
    ax.set_xlim(4, 30)
    ax.set_xticks([4, 5, 6, 8, 10, 15, 20, 30])
    ax.set_xticks([], minor=True)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlabel(r"Zellen pro Wellenlänge $N_\lambda$")
    ax.grid(True, which="both", alpha=0.3)


def dispersion():
    fig, (links, rechts) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for S, farbe in zip(COURANT, FARBEN):
        sim = FDTD2D(500, 400, D, courant=S, npml=NPML)
        linienquelle(sim, 34, 25e9)
        s1 = sim.add_probe(150, 200)
        s2 = sim.add_probe(350, 200)
        sim.run(schritte(sim, 450))
        f, v = phasengeschwindigkeit(sim, s1, s2, 200 * D)
        nl = C0 / (f * D)
        m = (nl >= 4) & (nl <= 30)
        links.plot(nl[m], v[m] / C0, color=farbe, lw=1.5, label=f"$S = {S}$")
        rechts.loglog(nl[m], np.abs(1 - v[m] / C0), color=farbe, lw=1.5, label=f"$S = {S}$")
        print(f"S = {S}: v_p/c bei N = 10: {np.interp(10, nl[::-1], v[::-1]) / C0:.5f}")

    n = np.array([5.0, 30.0])
    rechts.loglog(n, 4e-3 * (n / 10) ** -2, ":", color="k", lw=1.2, label=r"Steigung $-2$")
    links.axhline(1.0, color="0.5", lw=0.8, ls=":")
    links.set_xscale("log")
    links.set_ylabel(r"$v_\mathrm{p}/c$")
    rechts.set_ylabel(r"$|1 - v_\mathrm{p}/c|$")
    log_achse(links)
    log_achse(rechts)
    links.legend(fontsize=8, loc="lower right")
    rechts.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_a_dispersion.pdf")
    plt.close(fig)


def anisotropie():
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    i0, j0 = 50, 140
    for S, farbe in zip(COURANT, FARBEN):
        winkel = []
        werte = []
        for grad in np.arange(0.0, 45.1, 2.5):
            i1 = i0 + int(round(60 * np.cos(np.deg2rad(grad))))
            j1 = j0 - int(round(60 * np.sin(np.deg2rad(grad))))
            theta = np.arctan2(j0 - j1, i1 - i0)
            sim = FDTD2D(280, 280, D, courant=S, npml=NPML)
            linienquelle(sim, 22, 25e9, theta, rand=60)
            s1 = sim.add_probe(i0, j0)
            s2 = sim.add_probe(i1, j1)
            sim.run(schritte(sim, 250))
            f, v = phasengeschwindigkeit(sim, s1, s2, D * np.hypot(i1 - i0, j1 - j0))
            winkel.append(np.rad2deg(theta))
            werte.append(np.interp(C0 / (10 * D), f, v) / C0)
        ax.plot(winkel, werte, "o-", color=farbe, lw=1.4, ms=3, label=f"$S = {S}$")
        print(f"S = {S}: v_p/c bei 0 Grad {werte[0]:.5f}, bei 45 Grad {werte[-1]:.5f}")

    ax.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel(r"Ausbreitungsrichtung $\theta$ in Grad")
    ax.set_ylabel(r"$v_\mathrm{p}/c$  bei  $N_\lambda = 10$")
    ax.set_xlim(0, 45)
    ax.set_xticks(np.arange(0.0, 46.0, 5.0))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_a_anisotropie.pdf")
    plt.close(fig)


def pulsform():
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    for ax, (nl, S) in zip(axes, ((5, 0.5), (30, 0.99))):
        sim = FDTD2D(400, 300, D, courant=S, npml=NPML)
        linienquelle(sim, 30, C0 / (nl * D))
        frueh = schritte(sim, 90)
        spaet = schritte(sim, 280)
        sim.add_snapshot([frueh, spaet])
        sim.run(spaet)
        start = sim._snapshots[frueh]["Ez"][:, 150]
        ende = sim._snapshots[spaet]["Ez"][:, 150]
        ax.plot(ende / start[np.argmax(np.abs(start))], color=FARBEN[COURANT.index(S)], lw=0.9)
        ax.text(0.03, 0.9, rf"$N_\lambda = {nl}$,  $S = {S}$", transform=ax.transAxes, fontsize=10)
        ax.set_xlabel("Gitterindex $i$")
        ax.set_xlim(0, 400)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"$E_z$ (normiert)")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_a_pulszerfliessen.pdf")
    plt.close(fig)


if __name__ == "__main__":
    ABB.mkdir(parents=True, exist_ok=True)
    dispersion()
    anisotropie()
    pulsform()
