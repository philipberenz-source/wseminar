
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from fdtd_core import FDTD2D, C0

ABB = Path(__file__).resolve().parent.parent / "seminararbeit" / "abbildungen"
D, NPML = 1e-3, 14
COURANT = (0.1, 0.5, 0.99)          # 1 = Stabilitaetsgrenze in der Normierung des Cores
FARBEN = ("tab:blue", "tab:orange", "tab:green")


def ricker(fc):
    a, t0 = (np.pi * fc) ** 2, 1.2 / fc
    return lambda t: (1.0 - 2 * a * (t - t0) ** 2) * np.exp(-a * (t - t0) ** 2)


def linienquelle(sim, i_src, fc, theta=0.0, taper=0):
    puls, j = ricker(fc), np.arange(sim.ny)
    w = (np.ones(sim.ny) if taper == 0 else
         np.sin(0.5 * np.pi * np.clip(np.minimum(j, sim.ny - 1 - j) / taper, 0, 1)) ** 2)
    for jj in j[w > 0]:
        td = (sim.ny - 1 - jj) * sim.d * np.sin(theta) / C0
        sim.add_current_source(i_src, jj, lambda t, td=td: puls(t - td), amplitude=w[jj])


def schritte(sim, zellen):
    return int(round(zellen * sim.d / (C0 * sim.dt)))


def phasengeschwindigkeit(sim, ip1, ip2, abstand, pad=4):
    f, Y1 = sim.get_spectrum(ip1, pad=pad)
    _, Y2 = sim.get_spectrum(ip2, pad=pad)
    i = np.where(np.abs(Y1) > 0.05 * np.abs(Y1).max())[0]
    sl = slice(i[0], i[-1] + 1)
    fb, phi = f[sl], -np.unwrap(np.angle(Y2[sl] / Y1[sl]))
    phi += 2 * np.pi * np.round((2 * np.pi * fb[0] / C0 * abstand - phi[0]) / (2 * np.pi))
    return fb, 2 * np.pi * fb * abstand / phi


def log_achse(ax):
    ax.set_xlim(4, 30)
    ax.set_xticks([4, 5, 6, 8, 10, 15, 20, 30])
    ax.set_xticks([], minor=True)
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel(r"Zellen pro Wellenlänge $N_\lambda$")
    ax.grid(True, which="both", alpha=0.3)


# Teil 1: v_p/c ueber N_lambda bei theta = 0. ny gross, weil Ez[:,0] / Ez[:,-1]
# nie aktualisiert werden (PEC hinter der PML) und das Feld stoeren.

def teil1():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for S, farbe in zip(COURANT, FARBEN):
        sim = FDTD2D(500, 400, D, courant=S, npml=NPML)
        linienquelle(sim, 34, 25e9)
        p1, p2 = sim.add_probe(150, 200), sim.add_probe(350, 200)
        sim.run(schritte(sim, 450))     # Sonde 2 nach 316 Zellen, danach abklingen
        f, v = phasengeschwindigkeit(sim, p1, p2, 200 * D)
        nl = C0 / (f * D)
        m = (nl >= 4) & (nl <= 30)
        a1.plot(nl[m], v[m] / C0, color=farbe, lw=1.5, label=f"$S = {S}$")
        a2.loglog(nl[m], np.abs(1 - v[m] / C0), color=farbe, lw=1.5, label=f"$S = {S}$")
        print(f"  S = {S}: v_p/c bei N_lambda = 10: "
              f"{np.interp(10, nl[::-1], v[::-1]) / C0:.5f}")

    n = np.array([5.0, 30.0])
    a2.loglog(n, 4e-3 * (n / 10.0) ** -2, ":", color="k", lw=1.2, label=r"Steigung $-2$")
    a1.axhline(1.0, color="0.5", lw=0.8, ls=":")
    a1.set_xscale("log")
    a1.set_ylabel(r"$v_\mathrm{p}/c$")
    a2.set_ylabel(r"$|1 - v_\mathrm{p}/c|$")
    for ax, loc in ((a1, "lower right"), (a2, "upper right")):
        log_achse(ax)
        ax.legend(fontsize=8, loc=loc)
    fig.tight_layout()
    fig.savefig(ABB / "testfall_a_dispersion.pdf")
    plt.close(fig)


# Teil 2: Anisotropie v_p(theta) bei N_lambda = 10 -------------------------

def teil2():
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for S, farbe in zip(COURANT, FARBEN):
        winkel, werte = [], []
        for g in np.arange(0.0, 45.1, 2.5):
            i0, j0 = 50, 140
            i1 = i0 + int(round(60 * np.cos(np.deg2rad(g))))
            j1 = j0 - int(round(60 * np.sin(np.deg2rad(g))))
            # Welle in Richtung des GERUNDETEN Sondenpaares anregen und den
            # Abstand daraus bilden - sonst misst man k*cos(Fehlwinkel).
            th = np.arctan2(j0 - j1, i1 - i0)
            sim = FDTD2D(280, 280, D, courant=S, npml=NPML)
            linienquelle(sim, 22, 25e9, th, taper=60)
            p1, p2 = sim.add_probe(i0, j0), sim.add_probe(i1, j1)
            sim.run(schritte(sim, 250))
            f, v = phasengeschwindigkeit(sim, p1, p2, D * np.hypot(i1 - i0, j1 - j0))
            winkel.append(np.rad2deg(th))
            werte.append(np.interp(C0 / (10 * D), f, v) / C0)
        ax.plot(winkel, werte, "o-", color=farbe, lw=1.4, ms=3, label=f"$S = {S}$")
        print(f"  S = {S}: v_p/c von {werte[0]:.5f} (0 Grad) auf {werte[-1]:.5f} (45 Grad)")

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


# Teil 3: Pulsform nach `lauf` Zellen, normiert auf die Amplitude beim Start ---

def teil3(faelle=((5, 0.5), (30, 0.99)), lauf=280):
    fig, axes = plt.subplots(1, len(faelle), figsize=(9.5, 3.6), sharey=True)
    for ax, (nl, S) in zip(axes, faelle):
        sim = FDTD2D(400, 300, D, courant=S, npml=NPML)
        linienquelle(sim, 30, C0 / (nl * D))
        frueh, spaet = schritte(sim, 90), schritte(sim, lauf)
        sim.add_snapshot([frueh, spaet])
        sim.run(spaet)
        e0 = sim._snapshots[frueh]["Ez"][:, 150]
        ax.plot(sim._snapshots[spaet]["Ez"][:, 150] / e0[np.argmax(np.abs(e0))],
                color=FARBEN[COURANT.index(S)], lw=0.9)
        ax.text(0.03, 0.9, rf"$N_\lambda = {nl}$,  $S = {S}$",
                transform=ax.transAxes, fontsize=10)
        ax.set_xlabel("Gitterindex $i$")
        ax.set_xlim(0, 400)
        ax.grid(True, alpha=0.3)
    axes[0].set_ylabel(r"$E_z$ (normiert)")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_a_pulszerfliessen.pdf")
    plt.close(fig)


if __name__ == "__main__":
    ABB.mkdir(parents=True, exist_ok=True)
    print("Teil 1: Dispersion (theta = 0)")
    teil1()
    print("Teil 2: Anisotropie")
    teil2()
    print("Teil 3: Pulszerfliessen")
    teil3()
    print(f"Abbildungen in {ABB}")
