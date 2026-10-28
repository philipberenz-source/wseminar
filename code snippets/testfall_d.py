import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import FDTD2D, C0
from testfall_a import ABB, D, NPML, COURANT, FARBEN, linienquelle, schritte, log_achse

N1 = 1.0
N2 = 2.0
N2_VERGLEICH = (np.sqrt(2.0), 3.0)
NX, NY = 300, 400
I_QUELLE, I_GRENZE, J = 34, 150, 200
I_VOR, I_HINTER = 140, 160
FC = C0 / (15 * D)
LAUFWEG = 500


def lauf(S, eps_r, linie=None):
    sim = FDTD2D(NX, NY, D, courant=S, npml=NPML)
    sim.set_material_box(I_GRENZE, NX, 0, NY, eps_r=eps_r)
    linienquelle(sim, I_QUELLE, FC)
    sim.add_probe(I_VOR, J)
    sim.add_probe(I_HINTER, J)
    if linie is None:
        sim.run(schritte(sim, LAUFWEG))
    else:
        sim.run(schritte(sim, LAUFWEG), lambda s: linie.append(s.Ez[:, J].copy()))
    return sim


def k_halbe(nl, S, n):
    sc = S / np.sqrt(2.0)
    arg = n / sc * np.sin(np.pi * sc / nl)
    return np.where(arg < 0.95, np.arcsin(np.minimum(arg, 1.0)), np.nan)


def koeffizienten(vakuum, medium, S, n2):
    f, ein_vor = vakuum.get_spectrum(0, 4)
    ein_hinter = vakuum.get_spectrum(1, 4)[1]
    vor = medium.get_spectrum(0, 4)[1]
    hinter = medium.get_spectrum(1, 4)[1]
    m = np.abs(ein_vor) > 0.05 * np.abs(ein_vor).max()
    nl = C0 / (f[m] * D)
    k1 = k_halbe(nl, S, N1)
    k2 = k_halbe(nl, S, n2)
    nenner = N1 * np.cos(k2) + n2 * np.cos(k1)
    gamma = np.abs(vor - ein_vor)[m] / np.abs(ein_vor[m])
    T = np.abs(hinter[m] / ein_hinter[m])
    gamma_gitter = (N1 * np.cos(k2) - n2 * np.cos(k1)) / nenner
    T_gitter = 2 * N1 * np.cos(k1) / nenner
    return nl, gamma, T, gamma_gitter, T_gitter


def fresnel(n2):
    return (N1 - n2) / (N1 + n2), 2 * N1 / (N1 + n2)


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    erg = {}
    linie = []
    for S in COURANT:
        vakuum = lauf(S, 1.0)
        if S == 0.99:
            medium = lauf(S, N2 ** 2, linie)
        else:
            medium = lauf(S, N2 ** 2)
        erg[N2, S] = koeffizienten(vakuum, medium, S, N2)
    for n2 in N2_VERGLEICH:
        erg[n2, 0.99] = koeffizienten(vakuum, lauf(0.99, n2 ** 2), 0.99, n2)

    for (n2, S), (nl, gamma, T, gamma_gitter, T_gitter) in erg.items():
        k = (nl >= 7) & (nl <= 30)
        print(f"eps_r = {n2 ** 2:.0f}, S = {S}: "
              f"max ||Gamma| - |Gamma~|| = {np.nanmax(np.abs(gamma - abs(gamma_gitter))[k]):.1e}, "
              f"max |T - T~| = {np.nanmax(np.abs(T - T_gitter)[k]):.1e}, "
              f"|Gamma| bei N = 10: {np.interp(10, nl[::-1], gamma[::-1]):.4f}")

    g0, t0 = fresnel(N2)
    fig, (links, rechts) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for S, farbe in zip(COURANT, FARBEN):
        nl, gamma, T, _, _ = erg[N2, S]
        links.plot(nl, gamma, color=farbe, lw=1.5, label=f"$S = {S}$")
        links.plot(nl, T, color=farbe, lw=1.5)
        rechts.loglog(nl, np.abs(gamma - abs(g0)) / abs(g0), color=farbe, lw=1.5)
        rechts.loglog(nl, np.abs(T - t0) / t0, "--", color=farbe, lw=1.2)
    links.axhline(-g0, color="k", lw=0.9, ls=":")
    links.text(29, -g0 + 0.012, r"$|\Gamma|$ (Fresnel)", ha="right", fontsize=8)
    links.axhline(t0, color="k", lw=0.9, ls=":")
    links.text(29, t0 + 0.012, r"$T$ (Fresnel)", ha="right", fontsize=8)
    n = np.array([7.0, 30.0])
    rechts.loglog(n, 0.12 * (n / 10) ** -2, ":", color="k", lw=1.1)
    rechts.plot([], [], "-", color="k", lw=1.1, label=r"$||\Gamma_\mathrm{Sim}| - |\Gamma||\,/\,|\Gamma|$")
    rechts.plot([], [], "--", color="k", lw=1.1, label=r"$|T_\mathrm{Sim} - T|\,/\,T$")
    rechts.plot([], [], ":", color="k", lw=1.1, label=r"Steigung $-2$")
    links.set_ylim(0.3, 0.8)
    links.set_ylabel(r"$|\Gamma_\mathrm{Sim}|$ bzw. $T_\mathrm{Sim}$")
    rechts.set_ylabel("relative Abweichung von Fresnel")
    rechts.set_ylim(2e-3, 1.5)
    for ax, ort in ((links, "center right"), (rechts, "lower left")):
        ax.set_xscale("log")
        log_achse(ax)
        ax.set_xlim(7, 30)
        ax.legend(fontsize=7.5, loc=ort)
    fig.tight_layout()
    fig.savefig(ABB / "testfall_d_courant.pdf")
    plt.close(fig)

    fig, (links, rechts) = plt.subplots(1, 2, figsize=(9.5, 3.6))
    for n2, farbe in zip((N2_VERGLEICH[0], N2, N2_VERGLEICH[1]), FARBEN):
        nl, gamma, _, gamma_gitter, _ = erg[n2, 0.99]
        g = fresnel(n2)[0]
        links.plot(nl, np.where(np.isfinite(gamma_gitter), gamma / abs(g), np.nan),
                   color=farbe, lw=2.4, alpha=0.5, label=rf"$\varepsilon_r = {n2 ** 2:.0f}$")
        links.plot(nl, gamma_gitter / g, ":", color="k", lw=1.1)
    links.plot([], [], ":", color="k", lw=1.1, label=r"$\tilde\Gamma / \Gamma$ (Gitter)")
    links.axhline(1.0, color="0.5", lw=0.8, ls="--")
    links.set_ylim(0.95, 2.0)
    links.set_ylabel(r"$|\Gamma_\mathrm{Sim}| \,/\, |\Gamma|$")
    links.set_xscale("log")
    log_achse(links)
    links.set_xlim(7, 30)
    links.legend(fontsize=7.5, loc="upper right")

    ez = np.array(linie)
    ez /= np.abs(ez).max()
    ct = np.arange(1, len(ez) + 1) * 0.99 / np.sqrt(2.0)
    rechts.pcolormesh(np.arange(NX), ct, ez, cmap="RdBu_r", vmin=-0.5, vmax=0.5,
                      shading="auto", rasterized=True)
    rechts.axvline(I_GRENZE - 0.5, color="k", lw=0.9, ls="--")
    rechts.axvspan(I_GRENZE - 0.5, NX - NPML, color="0.5", alpha=0.08, lw=0)
    kasten = dict(fc="white", ec="none", alpha=0.8, pad=1)
    rechts.text(80, 95, "einfallend", fontsize=8, ha="center", bbox=kasten)
    rechts.text(85, 235, "reflektiert", fontsize=8, ha="center", bbox=kasten)
    rechts.text(205, 330, "transmittiert", fontsize=8, ha="center", bbox=kasten)
    rechts.text(I_GRENZE - 6, 20, r"$\varepsilon_r = 1$", ha="right", fontsize=8)
    rechts.text(I_GRENZE + 6, 20, r"$\varepsilon_r = 4$", ha="left", fontsize=8)
    rechts.set_xlim(NPML, NX - NPML)
    rechts.set_ylim(0, 420)
    rechts.set_xlabel("Gitterindex $i$")
    rechts.set_ylabel(r"Laufweg $c\,t/\Delta$ in Zellen")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_d_material.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
