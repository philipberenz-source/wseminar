
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import FDTD2D, C0
from testfall_a import ABB, D, NPML, COURANT, FARBEN, linienquelle, schritte, log_achse

# Senkrechter Einfall aus dem Vakuum auf einen dielektrischen Halbraum (Schneider Abschn. 5.8, 7.6).
# Zwei Laeufe (ohne/mit Dielektrikum), Betraege der Spektrenverhaeltnisse wie Schneider Abb. 5.6(a).
# Die Grenzflaeche liegt auf dem H_y-Knoten zwischen i = I_G - 1 und I_G.
N1, N_HAUPT, N_VERGL = 1.0, 2.0, (np.sqrt(2.0), 3.0)   # Brechzahlen: eps_r = 4, Vergleich eps_r = 2, 9
NX, NY, I_Q, I_G, J = 300, 400, 34, 150, 200
I_R, I_T = 140, 160                        # Sonden 9,5 bzw. 10,5 Zellen vor/hinter der Grenzflaeche
FC, LAUF = C0 / (15 * D), 500              # Ricker-Spitze bei N_lambda = 15; Laufweg in Zellen


def lauf(S, eps_r, linie=None):
    """Ein Lauf; optional wird E_z entlang der Mittellinie j = J in die Liste linie geschrieben."""
    sim = FDTD2D(NX, NY, D, courant=S, npml=NPML)
    sim.set_material_box(I_G, NX, 0, NY, eps_r=eps_r)
    linienquelle(sim, I_Q, FC)
    sim.add_probe(I_R, J), sim.add_probe(I_T, J)
    sim.run(schritte(sim, LAUF), None if linie is None else lambda s: linie.append(s.Ez[:, J].copy()))
    return sim


def kh(nl, S, n):
    """k~ Delta/2 aus der 1D-Dispersionsrelation (Schneider Gl. 7.36), S_c = S/sqrt(2);
    NaN, wo sich im Medium keine Welle mehr ausbreitet."""
    sc = S / np.sqrt(2.0)
    arg = n / sc * np.sin(np.pi * sc / nl)
    return np.where(arg < 0.95, np.arcsin(np.minimum(arg, 1.0)), np.nan)


def koeffizienten(a, b, S, n2, pad=4):
    """Lauf a: Vakuum (einfallend), Lauf b: mit Dielektrikum. Betraege der Spektren-
    verhaeltnisse (Schneider Gl. 5.98, Abb. 5.6a) und Gitterformeln (Gl. 7.71/7.72).
    Rueckgabe N_lambda, |Gamma|, |T|, Gamma~, T~."""
    f, Ar = a.get_spectrum(0, pad)
    At, Br, Bt = a.get_spectrum(1, pad)[1], b.get_spectrum(0, pad)[1], b.get_spectrum(1, pad)[1]
    m = np.abs(Ar) > 0.05 * np.abs(Ar).max()
    nl = C0 / (f[m] * D)
    k1, k2 = kh(nl, S, N1), kh(nl, S, n2)
    nenn = N1 * np.cos(k2) + n2 * np.cos(k1)
    return (nl, np.abs(Br - Ar)[m] / np.abs(Ar[m]), np.abs(Bt[m] / At[m]),
            (N1 * np.cos(k2) - n2 * np.cos(k1)) / nenn, 2 * N1 * np.cos(k1) / nenn)


def fresnel(n2):                                               # Schneider Gl. 5.94/5.95
    return (N1 - n2) / (N1 + n2), 2 * N1 / (N1 + n2)


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    erg = {}
    for S in COURANT:                                          # letzter Wert: S = 0.99
        a = lauf(S, 1.0)
        linie = [] if S == 0.99 else None                       # Weg-Zeit-Diagramm bei S = 0.99
        erg[N_HAUPT, S] = koeffizienten(a, lauf(S, N_HAUPT ** 2, linie), S, N_HAUPT)
    for n2 in N_VERGL:                                         # Vakuumlauf a (S = 0.99) wiederverwendet
        erg[n2, 0.99] = koeffizienten(a, lauf(0.99, n2 ** 2), 0.99, n2)
    for (n2, S), (nl, gam, T, gs, ts) in erg.items():
        k = (nl >= 7) & (nl <= 30)
        print(f"  eps_r = {n2 ** 2:.0f}, S = {S}: max||Gamma| - |Gamma~|| = {np.nanmax(np.abs(gam - abs(gs))[k]):.1e}, "
              f"max|T - T~| = {np.nanmax(np.abs(T - ts)[k]):.1e}, "
              f"|Gamma|(N=10) = {np.interp(10, nl[::-1], gam[::-1]):.4f}")

    # Abbildung 1: eps_r = 4, verschiedene Courant-Zahlen
    g0, t0 = fresnel(N_HAUPT)
    n = np.array([7.0, 30.0])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for S, farbe in zip(COURANT, FARBEN):
        nl, gam, T, gs, ts = erg[N_HAUPT, S]
        a1.plot(nl, gam, color=farbe, lw=1.5, label=f"$S = {S}$")
        a1.plot(nl, T, color=farbe, lw=1.5)
        a2.loglog(nl, np.abs(gam - abs(g0)) / abs(g0), color=farbe, lw=1.5)
        a2.loglog(nl, np.abs(T - t0) / t0, "--", color=farbe, lw=1.2)
    for y0, name in ((-g0, r"$|\Gamma|$"), (t0, r"$T$")):
        a1.axhline(y0, color="k", lw=0.9, ls=":")
        a1.text(29, y0 + 0.012, name + " (Fresnel)", ha="right", fontsize=8)
    a2.loglog(n, 0.12 * (n / 10.0) ** -2, ":", color="k", lw=1.1)
    for stil, name in (("-", r"$||\Gamma_\mathrm{Sim}| - |\Gamma||\,/\,|\Gamma|$"),
                       ("--", r"$|T_\mathrm{Sim} - T|\,/\,T$"), (":", r"Steigung $-2$")):
        a2.plot([], [], stil, color="k", lw=1.1, label=name)
    a1.set_ylim(0.3, 0.8)
    a1.set_ylabel(r"$|\Gamma_\mathrm{Sim}|$ bzw. $T_\mathrm{Sim}$")
    a2.set_ylabel("relative Abweichung von Fresnel")
    a2.set_ylim(2e-3, 1.5)
    for ax, loc in ((a1, "center right"), (a2, "lower left")):
        ax.set_xscale("log")
        log_achse(ax)
        ax.set_xlim(7, 30)
        ax.legend(fontsize=7.5, loc=loc)
    fig.tight_layout()
    fig.savefig(ABB / "testfall_d_courant.pdf")
    plt.close(fig)

    # Abbildung 2: links verschiedene eps_r, rechts Weg-Zeit-Diagramm (jeweils S = 0.99)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.6))
    for n2, farbe in zip((N_VERGL[0], N_HAUPT, N_VERGL[1]), FARBEN):
        nl, gam, _, gs, _ = erg[n2, 0.99]
        g0 = fresnel(n2)[0]
        a1.plot(nl, np.where(np.isfinite(gs), gam / abs(g0), np.nan), color=farbe, lw=2.4, alpha=0.5,
                label=rf"$\varepsilon_r = {n2 ** 2:.0f}$")
        a1.plot(nl, gs / g0, ":", color="k", lw=1.1)
    a1.plot([], [], ":", color="k", lw=1.1, label=r"$\tilde\Gamma / \Gamma$ (Gitter)")
    a1.axhline(1.0, color="0.5", lw=0.8, ls="--")
    a1.set_ylim(0.95, 2.0)
    a1.set_ylabel(r"$|\Gamma_\mathrm{Sim}| \,/\, |\Gamma|$")
    # rechts: Weg-Zeit-Diagramm von E_z auf der Mittellinie (eps_r = 4, S = 0.99)
    ez = np.array(linie)
    ez /= np.abs(ez).max()
    ct = np.arange(1, len(ez) + 1) * 0.99 / np.sqrt(2.0)             # Laufweg c t / Delta
    a2.pcolormesh(np.arange(NX), ct, ez, cmap="RdBu_r", vmin=-0.5, vmax=0.5,
                  shading="auto", rasterized=True)
    a2.axvline(I_G - 0.5, color="k", lw=0.9, ls="--")
    a2.axvspan(I_G - 0.5, NX - NPML, color="0.5", alpha=0.08, lw=0)
    for x, y, txt in ((80, 95, "einfallend"), (85, 235, "reflektiert"), (205, 330, "transmittiert")):
        a2.text(x, y, txt, fontsize=8, ha="center", bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
    a2.text(I_G - 6, 20, r"$\varepsilon_r = 1$", ha="right", fontsize=8)
    a2.text(I_G + 6, 20, r"$\varepsilon_r = 4$", ha="left", fontsize=8)
    a2.set_xlim(NPML, NX - NPML)
    a2.set_ylim(0, 420)
    a2.set_xlabel("Gitterindex $i$")
    a2.set_ylabel(r"Laufweg $c\,t/\Delta$ in Zellen")
    a1.set_xscale("log")
    log_achse(a1)
    a1.set_xlim(7, 30)
    a1.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_d_material.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
    print(f"Abbildungen in {ABB}")
