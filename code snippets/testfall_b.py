
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from fdtd_core import FDTD2D, C0, EPS0, MU0

ABB = Path(__file__).resolve().parent.parent / "seminararbeit" / "abbildungen"
D, N = 1e-3, 101                     # Zellbreite, Gittergroesse (PEC-Kasten)
FC = C0 / (20 * D)                   # Ricker-Spitzenfrequenz: N_lambda = 20
COURANT = (0.5, 0.99, 1.0)
FARBEN = ("tab:blue", "tab:orange", "tab:green")
INSTABIL = (1.01, 1.05)                  # knapp oberhalb der Stabilitaetsgrenze S = 1
FARBEN_INSTABIL = ("tab:red", "tab:purple")
LAUFWEG = 20000                      # Zellen Laufweg, rund 200 Durchquerungen


def ricker_puls(fc):
    """Ricker-Puls, der nach t = 2 t0 hart abgeschaltet wird: Die Quelle liefert
    nur diesen einen Puls, danach ist J_z = 0 (Restwert beim Abschalten ~1e-5)."""
    a, t0 = (np.pi * fc) ** 2, 1.2 / fc
    return lambda t: ((1.0 - 2 * a * (t - t0) ** 2) * np.exp(-a * (t - t0) ** 2)
                      if t <= 2 * t0 else 0.0)


def lauf(S, zellen, abbruch=None):
    """PEC-Kasten (npml = 0: die nie aktualisierten Randwerte E_z = 0 bilden
    den ideal leitenden Rand), Stromquelle im Zentrum.
    Rueckgabe je Zeitpunkt n: Laufweg c t / Delta, Feldenergie W und
    eingespeiste Energie W_Q.
      W   = Delta^2/2 * sum(eps0 Ez^n^2 + mu0 H^(n-1/2) H^(n+1/2))
      W_Q = -sum_k Delta^2 dt J^(k+1/2) (Ez^k + Ez^(k+1))/2   an der Quellzelle"""
    sim = FDTD2D(N, N, D, courant=S, npml=0)
    c = N // 2
    sim.add_current_source(c, c, ricker_puls(FC))
    schritte = int(round(zellen * D / (C0 * sim.dt)))
    W, WQ, WN = np.empty(schritte), np.empty(schritte), np.empty(schritte)
    wq = 0.0
    hx, hy = sim.Hx.copy(), sim.Hy.copy()
    for n in range(schritte):
        ez = sim.Ez.copy()                         # E^n
        sim.step()                                 # -> H^(n+1/2), E^(n+1)
        W[n] = (0.5 * EPS0 * np.sum(ez ** 2)
                + 0.5 * MU0 * (np.sum(hx * sim.Hx) + np.sum(hy * sim.Hy))) * D * D
        WN[n] = sim.total_energy()                 # E^(n+1), H^(n+1/2): immer positiv
        WQ[n] = wq                                 # eingespeist bis zum Zeitpunkt n
        wq -= D * D * sim.dt * sim.Jz[c, c] * 0.5 * (ez[c, c] + sim.Ez[c, c])
        hx, hy = sim.Hx.copy(), sim.Hy.copy()
        if abbruch is not None and WN[n] > abbruch:  # instabil: vor dem Ueberlauf aufhoeren
            W, WQ, WN = W[:n + 1], WQ[:n + 1], WN[:n + 1]
            break
    x = np.arange(len(W)) * sim.dt * C0 / D
    return x, W, WQ, WN


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    daten = {S: lauf(S, LAUFWEG) for S in COURANT}
    w_ref = daten[0.99][2][-1]                     # eingespeiste Energie bei S = 0.99
    instabil = {S: lauf(S, 400, abbruch=1e32 * w_ref) for S in INSTABIL}

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.6))

    # (a) Anregung und erste Durchlaeufe, logarithmisch: stabile und instabile S
    for S, farbe in zip(COURANT, FARBEN):
        x, W, WQ, _ = daten[S]
        m = x <= 400
        a1.semilogy(x[m], W[m] / w_ref, color=farbe, lw=1.4, label=f"$S = {S}$")
    x, W, WQ, _ = daten[0.99]
    k = np.flatnonzero((x > 1) & (x <= 400))[::60]
    a1.semilogy(x[k], WQ[k] / w_ref, "o", color="k", ms=3.5, mfc="none",
                label=r"$W_Q$ ($S = 0.99$)")
    for S, farbe in zip(INSTABIL, FARBEN_INSTABIL):
        x, _, _, WN = instabil[S]
        a1.semilogy(x, WN / w_ref, color=farbe, lw=1.4, label=f"$S = {S}$")
    a1.set_xlim(0, 400)
    a1.set_ylim(1e-3, 1e30)
    a1.set_xlabel(r"Laufweg $c\,t/\Delta$ in Zellen")
    a1.set_ylabel(r"$W \,/\, W_Q$")
    a1.grid(True, alpha=0.3)
    a1.legend(fontsize=8, loc="upper right", ncol=2)

    # (b) Langzeit bei S <= 1: groesste bis t aufgetretene Abweichung |W - W_Q| / W_Q
    for S, farbe in zip(COURANT, FARBEN):
        x, W, WQ, _ = daten[S]
        abw = np.maximum.accumulate(np.abs(W - WQ) / WQ[-1])
        a2.semilogy(x, abw + 1e-17, color=farbe, lw=1.3, label=f"$S = {S}$")
    a2.set_xlim(0, LAUFWEG)
    a2.set_ylim(1e-17, 1e-13)
    a2.set_xlabel(r"Laufweg $c\,t/\Delta$ in Zellen")
    a2.set_ylabel(r"$\max\,|W - W_Q|\,/\,W_Q$ bis $t$")
    a2.grid(True, alpha=0.3)
    a2.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(ABB / "testfall_b_energie.pdf")
    plt.close(fig)

    for S in COURANT:
        x, W, WQ, _ = daten[S]
        print(f"  S = {S}: W_Q = {WQ[-1]:.6e} J/m ({WQ[-1] / w_ref - 1:+.2%} ggü. S = 0.99), "
              f"max |W - W_Q|/W_Q = {np.abs(W - WQ).max() / WQ[-1]:.1e}")
    for S in INSTABIL:
        x, _, _, WN = instabil[S]
        m = (WN > 1e8 * w_ref) & (WN < 1e28 * w_ref)
        steig = np.polyfit(x[m], np.log10(WN[m]), 1)[0]
        print(f"  S = {S}: W > 1e32 W_Q nach {x[-1]:.0f} Zellen, Anstieg {steig:.2f} Dekaden je Zelle")


if __name__ == "__main__":
    main()
    print(f"Abbildung in {ABB}")
