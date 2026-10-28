import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import FDTD2D, C0, EPS0, MU0
from testfall_a import ABB, D, FARBEN

N = 101
FC = C0 / (20 * D)
STABIL = (0.5, 0.99, 1.0)
INSTABIL = (1.01, 1.05)
FARBEN_INSTABIL = ("tab:red", "tab:purple")
LAUFWEG = 20000


def kurzer_puls(fc):
    a = (np.pi * fc) ** 2
    t0 = 1.2 / fc

    def puls(t):
        if t > 2 * t0:
            return 0.0
        return (1 - 2 * a * (t - t0) ** 2) * np.exp(-a * (t - t0) ** 2)
    return puls


def lauf(S, zellen, grenze=None):
    sim = FDTD2D(N, N, D, courant=S, npml=0)
    m = N // 2
    sim.add_current_source(m, m, kurzer_puls(FC))
    anzahl = int(round(zellen * D / (C0 * sim.dt)))

    W = np.empty(anzahl)
    WQ = np.empty(anzahl)
    W_gesamt = np.empty(anzahl)
    eingespeist = 0.0
    hx, hy = sim.Hx.copy(), sim.Hy.copy()

    for n in range(anzahl):
        ez = sim.Ez.copy()
        sim.step()
        W[n] = (0.5 * EPS0 * np.sum(ez ** 2)
                + 0.5 * MU0 * (np.sum(hx * sim.Hx) + np.sum(hy * sim.Hy))) * D * D
        W_gesamt[n] = sim.total_energy()
        WQ[n] = eingespeist
        eingespeist -= D * D * sim.dt * sim.Jz[m, m] * 0.5 * (ez[m, m] + sim.Ez[m, m])
        hx, hy = sim.Hx.copy(), sim.Hy.copy()
        if grenze is not None and W_gesamt[n] > grenze:
            W, WQ, W_gesamt = W[:n + 1], WQ[:n + 1], W_gesamt[:n + 1]
            break

    x = np.arange(len(W)) * sim.dt * C0 / D
    return x, W, WQ, W_gesamt


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    stabil = {S: lauf(S, LAUFWEG) for S in STABIL}
    w_ref = stabil[0.99][2][-1]
    instabil = {S: lauf(S, 400, grenze=1e32 * w_ref) for S in INSTABIL}

    fig, (links, rechts) = plt.subplots(1, 2, figsize=(9.5, 3.6))

    for S, farbe in zip(STABIL, FARBEN):
        x, W, WQ, _ = stabil[S]
        m = x <= 400
        links.semilogy(x[m], W[m] / w_ref, color=farbe, lw=1.4, label=f"$S = {S}$")
    x, W, WQ, _ = stabil[0.99]
    k = np.flatnonzero((x > 1) & (x <= 400))[::60]
    links.semilogy(x[k], WQ[k] / w_ref, "o", color="k", ms=3.5, mfc="none",
                   label=r"$W_Q$ ($S = 0.99$)")
    for S, farbe in zip(INSTABIL, FARBEN_INSTABIL):
        x, _, _, W_gesamt = instabil[S]
        links.semilogy(x, W_gesamt / w_ref, color=farbe, lw=1.4, label=f"$S = {S}$")
    links.set_xlim(0, 400)
    links.set_ylim(1e-3, 1e30)
    links.set_xlabel(r"Laufweg $c\,t/\Delta$ in Zellen")
    links.set_ylabel(r"$W \,/\, W_Q$")
    links.grid(True, alpha=0.3)
    links.legend(fontsize=8, loc="upper right", ncol=2)

    for S, farbe in zip(STABIL, FARBEN):
        x, W, WQ, _ = stabil[S]
        abweichung = np.maximum.accumulate(np.abs(W - WQ) / WQ[-1])
        rechts.semilogy(x, abweichung + 1e-17, color=farbe, lw=1.3, label=f"$S = {S}$")
    rechts.set_xlim(0, LAUFWEG)
    rechts.set_ylim(1e-17, 1e-13)
    rechts.set_xlabel(r"Laufweg $c\,t/\Delta$ in Zellen")
    rechts.set_ylabel(r"$\max\,|W - W_Q|\,/\,W_Q$ bis $t$")
    rechts.grid(True, alpha=0.3)
    rechts.legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(ABB / "testfall_b_energie.pdf")
    plt.close(fig)

    for S in STABIL:
        x, W, WQ, _ = stabil[S]
        print(f"S = {S}: W_Q = {WQ[-1]:.6e} J/m ({WQ[-1] / w_ref - 1:+.2%} gegenüber S = 0.99), "
              f"max |W - W_Q|/W_Q = {np.abs(W - WQ).max() / WQ[-1]:.1e}")
    for S in INSTABIL:
        x, _, _, W_gesamt = instabil[S]
        m = (W_gesamt > 1e8 * w_ref) & (W_gesamt < 1e28 * w_ref)
        steigung = np.polyfit(x[m], np.log10(W_gesamt[m]), 1)[0]
        print(f"S = {S}: W > 1e32 W_Q nach {x[-1]:.0f} Zellen, {steigung:.2f} Dekaden pro Zelle")


if __name__ == "__main__":
    main()
