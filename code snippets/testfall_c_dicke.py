import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from testfall_a import ABB
from testfall_c import lauf, senkrecht, r0_text, L_GROSS, R0_WERTE, R0_DREI, FARBEN_R0

# Gemessene Reflexion bei senkrechtem Einfall ueber der PML-Dicke (N_lambda = 10, m = 3).
# Aufbau wie Testfall 4 (testfall_c.py); zusaetzlich das Minimum ueber alle Vorgaben |R(0)|.
DICKEN = np.array([2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 25, 30])


def main():
    ABB.mkdir(parents=True, exist_ok=True)
    ref = lauf(L_GROSS, 10, 1e-6)
    feste = {r0: np.array([senkrecht(ref, n, r0) for n in DICKEN]) for r0 in R0_DREI}
    opt, opt_r0 = [], []
    for n in DICKEN:
        w = np.array([senkrecht(ref, n, r0) for r0 in R0_WERTE])
        k = np.argmin(w)
        opt.append(w[k]); opt_r0.append(R0_WERTE[k])
    opt = np.array(opt)

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for (r0, r), farbe in zip(feste.items(), FARBEN_R0):
        ax.semilogy(DICKEN, r, "o-", ms=3.5, color=farbe, lw=1.5, label=r0_text(r0))
        if r0 > 1e-7:
            ax.axhline(r0, color=farbe, lw=0.9, ls=":")
    ax.semilogy(DICKEN, opt, "s--", ms=3.5, color="k", lw=1.2, label=r"Minimum über $|R(0)|$")
    ax.axvline(14, color="0.5", lw=0.8, ls="--")
    ax.text(14.4, 3e-2, "Standardwert", fontsize=8, color="0.35")
    ax.set_xlabel("Dicke der PML in Zellen")
    ax.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    ax.set_xlim(0, 31)
    ax.grid(True, which="major", alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower left")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_c_dicke.pdf")
    plt.close(fig)

    for i, n in enumerate(DICKEN):
        print(f"  {n:2d} Zellen: " + ", ".join(f"{r0:.0e}: {feste[r0][i]:.1e}" for r0 in R0_DREI)
              + f" | Minimum {opt[i]:.1e} bei {opt_r0[i]:.0e}")


if __name__ == "__main__":
    main()
    print(f"Abbildung in {ABB}")
