import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fdtd_core import FDTD2D, C0
from testfall_a import ABB, D, FARBEN, linienquelle, schritte

# Ebene Welle (Linienquelle J_z wie in Testfall 1, Ricker-Puls), senkrechter Einfall auf die rechte PML.
# Lauf A: Gebiet, dessen rechter Rand waehrend der Messdauer nicht erreicht wird (Referenz),
# Lauf B: Gebiet mit PML. Gemessen wird |E_B - E_A| / |E_A| im Frequenzbereich (wie Testfall 3).
# Die ebene Welle nimmt mit dem Laufweg nicht ab, eine Wegkorrektur entfaellt.
# Das Gebiet ist in y so hoch, dass Stoerungen von den oberen/unteren Raendern die Sonde
# innerhalb der Messdauer nicht erreichen.
L_KLEIN, L_GROSS = 200, 600          # Innengebiet in x (ohne PML)
H_INNEN = 500                        # Innengebiet in y (ohne PML)
ABSTAND = 95                         # Quelllinie -> Sonde, Sonde 5 Zellen vor der PML
LAUF = 210                           # Laufweg in Zellen; Reflexionen der anderen Raender kommen spaeter
FC = C0 / (15 * D)                   # Ricker-Spitze bei N_lambda = 15
S = 0.99
M = 3
N_BEZUG = 10                         # Auswertung bei N_lambda = 10 wie in den anderen Testfaellen

DICKEN = ((5, FARBEN[0]), (10, FARBEN[1]), (20, FARBEN[2]))
R0_WERTE = np.logspace(-1, -20, 20)


def lauf(innen, npml, r0):
    nx, ny = innen + 2 * npml, H_INNEN + 2 * npml
    xs = nx // 2
    sim = FDTD2D(nx, ny, D, courant=S, npml=npml, pml_m=M, r0=r0)
    linienquelle(sim, xs, FC)
    sim.add_probe(xs + ABSTAND, ny // 2)
    sim.run(schritte(sim, LAUF))
    return sim


def restreflexion(ref, sim, k=0, pad=8):
    """|R_Sim| ueber N_lambda aus Referenzlauf A und Lauf B (Sonde k)."""
    f, A = ref.get_spectrum(k, pad)
    B = sim.get_spectrum(k, pad)[1]
    m = (f > 0) & (np.abs(A) > 0.02 * np.abs(A).max())
    return C0 / (f[m] * D), np.abs(B - A)[m] / np.abs(A[m])


def bei(nl, r, n):
    return np.interp(n, nl[::-1], r[::-1])


_REF, _CACHE = None, {}


def spektrum(npml, r0):
    """|R_Sim| ueber N_lambda fuer eine Dicke und Vorgabe; der Referenzlauf A wird nur einmal gerechnet."""
    global _REF
    if _REF is None:
        _REF = lauf(L_GROSS, 10, 1e-6)
    key = (int(npml), round(float(np.log10(r0)), 6))
    if key not in _CACHE:
        _CACHE[key] = restreflexion(_REF, lauf(L_KLEIN, npml, r0))
    return _CACHE[key]


def senkrecht(npml, r0):
    return bei(*spektrum(npml, r0), N_BEZUG)


def main():
    ABB.mkdir(parents=True, exist_ok=True)

    # Abbildung: links Vorgabe |R(0)| fuer drei Dicken (N_lambda = 10), rechts Dicke fuer vier Aufloesungen
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for npml, farbe in DICKEN:
        werte = np.array([senkrecht(npml, r0) for r0 in R0_WERTE])
        a1.loglog(R0_WERTE, werte, color=farbe, lw=1.6, label=f"{npml} PML-Zellen")
        k = np.argmin(werte)
        print(f"  {npml} Zellen: Minimum {werte[k]:.1e} bei r0 = {R0_WERTE[k]:.0e}, "
              f"bei 1e-6: {werte[np.argmin(abs(np.log(R0_WERTE / 1e-6)))]:.1e}, bei 1e-20: {werte[-1]:.1e}")
    a1.loglog(R0_WERTE, R0_WERTE, ":", color="k", lw=1.1, label=r"$|R_\mathrm{Sim}| = |R(0)|$")
    a1.axvline(1e-6, color="0.5", lw=0.8, ls="--")
    a1.text(1.3e-6, 2.2e-1, "Standardwert", fontsize=8, color="0.35", ha="right")
    a1.set_xlim(1e-1, 1e-20)
    a1.set_ylim(1e-8, 1)
    a1.set_xlabel(r"vorgegebene Restreflexion $|R(0)|$")
    a1.set_ylabel(r"gemessene Reflexion $|R_\mathrm{Sim}|$")
    a1.grid(True, which="major", alpha=0.3)
    a1.legend(fontsize=7, loc="lower right", ncol=2)

    # rechts: kleinste Reflexion je Dicke (beste Vorgabe) fuer vier Aufloesungen, Daten aus testfall_c_karten.py
    from testfall_c_karten import laden, N_PML, NL
    best = laden().min(axis=1)
    # eigene, einfarbige Abstufung: N_lambda ist eine Groesse, keine Kategorie wie die Dicken links
    for nl, farbe in zip((7, 10, 20, 30), plt.cm.Purples(np.linspace(0.45, 1.0, 4))):
        werte = np.array([np.interp(nl, NL, best[i]) for i in range(len(N_PML))])
        a2.semilogy(N_PML, werte, "o-", ms=3.5, lw=1.5, color=farbe, label=rf"$N_\lambda = {nl}$")
        print(f"  N = {nl}: beste Reflexion bei 5, 10, 30 Zellen: "
              + ", ".join(f"{werte[list(N_PML).index(n)]:.1e}" for n in (5, 10, 30)))
    a2.set_xlim(0, 31)
    a2.set_ylim(1e-7, 1)
    a2.grid(True, which="major", alpha=0.3)
    a2.set_xlabel("Dicke der PML in Zellen")
    a2.set_ylabel(r"kleinste gemessene Reflexion $|R_\mathrm{Sim}|$")
    a2.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(ABB / "testfall_c_restreflexion.pdf")
    plt.close(fig)



if __name__ == "__main__":
    main()
    print(f"Abbildungen in {ABB}")
