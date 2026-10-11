import numpy as np
from pathlib import Path

from testfall_c import lauf, restreflexion, L_GROSS, L_KLEIN, ABSTAND

# Testfall 4, Messreihe: gemessene Reflexion |R_Sim| bei senkrechtem Einfall in Abhaengigkeit von
# PML-Dicke n (Zellen), Vorgabe |R(0)| und Aufloesung N_lambda. Ein Ricker-Lauf liefert das ganze
# Spektrum, also alle N_lambda zugleich; Messung und Wegkorrektur wie in testfall_c.py.
N_PML = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 22, 25, 28, 30])
R0 = np.logspace(-1, -20, 39)
NL = np.geomspace(7, 30, 60)
DATEI = Path(__file__).resolve().parent / "testfall_c_karten.npz"


def messen():
    ref = lauf(L_GROSS, 10, 1e-6)
    R = np.empty((len(N_PML), len(R0), len(NL)))
    for i, n in enumerate(N_PML):
        korr = np.sqrt((ABSTAND + 2 * (5 + n)) / ABSTAND)
        for j, r0 in enumerate(R0):
            nl, r = restreflexion(ref, lauf(L_KLEIN, n, r0))
            R[i, j] = np.interp(NL, nl[::-1], r[::-1]) * korr
        print(f"  {n} Zellen fertig", flush=True)
    np.savez(DATEI, R=R)
    return R


def laden():
    return np.load(DATEI)["R"] if DATEI.exists() else messen()


if __name__ == "__main__":
    laden()
    print(f"Messdaten in {DATEI}")
