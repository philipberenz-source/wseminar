import numpy as np
from pathlib import Path

from testfall_c import spektrum

# Testfall 4, Messreihe: gemessene Reflexion |R_Sim| bei senkrechtem Einfall in Abhaengigkeit von
# PML-Dicke n (Zellen), Vorgabe |R(0)| und Aufloesung N_lambda. Ein Ricker-Lauf liefert das ganze
# Spektrum, also alle N_lambda zugleich; Messung wie in testfall_c.py (ebene Welle).
N_PML = np.array([2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 25, 30])
R0 = np.logspace(-1, -10, 10)            # beste Vorgabe liegt stets in diesem Bereich
NL = np.geomspace(7, 30, 60)
DATEI = Path(__file__).resolve().parent / "testfall_c_karten.npz"


def messen():
    R = np.empty((len(N_PML), len(R0), len(NL)))
    for i, n in enumerate(N_PML):
        for j, r0 in enumerate(R0):
            nl, r = spektrum(n, r0)
            R[i, j] = np.interp(NL, nl[::-1], r[::-1])
        print(f"  {n} Zellen fertig", flush=True)
    np.savez(DATEI, R=R)
    return R


def laden():
    return np.load(DATEI)["R"] if DATEI.exists() else messen()


if __name__ == "__main__":
    laden()
    print(f"Messdaten in {DATEI}")
