import numpy as np
from pathlib import Path

from testfall_c import referenz, reflexion

N_PML = np.array([2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 25, 30])
R0 = np.logspace(-1, -10, 10)
NL = np.geomspace(7, 30, 60)
DATEI = Path(__file__).resolve().parent / "testfall_c_karten.npz"


def messen():
    ref = referenz()
    R = np.empty((len(N_PML), len(R0), len(NL)))
    for i, npml in enumerate(N_PML):
        for j, r0 in enumerate(R0):
            nl, r = reflexion(ref, npml, r0)
            R[i, j] = np.interp(NL, nl[::-1], r[::-1])
        print(f"{npml} Zellen fertig")
    np.savez(DATEI, R=R)
    return R


def laden():
    if DATEI.exists():
        return np.load(DATEI)["R"]
    return messen()


if __name__ == "__main__":
    laden()
