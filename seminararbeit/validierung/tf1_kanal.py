"""Testfall 1a --- ebene Welle im Kanal.

Eine Linienquelle regt eine ebene Welle an, die laengs der Gitterachse
laeuft. Weil das Feld nicht von y abhaengt, ist die Anordnung im Kern
eindimensional; die Kruemmung der Wellenfront, die bei der Kreiswelle
die Phasenmessung verfaelscht, entfaellt hier vollstaendig. Dadurch
laesst sich die Phasengeschwindigkeit ohne systematischen Beitrag der
Auswertung messen.
"""
import json, time
import numpy as np
from fdtd_core import FDTD2D, C0
import analytik as A

F0 = 1.0e9
LAM0 = C0 / F0


def kanal(N_lam, L_lam=14.0, H_lam=4.0, S=0.99, d_pml=1.5,
          n_perioden_dft=8, extra=6.0):
    dx = LAM0 / N_lam
    npml = max(8, int(round(d_pml * N_lam)))
    nx = int(round(L_lam * N_lam)) + 2 * npml
    ny = int(round(H_lam * N_lam)) + 2 * npml
    sim = FDTD2D(nx, ny, dx, dx, courant=S, npml=npml)
    i_src = npml + 2
    for j in range(1, ny - 1):
        sim.add_soft_source(i_src, j, lambda t: np.sin(2 * np.pi * F0 * t))

    spp = (1.0 / F0) / sim.dt
    n_settle = int(round((L_lam + extra) * spp))
    n_dft = int(round(n_perioden_dft * spp))

    t0 = time.time()
    sim.run(n_settle)
    jc = ny // 2
    phas = np.zeros(nx, dtype=complex)
    w = np.hanning(n_dft + 2)[1:-1]
    for k in range(n_dft):
        sim.step()
        phas += w[k] * sim.Ez[:, jc] * np.exp(-2j * np.pi * F0 * sim.t)
    phas *= 2.0 / w.sum()
    return sim, phas, i_src, dict(dx=dx, dt=sim.dt, nx=nx, ny=ny, npml=npml,
                                  i_src=i_src, sek=time.time() - t0,
                                  schritte=n_settle + n_dft)


def fit_beta(r, ph):
    Am = np.vstack([r, np.ones_like(r)]).T
    koef, *_ = np.linalg.lstsq(Am, ph, rcond=None)
    rest = ph - Am @ koef
    s2 = np.sum(rest ** 2) / max(len(r) - 2, 1)
    serr = np.sqrt(s2 / np.sum((r - r.mean()) ** 2))
    return -koef[0], serr, rest


def auswerten(phas, i_src, dx, r_von=3.0, r_bis=11.0):
    m = np.arange(phas.size) - i_src
    r_lam = m * dx / LAM0
    sel = (r_lam >= r_von) & (r_lam <= r_bis)
    ph = np.unwrap(np.angle(phas[sel]))
    r = r_lam[sel] * LAM0
    beta, serr, rest = fit_beta(r, ph)
    v = 2 * np.pi * F0 / beta / C0
    return v, v * serr / beta, r_lam[sel], ph, rest


if __name__ == "__main__":
    erg = {}

    print("=== Testfall 1a: ebene Welle, Aufloesungsstudie ===")
    print("  N_lam  Gitter        s     c~p/c gemessen        Theorie      Diff/ppm")
    aufl = []
    for N_lam in (5, 10, 20, 40, 80):
        sim, phas, i_src, meta = kanal(N_lam)
        S_eff = np.sqrt(2) * C0 * meta["dt"] / meta["dx"]
        v, verr, r_lam, ph, rest = auswerten(phas, i_src, meta["dx"])
        th = float(A.cp_achse(N_lam, S_eff))
        aufl.append(dict(N_lam=N_lam, nx=meta["nx"], ny=meta["ny"],
                         schritte=meta["schritte"], sek=meta["sek"], S=S_eff,
                         v=v, v_unsich=verr, v_theorie=th,
                         phasenrest_rad=float(np.std(rest))))
        print(f"  {N_lam:5}  {meta['nx']}x{meta['ny']:<5} {meta['sek']:6.1f}  "
              f"{v:.8f} +- {verr:.1e}  {th:.8f}  {1e6*(v-th):+8.2f}")
    erg["aufloesung"] = aufl

    # Konvergenzordnung
    N = np.array([a["N_lam"] for a in aufl], dtype=float)
    fehl = np.array([1 - a["v"] for a in aufl])
    lx, ly = np.log(1 / N), np.log(fehl)
    Am = np.vstack([lx, np.ones_like(lx)]).T
    koef, *_ = np.linalg.lstsq(Am, ly, rcond=None)
    rest = ly - Am @ koef
    s2 = np.sum(rest ** 2) / (len(lx) - 2)
    perr = np.sqrt(s2 / np.sum((lx - lx.mean()) ** 2))
    erg["ordnung"] = dict(p=float(koef[0]), p_unsich=float(perr),
                          verhaeltnisse=list(fehl[:-1] / fehl[1:]))
    print(f"\n  Konvergenzordnung p = {koef[0]:.3f} +- {perr:.3f}")
    print(f"  Fehlerverhaeltnisse bei Verdopplung: "
          f"{np.array2string(fehl[:-1]/fehl[1:], precision=2)}")

    # ---- Courant-Studie -------------------------------------------------
    print("\n=== Courant-Studie (N_lam = 20, ebene Welle) ===")
    cour = []
    for S in (0.35, 0.50, 0.70, 0.90, 0.99):
        sim, phas, i_src, meta = kanal(20, S=S)
        S_eff = np.sqrt(2) * C0 * meta["dt"] / meta["dx"]
        v, verr, *_ = auswerten(phas, i_src, meta["dx"])
        th = float(A.cp_achse(20, S_eff))
        cour.append(dict(S=S_eff, v=v, v_unsich=verr, v_theorie=th))
        print(f"  S={S_eff:.2f}   gemessen {v:.8f} +- {verr:.1e}   "
              f"Theorie {th:.8f}   Diff {1e6*(v-th):+7.2f} ppm")
    erg["courant"] = cour

    # ---- Verzug ueber lange Laufstrecke ---------------------------------
    print("\n=== Verzug ueber 50 Wellenlaengen ===")
    verz = {}
    for N_lam in (10, 20):
        sim, phas, i_src, meta = kanal(N_lam, L_lam=56.0, H_lam=3.0, extra=8.0)
        dx = meta["dx"]
        m = np.arange(phas.size) - i_src
        r_lam = m * dx / LAM0
        sel = (r_lam >= 2.0) & (r_lam <= 52.0)
        ph = np.unwrap(np.angle(phas[sel]))
        rr = r_lam[sel]
        ph_ex = -2 * np.pi * rr
        ph_ex += ph[0] - ph_ex[0]
        vz = (ph_ex - ph) / (2 * np.pi)
        stg = np.polyfit(rr, vz, 1)[0]
        S_eff = np.sqrt(2) * C0 * meta["dt"] / dx
        soll = 1 / float(A.cp_achse(N_lam, S_eff)) - 1
        verz[N_lam] = dict(r_lam=list(rr), verzug=list(vz), steigung=float(stg),
                           steigung_theorie=float(soll))
        print(f"  N_lam={N_lam}: Steigung {stg:.6f} lam/lam  (Theorie {soll:.6f})   "
              f"halbe Wellenlaenge nach {0.5/stg:.0f} lam")
    erg["verzug"] = verz

    json.dump(erg, open("data/tf1a.json", "w"), indent=1, default=float)
    print("\ngeschrieben: data/tf1a.json")
