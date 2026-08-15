"""Testfall 1b --- Kreiswelle im Vakuum.

Die Kreiswelle ist die einzige Anordnung, in der sich das Feld in alle
Richtungen zugleich ausbreitet. Geprueft werden daher (i) die
Richtungsabhaengigkeit der Phasengeschwindigkeit, (ii) das
Amplitudengesetz und (iii) das gesamte Feld gegen die geschlossene
Loesung mit der Hankel-Funktion.

Wichtig fuer die Auswertung: Auf einer gekruemmten Wellenfront ist die
oertliche Phasensteigung nicht exakt k, sondern enthaelt einen Zusatz der
Ordnung 1/(k r). Damit dieser Beitrag der Auswertung nicht als Fehler des
Verfahrens erscheint, wird dieselbe Auswertung zur Kontrolle auch auf die
exakte Loesung angewandt (Funktion `referenzlauf`).
"""
import json, time
import numpy as np
from fdtd_core import FDTD2D, C0
import analytik as A
from tf1_kanal import fit_beta

F0 = 1.0e9
LAM0 = C0 / F0
R_AUS = 8.5
D_PML = 1.5
R_PHASE = (3.0, 7.0)
R_AMP = (1.5, 7.5)
K0 = 2 * np.pi / LAM0


def lauf(N_lam, S=0.99, n_perioden_dft=8, extra=4.0):
    dx = LAM0 / N_lam
    npml = max(8, int(round(D_PML * N_lam)))
    n = 2 * int(round(R_AUS * N_lam)) + 1 + 2 * npml
    sim = FDTD2D(n, n, dx, dx, courant=S, npml=npml)
    c_ = n // 2
    sim.add_soft_source(c_, c_, lambda t: np.sin(2 * np.pi * F0 * t))
    spp = (1.0 / F0) / sim.dt
    n_settle = int(round((np.sqrt(2) * R_AUS + extra) * spp))
    n_dft = int(round(n_perioden_dft * spp))
    t0 = time.time()
    sim.run(n_settle)
    phas = np.zeros((n, n), dtype=complex)
    w = np.hanning(n_dft + 2)[1:-1]
    for k in range(n_dft):
        sim.step()
        phas += w[k] * sim.Ez * np.exp(-2j * np.pi * F0 * sim.t)
    phas *= 2.0 / w.sum()
    return phas, c_, dict(N_lam=N_lam, n=n, npml=npml, dx=dx, dt=sim.dt,
                          S=np.sqrt(2) * C0 * sim.dt / dx,
                          schritte=n_settle + n_dft, sek=time.time() - t0)


def _keys(t):
    t2, t3 = t * t, t * t * t
    return np.stack([-0.5 * t3 + t2 - 0.5 * t, 1.5 * t3 - 2.5 * t2 + 1.0,
                     -1.5 * t3 + 2.0 * t2 + 0.5 * t, 0.5 * t3 - 0.5 * t2])


def strahl(phas, c_, dx, theta, r_lam):
    xi = c_ + r_lam * LAM0 * np.cos(theta) / dx
    yj = c_ + r_lam * LAM0 * np.sin(theta) / dx
    i0, j0 = np.floor(xi).astype(int), np.floor(yj).astype(int)
    wx, wy = _keys(xi - i0), _keys(yj - j0)
    out = np.zeros(r_lam.shape, dtype=complex)
    for a in range(4):
        for b in range(4):
            out += wx[a] * wy[b] * phas[i0 - 1 + a, j0 - 1 + b]
    return out


# Ganzzahlige Richtungen (p, q): die Punkte (c + m p, c + m q) sind exakte
# Gitterpunkte, sodass die Winkelstudie ohne jede Interpolation auskommt.
RICHTUNGEN = [(1, 0), (8, 1), (5, 1), (4, 1), (3, 1), (5, 2), (2, 1),
              (5, 3), (3, 2), (4, 3), (5, 4), (7, 6), (1, 1)]


def abtast_m(p, q, dx):
    """Vielfache m, fuer die m*sqrt(p^2+q^2)*dx im Auswertebereich liegt."""
    schritt = np.hypot(p, q) * dx / LAM0
    m = np.arange(1, int(R_PHASE[1] / schritt) + 1)
    r = m * schritt
    sel = (r >= R_PHASE[0]) & (r <= R_PHASE[1])
    return m[sel], r[sel]


def v_gemessen(phas, c_, dx, p, q):
    """Phasengeschwindigkeit in Richtung (p, q).

    Die Anordnung ist achtfach symmetrisch; gemittelt wird deshalb ueber
    alle acht gleichwertigen Strahlen, was die Streuung entsprechend senkt.
    """
    m, r_lam = abtast_m(p, q, dx)
    z = np.zeros(m.shape, dtype=complex)
    aeq = {(p, q), (q, p), (-p, q), (q, -p), (p, -q), (-q, p), (-p, -q), (-q, -p)}
    for a_, b_ in aeq:
        z = z + phas[c_ + m * a_, c_ + m * b_]
    z /= len(aeq)
    beta, serr, _ = fit_beta(r_lam * LAM0, np.unwrap(np.angle(z)))
    v = 2 * np.pi * F0 / beta / C0
    return v, v * serr / beta


def v_referenz(dx, p, q):
    """Dieselbe Auswertung auf die exakte Loesung H_0^(2)(k r).

    Ohne Kruemmungseinfluss ergaebe sich exakt 1; der Rueckgabewert ist
    also der systematische Beitrag der Auswertung selbst.
    """
    _, r_lam = abtast_m(p, q, dx)
    z = np.conj(A.hankel1_0(K0 * r_lam * LAM0))
    beta, _, _ = fit_beta(r_lam * LAM0, np.unwrap(np.angle(z)))
    return 2 * np.pi * F0 / beta / C0


if __name__ == "__main__":
    erg = {}
    phasoren = {}
    for N_lam in (10, 20, 40):
        phas, c_, meta = lauf(N_lam)
        phasoren[N_lam] = (phas, c_, meta)
        print(f"  Lauf N_lam={N_lam}: Gitter {meta['n']}^2, "
              f"{meta['schritte']} Schritte, {meta['sek']:.1f} s")

    # ---- (a) Richtungsabhaengigkeit ------------------------------------
    print("\n=== Richtungsabhaengigkeit der Phasengeschwindigkeit (N_lam = 20) ===")
    print("  (p,q)   Winkel   gemessen   Kontrolle  korrigiert  Theorie   Diff/ppm")
    phas, c_, meta = phasoren[20]
    dx, S = meta["dx"], meta["S"]
    tab = []
    for p, q in RICHTUNGEN:
        th = np.arctan2(q, p)
        v, verr = v_gemessen(phas, c_, dx, p, q)
        vref = v_referenz(dx, p, q)
        vkorr = v / vref
        vth = float(A.cp_verhaeltnis(20.0, S, th))
        tab.append(dict(p=p, q=q, winkel=float(np.rad2deg(th)), roh=v, unsich=verr,
                        kontrolle=vref, korrigiert=vkorr, theorie=vth))
        print(f"  ({p},{q})   {np.rad2deg(th):5.2f}   {v:.6f}   {vref:.6f}   "
              f"{vkorr:.6f}   {vth:.6f}   {1e6*(vkorr-vth):+7.1f}")
    erg["richtung"] = tab
    d = np.array([abs(t["korrigiert"] - t["theorie"]) for t in tab])
    print(f"  groesste Abweichung von der Theorie: {1e6*d.max():.0f} ppm, "
          f"Streuung {1e6*d.std():.0f} ppm")
    erg["richtung_maxabw_ppm"] = float(1e6 * d.max())

    fa, fd = 1 - tab[0]["theorie"], 1 - tab[-1]["theorie"]
    erg["anisotropie"] = dict(fehler_achse=float(fa), fehler_diagonal=float(fd),
                              verhaeltnis=float(fa / fd))
    print(f"  Fehler Achse {100*fa:.4f} %, Diagonale {100*fd:.4f} %, "
          f"Verhaeltnis {fa/fd:.0f}")

    # ---- (b) Amplitudengesetz ------------------------------------------
    print("\n=== Amplitudengesetz ===")
    amp_tab = []
    for N_lam in (10, 20, 40):
        phas, c_, meta = phasoren[N_lam]
        dx, n = meta["dx"], meta["n"]
        ii = (np.arange(n) - c_) * dx
        X, Y = np.meshgrid(ii, ii, indexing="ij")
        R = np.hypot(X, Y) / LAM0
        # Ringmittel: alle Zellen eines Radiusintervalls, ohne Interpolation
        kanten = np.arange(R_AMP[0], R_AMP[1] + 1e-9, 0.1)
        r_mit, amp = [], []
        for a_, b_ in zip(kanten[:-1], kanten[1:]):
            sel = (R >= a_) & (R < b_)
            if sel.sum() > 8:
                r_mit.append(0.5 * (a_ + b_))
                amp.append(np.abs(phas[sel]).mean())
        r_mit, amp = np.array(r_mit), np.array(amp)
        lr = np.log(r_mit)
        Am = np.vstack([lr, np.ones_like(lr)]).T
        koef, *_ = np.linalg.lstsq(Am, np.log(amp), rcond=None)
        rest = np.log(amp) - Am @ koef
        s2 = np.sum(rest ** 2) / (len(lr) - 2)
        serr = np.sqrt(s2 / np.sum((lr - lr.mean()) ** 2))
        soll = -np.linalg.lstsq(Am, np.log(np.abs(A.hankel1_0(K0 * r_mit * LAM0))),
                                rcond=None)[0][0]
        amp_tab.append(dict(N_lam=N_lam, p=float(-koef[0]), p_unsich=float(serr),
                            p_soll=float(soll)))
        print(f"  N_lam={N_lam:3}: p = {-koef[0]:.4f} +- {serr:.4f}   "
              f"(exakte Loesung, gleicher Bereich: {soll:.4f};  "
              f"reines Potenzgesetz: 0.5)")
        if N_lam == 20:
            np.savez("data/tf1_amplitude.npz", r_lam=r_mit, amp=amp,
                     p=-koef[0], perr=serr, p_soll=soll)
    erg["amplitude"] = amp_tab

    # ---- (c) Differenzfeld ---------------------------------------------
    print("\n=== Differenzfeld gegen die exakte Loesung ===")
    R0 = 2.0                      # Radius, an dem die Quellstaerke angepasst wird
    dtab, profile = [], {}
    for N_lam in (10, 20, 40):
        phas, c_, meta = phasoren[N_lam]
        dx, n, S = meta["dx"], meta["n"], meta["S"]
        ii = (np.arange(n) - c_) * dx
        X, Y = np.meshgrid(ii, ii, indexing="ij")
        R = np.hypot(X, Y)
        Rl = R / LAM0
        gueltig = (Rl > 0.8) & (Rl < 7.5)
        H = np.zeros_like(phas)
        H[gueltig] = np.conj(A.hankel1_0(K0 * R[gueltig]))
        ring0 = (Rl > R0 - 0.05) & (Rl < R0 + 0.05)
        skal = np.vdot(H[ring0], phas[ring0]) / np.vdot(H[ring0], H[ring0])
        dfeld = np.zeros_like(phas)
        dfeld[gueltig] = phas[gueltig] - skal * H[gueltig]

        kanten = np.arange(2.0, 7.01, 0.25)
        rm, ef = [], []
        for a_, b_ in zip(kanten[:-1], kanten[1:]):
            sel = (Rl >= a_) & (Rl < b_)
            rm.append(0.5 * (a_ + b_))
            ef.append(np.sqrt((np.abs(dfeld[sel]) ** 2).mean())
                      / np.sqrt((np.abs(skal * H[sel]) ** 2).mean()))
        rm, ef = np.array(rm), np.array(ef)
        # Vorhersage: Phasenverzug seit R0, gemittelt ueber alle Richtungen
        thv = np.linspace(0, np.pi / 4, 46)
        dphi = np.array([2 * np.pi * (rm - R0) * (1 / float(A.cp_verhaeltnis(N_lam, S, t)) - 1)
                         for t in thv])
        vor = np.sqrt((np.abs(2 * np.sin(dphi / 2)) ** 2).mean(axis=0))
        profile[N_lam] = dict(r_lam=list(rm), gemessen=list(ef), vorhersage=list(vor))
        i6 = int(np.argmin(np.abs(rm - 6.0)))
        dtab.append(dict(N_lam=N_lam, ring6=float(ef[i6]), ring6_vorhersage=float(vor[i6]),
                         max_rel=float(np.abs(dfeld[gueltig]).max()
                                       / np.abs(skal * H[gueltig]).max())))
        print(f"  N_lam={N_lam:3}: bei r=6 lam  gemessen {100*ef[i6]:6.2f} %   "
              f"vorhergesagt {100*vor[i6]:6.2f} %")
        np.save(f"data/tf1_diff_N{N_lam}.npy",
                (np.abs(dfeld) / np.abs(skal * H[gueltig]).max() * gueltig).astype(np.float32))
    v = np.array([d["ring6"] for d in dtab])
    print(f"  Verhaeltnisse bei Verdopplung der Aufloesung: "
          f"{np.array2string(v[:-1]/v[1:], precision=2)}")
    erg["differenzfeld"] = dtab
    erg["differenzfeld_profil"] = profile

    # Rohdaten fuer die Abbildungen
    phas, c_, meta = phasoren[20]
    np.save("data/tf1_feld_N20.npy", np.real(phas).astype(np.float32))
    json.dump(erg, open("data/tf1b.json", "w"), indent=1, default=float)
    print("\ngeschrieben: data/tf1b.json")
