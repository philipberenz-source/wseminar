"""Testfall 2 --- PEC-Hohlraumresonator.

Geprueft werden drei voneinander unabhaengige Aussagen:
  (a) das Spektrum enthaelt genau die vorhergesagten Eigenfrequenzen und
      keine zusaetzlichen Linien,
  (b) die zum Leapfrog-Schema gehoerende Energie bleibt erhalten,
  (c) oberhalb der Stabilitaetsgrenze waechst sie exponentiell an.

Zur Energie: Die naheliegende Summe (eps/2)|E^n|^2 + (mu/2)|H^(n+1/2)|^2 ist
keine Erhaltungsgroesse des Schemas, weil elektrisches und magnetisches Feld
um einen halben Zeitschritt versetzt sind. Erhalten ist stattdessen
    W^n = (eps/2)|E^n|^2 + (mu/2) H^(n-1/2) . H^(n+1/2)
Beide Groessen werden gemessen, um den Unterschied sichtbar zu machen.
"""
import json, time
import numpy as np
from fdtd_core import FDTD2D, C0, EPS0, MU0
import analytik as A

LX, LY = 0.600, 0.400
W_PEC, N_PML = 3, 8
T_REC = 3.0e-6
F_MIN, F_MAX = 0.30e9, 1.20e9

MODEN = [(m, n) for m in range(1, 8) for n in range(1, 6)
         if F_MIN < A.f_kontinuum(m, n, LX, LY, C0) < F_MAX]
MODEN.sort(key=lambda mn: A.f_kontinuum(*mn, LX, LY, C0))

# Quell- und Sondenpositionen als Bruchteile der Kantenlaengen. Bewusst
# unsymmetrisch, damit keine Mode auf einer Knotenlinie liegt.
QUELLEN = [(0.31, 0.27), (0.58, 0.74), (0.19, 0.63)]
SONDEN = [(0.63, 0.41), (0.37, 0.83), (0.81, 0.22)]


def aufbau(M, S=0.99):
    dx = LX / M
    nx_frei, ny_frei = M - 1, int(round(LY / dx)) - 1
    nx = nx_frei + 2 * W_PEC + 2 * N_PML
    ny = ny_frei + 2 * W_PEC + 2 * N_PML
    sim = FDTD2D(nx, ny, dx, dx, courant=S, npml=N_PML)
    i0 = j0 = N_PML + W_PEC
    i1, j1 = i0 + nx_frei, j0 + ny_frei
    sim.set_pec_box(i0 - W_PEC, i0, 0, ny)
    sim.set_pec_box(i1, i1 + W_PEC, 0, ny)
    sim.set_pec_box(0, nx, j0 - W_PEC, j0)
    sim.set_pec_box(0, nx, j1, j1 + W_PEC)
    return sim, dx, (i0, i1, j0, j1), (nx_frei + 1) * dx, (ny_frei + 1) * dx


def energie_leapfrog(sim, H_alt):
    u_e = np.sum(0.5 * EPS0 * sim.eps_r * sim.Ez ** 2)
    u_h = 0.5 * MU0 * (np.sum(sim.Hx * H_alt[0]) + np.sum(sim.Hy * H_alt[1]))
    return (u_e + u_h) * sim.dx * sim.dy


def lauf_spektrum(M, S=0.99):
    sim, dx, (i0, i1, j0, j1), Lx_eff, Ly_eff = aufbau(M, S)
    tau = 2.0e-10
    t0 = 4 * tau
    for fx, fy in QUELLEN:
        iq = i0 + int(round(fx * (i1 - i0)))
        jq = j0 + int(round(fy * (j1 - j0)))
        sim.add_soft_source(iq, jq, lambda t: FDTD2D.gaussian_pulse(t, t0, tau))
    ks = [sim.add_probe(i0 + int(round(fx * (i1 - i0))),
                        j0 + int(round(fy * (j1 - j0)))) for fx, fy in SONDEN]
    n_steps = int(round(T_REC / sim.dt))
    t_ = time.time()
    sim.run(n_steps)

    w = np.hanning(n_steps)
    Y = None
    for k in ks:
        y = np.array(sim._probes[k]["Ez"])
        Yk = np.abs(np.fft.rfft(y * w, n_steps * 8))
        Y = Yk if Y is None else Y + Yk
    f = np.fft.rfftfreq(n_steps * 8, sim.dt)
    return dict(M=M, dx=dx, dt=sim.dt, n_steps=n_steps, nx=sim.nx, ny=sim.ny,
                Lx=Lx_eff, Ly=Ly_eff, f=f, Y=Y, sek=time.time() - t_,
                S=np.sqrt(2) * C0 * sim.dt / dx)


def linien_finden(f, Y, schwelle=0.01):
    """Unabhaengige Peaksuche: alle lokalen Maxima ueber der Schwelle."""
    sel = (f > F_MIN) & (f < F_MAX)
    fs, Ys = f[sel], Y[sel] / Y[sel].max()
    ist = (Ys[1:-1] >= Ys[:-2]) & (Ys[1:-1] > Ys[2:]) & (Ys[1:-1] > schwelle)
    idx = np.where(ist)[0] + 1
    linien = []
    for i in idx:
        a, b, c = Ys[i - 1], Ys[i], Ys[i + 1]
        nen = a - 2 * b + c
        d = 0.5 * (a - c) / nen if nen != 0 else 0.0
        linien.append((fs[i] + d * (fs[1] - fs[0]), b))
    # dicht beieinander liegende Kandidaten zusammenfassen
    linien.sort()
    zus = []
    for fk, ak in linien:
        if zus and fk - zus[-1][0] < 0.5e6:
            if ak > zus[-1][1]:
                zus[-1] = (fk, ak)
        else:
            zus.append((fk, ak))
    return zus


if __name__ == "__main__":
    erg = {}

    print("=== Testfall 2: Eigenfrequenzen ===")
    tab_aufl, spektren = [], {}
    for M in (9, 18, 36, 72):
        r = lauf_spektrum(M)
        linien = linien_finden(r["f"], r["Y"])
        zeilen = []
        for (m, n) in MODEN:
            f_k = A.f_kontinuum(m, n, LX, LY, C0)
            f_g = A.f_gitter(m, n, LX, LY, r["dx"], r["dt"], C0)
            # Zuordnung: naechstgelegene unabhaengig gefundene Linie
            if linien:
                fm, am = min(linien, key=lambda x: abs(x[0] - f_g))
            else:
                fm, am = np.nan, 0.0
            zeilen.append(dict(m=m, n=n, f_kont=float(f_k), f_gitter=float(f_g),
                               f_mess=float(fm), amp=float(am),
                               N_lam=float((C0 / f_k) / r["dx"])))
        d_k = np.array([abs(z["f_mess"] / z["f_kont"] - 1) for z in zeilen])
        d_g = np.array([abs(z["f_mess"] / z["f_gitter"] - 1) for z in zeilen])
        N_lam11 = (C0 / A.f_kontinuum(1, 1, LX, LY, C0)) / r["dx"]
        tab_aufl.append(dict(M=M, N_lam=float(N_lam11), dx=r["dx"], dt=r["dt"],
                             n_steps=r["n_steps"], nx=r["nx"], ny=r["ny"],
                             sek=r["sek"], moden=zeilen, n_linien=len(linien),
                             mittel_kont=float(d_k.mean()), max_kont=float(d_k.max()),
                             mittel_gitter=float(d_g.mean()), max_gitter=float(d_g.max()),
                             f_aufloesung=float(1 / T_REC)))
        spektren[M] = r
        print(f"  M={M:3} (N_lam,11={N_lam11:5.1f}, {r['nx']}x{r['ny']}, "
              f"{r['n_steps']} Schritte, {r['sek']:.0f}s): {len(linien)} Linien gefunden "
              f"(erwartet {len(MODEN)});  mittl. Abw. Kontinuum {100*d_k.mean():.4f} %, "
              f"Gitterformel {100*d_g.mean():.2e} %")
    erg["aufloesung"] = tab_aufl

    print("\n  Moden bei M = 18 (N_lam der Grundmode = 20):")
    print("   (m,n) N_lam    f_Kont/MHz  f_Gitter/MHz  f_gem/MHz    Abw.Kont     Abw.Gitter")
    for z in tab_aufl[1]["moden"]:
        print(f"   ({z['m']},{z['n']}) {z['N_lam']:5.1f}  {z['f_kont']/1e6:10.4f}  "
              f"{z['f_gitter']/1e6:10.4f}  {z['f_mess']/1e6:10.4f}  "
              f"{100*(z['f_mess']/z['f_kont']-1):+9.4f} %  "
              f"{1e6*(z['f_mess']/z['f_gitter']-1):+9.2f} ppm")

    Nl = np.array([t["N_lam"] for t in tab_aufl])
    dk = np.array([t["mittel_kont"] for t in tab_aufl])
    kk = np.polyfit(np.log(1 / Nl), np.log(dk), 1)
    erg["konvergenz_ordnung"] = float(kk[0])
    print(f"\n  Konvergenzordnung der Abweichung von der Kontinuumsformel: "
          f"p = {kk[0]:.2f}")
    print(f"  Verhaeltnisse bei Verdopplung: {np.array2string(dk[:-1]/dk[1:], precision=2)}")

    r = spektren[18]
    np.savez("data/tf2_spektrum.npz", f=r["f"][r["f"] < 1.35e9],
             Y=r["Y"][r["f"] < 1.35e9],
             f_kont=np.array([A.f_kontinuum(m, n, LX, LY, C0) for m, n in MODEN]),
             moden=np.array(MODEN))

    # ---- (b) Energieerhaltung ------------------------------------------
    print("\n=== Energieerhaltung (M = 18) ===")
    sim, dx, (i0, i1, j0, j1), *_ = aufbau(18)
    tau, t0 = 2.0e-10, 8.0e-10
    for fx, fy in QUELLEN:
        sim.add_soft_source(i0 + int(round(fx * (i1 - i0))),
                            j0 + int(round(fy * (j1 - j0))),
                            lambda t: FDTD2D.gaussian_pulse(t, t0, tau))
    sim.run(int(round(6 * t0 / sim.dt)))
    T11 = 1.0 / A.f_kontinuum(1, 1, LX, LY, C0)
    n_ges = 200_000
    W_naiv, W_lf = np.empty(n_ges), np.empty(n_ges)
    t_ = time.time()
    for i in range(n_ges):
        H_alt = (sim.Hx.copy(), sim.Hy.copy())
        sim.step()
        W_naiv[i] = sim.total_energy()
        W_lf[i] = energie_leapfrog(sim, H_alt)
    t_ges = n_ges * sim.dt
    stg = np.polyfit(np.arange(n_ges), W_lf / W_lf.mean(), 1)[0]
    abschn = (W_lf.reshape(10, -1).mean(axis=1) / W_lf.mean())
    erg["energie"] = dict(
        n_schritte=n_ges, perioden=float(t_ges / T11), sek=time.time() - t_,
        schwankung_naiv=float((W_naiv.max() - W_naiv.min()) / W_naiv.mean()),
        schwankung_leapfrog=float((W_lf.max() - W_lf.min()) / W_lf.mean()),
        drift_je_1000=float(stg * 1000),
        abschnittsmittel=[float(x) for x in abschn])
    e = erg["energie"]
    print(f"  {n_ges} Schritte = {t_ges/T11:.0f} Perioden der Grundmode "
          f"({e['sek']:.0f} s)")
    print(f"  Schwankung der naiven Summe E^2+H^2 : {100*e['schwankung_naiv']:.2f} %")
    print(f"  Schwankung der Erhaltungsgroesse    : {100*e['schwankung_leapfrog']:.2e} %")
    print(f"  Drift der Ausgleichsgeraden         : {100*e['drift_je_1000']:.2e} % je 1000 Schritte")
    print(f"  zehn Abschnittsmittel: {np.array2string(abschn, precision=10)}")
    np.savez("data/tf2_energie.npz", t=np.arange(0, n_ges, 25) * sim.dt,
             W_naiv=W_naiv[::25] / W_naiv.mean(), W_lf=W_lf[::25] / W_lf.mean(),
             T11=T11)

    # ---- (c) Stabilitaetsgrenze ----------------------------------------
    print("\n=== Stabilitaetsgrenze ===")
    stab, kurven = {}, {}
    for S in (0.9990, 0.9999, 1.0001, 1.0010, 1.0100):
        sim, dx, (i0, i1, j0, j1), *_ = aufbau(18, S=0.99)
        sim.dt = S / (C0 * np.sqrt(1 / dx ** 2 + 1 / dx ** 2))
        sim._update_material_coeffs()
        sim._setup_pml()
        sim.add_soft_source(i0 + 5, j0 + 3,
                            lambda t: FDTD2D.gaussian_pulse(t, 8e-10, 2e-10))
        n = 60000
        W = np.empty(n)
        with np.errstate(over="ignore", invalid="ignore"):
            for i in range(n):
                sim.step()
                W[i] = sim.total_energy()
                if not np.isfinite(W[i]) or W[i] > 1e60:
                    W[i:] = np.inf
                    break
        endlich = np.isfinite(W)
        W0 = W[2000]
        Wn = W / W0
        # Wachstumsrate aus dem Bereich, in dem W noch endlich und > 10 ist
        m = endlich & (Wn > 10)
        if m.sum() > 50:
            idx = np.where(m)[0]
            rate = np.polyfit(idx, np.log(Wn[idx]), 1)[0]
        else:
            idx = np.arange(2000, min(n, endlich.sum()))
            rate = np.polyfit(idx, np.log(np.maximum(Wn[idx], 1e-300)), 1)[0]
        stab[S] = dict(rate=float(rate), n_endlich=int(endlich.sum()),
                       ende=float(Wn[endlich][-1]))
        kurven[S] = Wn[::20]
        print(f"  S = {S:<7}: nach {int(endlich.sum())} Schritten W/W0 = "
              f"{Wn[endlich][-1]:.3e},  Wachstum {rate:+.5f} je Schritt")
    erg["stabilitaet"] = {str(k): v for k, v in stab.items()}
    np.savez("data/tf2_stabilitaet.npz",
             **{f"S_{str(k).replace('.', 'p')}": v for k, v in kurven.items()})

    json.dump(erg, open("data/tf2.json", "w"), indent=1, default=float)
    print("\ngeschrieben: data/tf2.json")
