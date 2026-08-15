"""Testfall 3 --- Reflexion an einer Materialgrenzflaeche.

Ausgewertet wird mit dem Differenzverfahren: Derselbe Aufbau wird zweimal
gerechnet, einmal mit und einmal ohne Materialsprung. Die Differenz beider
Aufzeichnungen an derselben Sonde ist exakt das reflektierte Feld, ohne
dass einfallende und reflektierte Welle im Zeitsignal getrennt werden
muessten. Aus dem Verhaeltnis der Fourier-Transformierten folgt der
Reflexionsgrad -- und zwar fuer jede im Puls enthaltene Frequenz zugleich.
Da eine Frequenz zugleich eine Aufloesung bedeutet (N_lambda = c/(f n dx)),
liefert ein einziger Lauf die vollstaendige Aufloesungsstudie.
"""
import json, time
import numpy as np
from fdtd_core import FDTD2D, C0, EPS0
import analytik as A

F0 = 1.0e9
LAM0 = C0 / F0


def lauf(eps2=1.0, sigma2=0.0, N_lam2=20.0, eps_ref=4.0, L_vor=12.0, L_nach=8.0,
         H_lam=3.0, d_pml=1.5, S=0.99):
    """Kanal mit Grenzflaeche. N_lam2 bezieht sich auf die Wellenlaenge im
    Referenzmaterial eps_ref bei F0 -- damit haben alle Laeufe einer Reihe
    dieselbe Zellweite."""
    lam2 = LAM0 / np.sqrt(eps_ref)
    dx = lam2 / N_lam2
    npml = max(8, int(round(d_pml * LAM0 / dx)))
    nx = int(round((L_vor + L_nach) * LAM0 / dx)) + 2 * npml
    ny = int(round(H_lam * LAM0 / dx)) + 2 * npml
    sim = FDTD2D(nx, ny, dx, dx, courant=S, npml=npml)
    i_src = npml + 2
    i_int = i_src + int(round(L_vor * LAM0 / dx))
    if eps2 != 1.0 or sigma2 != 0.0:
        sim.set_material_box(i_int, nx, 0, ny, eps_r=eps2, sigma=sigma2)
    tau, t0 = 1.6e-10, 6.4e-10
    for j in range(1, ny - 1):
        sim.add_soft_source(i_src, j, lambda t: FDTD2D.gaussian_pulse(t, t0, tau))
    jc = ny // 2
    i_p1 = i_src + int(round(0.45 * (i_int - i_src)))      # vor der Grenzflaeche
    i_p2 = i_int + 1                                        # dahinter
    k1 = sim.add_probe(i_p1, jc)
    k2 = sim.add_probe(i_p2, jc)
    n = int(round((L_vor + L_nach + 4) * LAM0 / C0 / sim.dt))
    sim.run(n)
    return (np.array(sim._probes[k1]["Ez"]), np.array(sim._probes[k2]["Ez"]),
            dict(dx=dx, dt=sim.dt, nx=nx, ny=ny, npml=npml, n=n,
                 i_src=i_src, i_int=i_int, i_p1=i_p1, i_p2=i_p2,
                 S=np.sqrt(2) * C0 * sim.dt / dx))


def spektren(y, dt, f_band):
    Y = np.fft.rfft(y)
    f = np.fft.rfftfreq(len(y), dt)
    sel = (f >= f_band[0]) & (f <= f_band[1])
    return f[sel], Y[sel]


if __name__ == "__main__":
    erg = {}
    F_BAND = (0.35e9, 2.4e9)
    N_LAM2 = 20.0
    EPS_REF = 4.0

    # Referenzlauf ohne Materialsprung
    y1_ref, y2_ref, meta = lauf(eps2=1.0, N_lam2=N_LAM2, eps_ref=EPS_REF)
    dt, dx = meta["dt"], meta["dx"]
    f, Y_inc = spektren(y1_ref, dt, F_BAND)
    print(f"Gitter {meta['nx']}x{meta['ny']}, dx = {dx*1000:.3f} mm, "
          f"{meta['n']} Zeitschritte")
    print(f"Auswertband {F_BAND[0]/1e9:.2f}-{F_BAND[1]/1e9:.2f} GHz  "
          f"entspricht N_lambda im Medium (eps_r=4) von "
          f"{(C0/F_BAND[1]/2)/dx:.1f} bis {(C0/F_BAND[0]/2)/dx:.1f}")

    # ---- Nachweisgrenze der Anordnung ----------------------------------
    # Im Referenzlauf ohne Materialsprung darf nach dem Durchgang des Pulses
    # kein Signal mehr an der Sonde ankommen. Was dennoch eintrifft, stammt
    # vom absorbierenden Rand und begrenzt die Aufloesung des Verfahrens.
    i_max = int(np.argmax(np.abs(y1_ref)))
    spaet = y1_ref[i_max + int(round(3 * LAM0 / C0 / dt)):]
    boden = np.abs(spaet).max() / np.abs(y1_ref).max()
    print(f"\nNachweisgrenze (Restsignal im Referenzlauf nach Durchgang des "
          f"Pulses): {boden:.2e}")
    erg["nachweisgrenze"] = float(boden)

    # ---- Verlustfreie Dielektrika --------------------------------------
    print("\n=== Reflexionsgrad, verlustfreie Dielektrika (bei 1 GHz) ===")
    print("  eps_r   N_lam,2   |r| gemessen   Fresnel     Abweichung   R+T-1")
    tab, kurven = [], {}
    i_f0 = int(np.argmin(np.abs(f - F0)))
    for eps2 in (1.5, 2.25, 4.0, 6.0, 9.0):
        y1, y2, _ = lauf(eps2=eps2, N_lam2=N_LAM2, eps_ref=EPS_REF)
        _, Y_ref = spektren(y1 - y1_ref, dt, F_BAND)
        _, Y_tr = spektren(y2, dt, F_BAND)
        r = np.abs(Y_ref / Y_inc)
        tt = np.abs(Y_tr / Y_inc)
        r_soll = abs(A.fresnel_r(1.0, eps2))
        n2 = np.sqrt(eps2)
        bilanz = r ** 2 + n2 * tt ** 2 - 1
        # Vorhersage aus der numerischen Dispersion: wirksame Brechungsindizes
        S_ = meta["S"]
        N1 = (C0 / F0) / dx
        N2 = (C0 / (F0 * n2)) / dx
        n1e = 1.0 / float(A.cp_achse(N1, S_))
        n2e = n2 / float(A.cp_achse(N2, S_ / n2))
        r_disp = abs((n1e - n2e) / (n1e + n2e))
        N_lam2_f = (C0 / (f * n2)) / dx
        kurven[eps2] = dict(f=list(f), r=list(r), N_lam2=list(N_lam2_f),
                            r_soll=float(r_soll))
        tab.append(dict(eps2=eps2, N_lam2=float(N_lam2_f[i_f0]),
                        r=float(r[i_f0]), r_soll=float(r_soll),
                        abw=float(r[i_f0] / r_soll - 1),
                        r_disp=float(r_disp), abw_disp=float(r_disp / r_soll - 1),
                        bilanz=float(bilanz[i_f0])))
        print(f"  {eps2:5.2f}   {N_lam2_f[i_f0]:6.1f}   {r[i_f0]:.6f}     "
              f"{r_soll:.6f}   {100*(r[i_f0]/r_soll-1):+8.3f} %  "
              f"{bilanz[i_f0]:+.2e}   (aus Dispersion allein: "
              f"{100*(r_disp/r_soll-1):+.3f} %)")
    erg["dielektrika"] = tab

    # ---- Aufloesungsstudie aus der Frequenzabhaengigkeit ----------------
    print("\n=== Aufloesungsstudie (eps_r = 4, aus demselben Lauf) ===")
    k = kurven[4.0]
    fk = np.array(k["f"]); rk = np.array(k["r"]); Nk = np.array(k["N_lam2"])
    print("  N_lam,2   |r|         Abweichung von Fresnel")
    aufl = []
    for Nz in (10, 15, 20, 30, 40):
        i = int(np.argmin(np.abs(Nk - Nz)))
        aufl.append(dict(N_lam2=float(Nk[i]), f=float(fk[i]), r=float(rk[i]),
                         abw=float(rk[i] / k["r_soll"] - 1)))
        print(f"  {Nk[i]:7.1f}   {rk[i]:.6f}    {100*(rk[i]/k['r_soll']-1):+8.4f} %")
    N_arr = np.array([a["N_lam2"] for a in aufl])
    d_arr = np.abs([a["abw"] for a in aufl])
    p = np.polyfit(np.log(1 / N_arr), np.log(d_arr), 1)[0]
    erg["aufloesung"] = dict(punkte=aufl, ordnung=float(p))
    print(f"  Konvergenzordnung: p = {p:.2f}")

    # ---- Verlustbehaftete Materialien ----------------------------------
    print("\n=== Verlustbehaftete Materialien (eps_r = 4, bei 1 GHz) ===")
    print("  sigma/(S/m)  sigma/(w eps)  |r| gemessen   |r| Fresnel   Abweichung")
    verl = []
    for sig in (0.01, 0.05, 0.20, 1.00):
        y1, y2, _ = lauf(eps2=4.0, sigma2=sig, N_lam2=N_LAM2, eps_ref=EPS_REF)
        _, Y_ref = spektren(y1 - y1_ref, dt, F_BAND)
        r = np.abs(Y_ref / Y_inc)
        r_soll = abs(A.fresnel_r(1.0, 4.0, sigma2=sig, f=F0))
        vt = sig / (2 * np.pi * F0 * EPS0 * 4.0)
        verl.append(dict(sigma=sig, verlusttangens=float(vt), r=float(r[i_f0]),
                         r_soll=float(r_soll), abw=float(r[i_f0] / r_soll - 1)))
        print(f"  {sig:10.2f}   {vt:11.3f}    {r[i_f0]:.6f}     {r_soll:.6f}   "
              f"{100*(r[i_f0]/r_soll-1):+8.3f} %")
    erg["verluste"] = verl

    # ---- Vergleich mit dem Stehwellenverfahren -------------------------
    print("\n=== Vergleichsverfahren: stehende Welle (harmonische Anregung) ===")
    swr = []
    for eps2 in (1.5, 2.25, 4.0, 6.0, 9.0):
        lam2 = LAM0 / np.sqrt(EPS_REF)
        dxs = lam2 / N_LAM2
        npml = max(8, int(round(1.5 * LAM0 / dxs)))
        nx = int(round(20 * LAM0 / dxs)) + 2 * npml
        ny = int(round(3 * LAM0 / dxs)) + 2 * npml
        sim = FDTD2D(nx, ny, dxs, dxs, courant=0.99, npml=npml)
        i_src = npml + 2
        i_int = i_src + int(round(12 * LAM0 / dxs))
        sim.set_material_box(i_int, nx, 0, ny, eps_r=eps2)
        rampe = 4.0 / F0
        for j in range(1, ny - 1):
            sim.add_soft_source(i_src, j, lambda t: np.sin(2 * np.pi * F0 * t)
                                * min(1.0, t / rampe))
        spp = (1 / F0) / sim.dt
        sim.run(int(round(30 * spp)))
        jc = ny // 2
        env = np.zeros(nx)
        sim.run(int(round(4 * spp)),
                callback=lambda s: np.maximum(env, np.abs(s.Ez[:, jc]), out=env))
        a = env[i_src + int(round(3 * LAM0 / dxs)): i_int - int(round(0.5 * LAM0 / dxs))]
        rr = (a.max() - a.min()) / (a.max() + a.min())
        r_soll = abs(A.fresnel_r(1.0, eps2))
        swr.append(dict(eps2=eps2, r=float(rr), r_soll=float(r_soll),
                        abw=float(rr / r_soll - 1)))
        print(f"  eps_r={eps2:5.2f}:  |r| = {rr:.5f}  (Fresnel {r_soll:.5f}, "
              f"{100*(rr/r_soll-1):+7.3f} %)")
        if eps2 == 4.0:
            np.savez("data/tf3_stehwelle.npz", x=(np.arange(nx) - i_int) * dxs / LAM0,
                     env=env / env.max(), i_int=i_int, dx=dxs, eps2=eps2,
                     r_soll=r_soll)
    erg["stehwelle"] = swr

    # ---- Rohdaten fuer Abbildungen -------------------------------------
    y1, y2, _ = lauf(eps2=4.0, N_lam2=N_LAM2, eps_ref=EPS_REF)
    np.savez("data/tf3_zeit.npz", t=np.arange(len(y1)) * dt,
             ges=y1, inc=y1_ref, ref=y1 - y1_ref, trans=y2, dt=dt)
    np.savez("data/tf3_kurven.npz",
             **{f"eps{str(e).replace('.','p')}_{k2}": np.array(v[k2])
                for e, v in kurven.items() for k2 in ("f", "r", "N_lam2")},
             eps_liste=np.array(list(kurven.keys())),
             r_soll=np.array([kurven[e]["r_soll"] for e in kurven]))

    json.dump(erg, open("data/tf3.json", "w"), indent=1, default=float)
    print("\ngeschrieben: data/tf3.json")
