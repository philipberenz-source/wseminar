"""Testfall 2b --- Energieerhaltung und Stabilitaetsgrenze.

Zur Energie: Die naheliegende Summe (eps/2)|E^n|^2 + (mu/2)|H^(n+1/2)|^2 ist
keine Erhaltungsgroesse des Leapfrog-Schemas, weil beide Felder um einen
halben Zeitschritt versetzt sind. Erhalten ist stattdessen

    W^n = (eps/2) |E^n|^2 + (mu/2) H^(n-1/2) . H^(n+1/2)

Beide Groessen werden gemessen, um den Unterschied sichtbar zu machen.

Zur Stabilitaet: Die CFL-Bedingung gilt fuer das unendliche Gitter, wo
Wellenzahlen bis zur Gitter-Nyquist-Grenze vorkommen. In einem endlichen
Resonator ist das Spektrum diskret und erreicht diese Grenze nicht ganz;
die tatsaechliche Grenze liegt deshalb geringfuegig ueber S = 1. Der
vorhergesagte Wert wird unten aus dem Modenspektrum berechnet und
gemessen.
"""
import json, time
import numpy as np
from fdtd_core import FDTD2D, C0, EPS0, MU0
import analytik as A
from tf2_resonator import LX, LY, QUELLEN, aufbau

M = 18


def s_krit(nx_frei, ny_frei):
    ux = np.arange(1, nx_frei + 1) * np.pi / (2 * (nx_frei + 1))
    uy = np.arange(1, ny_frei + 1) * np.pi / (2 * (ny_frei + 1))
    return float(np.sqrt(2.0 / (np.sin(ux) ** 2).max() + 0 * 1)) if False else \
        float(np.sqrt(2.0 / ((np.sin(ux) ** 2).max() + (np.sin(uy) ** 2).max())))


if __name__ == "__main__":
    erg = {}

    # =============================== Energieerhaltung ====================
    print("=== Energieerhaltung (M = 18) ===")
    sim, dx, (i0, i1, j0, j1), Lx_eff, Ly_eff = aufbau(M)
    tau, t0 = 2.0e-10, 8.0e-10
    for fx, fy in QUELLEN:
        sim.add_soft_source(i0 + int(round(fx * (i1 - i0))),
                            j0 + int(round(fy * (j1 - j0))),
                            lambda t: FDTD2D.gaussian_pulse(t, t0, tau))
    sim.run(int(round(8 * t0 / sim.dt)))          # Anregung vollstaendig abklingen lassen

    T11 = 1.0 / A.f_kontinuum(1, 1, LX, LY, C0)
    n_ges = 300_000
    W_naiv = np.empty(n_ges)
    W_lf = np.empty(n_ges)
    flaeche = sim.dx * sim.dy
    t_ = time.time()
    for i in range(n_ges):
        E_alt = sim.Ez.copy()                     # E^n
        Hx_alt, Hy_alt = sim.Hx.copy(), sim.Hy.copy()   # H^(n-1/2)
        sim.step()                                # -> H^(n+1/2), E^(n+1)
        W_naiv[i] = sim.total_energy()
        u_e = np.sum(0.5 * EPS0 * sim.eps_r * E_alt ** 2)
        u_h = 0.5 * MU0 * (np.sum(sim.Hx * Hx_alt) + np.sum(sim.Hy * Hy_alt))
        W_lf[i] = (u_e + u_h) * flaeche

    t_ges = n_ges * sim.dt
    mw = W_lf.mean()
    stg = np.polyfit(np.arange(n_ges), W_lf / mw, 1)[0]
    abschn = W_lf.reshape(10, -1).mean(axis=1) / mw
    erg["energie"] = dict(
        M=M, n_schritte=n_ges, perioden=float(t_ges / T11), sek=time.time() - t_,
        dt=float(sim.dt), t_ges=float(t_ges),
        schwankung_naiv=float((W_naiv.max() - W_naiv.min()) / W_naiv.mean()),
        schwankung_leapfrog=float((W_lf.max() - W_lf.min()) / mw),
        streuung_leapfrog=float(W_lf.std() / mw),
        drift_je_1000=float(stg * 1000),
        abschnittsmittel=[float(x) for x in abschn])
    e = erg["energie"]
    print(f"  {n_ges} Schritte = {t_ges/T11:.0f} Perioden der Grundmode ({e['sek']:.0f} s)")
    print(f"  Schwankung der naiven Summe E^2 + H^2 : {100*e['schwankung_naiv']:8.3f} %")
    print(f"  Schwankung der Erhaltungsgroesse      : {100*e['schwankung_leapfrog']:8.2e} %")
    print(f"  Standardabweichung derselben          : {100*e['streuung_leapfrog']:8.2e} %")
    print(f"  Drift der Ausgleichsgeraden           : {100*e['drift_je_1000']:8.2e} % je 1000 Schritte")
    print(f"  zehn Abschnittsmittel: {np.array2string(abschn, precision=12)}")
    np.savez("data/tf2_energie.npz", t=np.arange(0, n_ges, 40) * sim.dt,
             W_naiv=W_naiv[::40] / W_naiv.mean(), W_lf=W_lf[::40] / mw, T11=T11)
    # Fuer die Abbildung zusaetzlich die ersten Perioden ohne Unterabtastung:
    # bei jedem 40. Schritt wuerde die schnelle Schwingung sonst verfaelscht.
    n_d = 400
    np.savez("data/tf2_energie_dicht.npz", t=np.arange(n_d) * sim.dt,
             W_naiv=W_naiv[:n_d] / W_naiv[:n_d].mean(),
             W_lf=W_lf[:n_d] / W_lf[:n_d].mean(), T11=T11, dt=sim.dt)
    del W_naiv, W_lf

    # =============================== Stabilitaetsgrenze ==================
    print("\n=== Stabilitaetsgrenze ===")
    nx_frei, ny_frei = M - 1, int(round(LY / dx)) - 1
    S_vorher = s_krit(nx_frei, ny_frei)
    print(f"  Resonator {nx_frei} x {ny_frei} freie Zellen")
    print(f"  Vorhersage fuer das unendliche Gitter (CFL): S = 1")
    print(f"  Vorhersage fuer diesen Resonator:            S = {S_vorher:.6f}")

    BLOCK = 50

    def wachstum(S, n=40000):
        sim, dxl, (a0, a1, b0, b1), *_ = aufbau(M, S=0.99)
        sim.dt = S / (C0 * np.sqrt(2.0 / dxl ** 2))
        sim._update_material_coeffs()
        sim._setup_pml()
        # hoechste Gittermode gezielt anregen, damit die Messung nicht auf
        # Rundungsfehler als Startwert angewiesen ist
        ii = np.arange(a0, a1) - (a0 - 1)
        jj = np.arange(b0, b1) - (b0 - 1)
        MX, MY = np.meshgrid(np.sin(nx_frei * np.pi * ii / (nx_frei + 1)),
                             np.sin(ny_frei * np.pi * jj / (ny_frei + 1)), indexing="ij")
        sim.Ez[a0:a1, b0:b1] = MX * MY
        # Als Mass dient die groesste Feldamplitude; sie waechst bei
        # Instabilitaet exponentiell und schwankt nicht mit der Periode.
        Amax = np.empty(n)
        with np.errstate(over="ignore", invalid="ignore"):
            for i in range(n):
                sim.step()
                Amax[i] = np.abs(sim.Ez).max()
                if not np.isfinite(Amax[i]) or Amax[i] > 1e80:
                    Amax[i:] = np.inf
                    break
        endl = np.isfinite(Amax)
        A0 = np.abs(Amax[:200]).max()
        An = Amax / A0
        # Huellkurve ueber Bloecke, damit die Schwingung das Wachstum nicht
        # ueberdeckt; sie ist zugleich das, was die Abbildung zeigt
        n_e = n // BLOCK * BLOCK
        h = np.array([np.max(An[k:k + BLOCK]) for k in range(0, n_e, BLOCK)])
        x = np.arange(len(h)) * BLOCK
        g = np.isfinite(h) & (h > 0)
        rate = np.polyfit(x[g][2:], np.log(h[g][2:]), 1)[0]
        return h, float(rate), int(endl.sum()), float(np.nanmax(h[np.isfinite(h)]))

    stab, kurven = [], {}
    for S in (0.9900, 1.0000, 1.0040, 1.0060, 1.00620, 1.00625, 1.00630, 1.0070, 1.0100):
        h, rate, n_endl, hmax = wachstum(S)
        stab.append(dict(S=S, rate_je_schritt=rate, n_endlich=n_endl,
                         groesste_amplitude=hmax))
        kurven[S] = h
        print(f"  S = {S:.5f}:  Wachstum {rate:+.3e} je Schritt   "
              f"groesste Amplitude nach {n_endl} Schritten: {hmax:.3e}")
    erg["stabilitaet"] = dict(S_vorhersage_resonator=S_vorher, S_vorhersage_cfl=1.0,
                              messungen=stab)
    np.savez("data/tf2_stabilitaet.npz",
             schritt=np.arange(0, 40000 // BLOCK * BLOCK, BLOCK),
             **{f"S_{str(k).replace('.', 'p')}": v for k, v in kurven.items()})

    json.dump(erg, open("data/tf2b.json", "w"), indent=1, default=float)
    print("\ngeschrieben: data/tf2b.json")
