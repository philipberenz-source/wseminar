"""Testfall 2a --- Eigenfrequenzen des PEC-Resonators.

Die Peaksuche laeuft unabhaengig von der Vorhersage: erst werden alle
Linien im Spektrum gesucht, danach werden sie den Moden zugeordnet. So
kann sowohl eine fehlende als auch eine ueberzaehlige Linie auffallen.

Nebenmaxima des Hann-Fensters liegen 2,4 Bins neben jeder Linie und
haben 2,7 % ihrer Hoehe; sie werden verworfen, weil innerhalb der
Hauptkeulenbreite ein staerkerer Peak liegt.
"""
import json, time
import numpy as np
from fdtd_core import FDTD2D, C0
import analytik as A
from tf2_resonator import LX, LY, W_PEC, N_PML, QUELLEN, SONDEN, aufbau

T_REC_KURZ = 3.0e-6      # Aufloesungsstudie
T_REC_LANG = 30.0e-6     # Referenzlauf, trennt auch das enge Modenpaar
F_MIN = 0.30e9
PAD = 4

def modenliste(f_max):
    mn = [(m, n) for m in range(1, 8) for n in range(1, 6)
          if F_MIN < A.f_kontinuum(m, n, LX, LY, C0) < f_max]
    return sorted(mn, key=lambda x: A.f_kontinuum(*x, LX, LY, C0))


def lauf(M, T_REC):
    sim, dx, (i0, i1, j0, j1), Lx_eff, Ly_eff = aufbau(M)
    tau, t0 = 2.0e-10, 8.0e-10
    for fx, fy in QUELLEN:
        sim.add_soft_source(i0 + int(round(fx * (i1 - i0))),
                            j0 + int(round(fy * (j1 - j0))),
                            lambda t: FDTD2D.gaussian_pulse(t, t0, tau))
    ks = [sim.add_probe(i0 + int(round(fx * (i1 - i0))),
                        j0 + int(round(fy * (j1 - j0)))) for fx, fy in SONDEN]
    n = int(round(T_REC / sim.dt))
    t_ = time.time()
    sim.run(n)
    w = np.hanning(n)
    Y = None
    for k in ks:
        y = np.array(sim._probes[k]["Ez"], dtype=float)
        Yk = np.abs(np.fft.rfft(y * w, n * PAD))
        Y = Yk if Y is None else Y + Yk
        del y, Yk
    f = np.fft.rfftfreq(n * PAD, sim.dt)
    sel = (f > 0.25e9) & (f < 1.4e9)
    return dict(M=M, dx=dx, dt=sim.dt, n=n, nx=sim.nx, ny=sim.ny, Lx=Lx_eff,
                Ly=Ly_eff, f=f[sel], Y=Y[sel], sek=time.time() - t_,
                df=1.0 / T_REC)


def linien(f, Y, f_max, schwelle=0.01, breite_bins=4.0, df=None):
    sel = (f > F_MIN) & (f < f_max)
    fs, Ys = f[sel], Y[sel] / Y[sel].max()
    ist = (Ys[1:-1] >= Ys[:-2]) & (Ys[1:-1] > Ys[2:]) & (Ys[1:-1] > schwelle)
    idx = np.where(ist)[0] + 1
    roh = []
    for i in idx:
        a, b, c = Ys[i - 1], Ys[i], Ys[i + 1]
        nen = a - 2 * b + c
        d = 0.5 * (a - c) / nen if nen != 0 else 0.0
        roh.append((fs[i] + d * (fs[1] - fs[0]), b))
    # Nebenmaxima verwerfen: innerhalb einer Hauptkeulenbreite gewinnt der groesste
    roh.sort(key=lambda x: -x[1])
    behalten = []
    for fk, ak in roh:
        if all(abs(fk - fb) > breite_bins * df for fb, _ in behalten):
            behalten.append((fk, ak))
    return sorted(behalten)


def auswerten(r, moden):
    lin = linien(r["f"], r["Y"], max(A.f_kontinuum(*mn, LX, LY, C0) for mn in moden) * 1.02,
                 df=r["df"])
    zeilen = []
    for (m, n) in moden:
        f_k = A.f_kontinuum(m, n, LX, LY, C0)
        f_g = A.f_gitter(m, n, LX, LY, r["dx"], r["dt"], C0)
        fm, am = min(lin, key=lambda x: abs(x[0] - f_g))
        zeilen.append(dict(m=m, n=n, f_kont=float(f_k), f_gitter=float(f_g),
                           f_mess=float(fm), amp=float(am),
                           N_lam=float((C0 / f_k) / r["dx"])))
    return lin, zeilen


if __name__ == "__main__":
    erg = {}

    # ---------------- Aufloesungsstudie, kurzer Lauf ----------------------
    F_MAX_K = 1.0e9
    MK = modenliste(F_MAX_K)
    print(f"=== Aufloesungsstudie (T = {T_REC_KURZ*1e6:.0f} us, "
          f"{len(MK)} Moden bis {F_MAX_K/1e9:.1f} GHz) ===")
    tab = []
    for M in (9, 18, 36, 72):
        r = lauf(M, T_REC_KURZ)
        lin, zeilen = auswerten(r, MK)
        d_k = np.array([abs(z["f_mess"] / z["f_kont"] - 1) for z in zeilen])
        d_g = np.array([abs(z["f_mess"] / z["f_gitter"] - 1) for z in zeilen])
        N11 = (C0 / A.f_kontinuum(1, 1, LX, LY, C0)) / r["dx"]
        tab.append(dict(M=M, N_lam=float(N11), dx=r["dx"], dt=r["dt"], n=r["n"],
                        nx=r["nx"], ny=r["ny"], sek=r["sek"], n_linien=len(lin),
                        moden=zeilen, mittel_kont=float(d_k.mean()),
                        max_kont=float(d_k.max()), mittel_gitter=float(d_g.mean()),
                        max_gitter=float(d_g.max())))
        print(f"  M={M:3} (N_lam,11={N11:5.1f}, {r['nx']}x{r['ny']}, {r['n']:7} Schritte, "
              f"{r['sek']:5.1f}s): {len(lin)} Linien (erwartet {len(MK)});  "
              f"Abw. Kontinuum {100*d_k.mean():.4f} %,  Gitterformel {1e6*d_g.mean():.2f} ppm")
        del r
    erg["aufloesung"] = tab
    erg["T_REC_kurz"] = T_REC_KURZ

    Nl = np.array([t["N_lam"] for t in tab])
    dk = np.array([t["mittel_kont"] for t in tab])
    p = np.polyfit(np.log(1 / Nl), np.log(dk), 1)[0]
    erg["konvergenz_ordnung"] = float(p)
    erg["konvergenz_verhaeltnisse"] = list(dk[:-1] / dk[1:])
    print(f"\n  Konvergenzordnung gegen die Kontinuumsformel: p = {p:.2f}")
    print(f"  Verhaeltnisse bei Verdopplung: {np.array2string(dk[:-1]/dk[1:], precision=2)}")

    # ---------------- Referenzlauf, langer Lauf ---------------------------
    F_MAX_L = 1.2e9
    ML = modenliste(F_MAX_L)
    print(f"\n=== Referenzlauf M = 18, T = {T_REC_LANG*1e6:.0f} us "
          f"({len(ML)} Moden bis {F_MAX_L/1e9:.1f} GHz) ===")
    r = lauf(18, T_REC_LANG)
    lin, zeilen = auswerten(r, ML)
    d_k = np.array([abs(z["f_mess"] / z["f_kont"] - 1) for z in zeilen])
    d_g = np.array([abs(z["f_mess"] / z["f_gitter"] - 1) for z in zeilen])
    print(f"  {r['n']} Schritte, {r['sek']:.0f}s, Frequenzaufloesung {r['df']/1e3:.0f} kHz")
    print(f"  {len(lin)} Linien gefunden, {len(ML)} vorhergesagt")
    print("   (m,n) N_lam   f_Kont/MHz  f_Gitter/MHz  f_gem/MHz    Abw.Kont    Abw.Gitter")
    for z in zeilen:
        print(f"   ({z['m']},{z['n']}) {z['N_lam']:5.1f} {z['f_kont']/1e6:11.4f} "
              f"{z['f_gitter']/1e6:11.4f} {z['f_mess']/1e6:11.4f} "
              f"{100*(z['f_mess']/z['f_kont']-1):+9.4f} % "
              f"{1e6*(z['f_mess']/z['f_gitter']-1):+9.2f} ppm")
    print(f"  Mittel: Kontinuum {100*d_k.mean():.4f} %, max {100*d_k.max():.4f} %;  "
          f"Gitterformel {1e6*d_g.mean():.2f} ppm, max {1e6*d_g.max():.2f} ppm")
    erg["referenzlauf"] = dict(M=18, T_REC=T_REC_LANG, df=r["df"], n=r["n"],
                               n_linien=len(lin), n_moden=len(ML), moden=zeilen,
                               mittel_kont=float(d_k.mean()), max_kont=float(d_k.max()),
                               mittel_gitter=float(d_g.mean()), max_gitter=float(d_g.max()),
                               linien=[[float(a), float(b)] for a, b in lin])
    np.savez("data/tf2_spektrum.npz", f=r["f"], Y=r["Y"] / r["Y"].max(),
             f_kont=np.array([A.f_kontinuum(m, n, LX, LY, C0) for m, n in ML]),
             f_git=np.array([A.f_gitter(m, n, LX, LY, r["dx"], r["dt"], C0) for m, n in ML]),
             moden=np.array(ML))

    json.dump(erg, open("data/tf2a.json", "w"), indent=1, default=float)
    print("\ngeschrieben: data/tf2a.json")
