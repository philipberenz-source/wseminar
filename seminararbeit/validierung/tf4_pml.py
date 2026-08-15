"""Testfall 4 --- absorbierender Rand (PML).

Gemessen wird gegen eine Referenzloesung: Derselbe Puls laeuft einmal in
einem so grossen Gebiet, dass dessen Rand waehrend der Beobachtungsdauer
keine Rolle spielt, und einmal im kleinen Gebiet mit dem zu pruefenden
Rand. Die Differenz beider Felder an denselben Punkten ist genau der
Fehler, den der Rand verursacht.

Ausgewertet wird frequenzaufgeloest, denn die Dicke der Schicht bemisst
sich an der Wellenlaenge: ein und dieselbe Schicht ist fuer hohe
Frequenzen dick und fuer tiefe duenn. Eine einzige Zahl fuer einen
breitbandigen Puls wuerde diesen Zusammenhang gerade verdecken.

Angeregt wird mit einem Ricker-Puls, dessen Spektrum bei der Frequenz
null verschwindet; ein Gauss-Puls hinterliesse ein statisches Restfeld,
das kein Rand absorbieren kann.
"""
import json, time
import numpy as np
from fdtd_core import FDTD2D, C0

F0 = 1.0e9
LAM0 = C0 / F0
FP, TP = 1.0e9, 1.5e-9
F_BAND = (0.4e9, 2.0e9)


def lauf(R_lam, d_lam, N_lam, r_sonde, n_steps, S=0.99):
    dx = LAM0 / N_lam
    npml = max(4, int(round(d_lam * N_lam)))
    n = 2 * int(round(R_lam * N_lam)) + 1 + 2 * npml
    sim = FDTD2D(n, n, dx, dx, courant=S, npml=npml)
    c_ = n // 2
    sim.add_soft_source(c_, c_, lambda t: FDTD2D.ricker_wavelet(t, TP, FP))
    m = int(round(r_sonde * N_lam))
    pkt = sorted({p for k in range(-m, m + 1, max(1, m // 6))
                  for p in ((c_ + m, c_ + k), (c_ - m, c_ + k),
                            (c_ + k, c_ + m), (c_ + k, c_ - m))})
    ks = [sim.add_probe(i, j) for i, j in pkt]
    W = np.empty(n_steps)
    t_ = time.time()
    for i in range(n_steps):
        sim.step()
        W[i] = sim.total_energy()
    y = np.array([sim._probes[k]["Ez"] for k in ks])
    return y, W, dict(n=n, npml=npml, dx=dx, dt=sim.dt, sek=time.time() - t_)


def reflexion(y, y_ref, dt):
    """|R(f)| = |FFT(Differenz)| / |FFT(Referenz)|, groesster Wert ueber alle Sonden."""
    n = y.shape[1]
    f = np.fft.rfftfreq(n, dt)
    sel = (f >= F_BAND[0]) & (f <= F_BAND[1])
    D = np.abs(np.fft.rfft(y - y_ref, axis=1))[:, sel]
    Rf = np.abs(np.fft.rfft(y_ref, axis=1))[:, sel]
    return f[sel], (D / Rf).max(axis=0)


if __name__ == "__main__":
    erg = {}
    R_TEST, R_SONDE, N_LAM = 3.0, 2.6, 20
    dt_ = FDTD2D(64, 64, LAM0 / N_LAM, LAM0 / N_LAM, courant=0.99, npml=8).dt
    N_STEPS = int(round(30 / F0 / dt_))
    print(f"Beobachtung {N_STEPS} Schritte = {N_STEPS*dt_*F0:.0f} Perioden; "
          f"Sondenring bei {R_SONDE} lam, Gebietsrand bei {R_TEST} lam")

    y_ref, W_ref, m_ref = lauf(22.0, 3.0, N_LAM, R_SONDE, N_STEPS)
    print(f"Referenzlauf: Gitter {m_ref['n']}^2, {m_ref['sek']:.1f} s")

    # ---- Dicke der Schicht in Wellenlaengen ----------------------------
    print("\n=== Restreflexion ueber der Schichtdicke (N_lam = 20) ===")
    print("  d/lam  Zellen   |R| bei 1 GHz   |R| max im Band 0.4-2 GHz   in dB")
    tab, zeit, kurven = [], {}, {}
    for d in (0.2, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
        y, W, m = lauf(R_TEST, d, N_LAM, R_SONDE, N_STEPS)
        f, R = reflexion(y, y_ref, m["dt"])
        i0 = int(np.argmin(np.abs(f - F0)))
        tab.append(dict(d_lam=d, zellen=m["npml"], gitter=m["n"],
                        R_1GHz=float(R[i0]), R_max=float(R.max()),
                        dB=float(20 * np.log10(R[i0]))))
        zeit[d] = W / W.max()
        kurven[d] = (f, R)
        print(f"  {d:5.3f}  {m['npml']:5}    {R[i0]:.3e}       {R.max():.3e}"
              f"          {20*np.log10(R[i0]):7.1f}")
    erg["dicke"] = tab
    np.savez("data/tf4_energie.npz", t=np.arange(N_STEPS) * dt_ * 1e9,
             **{f"d{str(d).replace('.','p')}": v for d, v in zeit.items()})
    np.savez("data/tf4_frequenz.npz", f=kurven[0.5][0],
             **{f"d{str(d).replace('.','p')}": v[1] for d, v in kurven.items()})

    # ---- Dieselbe Zellenzahl bei verschiedener Aufloesung --------------
    # Gleiche Anordnung wie oben, damit die Zahlen vergleichbar bleiben;
    # nur das Gebiet ist kleiner, damit die feinsten Gitter rechenbar sind.
    print("\n=== Feste Zellenzahl (10 Zellen) bei verschiedener Aufloesung ===")
    print("  N_lam   d/lam   |R| bei 1 GHz     in dB")
    fest = []
    for N in (10, 20, 40):
        dtn = FDTD2D(64, 64, LAM0 / N, LAM0 / N, courant=0.99, npml=8).dt
        ns = int(round(12 / F0 / dtn))
        yr, _, _ = lauf(8.0, 3.0, N, R_SONDE, ns)
        y, _, m = lauf(R_TEST, 10.0 / N, N, R_SONDE, ns)
        f, R = reflexion(y, yr, m["dt"])
        i0 = int(np.argmin(np.abs(f - F0)))
        fest.append(dict(N_lam=N, d_lam=10.0 / N, zellen=m["npml"],
                         R_1GHz=float(R[i0]), dB=float(20 * np.log10(R[i0]))))
        print(f"  {N:5}   {10.0/N:5.2f}   {R[i0]:.3e}       {20*np.log10(R[i0]):7.1f}")
    erg["feste_zellenzahl"] = fest

    # ---- Senkrechter Einfall (ebene Welle im Kanal) ---------------------
    # Dieselbe Messung, aber die Welle trifft frontal auf die Schicht. Weil das
    # Feld nicht von y abhaengt, ist die Anordnung im Kern eindimensional.
    def kanal(L_lam, d_lam, N_lam, n_steps, r_sonde=1.0):
        dx = LAM0 / N_lam
        npml = max(4, int(round(d_lam * N_lam)))
        nx = int(round(L_lam * N_lam)) + 2 * npml
        ny = int(round(2.0 * N_lam)) + 2 * npml
        sim = FDTD2D(nx, ny, dx, dx, courant=0.99, npml=npml)
        i_src = npml + 2
        for j in range(1, ny - 1):
            sim.add_soft_source(i_src, j,
                                lambda t: FDTD2D.ricker_wavelet(t, TP, FP))
        k = sim.add_probe(i_src + int(round(r_sonde * N_lam)), ny // 2)
        sim.run(n_steps)
        return np.array(sim._probes[k]["Ez"])[None, :], sim.dt

    print("\n=== Senkrechter Einfall (ebene Welle), N_lam = 20 ===")
    print("  d/lam  Zellen   |R| bei 1 GHz     in dB")
    ns_k = int(round(40 / F0 / dt_))
    yr_k, _ = kanal(40.0, 3.0, N_LAM, ns_k)          # Referenz: sehr langer Kanal
    senk = []
    for d in (0.2, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
        y, dtk = kanal(4.0, d, N_LAM, ns_k)
        f, R = reflexion(y, yr_k, dtk)
        i0 = int(np.argmin(np.abs(f - F0)))
        senk.append(dict(d_lam=d, zellen=max(4, int(round(d * N_LAM))),
                         R_1GHz=float(R[i0]), dB=float(20 * np.log10(R[i0]))))
        print(f"  {d:5.3f}  {max(4,int(round(d*N_LAM))):5}    {R[i0]:.3e}    "
              f"{20*np.log10(R[i0]):7.1f}")
    erg["senkrecht"] = senk

    # ---- Ist der Daempfungsparameter guenstig gewaehlt? -----------------
    print("\n=== Daempfungsprofil: Faktor auf sigma_max (d = 0.5 lam) ===")
    from fdtd_core import EPS0
    orig = FDTD2D._pml_profile

    def profil(skal, m_ord):
        def prof(self, n):
            tiefe = self.npml
            s_max = skal * (m_ord + 1) / (150.0 * np.pi * min(self.dx, self.dy))
            idx = np.arange(n)
            d_ = np.maximum(np.maximum(tiefe - idx, idx - (n - 1 - tiefe)), 0)
            sig = s_max * np.minimum(d_ / tiefe, 1.0) ** m_ord
            x = sig * self.dt / (2.0 * EPS0)
            return sig, 2.0 * x / (1.0 + x), (1.0 - x) / (1.0 + x)
        return prof

    print("  Faktor   m   |R|(1 GHz)     dB")
    prof_tab = []
    for m_ord, skalen in ((3, (0.2, 0.3, 0.5, 1.0, 2.0)), (2, (0.5,)), (4, (0.5,))):
        for skal in skalen:
            FDTD2D._pml_profile = profil(skal, m_ord)
            y, W, mm = lauf(R_TEST, 0.5, N_LAM, R_SONDE, N_STEPS)
            f, R = reflexion(y, y_ref, mm["dt"])
            i0 = int(np.argmin(np.abs(f - F0)))
            prof_tab.append(dict(faktor=skal, m=m_ord, R_1GHz=float(R[i0]),
                                 dB=float(20 * np.log10(R[i0]))))
            print(f"  {skal:6.2f}  {m_ord:2}   {R[i0]:.3e}   {20*np.log10(R[i0]):7.1f}")
    FDTD2D._pml_profile = orig
    erg["daempfungsprofil"] = prof_tab

    # ---- Energieverfahren als einfacheres Vergleichsverfahren ----------
    print("\n=== Energieverfahren (Restenergie im Gebiet) ===")
    print("  d/lam   W_rest/W_max    rho = sqrt(W_rest/W_max)")
    en = []
    for d, W in zeit.items():
        en.append(dict(d_lam=d, W_rel=float(W[-1]), rho=float(np.sqrt(W[-1]))))
        print(f"  {d:5.3f}   {W[-1]:.3e}       {np.sqrt(W[-1]):.3e}")
    erg["energieverfahren"] = en

    json.dump(erg, open("data/tf4.json", "w"), indent=1, default=float)
    print("\ngeschrieben: data/tf4.json")
