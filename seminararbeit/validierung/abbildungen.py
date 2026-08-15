"""Erzeugt alle Abbildungen des Validierungskapitels als PDF."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import analytik as A

plt.rcParams.update({
    "font.family": "serif", "font.size": 9, "axes.labelsize": 9,
    "axes.titlesize": 9.5, "legend.fontsize": 8, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "axes.grid": True, "grid.alpha": 0.25,
    "grid.linewidth": 0.5, "lines.linewidth": 1.2, "axes.linewidth": 0.7,
    "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "axes.spines.top": False, "axes.spines.right": False,
})
SCHWARZ, BLAU, ROT, GRAU = "#1a1a1a", "#1f5fa8", "#c0392b", "#7f7f7f"
AUS = "/sessions/wonderful-vibrant-cannon/mnt/seminararbeit/abbildungen/"

t1a = json.load(open("data/tf1a.json"))
t1b = json.load(open("data/tf1b.json"))
t2a = json.load(open("data/tf2a.json"))
t2b = json.load(open("data/tf2b.json"))
t3 = json.load(open("data/tf3.json"))
t4 = json.load(open("data/tf4.json"))


def sichern(fig, name):
    fig.savefig(AUS + name, format="pdf")
    plt.close(fig)
    print("  ", name)


# ============================== 1  Dispersion ==============================
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
N = np.logspace(np.log10(4), np.log10(100), 250)
for S, ls in ((0.35, ":"), (0.7, "--"), (0.99, "-")):
    ax[0].plot(N, A.cp_achse(N, S), ls, color=SCHWARZ, lw=1.1)
    ax[0].plot(N, A.cp_diagonal(N, S), ls, color=BLAU, lw=1.1)
ax[0].plot([], [], "-", color=SCHWARZ, label="längs der Gitterachse")
ax[0].plot([], [], "-", color=BLAU, label="in der Diagonalen")
ax[0].plot([], [], "-", color=GRAU, label="$S=0{,}99$")
ax[0].plot([], [], "--", color=GRAU, label="$S=0{,}70$")
ax[0].plot([], [], ":", color=GRAU, label="$S=0{,}35$")
m = t1a["aufloesung"]
ax[0].plot([a["N_lam"] for a in m], [a["v"] for a in m], "o", ms=4.5,
           mfc="white", mec=ROT, mew=1.2, label="gemessen (ebene Welle)")
ax[0].set_xscale("log")
ax[0].set_xlabel(r"Zellen je Wellenlänge $N_\lambda$")
ax[0].set_ylabel(r"$\tilde{c}_p / c$")
ax[0].set_ylim(0.90, 1.005)
ax[0].axhline(1, color=GRAU, lw=0.6)
ax[0].legend(loc="lower right", frameon=False, ncol=1, fontsize=7,
             handlelength=1.6, borderpad=0.1, labelspacing=0.3)
ax[0].set_title("(a) Auflösungsabhängigkeit")

r = t1b["richtung"]
th = np.linspace(0, 45, 200)
ax[1].plot(th, [A.cp_verhaeltnis(20.0, 0.99, np.deg2rad(x)) for x in th],
           "-", color=SCHWARZ, label="Dispersionsrelation")
ax[1].plot([x["winkel"] for x in r], [x["korrigiert"] for x in r], "o", ms=4.5,
           mfc="white", mec=ROT, mew=1.2, label="gemessen (Kreiswelle)")
ax[1].axhline(1, color=GRAU, lw=0.6)
ax[1].set_xlabel(r"Ausbreitungsrichtung $\theta$ / Grad")
ax[1].set_ylabel(r"$\tilde{c}_p / c$")
ax[1].set_xticks([0, 15, 30, 45])
ax[1].legend(loc="lower right", frameon=False, fontsize=7.5,
             handlelength=1.6, borderpad=0.1)
ax[1].set_title(r"(b) Richtungsabhängigkeit, $N_\lambda = 20$")
fig.tight_layout()
sichern(fig, "abb_dispersion.pdf")

# ============================== 2  Konvergenz ==============================
fig, ax = plt.subplots(figsize=(4.0, 3.0))
m = t1a["aufloesung"]
ax.loglog([a["N_lam"] for a in m], [1 - a["v"] for a in m], "o-", color=SCHWARZ,
          ms=4.5, mfc="white", label="1: Phasengeschwindigkeit")
m2 = t2a["aufloesung"]
ax.loglog([a["N_lam"] for a in m2], [a["mittel_kont"] for a in m2], "s-",
          color=BLAU, ms=4.5, mfc="white", label="2: Eigenfrequenzen")
m3 = t3["aufloesung"]["punkte"]
ax.loglog([a["N_lam2"] for a in m3], [abs(a["abw"]) for a in m3], "^-",
          color=ROT, ms=4.5, mfc="white", label="3: Reflexionsgrad")
ax.loglog([5, 90], [2e-1, 2e-1 * (90 / 5) ** -2], "--", color=GRAU, lw=1.0)
ax.text(17, 1.1e-2, r"Steigung $-2$", color=GRAU, fontsize=8, rotation=-31)
ax.set_xlabel(r"Zellen je Wellenlänge $N_\lambda$")
ax.set_ylabel("relative Abweichung")
ax.legend(frameon=False, loc="lower left")
fig.tight_layout()
sichern(fig, "abb_konvergenz.pdf")

# ============================== 3  Verzug ==================================
# Nur die Sinusquelle; die frueher hier gezeigte Verformung eines
# Ricker-Pulses ist entfallen, damit im ganzen Testfall dieselbe Anregung
# verwendet wird.
fig, axL = plt.subplots(figsize=(4.2, 3.0))

for k, c in (("10", BLAU), ("20", SCHWARZ)):
    v = t1a["verzug"][k]
    r_ = np.array(v["r_lam"])
    axL.plot(r_, v["verzug"], color=c, label=rf"$N_\lambda = {k}$")
    axL.plot(r_, v["steigung_theorie"] * (r_ - r_[0]), "--", color="white", lw=1.6,
             zorder=3)
    axL.plot(r_, v["steigung_theorie"] * (r_ - r_[0]), "--", color=c, lw=0.9,
             zorder=4)
axL.plot([], [], "--", color=GRAU, lw=0.9, label="aus der Dispersionsrelation")
axL.set_xlabel(r"zurückgelegte Strecke / $\lambda$")
axL.set_ylabel(r"Verzug / $\lambda$")
axL.legend(frameon=False, loc="upper left")
fig.tight_layout()
sichern(fig, "abb_verzug.pdf")

# ============================== 4  Amplitude ===============================
d = np.load("data/tf1_amplitude.npz")
fig, ax = plt.subplots(figsize=(4.0, 2.9))
r_ = d["r_lam"]
amp = d["amp"] / d["amp"][0]
ax.loglog(r_, amp, "o", ms=3, mfc="white", mec=ROT, mew=0.9, label="gemessen")
ax.loglog(r_, (r_ / r_[0]) ** -0.5, "-", color=SCHWARZ,
          label=r"$A \propto r^{-1/2}$  (zwei Dimensionen)")
ax.loglog(r_, (r_ / r_[0]) ** -1.0, "--", color=GRAU,
          label=r"$A \propto r^{-1}$  (drei Dimensionen)")
ax.set_xlabel(r"Abstand von der Quelle / $\lambda$")
ax.set_ylabel("Amplitude (normiert)")
ax.set_xticks([2, 3, 4, 5, 6, 7]); ax.set_yticks([0.2, 0.3, 0.5, 0.7, 1.0])
ax.set_xticklabels(["2", "3", "4", "5", "6", "7"])
ax.set_yticklabels(["0,2", "0,3", "0,5", "0,7", "1,0"])
ax.legend(frameon=False, fontsize=7.5)
fig.tight_layout()
sichern(fig, "abb_amplitude.pdf")

# ============================== 5  Differenzfeld ===========================
fig, ax = plt.subplots(1, 3, figsize=(6.6, 2.5))
for k, N in enumerate((10, 20, 40)):
    D = np.load(f"data/tf1_diff_N{N}.npy")
    n = D.shape[0]
    e = (n // 2) / (N)          # halbe Kantenlaenge in Wellenlaengen
    im = ax[k].imshow(D.T * 100, origin="lower", cmap="magma",
                      extent=[-e, e, -e, e], vmin=0, vmax=18)
    ax[k].set_xlim(-7.6, 7.6); ax[k].set_ylim(-7.6, 7.6)
    ax[k].set_title(rf"$N_\lambda = {N}$")
    ax[k].set_xlabel(r"$x/\lambda$")
    ax[k].grid(False)
    if k == 0:
        ax[k].set_ylabel(r"$y/\lambda$")
    else:
        ax[k].set_yticklabels([])
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label("Abweichung / %", fontsize=8)
sichern(fig, "abb_differenzfeld.pdf")

# ============================== 6  Resonator ===============================
d = np.load("data/tf2_spektrum.npz")
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.8))
f, Y = d["f"] / 1e9, d["Y"]
ax[0].semilogy(f, np.maximum(Y, 1e-5), color=SCHWARZ, lw=0.8)
# das enge Paar (3,2)/(4,1) gemeinsam beschriften
gruppen, akt = [], []
for fk, mn in zip(d["f_kont"], d["moden"]):
    if akt and fk / 1e9 - akt[-1][0] < 0.02:
        akt.append((fk / 1e9, mn))
    else:
        if akt: gruppen.append(akt)
        akt = [(fk / 1e9, mn)]
gruppen.append(akt)
for g in gruppen:
    fm = sum(x[0] for x in g) / len(g)
    txt = ",".join(f"({m[0]},{m[1]})" for _, m in g)
    for fk, _ in g:
        ax[0].axvline(fk, color=ROT, ls=":", lw=0.8)
    ax[0].annotate(txt, (fm, 4e1), ha="center", va="bottom", fontsize=6.5,
                   color=ROT, rotation=90)
ax[0].set_xlim(0.35, 1.21); ax[0].set_ylim(1e-4, 1e7)
ax[0].set_yticks([1e-4, 1e-2, 1e0])
ax[0].set_xlabel("Frequenz / GHz")
ax[0].set_ylabel("Amplitude (normiert)")
ax[0].set_title("(a) Spektrum, punktiert die Vorhersage")

mo = t2a["referenzlauf"]["moden"]
Nl = [z["N_lam"] for z in mo]
dk = [abs(z["f_mess"] / z["f_kont"] - 1) for z in mo]
dg = [max(abs(z["f_mess"] / z["f_gitter"] - 1), 3e-9) for z in mo]
ax[1].semilogy(Nl, dk, "o", ms=5, color=ROT, label="gegen die Kontinuumsformel")
ax[1].semilogy(Nl, dg, "s", ms=5, mfc="white", mec=BLAU, mew=1.2,
               label="gegen die Gitterformel")
ax[1].axhline(3.3e4 / 4.5e8, color=GRAU, ls="--", lw=0.8)
ax[1].text(7.0, 5.5e-5, "Grenze der Frequenzmessung", fontsize=6.5,
           color=GRAU, ha="left", va="top")
ax[1].set_yscale("log"); ax[1].set_ylim(1e-9, 3e2)
ax[1].set_yticks([1e-8, 1e-6, 1e-4, 1e-2])
ax[1].set_xlim(6.5, 21.5)
ax[1].set_xlabel(r"Zellen je Wellenlänge der Mode $N_\lambda$")
ax[1].set_ylabel("Betrag der relativen Abweichung")
ax[1].legend(frameon=False, loc="upper center", fontsize=7, ncol=2,
             handlelength=1.4, borderpad=0.1, columnspacing=1.0)
ax[1].set_title("(b) Abweichung je Mode")
fig.tight_layout()
sichern(fig, "abb_resonator.pdf")

# ============================== 7  Energie + Stabilitaet ===================
d = np.load("data/tf2_energie_dicht.npz")
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
t_ = d["t"] / d["T11"]
sel = t_ < 6
ax[0].plot(t_[sel], d["W_naiv"][sel], color=GRAU, lw=0.9,
           label=r"$\frac{\varepsilon}{2}|E^n|^2 + \frac{\mu}{2}|H^{n+1/2}|^2$")
ax[0].plot(t_[sel], d["W_lf"][sel], color=SCHWARZ,
           label=r"$\frac{\varepsilon}{2}|E^n|^2 + \frac{\mu}{2}H^{n-1/2}\!\cdot\!H^{n+1/2}$")
ax[0].set_xlabel("Zeit / Periode der Grundmode")
ax[0].set_ylabel(r"$W / \bar{W}$")
ax[0].set_ylim(0.90, 1.42)
ax[0].legend(frameon=False, loc="upper center", fontsize=6.8, ncol=1,
             handlelength=1.6, borderpad=0.1)
ax[0].set_title("(a) Energie über der Zeit")

d = np.load("data/tf2_stabilitaet.npz")
xs = d["schritt"]
for key, lab, c in (("S_0p99", "0{,}9900", GRAU), ("S_1p006", "1{,}00600", BLAU),
                    ("S_1p0062", "1{,}00620", SCHWARZ),
                    ("S_1p00625", "1{,}00625", ROT),
                    ("S_1p007", "1{,}00700", "#e07b39")):
    y = d[key]
    g = np.isfinite(y) & (y > 0)
    ax[1].semilogy(xs[g], y[g], color=c, label=f"$S = {lab}$")
ax[1].set_ylim(3e-1, 1e45)
ax[1].set_xlim(0, 12000)
ax[1].set_yticks([1e0, 1e10, 1e20, 1e30, 1e40])
ax[1].set_xlabel("Zeitschritt")
ax[1].set_ylabel(r"Hüllkurve von $\max|E_z|$")
ax[1].legend(frameon=False, loc="lower right", fontsize=7,
             handlelength=1.6, borderpad=0.1)
ax[1].set_title("(b) Verhalten an der Stabilitätsgrenze")
fig.tight_layout()
sichern(fig, "abb_stabilitaet.pdf")

# ============================== 8  Grenzflaeche ============================
dz = np.load("data/tf3_zeit.npz")
dk = np.load("data/tf3_kurven.npz")
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
t_ = dz["t"] * 1e9
ax[0].plot(t_, dz["ges"], color=SCHWARZ, lw=1.0, label="Gesamtfeld an der Sonde")
ax[0].plot(t_, dz["inc"], color=GRAU, lw=0.9, ls="--", label="Referenzlauf ohne Sprung")
ax[0].plot(t_, dz["ref"], color=ROT, lw=1.0, label="Differenz = reflektiertes Feld")
ax[0].set_xlim(0, 23)
ax[0].set_xlabel("Zeit / ns")
ax[0].set_ylabel(r"$E_z$ (normiert)")
ax[0].set_ylim(-0.45, 1.05)
ax[0].legend(frameon=False, fontsize=7, loc="upper right")
ax[0].set_title("(a) Trennung durch Differenzbildung")

for e, c in ((1.5, "#8fb8de"), (2.25, BLAU), (4.0, SCHWARZ), (6.0, "#e07b39"),
             (9.0, ROT)):
    k = f"eps{str(e).replace('.','p')}"
    N_ = dk[k + "_N_lam2"]
    r_ = dk[k + "_r"]
    s = (N_ > 8) & (N_ < 55)
    ax[1].plot(N_[s], r_[s], color=c, label=rf"$\varepsilon_r = {e}$".replace(".", "{,}"))
    ax[1].axhline(abs(A.fresnel_r(1, e)), color=c, ls=":", lw=0.8)
ax[1].set_xlim(8, 55)
ax[1].set_xlabel(r"Zellen je Wellenlänge im Medium $N_{\lambda,2}$")
ax[1].set_ylabel(r"$|r|$")
ax[1].set_ylim(0.05, 0.72)
ax[1].legend(frameon=False, ncol=3, fontsize=6.8, loc="upper center",
             columnspacing=1.0, handlelength=1.4, borderpad=0.1)
ax[1].set_title("(b) Reflexionsgrad; punktiert: Fresnel")
fig.tight_layout()
sichern(fig, "abb_grenzflaeche.pdf")

# ============================== 9  PML =====================================
d4e = np.load("data/tf4_energie.npz")
fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.7))
for dd, c in ((0.2, "#8fb8de"), (0.5, BLAU), (1.0, SCHWARZ), (2.0, ROT)):
    k = "d" + str(dd).replace(".", "p")
    ax[0].semilogy(d4e["t"], np.maximum(d4e[k], 1e-16), color=c,
                   label=rf"$d = {dd}\,\lambda$".replace(".", "{,}"))
ax[0].set_xlabel("Zeit / ns")
ax[0].set_ylabel(r"$W(t) / W_{\max}$")
ax[0].set_ylim(1e-14, 3)
ax[0].legend(frameon=False, fontsize=7.5)
ax[0].set_title("(a) Energie im Gebiet über der Zeit")

t = t4["dicke"]
sk = {a["d_lam"]: a["R_1GHz"] for a in t4["senkrecht"]}
ax[1].semilogy([a["d_lam"] for a in t], [a["R_1GHz"] for a in t], "o-",
               color=SCHWARZ, ms=4.5, mfc="white", label="Punktquelle, alle Winkel")
ax[1].semilogy(list(sk.keys()), list(sk.values()), "s--", color=BLAU, ms=4,
               mfc="white", label="ebene Welle, senkrecht")
ax[1].set_xlabel(r"Dicke der Schicht $d / \lambda$")
ax[1].set_ylabel(r"$|R|$ bei 1 GHz")
ax[1].set_ylim(5e-5, 3)
ax[1].legend(frameon=False, fontsize=7, loc="upper right",
             handlelength=1.6, borderpad=0.1)
ax[1].set_title("(b) Restreflexion über der Dicke")
fig.tight_layout()
sichern(fig, "abb_pml.pdf")

# ============================== 10  Kreiswelle Feld ========================
F = np.load("data/tf1_feld_N20.npy")
n = F.shape[0]
e = (n // 2) / 20
fig, ax = plt.subplots(figsize=(3.4, 3.2))
v = np.percentile(np.abs(F), 99.5)
ax.imshow(F.T, origin="lower", cmap="RdBu_r", vmin=-v, vmax=v,
          extent=[-e, e, -e, e])
for a_ in (0, 45):
    ax.plot([0, 8 * np.cos(np.deg2rad(a_))], [0, 8 * np.sin(np.deg2rad(a_))],
            color="k", lw=0.8, ls="--")
ax.text(7.0, 0.5, "Achse", fontsize=7.5)
ax.text(5.0, 5.6, "Diagonale", fontsize=7.5)
ax.set_xlim(-8.5, 8.5); ax.set_ylim(-8.5, 8.5)
ax.set_xlabel(r"$x/\lambda$"); ax.set_ylabel(r"$y/\lambda$")
ax.grid(False)
fig.tight_layout()
sichern(fig, "abb_kreiswelle.pdf")

print("fertig.")
