"""Analytische Vergleichsgroessen fuer die Validierung.

Enthaelt: Bessel-/Hankel-Funktionen (ohne scipy), die numerische
Dispersionsrelation des 2D-Yee-Gitters und die Fresnel-Koeffizienten.
Jede Funktion ist am Ende der Datei gegen eine unabhaengige Methode geprueft.
"""
import numpy as np

EPS0 = 8.8541878128e-12
EULER_GAMMA = 0.5772156649015328606

# ---------------------------------------------------------------- Bessel ----

def j0(z):
    """J_0(z) ueber die periodische Trapezregel:  J0 = 1/(2pi) int_0^2pi cos(z sin t) dt.

    Der Integrand ist glatt und periodisch; die Trapezregel konvergiert
    dann exponentiell und erreicht fuer |z| < 100 Maschinengenauigkeit.
    """
    z = np.asarray(z, dtype=float)
    M = 512
    t = np.arange(M) * (2 * np.pi / M)
    return np.cos(z[..., None] * np.sin(t)).mean(axis=-1)


# Koeffizienten nach Abramowitz & Stegun, 9.4.2 / 9.4.3
_Y0_SMALL = [0.36746691, 0.60559366, -0.74350384, 0.25300117,
             -0.04261214, 0.00427916, -0.00024846]
_F0_LARGE = [0.79788456, -0.00000077, -0.00552740, -0.00009512,
             0.00137237, -0.00072805, 0.00014476]
_T0_LARGE = [-0.78539816, -0.04166397, -0.00003954, 0.00262573,
             -0.00054125, -0.00029333, 0.00013558]


def y0(z):
    """Y_0(z), Polynomnaeherung nach Abramowitz & Stegun (Fehler < 1.4e-8)."""
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)

    klein = z <= 3.0
    if np.any(klein):
        zs = z[klein]
        t2 = (zs / 3.0) ** 2
        s = np.zeros_like(zs)
        for k, c in enumerate(_Y0_SMALL):
            s += c * t2 ** k
        out[klein] = (2 / np.pi) * np.log(zs / 2) * j0(zs) + s

    gross = ~klein
    if np.any(gross):
        zg = z[gross]
        t = 3.0 / zg
        f0 = np.zeros_like(zg)
        th = np.zeros_like(zg)
        for k, c in enumerate(_F0_LARGE):
            f0 += c * t ** k
        th += zg
        for k, c in enumerate(_T0_LARGE):
            th += c * t ** k
        out[gross] = f0 / np.sqrt(zg) * np.sin(th)
    return out


def hankel1_0(z):
    """H_0^(1)(z) = J_0(z) + i Y_0(z)."""
    return j0(z) + 1j * y0(z)


def hankel1_0_asym(z, n_terms=4):
    """Asymptotische Entwicklung von H_0^(1)(z) fuer grosse Argumente."""
    z = np.asarray(z, dtype=float)
    s = np.ones(z.shape, dtype=complex)
    a = 1.0 + 0j
    for k in range(1, n_terms + 1):
        a *= -((2 * k - 1) ** 2) / (k * 8.0 * 1j)
        s += a / z ** k
    return np.sqrt(2 / (np.pi * z)) * np.exp(1j * (z - np.pi / 4)) * s


# ------------------------------------------------- numerische Dispersion ----

def _disp_residual(u, N_lam, S, theta):
    lhs = np.sin(np.pi * S / (np.sqrt(2) * N_lam)) ** 2
    rhs = S ** 2 / 2 * (np.sin(u * np.cos(theta)) ** 2 + np.sin(u * np.sin(theta)) ** 2)
    return rhs - lhs


def u_numerisch(N_lam, S, theta):
    """Loest die 2D-Dispersionsrelation nach u = beta~ * Delta / 2 (Bisektion)."""
    skalar = np.isscalar(N_lam)
    N_lam = np.atleast_1d(np.asarray(N_lam, dtype=float))
    u_max = np.pi / (2 * max(abs(np.cos(theta)), abs(np.sin(theta)))) * (1 - 1e-12)
    out = np.empty_like(N_lam)
    for i, N in enumerate(N_lam):
        lo, hi = 1e-14, u_max
        f_lo = _disp_residual(lo, N, S, theta)
        f_hi = _disp_residual(hi, N, S, theta)
        if f_lo * f_hi > 0:
            out[i] = np.nan
            continue
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f_lo * _disp_residual(mid, N, S, theta) <= 0:
                hi = mid
            else:
                lo, f_lo = mid, _disp_residual(mid, N, S, theta)
        out[i] = 0.5 * (lo + hi)
    return out[0] if skalar else out


def cp_verhaeltnis(N_lam, S, theta):
    """c~_p / c fuer das 2D-Yee-Gitter, beliebiger Winkel."""
    u = u_numerisch(N_lam, S, theta)
    return np.pi / (np.asarray(N_lam, dtype=float) * u)


def cp_achse(N_lam, S):
    """Geschlossene Loesung fuer theta = 0."""
    N_lam = np.asarray(N_lam, dtype=float)
    u = np.arcsin(np.sqrt(2) / S * np.sin(np.pi * S / (np.sqrt(2) * N_lam)))
    return np.pi / (N_lam * u)


def cp_diagonal(N_lam, S):
    """Geschlossene Loesung fuer theta = 45 Grad."""
    N_lam = np.asarray(N_lam, dtype=float)
    u = np.sqrt(2) * np.arcsin(1.0 / S * np.sin(np.pi * S / (np.sqrt(2) * N_lam)))
    return np.pi / (N_lam * u)


def f_gitter(m, n, Lx, Ly, dx, dt, c):
    """Eigenfrequenz einer PEC-Mode, wie sie das Gitter zeigt.

    Die Wellenzahlen sind durch die Waende exakt vorgegeben; die
    Dispersionsrelation des Yee-Schemas liefert daraus die zugehoerige
    Frequenz.
    """
    bx, by = m * np.pi / Lx, n * np.pi / Ly
    arg = (c * dt) ** 2 * ((np.sin(bx * dx / 2) / (dx / 2)) ** 2
                           + (np.sin(by * dx / 2) / (dx / 2)) ** 2)
    return np.arcsin(0.5 * np.sqrt(arg)) / (np.pi * dt)


def f_kontinuum(m, n, Lx, Ly, c):
    return c / 2 * np.sqrt((m / Lx) ** 2 + (n / Ly) ** 2)


# ------------------------------------------------------------- Fresnel ------

def fresnel_r(eps1, eps2, sigma2=0.0, f=None):
    """Reflexionskoeffizient bei senkrechtem Einfall, optional mit Verlusten."""
    n1 = np.sqrt(complex(eps1))
    eps2c = complex(eps2) if sigma2 == 0 else eps2 - 1j * sigma2 / (2 * np.pi * f * EPS0)
    n2 = np.sqrt(eps2c)
    return (n1 - n2) / (n1 + n2)


# ----------------------------------------------------------- Selbsttest -----

if __name__ == "__main__":
    print("== Selbsttest der analytischen Hilfsfunktionen ==\n")

    # 1) Tabellierte Werte gegen die exakte Reihe in Bruchrechnung
    from fractions import Fraction as F
    import math

    def j0_exakt(z, terms=400):
        zz, s = F(z), F(0)
        for m in range(terms):
            s += F((-1) ** m) * (zz / 2) ** (2 * m) / F(math.factorial(m)) ** 2
        return float(s)

    def y0_exakt(z, terms=400):
        zz, s, H = F(z), F(0), F(0)
        for m in range(1, terms):
            H += F(1, m)
            s += F((-1) ** (m + 1)) * H * (zz / 2) ** (2 * m) / F(math.factorial(m)) ** 2
        return (2 / np.pi) * (np.log(z / 2) + EULER_GAMMA) * j0_exakt(z) + (2 / np.pi) * float(s)

    print("   z      J0 Fehler     Y0 Fehler   (gegen exakte Reihe in Bruchrechnung)")
    for z in (1, 2, 5, 10, 20, 30, 50, 80):
        print(f"  {z:5}   {abs(j0(np.array(float(z)))-j0_exakt(z)):.3e}    "
              f"{abs(y0(np.array(float(z)))-y0_exakt(z)):.3e}")

    # 2) Bessel-Differentialgleichung als unabhaengige Probe
    z = np.linspace(2.0, 60.0, 500)
    h = 1e-4
    for fun, name in ((j0, "J0"), (y0, "Y0")):
        f2 = (fun(z + h) - 2 * fun(z) + fun(z - h)) / h ** 2
        f1 = (fun(z + h) - fun(z - h)) / (2 * h)
        res = z ** 2 * f2 + z * f1 + z ** 2 * fun(z)
        print(f"\n  {name}: max |DGL-Residuum| / max|z^2 f| = "
              f"{np.max(np.abs(res))/np.max(np.abs(z**2*fun(z))):.2e}")

    # 3) Reihe gegen asymptotische Entwicklung
    z = np.linspace(2 * np.pi, 80.0, 300)
    rel = np.abs(hankel1_0(z) - hankel1_0_asym(z)) / np.abs(hankel1_0(z))
    print(f"\n  H0^(1): max. rel. Abweichung exakt <-> Asymptotik (z >= 2pi): {rel.max():.2e}")

    # 4) Dispersionsrelation
    print("\n  N_lam   Achse         Diagonale     Bisektion(0 Grad)")
    for N in (5, 10, 20, 40, 80, 1000):
        print(f"  {N:5}   {cp_achse(N,0.99):.9f}   {cp_diagonal(N,0.99):.9f}   "
              f"{cp_verhaeltnis(N,0.99,0.0):.9f}")
    print(f"\n  Bisektion(45 Grad) vs. geschlossen, N=20: "
          f"{cp_verhaeltnis(20,0.99,np.pi/4):.9f} / {cp_diagonal(20,0.99):.9f}")
    print(f"  magic time-step (S=1, Diagonale, N=20): {cp_diagonal(20,1.0):.12f}  (soll 1)")

    # 5) Konvergenzordnung der Dispersionsrelation
    Ns = np.array([10., 20., 40., 80.])
    err = 1 - cp_achse(Ns, 0.99)
    print(f"\n  Fehlerverhaeltnisse bei Verdopplung von N_lam (soll 4): "
          f"{err[:-1]/err[1:]}")
