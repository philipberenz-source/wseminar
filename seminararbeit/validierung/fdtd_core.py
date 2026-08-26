import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 4e-7 * np.pi
C0 = 1.0 / np.sqrt(EPS0 * MU0)


class FDTD2D:
    """2D-FDTD-Loeser fuer den TM_z-Fall (Ez, Hx, Hy) auf einem quadratischen
    Yee-Gitter (dx = dy = d). Es gilt mu = mu0 ueberall; eps und sigma sind
    ortsabhaengig und werden ueber set_material_box() gesetzt. Der Rand wird
    durch eine PML in gestreckten Koordinaten (rekursive Faltung) absorbiert."""

    def __init__(self, nx, ny, dx, dy, courant=0.99,
                 eps_r=None, sigma=None, npml=10, pml_m=3, r0=1e-6):
        if dx <= 0 or dy <= 0:
            raise ValueError("dx und dy muessen positiv sein.")
        if dx != dy:
            raise ValueError("Das Gitter muss quadratisch sein: dx muss gleich dy sein.")
        if not (0.0 < courant < 1.0):
            raise ValueError("courant (Courant-Zahl S) muss im Intervall (0, 1) liegen.")
        if nx <= 2 * npml + 2 or ny <= 2 * npml + 2:
            raise ValueError(
                "nx/ny sind zu klein fuer die gewaehlte PML-Dicke npml "
                "(mindestens 2*npml + 3 Gitterpunkte je Richtung noetig)."
            )

        self.nx, self.ny = nx, ny
        self.d = dx
        self.dx, self.dy = dx, dy  # Aliase, weil Auswertungsskripte sie lesen

        self.dt = courant * self.d / (C0 * np.sqrt(2.0))

        self.eps_r = np.ones((nx, ny)) if eps_r is None else np.asarray(eps_r, dtype=float).copy()
        self.sigma = np.zeros((nx, ny)) if sigma is None else np.asarray(sigma, dtype=float).copy()
        if np.any(self.eps_r < 1.0) or np.any(self.sigma < 0.0):
            raise ValueError("eps_r muss groesser als eins sein und sigma muss groesser als null sein.")

        self.Ez = np.zeros((nx, ny))
        self.Hx = np.zeros((nx, ny - 1))
        self.Hy = np.zeros((nx - 1, ny))
        self.Jz = np.zeros((nx, ny))
        self._pec_mask = np.zeros((nx, ny), dtype=bool)

        self._update_material_coeffs()

        self.npml = npml
        self.pml_m = pml_m
        self.r0 = r0
        self._setup_pml()

        self._soft_sources = []
        self._current_sources = []
        self._probes = []
        self._snapshot_steps = set()
        self._snapshots = {}

        self.time_step = 0
        self.t = 0.0

    def get_grid(self):
        """Liefert die (x, y)-Koordinaten aller Ez-Gitterpunkte (indexing='ij')."""
        x = np.arange(self.nx) * self.d
        y = np.arange(self.ny) * self.d
        return np.meshgrid(x, y, indexing="ij")

    def _update_material_coeffs(self):
        """Berechnet die Koeffizienten Ca, Cb des Ez-Updates aus eps_r und sigma.
        Cb enthaelt bereits den Faktor 1/d, siehe _update_E."""
        loss = self.sigma * self.dt / (2.0 * EPS0 * self.eps_r)
        self.Ca = (1.0 - loss) / (1.0 + loss)
        self.Cb = (self.dt / (EPS0 * self.eps_r * self.d)) / (1.0 + loss)

    def set_material_box(self, i0, i1, j0, j1, eps_r=1.0, sigma=0.0):
        """Setzt eps_r, sigma auf dem Indexbereich [i0,i1) x [j0,j1) (z. B. Wand, Dielektrikum)."""
        self.eps_r[i0:i1, j0:j1] = eps_r
        self.sigma[i0:i1, j0:j1] = sigma
        self._update_material_coeffs()

    def set_pec_box(self, i0, i1, j0, j1):
        """Markiert [i0,i1) x [j0,j1) als ideal leitend (PEC): Ez wird dort auf 0 gehalten."""
        self._pec_mask[i0:i1, j0:j1] = True
        self.Ez[self._pec_mask] = 0.0

    def _pml_koeffizienten(self, positionen, n_zellen):
        """Rekursionskoeffizienten (b, C) der PML an den gegebenen Positionen
        (in Zelleinheiten, auch halbzahlig). n_zellen ist die Anzahl der
        Gitterpunkte des zugehoerigen ganzzahligen Gitters. tiefe misst den
        Abstand zum inneren Rand der Schicht, rho normiert ihn auf [0, 1]."""
        if self.npml == 0:
            return np.ones_like(positionen), np.zeros_like(positionen)
        d_links = self.npml - positionen
        d_rechts = positionen - (n_zellen - 1 - self.npml)
        tiefe = np.maximum(np.maximum(d_links, d_rechts), 0.0)
        rho = np.minimum(tiefe / self.npml, 1.0)
        sig = self.sigma_max * rho ** self.pml_m
        b = np.exp(-sig * self.dt / EPS0)
        return b, b - 1.0

    def _setup_pml(self):
        """Legt sigma_max aus der Zielreflexion r0 fest und berechnet daraus die
        vier PML-Profile (eines je Feldversatz) sowie die Hilfsfelder der
        rekursiven Faltung."""
        if self.npml > 0:
            eta0 = np.sqrt(MU0 / EPS0)
            D = self.npml * self.d
            self.sigma_max = -(self.pml_m + 1) * np.log(self.r0) / (2.0 * eta0 * D)
        else:
            self.sigma_max = 0.0

        nx, ny = self.nx, self.ny
        self.bx_i, self.Cx_i = self._pml_koeffizienten(np.arange(nx), nx)
        self.by_j, self.Cy_j = self._pml_koeffizienten(np.arange(ny), ny)
        self.bx_h, self.Cx_h = self._pml_koeffizienten(np.arange(nx - 1) + 0.5, nx)
        self.by_h, self.Cy_h = self._pml_koeffizienten(np.arange(ny - 1) + 0.5, ny)

        self.psi_hx = np.zeros_like(self.Hx)
        self.psi_hy = np.zeros_like(self.Hy)
        self.psi_ez_x = np.zeros((nx - 2, ny - 2))
        self.psi_ez_y = np.zeros((nx - 2, ny - 2))

    def add_soft_source(self, i, j, waveform):
        """Additive Quelle: Ez(i,j) += waveform(t) nach jedem Zeitschritt."""
        self._soft_sources.append((i, j, waveform))

    def add_current_source(self, i, j, waveform, amplitude=1.0):
        """Reale Stromquelle: Jz(i,j) = amplitude * waveform(t) im Ampere-Maxwell-Gesetz."""
        self._current_sources.append((i, j, waveform, amplitude))

    def add_probe(self, i, j, freqs=None):
        """Registriert eine Sonde bei (i,j). Optional: freqs (Hz) fuer eine laufende DFT,
        z. B. sim.add_probe(i, j, freqs=[2.4e9, 5.0e9]) fuer WLAN-Frequenzen."""
        probe = {"i": i, "j": j, "t": [], "Ez": []}
        if freqs is not None:
            probe["freqs"] = np.asarray(freqs, dtype=float)
            probe["dft"] = np.zeros(probe["freqs"].shape, dtype=complex)
        self._probes.append(probe)
        return len(self._probes) - 1

    def add_snapshot(self, steps):
        """Speichert das volle Ez-Feld bei den angegebenen Zeitschritt-Nummern."""
        self._snapshot_steps.update(steps)

    def get_spectrum(self, probe_index):
        """Gibt (freqs, dft) der Sonde probe_index zurueck; die Sonde muss mit
        freqs angelegt worden sein."""
        p = self._probes[probe_index]
        if "freqs" not in p:
            raise ValueError("Diese Sonde wurde ohne freqs angelegt.")
        return p["freqs"], p["dft"]

    @staticmethod
    def gaussian_pulse(t, t0, tau):
        return np.exp(-((t - t0) / tau) ** 2)

    @staticmethod
    def ricker_wavelet(t, t0, fp):
        arg = (np.pi * fp * (t - t0)) ** 2
        return (1.0 - 2.0 * arg) * np.exp(-arg)

    def _update_H(self):
        """Aktualisiert Hx und Hy aus Ez (halben Zeitschritt vor Ez)."""
        d = self.d

        dEz_dy = (self.Ez[:, 1:] - self.Ez[:, :-1]) / d
        self.psi_hx = self.by_h[None, :] * self.psi_hx + self.Cy_h[None, :] * dEz_dy
        self.Hx -= (self.dt / MU0) * (dEz_dy + self.psi_hx)

        dEz_dx = (self.Ez[1:, :] - self.Ez[:-1, :]) / d
        self.psi_hy = self.bx_h[:, None] * self.psi_hy + self.Cx_h[:, None] * dEz_dx
        self.Hy += (self.dt / MU0) * (dEz_dx + self.psi_hy)

    def _update_E(self):
        """Aktualisiert Ez auf dem Inneren [1:-1, 1:-1] aus Hx, Hy und Jz und
        setzt danach die PEC-Zellen zurueck. Die aeusserste Zeile/Spalte von Ez
        wird nie beschrieben und bleibt damit dauerhaft null (PEC-Abschluss
        hinter der PML)."""
        d = self.d

        dHy = self.Hy[1:, 1:-1] - self.Hy[:-1, 1:-1]      # nackte Differenz
        dHx = self.Hx[1:-1, 1:] - self.Hx[1:-1, :-1]      # nackte Differenz

        self.psi_ez_x = self.bx_i[1:-1, None] * self.psi_ez_x + self.Cx_i[1:-1, None] * (dHy / d)
        self.psi_ez_y = self.by_j[None, 1:-1] * self.psi_ez_y + self.Cy_j[None, 1:-1] * (dHx / d)

        self.Ez[1:-1, 1:-1] = (self.Ca[1:-1, 1:-1] * self.Ez[1:-1, 1:-1]
                                + self.Cb[1:-1, 1:-1] * (dHy - dHx
                                                          + d * (self.psi_ez_x - self.psi_ez_y)
                                                          - self.Jz[1:-1, 1:-1] * d))

        if self._pec_mask.any():
            self.Ez[self._pec_mask] = 0.0

    def _inject_current_sources(self):
        self.Jz[:, :] = 0.0
        for i, j, waveform, amplitude in self._current_sources:
            self.Jz[i, j] = amplitude * waveform(self.t)

    def _inject_soft_sources(self):
        for i, j, waveform in self._soft_sources:
            self.Ez[i, j] += waveform(self.t)

    def _record_probes(self):
        for p in self._probes:
            ez = self.Ez[p["i"], p["j"]]
            p["t"].append(self.t)
            p["Ez"].append(ez)
            if "freqs" in p:
                p["dft"] += ez * np.exp(-1j * 2.0 * np.pi * p["freqs"] * self.t) * self.dt

    def _record_snapshots(self):
        if self.time_step in self._snapshot_steps:
            self._snapshots[self.time_step] = {"t": self.t, "Ez": self.Ez.copy()}

    def step(self):
        self._update_H()
        self._inject_current_sources()
        self._update_E()
        self._inject_soft_sources()
        self.time_step += 1
        self.t = self.time_step * self.dt
        self._record_probes()
        self._record_snapshots()

    def run(self, n_steps, callback=None):
        for _ in range(n_steps):
            self.step()
            if callback is not None:
                callback(self)

    def total_energy(self):
        u_e = np.sum(0.5 * EPS0 * self.eps_r * self.Ez**2)
        u_h = np.sum(0.5 * MU0 * self.Hx**2) + np.sum(0.5 * MU0 * self.Hy**2)
        return (u_e + u_h) * self.d**2

    def reset_fields(self):
        """Setzt Felder, Zeit und aufgezeichnete Daten zurueck; Gitter, Material, PML und
        registrierte Quellen/Sonden bleiben erhalten (nuetzlich fuer Parameterstudien)."""
        for arr in (self.Ez, self.Hx, self.Hy, self.Jz,
                    self.psi_hx, self.psi_hy, self.psi_ez_x, self.psi_ez_y):
            arr[...] = 0.0
        self.time_step = 0
        self.t = 0.0
        for p in self._probes:
            p["t"].clear()
            p["Ez"].clear()
            if "freqs" in p:
                p["dft"][:] = 0.0
        self._snapshots.clear()


if __name__ == "__main__":
    # Drei einfache Selbsttests: Energieerhaltung in der geschlossenen,
    # verlustlosen Box, Laufzeit im Vakuum und Restreflexion an der PML.
    # Kein Testframework, nur Zahlen und ein bestanden/nicht bestanden.

    def test_energieerhaltung():
        nx = ny = 80
        d = 1.0e-3
        sim = FDTD2D(nx, ny, d, d, courant=0.99, npml=0)  # ohne PML: Rand ist
        # bereits PEC, da die aeusserste Zeile/Spalte von Ez nie aktualisiert wird.

        tau = 40 * sim.dt
        t0 = 8 * tau
        ic, jc = nx // 2, ny // 2
        sim.add_soft_source(ic, jc, lambda t: FDTD2D.gaussian_pulse(t, t0, tau))

        n_start = int(round((t0 + 10 * tau) / sim.dt))  # Puls dann abgeklungen
        n_nach = 2400
        energie = np.empty(n_start + n_nach)
        for n in range(len(energie)):
            sim.step()
            energie[n] = sim.total_energy()

        # total_energy() mischt Ez^n mit H^(n-1/2): dieser halbe Zeitschritt
        # Versatz erzeugt eine kleine, aber nicht abklingende Welligkeit der
        # Momentanenergie (mehrere Prozent), die keine echte Drift ist. Um
        # trotzdem eine echte Drift (z. B. durch einen Fehler im Update) zu
        # erkennen, wird ueber je eine Haelfte des Beobachtungsfensters
        # gemittelt; die Welligkeit mittelt sich dabei heraus.
        fenster = energie[n_start:]
        haelfte = len(fenster) // 2
        e_anfang = fenster[:haelfte].mean()
        e_ende = fenster[haelfte:].mean()
        drift = abs(e_ende - e_anfang) / ((e_anfang + e_ende) / 2)
        bestanden = drift < 1.0e-3
        print(f"Testfall 1 - Energieerhaltung: Drift zwischen erster und zweiter "
              f"Haelfte des Beobachtungsfensters {drift * 100:.4f} % "
              f"({'bestanden' if bestanden else 'NICHT bestanden'})")
        return bestanden

    def test_ausbreitung():
        nx = ny = 220
        d = 1.0e-3
        sim = FDTD2D(nx, ny, d, d, courant=0.99, npml=10)

        tau = 5 * sim.dt
        t0 = 8 * tau
        ic, jc = nx // 2, ny // 2
        sim.add_soft_source(ic, jc, lambda t: FDTD2D.gaussian_pulse(t, t0, tau))

        m = 90
        probe = sim.add_probe(ic + m, jc)
        r = m * d

        n_steps = int(round((t0 + r / C0 + 30 * tau) / sim.dt))
        sim.run(n_steps)

        # Die Punktquelle strahlt in 2D kreisfoermig ab; Laufzeit wird deshalb
        # am Maximum des Pulses gemessen (nicht an einer Amplitudenschwelle,
        # die durch die Abklingform der Zylinderwelle verzerrt waere).
        ez = np.abs(np.array(sim._probes[probe]["Ez"]))
        t_arr = np.array(sim._probes[probe]["t"])
        t_maximum = t_arr[np.argmax(ez)]

        laufzeit_gemessen = t_maximum - t0
        laufzeit_theorie = r / C0
        abweichung = abs(laufzeit_gemessen - laufzeit_theorie) / laufzeit_theorie
        bestanden = abweichung < 0.02
        print(f"Testfall 2 - Vakuum-Ausbreitung: Laufzeit gemessen "
              f"{laufzeit_gemessen * 1e9:.3f} ns, Theorie {laufzeit_theorie * 1e9:.3f} ns, "
              f"Abweichung {abweichung * 100:.2f} % "
              f"({'bestanden' if bestanden else 'NICHT bestanden'})")
        return bestanden

    def test_pml_reflexion():
        d = 1.0e-3
        fp = C0 / (40 * d)  # Wellenlaenge = 40 d, damit die Schicht (10 Zellen) gut aufgeloest ist
        t0 = 3.0 / fp
        npml = 10

        def bauen(n):
            sim = FDTD2D(n, n, d, d, courant=0.99, npml=npml)
            c_ = n // 2
            sim.add_soft_source(c_, c_, lambda t: FDTD2D.ricker_wavelet(t, t0, fp))
            probe = sim.add_probe(c_ + 15, c_ + 15)  # diagonal versetzt
            return sim, probe

        n_steps = 900
        sim_klein, p_klein = bauen(60)
        sim_gross, p_gross = bauen(260)
        sim_klein.run(n_steps)
        sim_gross.run(n_steps)

        ez_klein = np.array(sim_klein._probes[p_klein]["Ez"])
        ez_gross = np.array(sim_gross._probes[p_gross]["Ez"])
        amp_referenz = np.abs(ez_gross).max()
        abweichung = np.abs(ez_klein - ez_gross).max() / amp_referenz
        bestanden = abweichung < 0.10
        print(f"Testfall 3 - PML-Restreflexion: max. Abweichung {abweichung * 100:.2f} % "
              f"der Referenzamplitude ({'bestanden' if bestanden else 'NICHT bestanden'})")
        return bestanden

    print("Selbsttest fuer fdtd_core.py\n")
    ok1 = test_energieerhaltung()
    ok2 = test_ausbreitung()
    ok3 = test_pml_reflexion()
    print()
    print("Alle Tests bestanden." if (ok1 and ok2 and ok3) else "Mindestens ein Test NICHT bestanden.")
