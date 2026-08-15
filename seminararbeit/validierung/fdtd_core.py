import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 4e-7 * np.pi
C0 = 1.0 / np.sqrt(EPS0 * MU0)

class FDTD2D:

    def __init__(self, nx, ny, dx, dy, courant=0.99,
                 eps_r=None, sigma=None, npml=10, pml_m=3):
        if dx <= 0 or dy <= 0:
            raise ValueError("dx und dy muessen positiv sein.")
        if not (0.0 < courant < 1.0):
            raise ValueError("courant (Courant-Zahl S) muss im Intervall (0, 1) liegen.")
        if nx <= 2 * npml + 2 or ny <= 2 * npml + 2:
            raise ValueError(
                "nx/ny sind zu klein fuer die gewaehlte PML-Dicke npml "
                "(mindestens 2*npml + 3 Gitterpunkte je Richtung noetig)."
            )

        self.nx, self.ny = nx, ny 
        self.dx, self.dy = dx, dy 

        self.dt = courant / (C0 * np.sqrt(1.0 / dx**2 + 1.0 / dy**2))

        self.eps_r = np.ones((nx, ny)) if eps_r is None else np.asarray(eps_r, dtype=float).copy()
        self.sigma = np.zeros((nx, ny)) if sigma is None else np.asarray(sigma, dtype=float).copy()
        if np.any(self.eps_r < 1.0) or np.any(self.sigma < 0.0):
            raise ValueError("eps_r muss größer als eins sein und sigma muss größer als null sein")

        self.Ez = np.zeros((nx, ny))
        self.Hx = np.zeros((nx, ny - 1))
        self.Hy = np.zeros((nx - 1, ny))
        self.Jz = np.zeros((nx, ny))
        self._pec_mask = np.zeros((nx, ny), dtype=bool)

        self._update_material_coeffs()

        self.npml = npml
        self.pml_m = pml_m
        self._setup_pml()

        self._soft_sources = []
        self._current_sources = []
        self._probes = []
        self._snapshot_steps = set()
        self._snapshots = {}

        self.time_step = 0
        self.t = 0.0

    def get_grid(self):
        x = np.arange(self.nx) * self.dx
        y = np.arange(self.ny) * self.dy
        return np.meshgrid(x, y, indexing="ij")

    
    def _update_material_coeffs(self):
        loss_term = self.sigma * self.dt / (2.0 * EPS0 * self.eps_r)
        self.Ca = (1.0 - loss_term) / (1.0 + loss_term)
        self.Cb = (self.dt / (EPS0 * self.eps_r)) / (1.0 + loss_term)

    def set_material_box(self, i0, i1, j0, j1, eps_r=1.0, sigma=0.0):
        """Setzt eps_r, sigma auf dem Indexbereich [i0,i1) x [j0,j1) (z. B. Wand, Dielektrikum)."""
        self.eps_r[i0:i1, j0:j1] = eps_r
        self.sigma[i0:i1, j0:j1] = sigma
        self._update_material_coeffs()
    
    def set_pec_box(self, i0, i1, j0, j1):
        """Markiert [i0,i1) x [j0,j1) als ideal leitend (PEC): Ez wird dort auf 0 gehalten."""
        self._pec_mask[i0:i1, j0:j1] = True
        self.Ez[self._pec_mask] = 0.0


    def _pml_profile(self, n):
        depth = self.npml
        sigma_max = (self.pml_m + 1) / (150.0 * np.pi * min(self.dx, self.dy))
        idx = np.arange(n)
        d_left = depth - idx
        d_right = idx - (n - 1 - depth)
        d = np.maximum(np.maximum(d_left, d_right), 0)
        rho = np.minimum(d / depth, 1.0)
        sigma = sigma_max * rho**self.pml_m
        x = sigma * self.dt / (2.0 * EPS0)
        g3 = (1.0 - x) / (1.0 + x)
        g2 = 2.0 * x / (1.0 + x)
        return sigma, g2, g3


    def _setup_pml(self):
        _, self.gi2, self.gi3 = self._pml_profile(self.nx)
        _, self.gj2, self.gj3 = self._pml_profile(self.ny)
        _, self.fi2, self.fi3 = self._pml_profile(self.nx - 1)
        _, self.fj2, self.fj3 = self._pml_profile(self.ny - 1)

        self.psi_hx = np.zeros_like(self.Hx)
        self.psi_hy = np.zeros_like(self.Hy)
        self.psi_ez_x = np.zeros((self.nx - 2, self.ny - 2))
        self.psi_ez_y = np.zeros((self.nx - 2, self.ny - 2))


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
        p=self._probes[probe_index]
        if "freqs" not in p:
            raise ValueError("Diese sonde wurde ohne f angelegt")
        return p["freqs"], p["dft"]

    @staticmethod
    def gaussian_pulse(t, t0, tau):
        return np.exp(-((t - t0) / tau) ** 2)

    @staticmethod
    def ricker_wavelet(t, t0, fp):
        arg = (np.pi * fp * (t - t0)) ** 2
        return (1.0 - 2.0 * arg) * np.exp(-arg)


    def _update_H(self):
        dEz_dy = (self.Ez[:, 1:] - self.Ez[:, :-1]) / self.dy
        self.psi_hx = self.fj3[None, :] * self.psi_hx - self.fj2[None, :] * dEz_dy
        self.Hx -= (self.dt / MU0) * (dEz_dy + self.psi_hx)

        dEz_dx = (self.Ez[1:, :] - self.Ez[:-1, :]) / self.dx
        self.psi_hy = self.fi3[:, None] * self.psi_hy - self.fi2[:, None] * dEz_dx
        self.Hy += (self.dt / MU0) * (dEz_dx + self.psi_hy)


    def _update_E(self):
        curl_x = (self.Hy[1:, 1:-1] - self.Hy[:-1, 1:-1]) / self.dx
        curl_y = (self.Hx[1:-1, 1:] - self.Hx[1:-1, :-1]) / self.dy

        gi3 = self.gi3[1:-1, None]
        gi2 = self.gi2[1:-1, None]
        gj3 = self.gj3[None, 1:-1]
        gj2 = self.gj2[None, 1:-1]

        self.psi_ez_x = gi3 * self.psi_ez_x - gi2 * curl_x
        self.psi_ez_y = gj3 * self.psi_ez_y - gj2 * curl_y

        curlH = (curl_x + self.psi_ez_x) - (curl_y + self.psi_ez_y)

        self.Ez[1:-1, 1:-1] = (self.Ca[1:-1, 1:-1] * self.Ez[1:-1, 1:-1]
                                + self.Cb[1:-1, 1:-1] * (curlH - self.Jz[1:-1, 1:-1]))

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
        return (u_e + u_h) * self.dx * self.dy

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