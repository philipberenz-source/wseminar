// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Quelltext des Moduls `fdtd_core.py` <anhang-quelltext>

Nachstehend das Modul in der Fassung, mit der sämtliche Rechnungen der
Kapitel 5 und 6 durchgeführt wurden. Die einzelnen Stellen sind in Kapitel 4
erläutert. Nicht abgedruckt ist der angehängte
`if __name__ == "__main__":`-Block mit den drei Selbsttests
(Energieerhaltung, Vakuum-Laufzeit, Restreflexion der Schicht); er gehört
nicht zum Modul selbst und wird in Kapitel 5 aufgegriffen.

```python
import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 4e-7 * np.pi
C0 = 1.0 / np.sqrt(EPS0 * MU0)


class FDTD2D:
    """2D-FDTD-Loeser fuer den TM_z-Fall (Ez, Hx, Hy) auf einem quadratischen
    Yee-Gitter mit der einen Zellweite d (Kapitel 3 setzt dx = dy = d voraus).
    Es gilt mu = mu0 ueberall; eps und sigma sind
    ortsabhaengig und werden ueber set_material_box() gesetzt. Der Rand wird
    durch eine PML in gestreckten Koordinaten (rekursive Faltung) absorbiert."""

    def __init__(self, nx, ny, d, courant=0.99,
                 eps_r=None, sigma=None, npml=10, pml_m=3, r0=1e-6):

        self.nx, self.ny = nx, ny
        self.d = d

        self.dt = courant * self.d / (C0 * np.sqrt(2.0))

        self.eps_r = np.ones((nx, ny)) if eps_r is None else np.asarray(eps_r, dtype=float).copy()
        self.sigma = np.zeros((nx, ny)) if sigma is None else np.asarray(sigma, dtype=float).copy()
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
        """Additive Quelle: Ez(i,j) += waveform(t) nach jedem Zeitschritt.

        Der Zuschlag wird auf den bereits aktualisierten Wert Ez^(n+1) addiert
        und die Wellenform deshalb bei t^(n+1) ausgewertet."""
        self._soft_sources.append((i, j, waveform))

    def add_current_source(self, i, j, waveform, amplitude=1.0):
        """Reale Stromquelle: Jz(i,j) = amplitude * waveform(t) im Ampere-Maxwell-Gesetz.

        Jz geht als Jz^(n+1/2) in das Ez-Update ein; die Wellenform wird deshalb
        bei t^(n+1/2) ausgewertet."""
        self._current_sources.append((i, j, waveform, amplitude))

    def add_probe(self, i, j):
        """Registriert eine Sonde bei (i,j), die Zeit und Ez aufzeichnet."""
        self._probes.append({"i": i, "j": j, "t": [], "Ez": []})
        return len(self._probes) - 1

    def add_snapshot(self, steps):
        """Speichert das volle Ez-Feld bei den angegebenen Zeitschritt-Nummern."""
        self._snapshot_steps.update(steps)

    def get_spectrum(self, probe_index, pad=1):
        """Spektrum einer Sonde ueber die schnelle Fourier-Transformation.

        Die aufgezeichnete Zeitreihe wird mit numpy.fft.rfft transformiert; das
        Ergebnis liegt auf dem festen Raster f_k = k / (pad*N*dt), wobei pad >= 1
        Nullen anhaengt und das Raster verfeinert, ohne die Aufloesung 1/(N*dt)
        zu verbessern. Der Faktor dt normiert die Summe auf das Zeitintegral,
        der Exponentialfaktor beruecksichtigt, dass der erste Abtastwert zu
        t = dt und nicht zu t = 0 gehoert. Rueckgabe: (f, Y)."""
        p = self._probes[probe_index]
        y = np.asarray(p["Ez"], dtype=float)
        if y.size == 0:
            raise ValueError("Diese Sonde hat noch keinen Wert aufgezeichnet.")
        if pad < 1:
            raise ValueError("pad muss mindestens 1 sein.")
        n = int(pad) * y.size
        f = np.fft.rfftfreq(n, self.dt)
        Y = np.fft.rfft(y, n) * self.dt
        return f, Y * np.exp(-2j * np.pi * f * p["t"][0])

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

    def _inject_current_sources(self, t):
        """Setzt Jz zum Zeitpunkt t (aufgerufen mit t^(n+1/2))."""
        self.Jz[:, :] = 0.0
        for i, j, waveform, amplitude in self._current_sources:
            self.Jz[i, j] += amplitude * waveform(t)

    def _inject_soft_sources(self, t):
        """Addiert die Wellenform zum Zeitpunkt t auf Ez (aufgerufen mit t^(n+1))."""
        for i, j, waveform in self._soft_sources:
            self.Ez[i, j] += waveform(t)

    def _record_probes(self):
        for p in self._probes:
            ez = self.Ez[p["i"], p["j"]]
            p["t"].append(self.t)
            p["Ez"].append(ez)

    def _record_snapshots(self):
        if self.time_step in self._snapshot_steps:
            self._snapshots[self.time_step] = {"t": self.t, "Ez": self.Ez.copy()}

    def step(self):

        t_halb = self.t + 0.5 * self.dt
        t_neu = self.t + self.dt

        self._update_H()
        self._inject_current_sources(t_halb)
        self._update_E()
        self._inject_soft_sources(t_neu)
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
        self._snapshots.clear()
```
