// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Implementierung in Python

Dieses Kapitel geht den Quelltext des Moduls `fdtd_core.py` vollständig durch
und beschreibt, was die einzelnen Abschnitte tun. Die zugrunde liegenden
Gleichungen stammen aus Kapitel 3 und werden hier nicht erneut behandelt.

== Weshalb Python

Python ist in den Naturwissenschaften weit verbreitet, vor allem wegen seiner
Bibliotheken: NumPy rechnet auf ganzen Zahlenfeldern, matplotlib erzeugt daraus
die Abbildungen. Beides geschieht in derselben Sprache, ein Übersetzungsschritt
entfällt.

NumPy ist hier nicht nur bequem, sondern notwendig. Eine Schleife über alle
Gitterzellen wäre in Python rund zweihundertmal langsamer als dieselbe Rechnung
als Array-Operation. Das Modul enthält deshalb an keiner Stelle eine Schleife
über Gitterzellen; jede Rechnung wirkt auf ganze Felder zugleich.

== Aufbau des Moduls und Konstruktor

Zu Beginn stehen die Naturkonstanten:

```python
import numpy as np

EPS0 = 8.8541878128e-12
MU0 = 4e-7 * np.pi
C0 = 1.0 / np.sqrt(EPS0 * MU0)
```

Der übrige Quelltext steht in der Klasse `FDTD2D`. Alle Daten einer Simulation
--- Felder, Material, Quellen, Sonden, Zeitzähler --- gehören damit zur
jeweiligen Instanz, sodass mehrere Simulationen nebeneinander laufen können.
Methoden mit führendem Unterstrich sind zur internen Verwendung gedacht.

Der Konstruktor prüft zuerst die übergebenen Werte und bricht bei unzulässigen
Eingaben mit einer Meldung ab:

```python
def __init__(self, nx, ny, dx, dy, courant=0.99,
             eps_r=None, sigma=None, npml=10, pml_m=3):
    if dx <= 0 or dy <= 0:
        raise ValueError("dx und dy muessen positiv sein.")
    if not (0.0 < courant < 1.0):
        raise ValueError("courant (Courant-Zahl S) muss im Intervall (0, 1) liegen.")
    if nx <= 2 * npml + 2 or ny <= 2 * npml + 2:
        raise ValueError("nx/ny sind zu klein fuer die gewaehlte PML-Dicke npml ...")
```

Danach werden der Zeitschritt aus der Stabilitätsbedingung und die
Materialfelder angelegt:

```python
self.dt = courant / (C0 * np.sqrt(1.0 / dx**2 + 1.0 / dy**2))

self.eps_r = np.ones((nx, ny)) if eps_r is None else np.asarray(eps_r, dtype=float).copy()
self.sigma = np.zeros((nx, ny)) if sigma is None else np.asarray(sigma, dtype=float).copy()
```

Ohne Angabe wird Vakuum angenommen. Übergibt der Aufrufer eigene Werte, legt
`.copy()` eine eigene Kopie an, damit Modul und Aufrufer nicht dasselbe Array
verändern.

== Die Feldgrößen als Arrays

```python
self.Ez = np.zeros((nx, ny))
self.Hx = np.zeros((nx, ny - 1))
self.Hy = np.zeros((nx - 1, ny))
self.Jz = np.zeros((nx, ny))
self._pec_mask = np.zeros((nx, ny), dtype=bool)
```

`Ez` ist das elektrische Feld, `Hx` und `Hy` sind die beiden magnetischen
Komponenten, `Jz` die eingeprägte Stromdichte. Dass die magnetischen Felder je
einen Eintrag weniger besitzen, folgt aus dem Versatz des Gitters (Abschnitt
3.2). `_pec_mask` merkt sich, welche Zellen ideal leitend sind.

In Kapitel 3 tragen die Feldgrößen halbzahlige Indizes wie
$H_x^(n+1\/2)(i, j+1\/2)$. Solche Indizes verschwinden im Programm vollständig,
denn ein Array lässt sich ausschließlich ganzzahlig ansprechen: Es gibt keinen
Platz `Hx[i, j+0.5]`. `Hx[i, j]` bezeichnet deshalb denjenigen Wert, der in der
Schreibweise des vorigen Kapitels bei $(i, j+1\/2)$ liegt. Der räumliche Versatz
ist damit nicht verschwunden, sondern nur unsichtbar geworden: Er steckt nicht
mehr im Index, sondern darin, in _welchem_ Array ein Wert steht.

Für den Zeitindex gilt dasselbe, und dort hat es eine Folge. Ein halber
Zeitschritt lässt sich so wenig zählen wie ein halber Arrayplatz, weshalb
`self.Hx` schlicht den zuletzt berechneten Stand des Magnetfelds enthält ---
gleichgültig, ob dieser in der Schreibweise aus Kapitel 3 auf $n - 1\/2$ oder
auf $n + 1\/2$ liegt. Der zeitliche Versatz steckt folglich weder im Index noch
im Array, sondern allein in der Reihenfolge, in der die beiden Felder
fortgeschrieben werden. Eben deshalb ist diese Reihenfolge nicht beliebig,
worauf Abschnitt 4.8 zurückkommt.

== Ableitungen als Array-Ausschnitte

Ein Differenzenquotient ist die Differenz zweier gegeneinander verschobener
Ausschnitte desselben Arrays: `Ez[:, 1:]` sind alle Werte ab dem zweiten
Eintrag, `Ez[:, :-1]` alle bis auf den letzten.

```python
def _update_H(self):
    dEz_dy = (self.Ez[:, 1:] - self.Ez[:, :-1]) / self.dy
    self.psi_hx = self.fj3[None, :] * self.psi_hx - self.fj2[None, :] * dEz_dy
    self.Hx -= (self.dt / MU0) * (dEz_dy + self.psi_hx)

    dEz_dx = (self.Ez[1:, :] - self.Ez[:-1, :]) / self.dx
    self.psi_hy = self.fi3[:, None] * self.psi_hy - self.fi2[:, None] * dEz_dx
    self.Hy += (self.dt / MU0) * (dEz_dx + self.psi_hy)
```

Beide Magnetfeldkomponenten werden aus der jeweiligen Ableitung von `Ez`
fortgeschrieben. `psi_hx` und `psi_hy` sind die Hilfsgrößen der absorbierenden
Randschicht; sie bleiben im Gebietsinneren null und wirken nur am Rand. Das
Kürzel `[None, :]` verteilt ein eindimensionales Profil auf das ganze Feld.

```python
def _update_E(self):
    curl_x = (self.Hy[1:, 1:-1] - self.Hy[:-1, 1:-1]) / self.dx
    curl_y = (self.Hx[1:-1, 1:] - self.Hx[1:-1, :-1]) / self.dy
    ...
    self.Ez[1:-1, 1:-1] = (self.Ca[1:-1, 1:-1] * self.Ez[1:-1, 1:-1]
                            + self.Cb[1:-1, 1:-1] * (curlH - self.Jz[1:-1, 1:-1]))

    if self._pec_mask.any():
        self.Ez[self._pec_mask] = 0.0
```

`curl_x` und `curl_y` sind die beiden Anteile der Rotation von $bold(H)$, aus
denen zusammen mit der Stromdichte das neue $E_z$ folgt. Der Ausschnitt
`[1:-1, 1:-1]` lässt die äußerste Zeile und Spalte aus, weil dort ein Nachbar
fehlt. Zuletzt werden alle als leitend markierten Zellen auf null gesetzt.

== Vorberechnete Koeffizienten und Randprofile

```python
def _update_material_coeffs(self):
    loss_term = self.sigma * self.dt / (2.0 * EPS0 * self.eps_r)
    self.Ca = (1.0 - loss_term) / (1.0 + loss_term)
    self.Cb = (self.dt / (EPS0 * self.eps_r)) / (1.0 + loss_term)
```

`Ca` und `Cb` sind die Faktoren des $E$-Updates. Sie hängen nur von Material
und Schrittweite ab und werden deshalb vorab für jede Zelle berechnet; die
Zeitschleife braucht dann keine Fallunterscheidung nach Material mehr.

```python
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
```

Die Methode berechnet den Verlauf der Dämpfung über die Breite des Gitters.
`d` ist der Abstand zum näheren der beiden Ränder --- `np.maximum` wählt
elementweise den größeren Wert und setzt alles im Inneren auf null. Daraus
folgen die Koeffizienten `g2` und `g3`, die im Inneren null beziehungsweise
eins betragen.

```python
def _setup_pml(self):
    _, self.gi2, self.gi3 = self._pml_profile(self.nx)
    _, self.gj2, self.gj3 = self._pml_profile(self.ny)
    _, self.fi2, self.fi3 = self._pml_profile(self.nx - 1)
    _, self.fj2, self.fj3 = self._pml_profile(self.ny - 1)

    self.psi_hx = np.zeros_like(self.Hx)
    self.psi_hy = np.zeros_like(self.Hy)
    self.psi_ez_x = np.zeros((self.nx - 2, self.ny - 2))
    self.psi_ez_y = np.zeros((self.nx - 2, self.ny - 2))
```

Vier Profile sind nötig, weil die versetzten Felder unterschiedlich lang sind.
Anschließend werden die vier Hilfsgrößen mit Nullen angelegt.

== Material und Geometrie

```python
def set_material_box(self, i0, i1, j0, j1, eps_r=1.0, sigma=0.0):
    self.eps_r[i0:i1, j0:j1] = eps_r
    self.sigma[i0:i1, j0:j1] = sigma
    self._update_material_coeffs()

def set_pec_box(self, i0, i1, j0, j1):
    self._pec_mask[i0:i1, j0:j1] = True
    self.Ez[self._pec_mask] = 0.0
```

Die erste Methode belegt einen rechteckigen Bereich mit einem Material und
berechnet danach die Koeffizienten neu. Die zweite markiert einen Bereich als
ideal leitend.

```python
def get_grid(self):
    x = np.arange(self.nx) * self.dx
    y = np.arange(self.ny) * self.dy
    return np.meshgrid(x, y, indexing="ij")
```

Liefert die Koordinaten aller Gitterpunkte. `indexing="ij"` legt fest, dass der
erste Index zur $x$-Richtung gehört.

== Quellen

```python
def add_soft_source(self, i, j, waveform):
    self._soft_sources.append((i, j, waveform))

def add_current_source(self, i, j, waveform, amplitude=1.0):
    self._current_sources.append((i, j, waveform, amplitude))

def _inject_soft_sources(self):
    for i, j, waveform in self._soft_sources:
        self.Ez[i, j] += waveform(self.t)

def _inject_current_sources(self):
    self.Jz[:, :] = 0.0
    for i, j, waveform, amplitude in self._current_sources:
        self.Jz[i, j] = amplitude * waveform(self.t)
```

Eine Quelle wird durch ihren Ort und eine Wellenform festgelegt. Die Wellenform
ist eine Funktion, die in jedem Zeitschritt mit der aktuellen Zeit aufgerufen
wird; dadurch lassen sich beliebige Anregungen übergeben, ohne das Modul zu
ändern. Die additive Quelle addiert ihren Wert auf das Feld, die Stromquelle
setzt die Stromdichte, die zuvor überall auf null gesetzt wird.

```python
@staticmethod
def gaussian_pulse(t, t0, tau):
    return np.exp(-((t - t0) / tau) ** 2)

@staticmethod
def ricker_wavelet(t, t0, fp):
    arg = (np.pi * fp * (t - t0)) ** 2
    return (1.0 - 2.0 * arg) * np.exp(-arg)
```

Zwei fertige Pulsformen, die als Wellenform übergeben werden können.

== Zeitschleife

```python
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
```

`step` führt einen vollständigen Zeitschritt aus: erst das Magnetfeld, dann die
Stromquelle, dann das elektrische Feld, dann die additive Quelle. Diese Abfolge
ist keine Setzung des Programms, sondern liegt bereits in der Zuordnung aus
Abschnitt 3.3 fest, die $E$ auf die ganzzahligen und $bold(H)$ auf die
halbzahligen Zeitpunkte legt. Da diese Halbschritte im Programm nur noch als
Reihenfolge vorliegen, muss sie genau das leisten, was dort die halbzahligen
Zeitindizes leisten: Das $H$-Update benötigt $E$ zum ganzzahligen Zeitpunkt $n$, das
anschließende $E$-Update benötigt $H$ zum dazwischenliegenden Zeitpunkt
$n + 1\/2$. Jeder der beiden Schritte greift also auf denjenigen Stand des
anderen Feldes zu, den der unmittelbar vorangegangene Schritt gerade erzeugt
hat. Vertauscht man beide, so rechnet das $E$-Update mit einem $H$, das einen
vollen Zeitschritt alt ist; die Differenz wäre dann nicht mehr um den
auszuwertenden Zeitpunkt zentriert, und mit der Zentrierung entfiele gerade die
Eigenschaft, aus der in Abschnitt 3.1 die zweite Fehlerordnung folgte.

Aus demselben Grund stehen die beiden Quellen an verschiedenen Stellen. Die
Stromdichte `Jz` wird innerhalb von `_update_E` gelesen und muss deshalb vorher
gesetzt sein. Die additive Quelle addiert dagegen unmittelbar auf `Ez` und steht
folglich hinter dem Update --- käme sie davor, ginge ihr Beitrag als _alter_
Feldwert in die Rechnung ein und würde dabei zusätzlich mit `Ca` multipliziert.
Erst danach werden Zeit und Zähler erhöht und die Messgrößen aufgezeichnet, die
somit stets ein vollständig fortgeschriebenes Feld sehen. Die Zeit wird aus dem
Zähler berechnet und nicht aufaddiert, damit sich keine Rundungsfehler
ansammeln. `run` wiederholt den Schritt und ruft dabei eine optionale Funktion
`callback` auf, über die sich von außen auf den Zustand zugreifen lässt.

== Aufzeichnung

```python
def add_probe(self, i, j, freqs=None):
    probe = {"i": i, "j": j, "t": [], "Ez": []}
    if freqs is not None:
        probe["freqs"] = np.asarray(freqs, dtype=float)
        probe["dft"] = np.zeros(probe["freqs"].shape, dtype=complex)
    self._probes.append(probe)
    return len(self._probes) - 1

def _record_probes(self):
    for p in self._probes:
        ez = self.Ez[p["i"], p["j"]]
        p["t"].append(self.t)
        p["Ez"].append(ez)
        if "freqs" in p:
            p["dft"] += ez * np.exp(-1j * 2.0 * np.pi * p["freqs"] * self.t) * self.dt
```

Eine Sonde zeichnet an einem festen Gitterpunkt den Verlauf des Feldes auf.
Werden zusätzlich Frequenzen angegeben, wird in jedem Zeitschritt auch die
Fourier-Transformierte an diesen Frequenzen fortgeschrieben, ohne dass am Ende
die ganze Zeitreihe transformiert werden müsste. Zurückgegeben wird der Index,
über den die Sonde später ansprechbar ist.

```python
def get_spectrum(self, probe_index):
    p=self._probes[probe_index]
    if "freqs" not in p:
        raise ValueError("Diese sonde wurde ohne f angelegt")
    return p["freqs"], p["dft"]

def add_snapshot(self, steps):
    self._snapshot_steps.update(steps)

def _record_snapshots(self):
    if self.time_step in self._snapshot_steps:
        self._snapshots[self.time_step] = {"t": self.t, "Ez": self.Ez.copy()}
```

`get_spectrum` gibt Frequenzen und zugehörige Amplituden einer Sonde zurück.
`add_snapshot` merkt sich Zeitschritte, zu denen das vollständige Feld
gespeichert werden soll; `_record_snapshots` legt dort eine Kopie ab --- ohne
`.copy()` würde nur ein Verweis gespeichert und im nächsten Schritt
überschrieben.

```python
def total_energy(self):
    u_e = np.sum(0.5 * EPS0 * self.eps_r * self.Ez**2)
    u_h = np.sum(0.5 * MU0 * self.Hx**2) + np.sum(0.5 * MU0 * self.Hy**2)
    return (u_e + u_h) * self.dx * self.dy

def reset_fields(self):
    for arr in (self.Ez, self.Hx, self.Hy, self.Jz,
                self.psi_hx, self.psi_hy, self.psi_ez_x, self.psi_ez_y):
        arr[...] = 0.0
    self.time_step = 0
    self.t = 0.0
    for p in self._probes:
        p["t"].clear()
        p["Ez"].clear()
```

`total_energy` summiert die im Gitter gespeicherte Energie; elektrischer und
magnetischer Anteil werden getrennt aufsummiert, weil die Felder
unterschiedliche Größen haben. `reset_fields` setzt Felder, Zeit und
aufgezeichnete Daten zurück, behält aber Gitter, Material und Sonden. Die
Arrays werden dabei mit Nullen überschrieben statt neu angelegt, sonst zeigte
die Klasse weiterhin auf die alten.

Damit ist der Quelltext vollständig beschrieben. Das Modul umfasst rund
zweihundertdreißig Zeilen, erzeugt selbst keine Abbildungen und stellt nur die
Rechnung und die Daten bereit, aus denen die Auswertungen des folgenden
Kapitels entstehen.
