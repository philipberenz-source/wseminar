// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

// --- Bausteine für das Klassendiagramm in Abschnitt 4.1 ------------------
// Absolute Koordinaten in cm; wird ein Kasten höher, sind die y-Werte der
// folgenden Elemente mit anzupassen.
#let kl-strich = 0.6pt + luma(70)

#let kl-kasten(x, y, w, h, fuellung: white) = place(top + left, dx: x, dy: y,
  rect(width: w, height: h, stroke: kl-strich, fill: fuellung, radius: 2pt))

#let kl-inhalt(x, y, w, h, body) = place(top + left, dx: x, dy: y,
  box(width: w, height: h, inset: (x: 7pt, y: 6pt))[#align(left + top)[#body]])

#let kl-mitte(x, y, w, h, body) = place(top + left, dx: x, dy: y,
  box(width: w, height: h)[#align(center + horizon)[#body]])

#let kl-linie(x1, y1, x2, y2) = place(top + left,
  line(start: (x1, y1), end: (x2, y2), stroke: kl-strich))

#let kl-kopf-r(x, y) = place(top + left, dx: x - 0.17cm, dy: y - 0.105cm,
  polygon(fill: luma(70), stroke: none, (0cm, 0cm), (0cm, 0.21cm), (0.17cm, 0.105cm)))
#let kl-kopf-l(x, y) = place(top + left, dx: x, dy: y - 0.105cm,
  polygon(fill: luma(70), stroke: none, (0.17cm, 0cm), (0.17cm, 0.21cm), (0cm, 0.105cm)))
#let kl-kopf-ab(x, y) = place(top + left, dx: x - 0.105cm, dy: y - 0.16cm,
  polygon(fill: luma(70), stroke: none, (0cm, 0cm), (0.21cm, 0cm), (0.105cm, 0.16cm)))

#let kl-pfeil-ab(x, y0, y1) = { kl-linie(x, y0, x, y1 - 0.13cm); kl-kopf-ab(x, y1) }
#let kl-pfeil-r(x0, x1, y) = { kl-linie(x0, y, x1 - 0.14cm, y); kl-kopf-r(x1, y) }
#let kl-pfeil-l(x0, x1, y) = { kl-linie(x0, y, x1 + 0.14cm, y); kl-kopf-l(x1, y) }
#let kl-pfeil-lr(x0, x1, y) = {
  kl-linie(x0 + 0.14cm, y, x1 - 0.14cm, y); kl-kopf-r(x1, y); kl-kopf-l(x0, y)
}
#let kl-label(x, y, body) = place(top + left, dx: x, dy: y,
  text(size: 7pt, fill: luma(60), style: "italic")[#body])
= Implementierung in Python

 Ziel dieses Kapitels ist es, die Update-Gleichungen und die Randbedingungen in ein Computerprogramm zu überführen, welches objektorientiert im Rahmen einer _Solver-Klasse_ Methoden zur numerischen Lösung der Maxwell-Gleichungen bereitstellt. 
Python wurde gewählt, da die Sprache eine vergleichsweise kompakte Implementierung mathematischer Operationen ermöglicht und mit NumPy eine Bibliothek zur effizienten Verarbeitung mehrdimensionaler Arrays bereitstellt.
== Aufbau der Klasse <sec-aufbau>



Der gesamte Löser wird als Klasse `FDTD2D` aufgefasst. Eine Simulation kann dementsprechend als ein Objekt betrachtet werden, welches sämtliche Feldwerte, Materialkoeffizienten, PML-Profile und Auflösungen hält. Ändern lassen sich diese Zustände durch verschiedene Methoden, die die Klasse ebenfalls bereitstellt. Der entscheidende Vorteil einer solchen Struktur ist die Möglichkeit, mehrere voneinander unabhängige Simulationen als separate Objektinstanzen zu erzeugen, was insbesondere in Kapitel 5 bei der Durchführung mehrerer Testfälle mit unterschiedlichen Parametern von Nutzen sein wird.

Das untenstehende _Klassendiagramm_ beschreibt den Aufbau der Klasse:



#figure(
  placement: auto,
  block(width: 15.0cm, height: 10.3cm)[
    #set text(size: 7.5pt)
    #set par(justify: false, leading: 0.45em)
    #set align(left)

    // ---------- Zustand (rechts) ----------
    #kl-kasten(8.9cm, 0.0cm, 6.1cm, 10.3cm, fuellung: luma(249))
    #kl-inhalt(8.9cm, 0.0cm, 6.1cm, 10.3cm)[
      #text(size: 8.4pt, weight: "bold")[Zustand des Objekts] \
      #v(0.45em)
      `Ez`, `Hx`, `Hy` \
      #text(fill: luma(80))[die drei Feldkomponenten] \
      #v(0.3em)
      `psi_hx`, `psi_hy`, `psi_ez_x`, `psi_ez_y` \
      #text(fill: luma(80))[Gedächtnis der absorbierenden Schicht] \
      #v(0.3em)
      `Ca`, `Cb` \
      #text(fill: luma(80))[Materialkoeffizienten, je Gitterpunkt ein Paar] \
      #v(0.3em)
      `bx_i`, `by_j`, `bx_h`, `by_h` \
      #text(fill: luma(80))[Dämpfungsprofile der Schicht] \
      #v(0.3em)
      `eps_r`, `sigma`, `_pec_mask` \
      #text(fill: luma(80))[Material und Geometrie] \
      #v(0.3em)
      `_soft_sources`, `_probes`, ... \
      #text(fill: luma(80))[angemeldete Quellen und Sonden] \
      #v(0.3em)
      `time_step`, `t` \
      #text(fill: luma(80))[Zähler und Zeit des aktuellen $E_z$] \
      #v(0.7em)
      #text(size: 7pt, style: "italic", fill: luma(60))[
        Bleibt zwischen allen Aufrufen erhalten. Nur die Methoden
        links verändern ihn.]
    ]

    // ---------- Aufbau ----------
    #kl-kasten(0.0cm, 0.0cm, 5.6cm, 2.75cm)
    #kl-inhalt(0.0cm, 0.0cm, 5.6cm, 2.75cm)[
      #text(size: 8.4pt, weight: "bold")[Aufbau] #h(0.3em)
      #text(size: 7pt, fill: luma(90))[einmal vor der Rechnung] \
      #v(0.4em)
      `__init__()` \
      `_update_material_coeffs()` \
      `_setup_pml()` \
      `set_material_box()`, `set_pec_box()` \
      `add_soft_source()`, `add_probe()`, ...
    ]

    // ---------- Zeitschritt ----------
    #kl-kasten(0.0cm, 3.55cm, 5.6cm, 4.55cm, fuellung: luma(252))
    #kl-inhalt(0.0cm, 3.55cm, 5.6cm, 0.6cm)[
      #text(size: 8.4pt, weight: "bold")[Zeitschritt] #h(0.3em)
      #text(size: 7pt, fill: luma(90))[`run(n)` ruft `step()` $n$-mal]
    ]
    #kl-kasten(0.5cm, 4.25cm, 4.6cm, 0.5cm)
    #kl-mitte(0.5cm, 4.25cm, 4.6cm, 0.5cm)[`_update_H()`]
    #kl-pfeil-ab(2.8cm, 4.75cm, 5.03cm)
    #kl-kasten(0.5cm, 5.03cm, 4.6cm, 0.5cm)
    #kl-mitte(0.5cm, 5.03cm, 4.6cm, 0.5cm)[`_inject_current_sources()`]
    #kl-pfeil-ab(2.8cm, 5.53cm, 5.81cm)
    #kl-kasten(0.5cm, 5.81cm, 4.6cm, 0.5cm)
    #kl-mitte(0.5cm, 5.81cm, 4.6cm, 0.5cm)[`_update_E()`]
    #kl-pfeil-ab(2.8cm, 6.31cm, 6.59cm)
    #kl-kasten(0.5cm, 6.59cm, 4.6cm, 0.5cm)
    #kl-mitte(0.5cm, 6.59cm, 4.6cm, 0.5cm)[`_inject_soft_sources()`]
    #kl-pfeil-ab(2.8cm, 7.09cm, 7.37cm)
    #kl-kasten(0.5cm, 7.37cm, 4.6cm, 0.5cm)
    #kl-mitte(0.5cm, 7.37cm, 4.6cm, 0.5cm)[`_record_probes()`, `_record_snapshots()`]

    // ---------- Auswertung ----------
    #kl-kasten(0.0cm, 8.90cm, 5.6cm, 1.40cm)
    #kl-inhalt(0.0cm, 8.90cm, 5.6cm, 1.40cm)[
      #text(size: 8.4pt, weight: "bold")[Auswertung] #h(0.3em)
      #text(size: 7pt, fill: luma(90))[nach der Rechnung] \
      #v(0.4em)
      `get_spectrum()`, `total_energy()`, `get_grid()` \
      `reset_fields()`
    ]

    // ---------- Verbindungen ----------
    #kl-pfeil-r(5.6cm, 8.9cm, 1.38cm)
    #kl-label(5.95cm, 1.02cm)[legt an und trägt ein]
    #kl-pfeil-lr(5.6cm, 8.9cm, 5.83cm)
    #kl-label(5.95cm, 5.47cm)[liest und schreibt]
    #kl-pfeil-l(8.9cm, 5.6cm, 9.60cm)
    #kl-label(6.55cm, 9.24cm)[liest]
  ],
  kind: image,
  caption: [Aufbau der Klasse `FDTD2D`. Links die Methoden in drei Gruppen,
    rechts der Zustand, den ein Objekt zwischen den Aufrufen hält. Innerhalb
    des Zeitschritts ist die Reihenfolge der Aufrufe festgelegt; sie folgt aus
    der Zuordnung der Felder zu ganz- und halbzahligen Zeitpunkten.],
) <fig-klasse>


== Die absorbierende Randschicht <sec-randschicht>
Um die in @sec-pml hergeleiteten Randbedingungen in Python zu implementieren, müssen zunächst die Koeffizienten berechnet werden:
```python
    def _pml_koeffizienten(self, positionen, n_zellen):
        if self.npml == 0:
            return np.ones_like(positionen), np.zeros_like(positionen)
        d_links = self.npml - positionen
        d_rechts = positionen - (n_zellen - 1 - self.npml)
        tiefe = np.maximum(np.maximum(d_links, d_rechts), 0.0)
        rho = np.minimum(tiefe / self.npml, 1.0)
        sig = self.sigma_max * rho ** self.pml_m
        b = np.exp(-sig * self.dt / EPS0)
        return b, b - 1.0
```
Die Methode `_pml_koeffizienten` berechnet die beiden Rekursionskoeffizienten nach @pml-bc an
beliebigen Positionen. 
Dafür nimmt die Methode als Argument ein eindimensionales _Array_ mit allen Positionen der Gitterpunkte in einer Richtung (bspw. Liste von ganzzahligen x-Koordinaten für $E_z$ oder Liste von  halben Zahlen für die x-Koordinaten von $H_y$).
`d_links` und `d_rechts` messen anschließend, wie tief welche Position in der PML ist. Dabei entstehen jeweils zwei _Arrays_, sodass jede einzelne Position behandelt wird. Da bei einer Position stets mindestens einer der Werte `d_links` und `d_rechts` negativ ist, wählt die Variable `tiefe` die größere Zahl aus. Sind beide Zahlen negativ, so befindet sich die Position nicht in der PML und wird damit auf null gesetzt. Damit ist `tiefe` ebenfalls ein Array, welches für jede Position die Tiefe in einer PML speichert. Darauf aufbauend wird die Tiefe nun auf ein Intervall $[0;1]$ normiert und daraus eine Funktion für die Dämpfung $sigma$ festgelegt. Da außerhalb der PML $rho$ stets null beträgt, ist an diesen Stellen die Dämpfung ebenfalls null. Nach der Berechnung des Koeffizienten _b_ nach @pml-bc wird dieser gemeinsam mit dem Koeffizienten _C_ jeweils als _Array_ zurückgegeben.


```python
    def _setup_pml(self):
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
```
Aus der Restreflexion `r0` folgt $sigma_max$ aus @pml-sigmamax. Erwähnenswert ist ein Namenskonflikt: Während dort $d$ die Dicke der Schicht darstellt, bezeichnet $d$ im Code die Zellweite --
die Variable `D` hält die beiden auseinander.

Insgesamt werden vier Dämpfungsprofile berechnet, da die Feldkomponenten aufgrund des Yee-Gitters an verschiedenen Positionen sind. 
Aus der ganzzahligen Position von $E_z$ folgt `np.arange(nx)` und aus dem Versatz von $H_x$ und $H_y$ folgt `np.arange(nx-1)+0.5`, also eine Position weniger und eine Verschiebung aller
Positionen um $0.5$.  

Abschließend werden die Hilfsfelder $Psi$ auf null initialisiert. `psi_hx` und `psi_hy` haben dabei die gleichen Dimensionen wie 
die beiden H-Felder, während die beiden Hilfsgrößen des E-Feldes mit  $(n_x - 2) times (n_y - 2)$ kleiner sind, da $E_z$  nur im Inneren aktualisiert wird und die beiden äußeren Spalten und Zeilen 
als PEC-Abschluss null bleiben.



== Die Update-Gleichungen im Code <sec-updateE>
```python
    def _update_H(self):
        d = self.d

        dEz_dy = (self.Ez[:, 1:] - self.Ez[:, :-1]) / d
        self.psi_hx = self.by_h[None, :] * self.psi_hx + self.Cy_h[None, :] * dEz_dy
        self.Hx -= (self.dt / MU0) * (dEz_dy + self.psi_hx)

        dEz_dx = (self.Ez[1:, :] - self.Ez[:-1, :]) / d
        self.psi_hy = self.bx_h[:, None] * self.psi_hy + self.Cx_h[:, None] * dEz_dx
        self.Hy += (self.dt / MU0) * (dEz_dx + self.psi_hy)
```
Der obenstehende Code setzt @pml-psi-hx und @pml-update-hx sowie das analoge Paar für $H_y$ um. Für jede Komponente werden jeweils zwei Zeilen Code benötigt:
Zunächst wird die Hilfsgröße nach @pml-bc fortgeschrieben, woraus danach das Feld nach @pml-update-hx berechnet werden kann.

Auffällig ist die Abwesenheit jeglicher _for-Schleifen_, die in der Programmierung klassischerweise für Array-Iterationen verwendet werden. 
Stattdessen greift der Code auf das Konzept der _Vektorisierung_ zurück, bei der Operationen elementweise auf gesamte Arrays wirken. Verwendet wird dabei die Ausschnitts-Schreibweise
(_slicing_) von NumPy: `a:b` wählt auf einer Achse alle Indizes von `a` (eingeschlossen) bis `b` (ausgeschlossen), wobei das Fehlen einer Zahl immer die jeweilige Grenze impliziert.  `self.Ez[:, 1:]` wählt folglich alle
$x$-Werte und alle $y$-Werte ab dem zweiten Gitterpunkt, `self.Ez[:, :-1]`
alle $x$-Werte und alle $y$-Werte bis zum vorletzten Gitterpunkt; beide
Ausschnitte sind damit um einen Gitterpunkt in $y$-Richtung
zueinander verschoben. `self.Ez[:, 1:] - self.Ez[:, :-1]` liefert dementsprechend in
einem einzigen Aufruf die Differenz zwischen benachbarten Gitterpunkten in $y$-Richtung an
allen Gitterpunkten zugleich --- also die diskrete Ableitung
$partial E_z\/partial y$ als vollständiges Array.
Das Verfahren der _Vektorisierung_ ist aufgrund der Ausführung in kompiliertem, für Array-Operationen optimiertem Code deutlich effizienter als eine elementweise Schleife in reinem Python.

Zu erläutern ist außerdem die unterschiedliche Schreibweise von `self.by_h[None, :]` und `self.bx_h[:, None]`. Da die Profile zunächst eindimensional vorliegen, müssen sie für die Verwendung in den zweidimensionalen Feldern entsprechend erweitert werden. Die eingefügte Achse bestimmt dabei, entlang welcher Raumrichtung das jeweilige Profil aufgespannt wird. `[None, :]` formt das $y$-Profil zu einer Zeile, die über alle $x$-Werte ausgedehnt wird, während `[:, None]` das $x$-Profil zu einer Spalte formt, die sich über alle $y$-Werte erstreckt. Würde ein $sigma_y$ versehentlich in $x$-Richtung aufgespannt, entstünde die entsprechende Schicht an den falschen Rändern.


```python
    def _update_E(self):
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
```


Der Ausschnitt `[1:-1, 1:-1]` lässt die äußerste Zeile und Spalte aus, sodass diese immer null bleiben; das ist die PEC-Bedingung @pec als
Abschluss _hinter_ der Schicht, wie in @sec-randschicht bereits angedeutet. Im Gebietsinneren bleiben zudem die Hilfsgrößen $Psi$ wegen $b = 1$, $C = 0$ dauerhaft null, sodass der Ausdruck dort exakt in @update-ez übergeht – die PML ist damit im Gebietsinneren wirkungslos.


// --- Bausteine für den Ablaufplan in Abschnitt 4.8 ------------------------
// Ergänzt die kl-*-Bausteine aus Abschnitt 4.1 um die DIN-Symbole eines
// Flussdiagramms: Raute (Verzweigung), Parallelogramm (Ein-/Ausgabe) und
// Terminator (Start/Ende). Wie dort absolute Koordinaten in cm; wird ein
// Kasten höher oder breiter, verschieben sich die folgenden Werte mit.
#let ad-raute(x, y, w, h, fuellung: white) = place(top + left, dx: x, dy: y,
  polygon(fill: fuellung, stroke: kl-strich,
    (w / 2, 0cm), (w, h / 2), (w / 2, h), (0cm, h / 2)))

#let ad-parallelo(x, y, w, h, skew: 0.4cm, fuellung: white) = place(top + left, dx: x, dy: y,
  polygon(fill: fuellung, stroke: kl-strich,
    (skew, 0cm), (w, 0cm), (w - skew, h), (0cm, h)))

#let ad-terminator(x, y, w, h, fuellung: luma(249)) = place(top + left, dx: x, dy: y,
  rect(width: w, height: h, stroke: kl-strich, fill: fuellung, radius: h / 2))

== Der Gesamtalgorithmus <sec-gesamtablauf>

Unter Berücksichtigung des zeitlichen Versatzes im Rahmen des Leapfrog-Schemas ergibt sich für `n` simulierte Zeitschritte der in @fig-gesamtablauf dargestellte Gesamtalgorithmus. Die abgebildete Reihenfolge folgt unmittelbar aus `__init__`, `step` und `run`; nach DIN-Konvention steht die Raute für die Abbruchbedingung der Zeitschleife und das Parallelogramm für die Eingabe der Parameter. Die Methoden `_record_probes()` und `_record_snapshots()` dienen hierbei der Sicherung der Simulationsergebnisse, welche dann nach der Simulation für weitere Untersuchungen aufgerufen werden können.

#figure(
  block(width: 12.8cm, height: 12.6cm)[
    #set text(size: 7.5pt)
    #set par(justify: false, leading: 0.42em)
    #set align(left)

    // ---------- Aufbau: Konstruktor (senkrecht, Mitte) ----------
    #ad-terminator(5.9cm, 0.0cm, 2.2cm, 0.6cm)
    #kl-mitte(5.9cm, 0.0cm, 2.2cm, 0.6cm)[Start]
    #kl-pfeil-ab(7.0cm, 0.6cm, 0.9cm)

    #ad-parallelo(3.0cm, 0.9cm, 8.0cm, 1.0cm)
    #kl-mitte(3.7cm, 0.9cm, 6.6cm, 1.0cm)[Parameter einlesen: `nx`, `ny`, `d`, `courant`, `eps_r`, `sigma`, `npml`, `pml_m`, `r0`]
    #kl-pfeil-ab(7.0cm, 1.9cm, 2.2cm)

    #kl-kasten(3.0cm, 2.2cm, 8.0cm, 0.6cm)
    #kl-mitte(3.0cm, 2.2cm, 8.0cm, 0.6cm)[Felder auf null setzen: `Ez`, `Hx`, `Hy`, `Jz`, `_pec_mask`]
    #kl-pfeil-ab(7.0cm, 2.8cm, 3.1cm)

    #kl-kasten(3.0cm, 3.1cm, 8.0cm, 0.6cm)
    #kl-mitte(3.0cm, 3.1cm, 8.0cm, 0.6cm)[Materialkoeffizienten berechnen: `Ca`, `Cb` ]
    #kl-pfeil-ab(7.0cm, 3.7cm, 4.0cm)

    #kl-kasten(3.0cm, 4.0cm, 8.0cm, 0.6cm)
    #kl-mitte(3.0cm, 4.0cm, 8.0cm, 0.6cm)[PML-Profile berechnen: `b`- und `C`-Koeffizienten, $Psi$-Felder auf null (@sec-randschicht)]
    #kl-pfeil-ab(7.0cm, 4.6cm, 4.9cm)

    #kl-kasten(3.0cm, 4.9cm, 8.0cm, 0.6cm)
    #kl-mitte(3.0cm, 4.9cm, 8.0cm, 0.6cm)[Quellen festlegen]
    #kl-pfeil-ab(7.0cm, 5.5cm, 5.85cm)

    // ---------- Zeitschleife ----------
    #ad-raute(4.8cm, 5.85cm, 4.4cm, 1.6cm)
    #kl-mitte(4.8cm, 5.85cm, 4.4cm, 1.6cm)[$n < n_"steps"$?]
    #kl-label(7.15cm, 7.48cm)[ja]
    #kl-label(9.35cm, 6.30cm)[nein]

    #kl-pfeil-ab(7.0cm, 7.45cm, 7.8cm)
    #kl-kasten(3.0cm, 7.8cm, 8.0cm, 0.7cm)
    #kl-mitte(3.0cm, 7.8cm, 8.0cm, 0.7cm)[`_update_H()`: $H_x$ und $H_y$ aktualisieren (inkl. $Psi_(H_x)$, $Psi_(H_y)$)]
    #kl-pfeil-ab(7.0cm, 8.5cm, 8.8cm)
    #kl-kasten(3.0cm, 8.8cm, 8.0cm, 0.6cm)
    #kl-mitte(3.0cm, 8.8cm, 8.0cm, 0.6cm)[`_inject_current_sources()`: $J_z$ setzen]
    #kl-pfeil-ab(7.0cm, 9.4cm, 9.7cm)

    #kl-kasten(3.0cm, 9.7cm, 8.0cm, 0.7cm)
    #kl-mitte(3.0cm, 9.7cm, 8.0cm, 0.7cm)[`_update_E()`: $E_z$ aktualisieren (inkl. $Psi_(E\,x)$, $Psi_(E\,y)$, PEC-Reset)]
    #kl-pfeil-ab(7.0cm, 10.4cm, 10.7cm)

    #kl-kasten(3.0cm, 10.7cm, 8.0cm, 0.6cm)
    #kl-mitte(3.0cm, 10.7cm, 8.0cm, 0.6cm)[`_inject_soft_sources()`: additive Quellen]
    #kl-pfeil-ab(7.0cm, 11.3cm, 11.6cm)

    #kl-kasten(3.0cm, 11.6cm, 8.0cm, 0.7cm)
    #kl-mitte(3.0cm, 11.6cm, 8.0cm, 0.7cm)[`_record_probes()`, `_record_snapshots()`; Zähler und Zeit fortschreiben]

    // Rücksprung: unten aus der Schleife heraus, links vorbei, oben wieder in die Raute hinein
    #kl-linie(7.0cm, 12.3cm, 0.6cm, 12.3cm)
    #kl-linie(0.6cm, 12.3cm, 0.6cm, 6.65cm)
    #kl-pfeil-r(0.6cm, 4.8cm, 6.65cm)

    // ---------- Ende ----------
    #kl-pfeil-r(9.2cm, 10.0cm, 6.65cm)
    #ad-terminator(10.0cm, 6.35cm, 2.4cm, 0.6cm)
    #kl-mitte(10.0cm, 6.35cm, 2.4cm, 0.6cm)[Ende]
  ],
  kind: image,
  caption: [Gesamtalgorithmus des Lösers, zusammengesetzt aus Konstruktor
     und Zeitschleife. Die Raute stellt
    die Abbruchbedingung der Zeitschleife dar; bleibt sie erfüllt (ja), wird
    der Zeitschritt erneut durchlaufen, andernfalls (nein) endet die
    Simulation. Die Auswertung  erfolgt danach separat und
    ist hier nicht mit dargestellt.],
) <fig-gesamtablauf>



Damit ist die aus den Maxwell-Gleichungen entwickelte FDTD-Methode vollständig in Python implementiert. Die konkrete Umsetzung bildet die Grundlage für die im folgenden Kapitel durchgeführten Testfälle, anhand derer die numerischen Ergebnisse mit analytisch bekannten Eigenschaften elektromagnetischer Wellen verglichen und damit die Zuverlässigkeit des implementierten Solvers untersucht wird. Der gesamte Code ist in dem im Anhang verlinkten _GitHub_-Archiv aufrufbar.
