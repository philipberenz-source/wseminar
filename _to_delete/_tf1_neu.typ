
Der erste Testfall betrachtet allein die Ausbreitung im Vakuum. Ohne Material
und leitende Körper werden so ausschließlich die Update-Gleichungen @update-hx
bis @update-ez mit $C_a = 1$, $C_b = Delta t\/(epsilon_0 Delta)$ und die Wahl
des Zeitschritts geprüft; die PML begrenzt lediglich das Gebiet.

=== Aufbau und Messgröße

Eine Gitterlinie wird nach @sec-quellen mit Stromquellen $J_z$ besetzt
(@abb:tf1-aufbau). Speisen alle Quellen gleichphasig ein, entsteht eine ebene
Welle in $x$-Richtung ($theta = 0$). Für eine schräge Ausbreitung unter dem
Winkel $theta$ wird jede Quelle um

$ t_d = (y sin theta)/c $ <tf1-verzoegerung>

verzögert, sodass sich die Elementarwellen nach dem Huygensschen Prinzip zu
einer um $theta$ gekippten Front überlagern. Die Enden der Quellenlinie werden
mit einer $sin^2$-Gewichtung ausgeblendet, um Beugung zu unterdrücken.

Als Zeitverlauf dient ein _Ricker-Puls_, die zweite Ableitung einer Gaußkurve
@Schneider2010[S. 118, Gl. (5.20)]. Er ist gleichanteilfrei und
breitbandig, sodass ein Lauf viele Frequenzen zugleich liefert. Jeder Frequenz
entspricht eine _Auflösung_ $N_lambda = lambda_0\/Delta = c\/(f Delta)$,
also die Zahl der Zellen je Wellenlänge.

#grid(
  columns: (1fr, 1fr), column-gutter: 1.2em, align: bottom,
  [#figure(
    image("../abbildungen/testfall_a_aufbau.pdf", width: 100%),
    caption: [Versuchsaufbau (schematisch). Die um $t_d$ verzögerten
      Stromquellen erzeugen eine unter $theta$ laufende ebene Welle.],
  ) <abb:tf1-aufbau>],
  [#figure(
    image("../abbildungen/testfall_a_anisotropie.pdf", width: 100%),
    caption: [$tilde(v)_p$ über der Ausbreitungsrichtung bei
      $N_lambda = 10$. Punkte: Simulation, gestrichelt: @num-dispersion-n.],
  ) <abb:tf1-anisotropie>],
)

An zwei Sonden $P_1$ und $P_2$ im Abstand $L$ in Ausbreitungsrichtung wird
$E_z$ aufgezeichnet; die FFT aus @sec-auswertung liefert für jede Frequenz die
Phase beider Signale. Da eine Fläche konstanter Phase die
Strecke $L$ in der Zeit $L\/tilde(v)_p$ zurücklegt, beträgt die Phasendifferenz
$Delta phi(f) = 2 pi f L\/tilde(v)_p$, und die Phasengeschwindigkeit folgt zu

$ tilde(v)_p (f) = (2 pi f L)/(Delta phi(f)). $ <tf1-vp>

Die Tilde kennzeichnet Größen im Gitter @Taflove2005[S. 109]; die nur bis auf
Vielfache von $2 pi$ bestimmte Phase wird über die Frequenz stetig fortgesetzt.
Sollgröße ist nach @phasengeschwindigkeit
die Vakuumlichtgeschwindigkeit $c = 1\/sqrt(mu_0 epsilon_0)$, unabhängig von
Frequenz und Richtung. Variiert werden $N_lambda$, die Courant-Zahl $S$
(Normierung aus @sec-cfl) und $theta$.

=== Ergebnis

#figure(
  image("../abbildungen/testfall_a_dispersion.pdf", width: 92%),
  caption: [Phasengeschwindigkeit bei achsenparalleler Ausbreitung
    ($theta = 0$) über der Auflösung, links linear, rechts als relativer Fehler
    doppelt-logarithmisch. Gestrichelt: numerische Dispersionsrelation
    @num-dispersion-n.],
) <abb:tf1-dispersion>

@abb:tf1-dispersion und @abb:tf1-anisotropie zeigen drei systematische Befunde.
Erstens ist die Welle im Gitter stets langsamer als $c$, und die Abweichung
wächst mit abnehmender Auflösung; im doppelt-logarithmischen Diagramm verläuft
sie parallel zur Geraden der Steigung $-2$, fällt also mit $N_lambda^(-2)$.
Zweitens ist sie bei gleicher Auflösung umso kleiner, je größer $S$ ist.
Drittens hängt $tilde(v)_p$ von der Richtung ab: Längs der Gitterachsen ist die
Welle am langsamsten, längs der Diagonalen am schnellsten; für $S = 0,99$
erreicht sie dort nahezu $c$. Zudem verändert sich die Pulsform (@abb:tf1-puls):
Bei grober Auflösung zerfällt der Puls in einen langen Nachlauf, bei feiner
bleibt er erhalten.

#figure(
  image("../abbildungen/testfall_a_pulszerfliessen.pdf", width: 80%),
  caption: [$E_z$ längs der Ausbreitungsrichtung nach $280$ Zellen Laufweg,
    normiert auf die Anfangsamplitude, für grobe (links) und feine (rechts)
    Auflösung.],
) <abb:tf1-puls>

=== Deutung: numerische Dispersion

Im Kontinuum folgt $c_p = c$ aus der Dispersionsrelation @dispersion, die sich in
@sec-dispersionsrelation durch Einsetzen einer ebenen Welle in die
Wellengleichung ergab. Dasselbe Vorgehen lässt sich auf das Gitter übertragen
@Taflove2005[S. 108 f.]: In die verlustfreien Update-Gleichungen wird die an den
Gitterpunkten abgetastete ebene Welle

$ E_z bar.v_(i\,j)^n = E_0 thin
  e^(upright(i)(omega n Delta t - tilde(k)_x i Delta - tilde(k)_y j Delta)) $ <tf1-ansatz>

sowie ein entsprechender Ansatz für $H_x$ und $H_y$ eingesetzt. Jede zentrale
Differenz liefert dann die Differenz zweier Exponentialfaktoren im Abstand einer
halben Schrittweite, in der Zeit etwa
$e^(upright(i) omega Delta t\/2) - e^(-upright(i) omega Delta t\/2)
= 2 upright(i) sin(omega Delta t\/2)$. An die Stelle des Faktors
$upright(i) omega$ der Zeitableitung tritt also
$(2 upright(i)\/Delta t) sin(omega Delta t\/2)$, im Ort entsprechend.
Eliminiert man die Amplituden, ergibt sich für das quadratische Gitter die
_numerische Dispersionsrelation_ @Taflove2005[S. 109, Gl. (4.5)]

$ [1/(c Delta t) sin((omega Delta t)/2)]^2
  = [1/Delta sin((tilde(k)_x Delta)/2)]^2
  + [1/Delta sin((tilde(k)_y Delta)/2)]^2 . $ <num-dispersion>

Für $Delta, Delta t -> 0$ gilt $sin x approx x$, und @num-dispersion geht in
$omega^2\/c^2 = k_x^2 + k_y^2$, also in @dispersion über @Taflove2005[S. 111].
Bei endlicher Schrittweite dagegen hängt $tilde(k)$ nichtlinear von $omega$ ab:
$tilde(v)_p = omega\/tilde(k)$ wird frequenzabhängig, und das Gitter verhält sich
wie ein dispersives Medium, obwohl im Vakuum gerechnet wird. Über
$tilde(k)_x = tilde(k) cos theta$ und $tilde(k)_y = tilde(k) sin theta$ geht
zudem die Richtung ein; das Gitter ist also auch anisotrop. Taflove spricht
von einem „numerischen Äther“, der dem Vakuum sehr nahe kommt, ohne ihm zu
gleichen @Taflove2005[S. 107].

Mit $S_c = S\/sqrt(2)$ aus @sec-cfl und $N_lambda$ lässt sich @num-dispersion
umformen zu @Taflove2005[S. 109, Gl. (4.6)]

$ 1/S_c^2 sin^2((pi S_c)/N_lambda)
  = sin^2((tilde(k) Delta cos theta)/2) + sin^2((tilde(k) Delta sin theta)/2) . $ <num-dispersion-n>

Diese Gleichung ist nicht geschlossen nach $tilde(k)$ auflösbar und wird daher
numerisch gelöst @Taflove2005[vgl. S. 114]; das Ergebnis ist in
@abb:tf1-dispersion und @abb:tf1-anisotropie gestrichelt eingetragen. Die
Messwerte liegen auf diesen Kurven. Nur nahe $theta = 45 degree$ streuen sie
geringfügig, vermutlich weil die endliche Quellenlinie keine ideal ebene Welle
erzeugt. Der Verlauf entspricht qualitativ @Taflove2005[S. 115, Abb. 4.2]. Die
Abweichung von $c$ ist damit vollständig durch die numerische Dispersion
erklärt; das Programm setzt das Yee-Verfahren korrekt um.

Woher die drei Befunde stammen, zeigt eine Taylor-Entwicklung von
@num-dispersion-n für $N_lambda >> 1$ mit $sin x approx x - x^3\/6$ (eigene
Rechnung):

$ tilde(v)_p/c approx 1 - pi^2/(6 N_lambda^2)
  (cos^4 theta + sin^4 theta - S_c^2) . $ <tf1-taylor>

Der Fehler fällt mit $N_lambda^(-2)$ und sinkt mit wachsendem $S_c$; der
Klammerausdruck ist längs der Achsen am größten und längs der Diagonalen am
kleinsten. Für $S = 1$, also $S_c = 1\/sqrt(2)$, verschwindet er dort sogar:
Längs der Gitterdiagonalen ist die Ausbreitung dann dispersionsfrei
@Taflove2005[S. 111]. Das Zerfließen des Pulses in @abb:tf1-puls ist die direkte
Folge: Seine hochfrequenten, schlechter aufgelösten Anteile laufen langsamer und
bleiben als Nachlauf zurück @Taflove2005[S. 107].

Für die folgenden Testfälle heißt das: $S$ nahe der Stabilitätsgrenze wählen und
die kleinste Wellenlänge fein auflösen --- umso feiner, je mehr Wellenlängen das
Gebiet umfasst, da sich der Fehler längs des Laufwegs aufsummiert
@Taflove2005[S. 107].
