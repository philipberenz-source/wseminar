
Der erste Testfall betrachtet allein die Ausbreitung im Vakuum. Ohne Material
und leitende Körper werden so ausschließlich die Update-Gleichungen @update-hx
bis @update-ez mit $C_a = 1$, $C_b = Delta t\/(epsilon_0 Delta)$ und die Wahl
des Zeitschritts geprüft; die PML begrenzt lediglich das Gebiet.

=== Aufbau und Messgröße

@abb:tf1-aufbau zeigt die beiden verwendeten Anordnungen. In (a) wird in einem
Gebiet von $520 times 400$ Zellen mit $Delta = 1 "mm"$ die Gitterlinie $i = 60$
nach @sec-quellen mit gleichphasigen Stromquellen besetzt, sodass eine ebene
Welle in $x$-Richtung ($theta = 0$) entsteht. Als Zeitverlauf dient ein
_Ricker-Puls_, die zweite Ableitung einer Gaußkurve, mit der Spitzenfrequenz
$f_P = 30 "GHz"$ @Schneider2010[S. 118, Gl. (5.20)]. Er ist frei von einem
Gleichanteil und breitbandig, sodass ein einziger Lauf alle Frequenzen von etwa
$9$ bis $75 "GHz"$ zugleich liefert. Jeder Frequenz entspricht eine
_Auflösung_ $N_lambda = lambda_0\/Delta = c\/(f Delta)$, also die Zahl der
Zellen je Wellenlänge; ausgewertet wird $4 <= N_lambda <= 32$. Um auch schräge Richtungen zu erfassen, wird in (b)
--- wie bei @Taflove2005[S. 112] --- eine einzelne Punktquelle im Zentrum eines
Gebiets von $500 times 500$ Zellen angeregt. Die Sondenpaare liegen auf Strahlen
unter $0 <= theta <= 45 degree$ zur $x$-Achse; gewählt werden nur Richtungen mit
ganzzahligem Schrittverhältnis (etwa $2 : 1$), damit jede Sonde exakt auf einem
$E_z$-Punkt liegt. Wegen der Symmetrie des quadratischen Gitters genügt dieser
Winkelbereich.

#figure(
  image("../abbildungen/tf1_aufbau.pdf", width: 82%),
  caption: [Versuchsaufbau von Testfall 1 (schematisch, nicht maßstäblich).
    (a) Eine mit Stromquellen besetzte Gitterlinie erzeugt eine ebene Welle in
    $x$-Richtung, gemessen wird an $P_1$ und $P_2$ im Abstand $L = 200 Delta$.
    (b) Eine Punktquelle erzeugt eine Zylinderwelle; Sondenpaare auf Strahlen
    unter dem Winkel $theta$ erfassen die Richtungsabhängigkeit.],
) <abb:tf1-aufbau>

An den Sonden $P_1$ und $P_2$ im Abstand $L$ (in (a) $L = 200 Delta$, in (b)
$L approx 120 Delta$) wird der Zeitverlauf von $E_z$ aufgezeichnet; die FFT aus
@sec-auswertung liefert daraus für jede Frequenz die Phase beider Signale. Da
eine Fläche konstanter Phase die Strecke $L$ in der Zeit $L\/tilde(v)_p$
zurücklegt, beträgt die Phasendifferenz $Delta phi(f) = 2 pi f L\/tilde(v)_p$,
und die Phasengeschwindigkeit folgt zu

$ tilde(v)_p (f) = (2 pi f L)/(Delta phi(f)). $ <tf1-vp>

Die Tilde kennzeichnet hier und im Folgenden Größen im Gitter
@Taflove2005[S. 109]. Die nur bis auf Vielfache von $2 pi$
bestimmte Phase wird über die Frequenz stetig fortgesetzt und bei der
niedrigsten Frequenz an $2 pi f L\/c$ angeschlossen. In (b) wird zudem der rein
geometrische, analytisch bekannte Phasenanteil der Zylinderwelle
(Hankel-Funktion $H_0^((2))$) abgezogen.

Sollgröße ist nach @phasengeschwindigkeit die Vakuumlichtgeschwindigkeit
$c = 1\/sqrt(mu_0 epsilon_0)$, unabhängig von Frequenz und Richtung. Variiert
werden $N_lambda$, die Courant-Zahl $S$ (Normierung aus @sec-cfl) und $theta$.

=== Ergebnis

#figure(
  image("../abbildungen/tf1_dispersion.pdf", width: 100%),
  caption: [Numerische Phasengeschwindigkeit $tilde(v)_p$: Simulation (Symbole)
    und numerische Dispersionsrelation @num-dispersion-n (Linien).
    (a) Relativer Fehler über der Auflösung bei $theta = 0$ für drei
    Courant-Zahlen, doppelt-logarithmisch. (b) $tilde(v)_p$ über der
    Ausbreitungsrichtung; $S = 0,71$ entspricht $S_c = 0,5$ wie in
    @Taflove2005[Abb. 4.2].],
) <abb:tf1-dispersion>

@abb:tf1-dispersion zeigt drei systematische Befunde. Erstens ist die Welle im
Gitter stets langsamer als $c$, und der Fehler $1 - tilde(v)_p\/c$ wächst mit
abnehmender Auflösung: Halbiert man $N_lambda$, vervierfacht er sich ungefähr
(bei $S = 0,71$ von $3,1 dot 10^(-3)$ für $N_lambda = 20$ auf
$1,27 dot 10^(-2)$ für $N_lambda = 10$), was in (a) der Steigung $-2$
entspricht. Zweitens sinkt der Fehler bei gleicher Auflösung mit wachsendem $S$.
Drittens hängt $tilde(v)_p$ von der Richtung ab: Längs der Gitterachsen ist die
Welle am langsamsten, längs der Diagonalen am schnellsten (bei $S = 0,71$ und
$N_lambda = 10$ etwa $0,9873 c$ gegenüber $0,9958 c$). Die Abweichung von der
Sollgröße ist damit keine zufällige Schwankung, sondern hängt reproduzierbar von
$N_lambda$, $S$ und $theta$ ab.

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
$(2 upright(i)\/Delta t) sin(omega Delta t\/2)$, im Ort entsprechend. Eliminiert man die Amplituden, ergibt sich für das quadratische
Gitter die _numerische Dispersionsrelation_ @Taflove2005[S. 109, Gl. (4.5)]

$ [1/(c Delta t) sin((omega Delta t)/2)]^2
  = [1/Delta sin((tilde(k)_x Delta)/2)]^2
  + [1/Delta sin((tilde(k)_y Delta)/2)]^2 . $ <num-dispersion>

Für $Delta, Delta t -> 0$ gilt $sin x approx x$, und @num-dispersion geht in
$omega^2\/c^2 = k_x^2 + k_y^2$, also in @dispersion über @Taflove2005[S. 111].
Bei endlicher Schrittweite dagegen hängt $tilde(k)$ nichtlinear von $omega$ ab:
$tilde(v)_p = omega\/tilde(k)$ wird frequenzabhängig, und das Gitter verhält sich
wie ein dispersives Medium, obwohl im Vakuum gerechnet wird; ein Puls zerfließt
@Taflove2005[S. 107]. Über $tilde(k)_x = tilde(k) cos theta$ und
$tilde(k)_y = tilde(k) sin theta$ geht zudem die Richtung ein; das Gitter ist
also auch anisotrop. Taflove beschreibt es deshalb als „numerischen Äther“,
dessen Eigenschaften denen des Vakuums sehr nahe kommen, ohne sie ganz zu
erreichen @Taflove2005[S. 107].

Mit $S_c = S\/sqrt(2)$ aus @sec-cfl und $N_lambda$ lässt sich @num-dispersion
umformen zu @Taflove2005[S. 109, Gl. (4.6)]

$ 1/S_c^2 sin^2((pi S_c)/N_lambda)
  = sin^2((tilde(k) Delta cos theta)/2) + sin^2((tilde(k) Delta sin theta)/2) . $ <num-dispersion-n>

Für allgemeines $theta$ ist diese Gleichung nicht geschlossen nach $tilde(k)$
auflösbar und wird daher numerisch gelöst @Taflove2005[vgl. S. 114]; das
Ergebnis ist in @abb:tf1-dispersion als Linie eingetragen. Die Messwerte liegen
auf diesen Kurven, in (a) höchstens $1,1 dot 10^(-4)$ und in (b) höchstens
$1,5 dot 10^(-4)$ davon entfernt.#footnote[Bei $N_lambda = 5$ (nicht dargestellt)
steigt sie unter Zwischenwinkeln auf bis zu $3 dot 10^(-3)$: Im anisotropen
Gitter läuft die Energie nicht exakt in Richtung von $tilde(bold(k))$, sodass
die Sonden auf dem Strahl unter $theta$ Wellenflächen leicht abweichender
Richtung messen (eigene Überlegung, numerisch nachgeprüft).] Für $S_c = 0,5$ und $N_lambda = 20$ gibt Taflove
$tilde(v)_p = 0,996892 c$ längs der Achsen und $0,998968 c$ längs der
Diagonalen an, also eine Anisotropie von $0,208 %$ @Taflove2005[S. 112]; die
Simulation liefert $0,99688 c$, $0,99898 c$ und $0,210 %$. Die gemessene
Abweichung von $c$ ist damit vollständig durch die numerische Dispersion
erklärt, und das Programm setzt das Yee-Verfahren korrekt um.

Woher die drei Befunde stammen, zeigt eine Taylor-Entwicklung von
@num-dispersion-n für $N_lambda >> 1$ mit $sin x approx x - x^3\/6$ (eigene
Rechnung):

$ tilde(v)_p/c approx 1 - pi^2/(6 N_lambda^2)
  (cos^4 theta + sin^4 theta - S_c^2) . $ <tf1-taylor>

Der Fehler fällt mit $N_lambda^(-2)$ und sinkt mit wachsendem $S_c$; der
Klammerausdruck ist längs der Achsen ($cos^4 theta + sin^4 theta = 1$) am
größten und längs der Diagonalen ($= 1\/2$) am kleinsten. Für $S = 1$, also
$S_c = 1\/sqrt(2)$, verschwindet er dort sogar: Längs der Gitterdiagonalen ist
die Ausbreitung dann dispersionsfrei @Taflove2005[S. 111]. Die Simulation
bestätigt dies mit $tilde(v)_p = 0,99984 c$ bei $S = 0,99$, $N_lambda = 10$ und
$theta = 45 degree$.

Für die folgenden Testfälle folgt daraus: Mit $S = 0,99$ und $N_lambda >= 20$
für die kleinste auftretende Wellenlänge bleibt der Fehler der
Phasengeschwindigkeit unter $0,25 %$. Da er sich längs des Laufwegs zu einem
Phasenfehler aufsummiert, wiegt er umso schwerer, je mehr Wellenlängen das
Gebiet umfasst @Taflove2005[S. 107].
