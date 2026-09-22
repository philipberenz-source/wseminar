// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Validierung an analytisch lösbaren Anordnungen

In Kapitel 4 wurden die zuvor hergeleiteten Update-Gleichungen in ein Python Programm überführt. Offen bleibt jedoch weiterhin, inwiefern daraus gewonnene 
Lösungen tatsächlich ausagekräftig sind und vor allem, ob diese die Eigenschaften elektromagnetischer Wellen korrekt reproduzieren - insbesondere verglichen mit 
analytischen Lösungsansätzen. Vor diesem Hintergrund werden im Folgenden verschiedene Testfälle mittels fdtd_core.py simuliert und anschließend mit analytischen Lösungen verglichen. 
Ziel ist es zu beurteiöen, ob das Verfahren genau rechnet und wovon mögliche Abweichungen abhängen.

== Fragestellung und Vorgehen <sec-vorgehen>

Ein Programm kann an mehreren Stellen gleichzeitig fehelerhaft sein; eine Anordnung, die zu viele Bestandteile 
der Simulation (PEC,PML,Dielektrikum etc.) könnte somit im Zweifel keine verwertbare Aussage über den Grund 
des Fehlers liefern. Aufgrund dessen rückt in den verschiedenen Testfällen jeweils eine andere Testkonfiguration in den Vordergrund, 
bei der dann jeweils spezifische Testgrößen erfasst und ausgwertet werden:
#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: (left, left, left),
    table.header([Testfall], [Geprüfter Bestandteil], [Analytische Sollgröße]),
    [1 --- Ausbreitung im Vakuum], [Update-Gleichungen, Zeitschritt],
      [Phasengeschwindigkeit],
    [2 --- PEC-Hohlraumresonator], [Randbedingung, Langzeitverhalten],
      [Eigenfrequenzen, Energieerhaltung],
    [3 --- Materialgrenzfläche], [Materialkoeffizienten],
      [Fresnel-Koeffizienten],
    [4 --- Absorbierender Rand], [Randschicht],
      [Reflexionsfreiheit],
  ),
  caption: [Die vier Testfälle und die jeweils geprüften Bestandteile.],
)

Von einem Laborversuch unterscheidet sich diese Prüfung in einem wesentlichen
Punkt: Es gibt keine Zufallsfehler. Dasselbe Programm liefert bei jedem Lauf
denselben Wert, und die Sollseite ist eine Formel und keine Messung. Eine
Fehlerrechnung im üblichen Sinn entfällt damit. Unsicherheiten gibt es dennoch
--- unter anderem durch Rundungsfehler bei der Umwandlung aus dem Binärsystem. Diese Unsicherheiten
werden jedoch als sehr gering eingeschätzt und werden somit im weiteren Verlauf nicht mehr explizit bei jedem Testfall 
berücksichtigt.

Von besonderer Bedeutung für das Kapitel wird ohnehin nicht der absolute Fehler sein, sondern vielmehr, wie sich dieser verhält -- beispielsweise 
in Abhängigkeit von der Gitterbreite. Dadurch gewinnt das Verfahren an Aussagekraft und durch die Erkentniss über das Fehlerverhalten können  
tatsächlichen Anwendungen auf einen möglichst kleinen Fehler optimiert werden. 

Aussagekräftig ist das Fehlerverhalten vor allem in Abhängigkeit von zwei
Größen. Nach @Schneider2010[S. 113] bestimmen nicht die absoluten Schrittweiten
$Delta$ und $Delta t$ das Verhalten der Felder im Gitter, sondern allein die
Courant-Zahl $S$, die nur ihr Verhältnis enthält, und die _Auflösung_

$ N_lambda = lambda_0/Delta = c/(f Delta), $ <nlambda>

also die Zahl der Zellen je Wellenlänge. Ein in diesen beiden Größen angegebenes
Ergebnis gilt somit für jede Frequenz und jede Gittergröße gleichermaßen; das
Fehlerverhalten wird deshalb in allen Testfällen über $N_lambda$ und $S$
untersucht.

Damit dafür nicht für jedes $N_lambda$ ein eigener Lauf nötig ist, wird in allen
Testfällen mit einem Puls statt mit einer harmonischen Quelle angeregt
@Schneider2010[S. 118]. Verwendet wird der _Ricker-Puls_, die zweite Ableitung
einer Gaußkurve @Schneider2010[S. 118, Gl. (5.20)]. Er enthält ein breites
Frequenzspektrum, im Gegensatz zum Gaußpuls aber keinen Gleichanteil, der im
Gitter unphysikalische Artefakte erzeugen kann @Schneider2010[S. 118]. Da nach
@nlambda jeder Frequenz eine Auflösung entspricht, liefert ein einziger Lauf
nach der Fourier-Transformation aus @sec-auswertung die Ergebnisse für einen
ganzen Bereich von $N_lambda$ zugleich.

Der Code ist in dem Folgenden Kapitel nicht einsehbar, jedoch befindet sich ein Link zu dem _Github_ Archiv mit allen Dateien im Anhang. 
Jegliche Diagramme wurden mithilfe der Bibliothek `matplotlib` erstellt.

== Testfall 1 --- Ebene Welle im Vakuum <sec-tf1>

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

#figure(
  image("../abbildungen/testfall_a_aufbau.pdf", width: 55%),
  caption: [Versuchsaufbau (schematisch). Die um $t_d$ verzögerten
    Stromquellen erzeugen eine unter $theta$ laufende ebene Welle; gemessen
    wird an $P_1$ und $P_2$ im Abstand $L$.],
) <abb:tf1-aufbau>

An zwei Sonden $P_1$ und $P_2$ im Abstand $L$ in Ausbreitungsrichtung wird
$E_z$ aufgezeichnet; die FFT aus @sec-auswertung liefert für jede Frequenz die
Phase beider Signale, also nach @sec-vorgehen für jedes $N_lambda$. Da eine Fläche konstanter Phase die
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
    doppelt-logarithmisch mit Vergleichsgerade der Steigung $-2$.],
) <abb:tf1-dispersion>

#figure(
  image("../abbildungen/testfall_a_anisotropie.pdf", width: 62%),
  caption: [Phasengeschwindigkeit über der Ausbreitungsrichtung $theta$ bei
    $N_lambda = 10$.],
) <abb:tf1-anisotropie>

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
numerisch gelöst @Taflove2005[vgl. S. 114]. Die so berechneten Werte stimmen
mit den Messwerten aus @abb:tf1-dispersion und @abb:tf1-anisotropie überein.
Kleine Abweichungen zeigen sich nur bei großem $N_lambda$, wo die begrenzte
Laufzeit den Nachlauf des Pulses abschneidet, und nahe $theta = 45 degree$,
vermutlich weil die endliche Quellenlinie keine ideal ebene Welle erzeugt. Der Verlauf entspricht qualitativ @Taflove2005[S. 115, Abb. 4.2]. Die
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
