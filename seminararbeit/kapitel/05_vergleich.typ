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
--- unter anderem durch Rundungsfehler bei der Umwandlung aus dem Binärsystem oder Frequenz- und Phasenfehler bei der diskreten Fourier-Transformation. Diese Unsicherheiten
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

Um nicht für jedes $N_lambda$ eigene Testläufe mit jeweils einer anderen Frequenz durchführen zu müssen, tritt in den Testfällen an die Stelle der harmonischen Quelle ein Puls. Konkret wird im Folgenden der _Ricker-Puls_ verwendet, also die zweite Ableitung der Gaußkurve @Schneider2010[S.118, eq.(5.20)]. Ähnlich wie der Gaußpuls enthält dieser ein breites Frequenzspektrum, besitzt jedoch keinen Gleichanteil, welcher im Gitter teils unphysikalische Artefakte erzeugen kann. Da nach @nlambda jede Frequenz bei gleicher Gitterbreite einer anderen Auflösung entspricht, liefert ein einziger Lauf nach der Fourier-Transformation aus @sec-auswertung die Ergebnisse für verschiedene $N_lambda$ zugleich. 


Der Code ist in dem Folgenden Kapitel nicht einsehbar, jedoch befindet sich ein Link zu dem _Github_ Archiv mit allen Dateien im Anhang. 
Jegliche Diagramme wurden mithilfe der Bibliothek `matplotlib` erstellt.

== Testfall 1 --- Ebene Welle im Vakuum <sec-tf1>

Im ersten Testfall steht alleine die Aubreitung im Vakuum im Zentrum. Dementsprechend werden 
lediglich die Update-Gleichungen @update-hx bis @update-ez mit  $C_a = 1$, $C_b = Delta t\/(epsilon_0 Delta)$ und die Wahl 
des Zeitschrittes geprüft; die PML begrenzt lediglich das Gebiet. 

=== Aufbau und Messgröße

Eine Gitterlinie wird nach @sec-quellen mit Stromquellen $J_z$ besetzt
(@abb:tf1-aufbau). Speisen alle Quellen gleichphasig ein, entsteht eine ebene
Welle in $x$-Richtung ($theta = 0$). Für eine schräge Ausbreitung unter dem
Winkel $theta$ wird jede Quelle um

Nach @sec-quellen wird eine Gitterlinie mit mehreren Stromquellen $J_z$ besetzt @abb:tf1-aufbau, welche durch gleichphasige Einspeisung eine ebene Welle in $x$-Richtung erzeugen ($theta=0$). Für eine schräge Ausbreitung  um $theta$ wird die jeweilige Einspeisung der Punktquellen um 

$ t_d = (y sin theta)/c $ <tf1-verzoegerung>

verzögert, sodass sich die Elementarwellen zu
einer um $theta$ gekippten Front überlagern.

#figure(
  image("../abbildungen/testfall_a_aufbau.pdf", width: 55%),
  caption: [Versuchsaufbau (schematisch). Die um $t_d$ verzögerten
    Stromquellen erzeugen eine unter $theta$ laufende ebene Welle; gemessen
    wird an $P_1$ und $P_2$ im Abstand $L$. (KI-generiert durch Opus 5 / Anthropic)],
) <abb:tf1-aufbau>



Parallel zur Ausbreitungsrichtung werden zwei Sonden $P_1$ und $P_2$ mit dem Abstand L registriert ; die FFT liefert für jegliche Frequenzen jeweils Frequenz und Phase. Da eine Fläche konstanter Phase die Strecke $L$ in der Zeit $L\/tilde(v)_p$ zurücklegt ergibt sich die Phasendifferenz $Delta phi(f) = 2 pi f L\/tilde(v)_p$, woraus für die Phasengeschwindigkeit folgt:

$ tilde(v)_p (f) = (2 pi f L)/(Delta phi(f)). $ <tf1-vp>

Die Tilde kennzeichnet Größen im Gitter @Taflove2005[S. 109].
Sollgröße ist nach @phasengeschwindigkeit
die Vakuumlichtgeschwindigkeit $c$ aus @sec-wellengleichung, unabhängig von
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


Aus @abb:tf1-dispersion und @abb:tf1-anisotropie folgen drei systematische Erkentnisse.
So ist die Phasengeschwindigkeit im Gitter stets langsamer als im Kontinuum, wobei die Abweichung mit abnehmender Auflösung zunimmt;
im doppelt-logarithmischen Diagramm verläuft sie annähernd parallel mit einer Geraden der Steigung $-2$, sodass die Abweichung mit $N_lambda^(-2)$ fällt. Außerdem  zeigt sich eine Anisotropie von $tilde(v)_p$: Längs der Gitterachsen ist die Welle am langsamsten, längs der Diagonalen schneller; Für $S=0.99$ beinahe $c$.
Visuell lässt sich das Phänomen der Frequenzunabhängigkeit besonders bei der Betrachtung eines Pulses beobachten (@abb:tf1-dispersion): 
Bei grober Auflösung zerfällt der Puls in einen langen Nachlauf, bei feiner bleibt er erhalten. 

 Diese Ergebnisse widersprechen zunächst den in @sec-material getroffenen Vereinfachungen der Richtungs- und Frequenzunabhängigkeit. Zurückzuführen sind sie auf die numerische Dispersion @Taflove2005[S.107].

#figure(
  image("../abbildungen/testfall_a_pulszerfliessen.pdf", width: 80%),
  caption: [$E_z$ längs der Ausbreitungsrichtung nach $280$ Zellen Laufweg,
    normiert auf die Anfangsamplitude, für grobe (links) und feine (rechts)
    Auflösung.],
) <abb:tf1-puls>

=== Deutung: numerische Dispersion


Im Kontinuum folgt die Aubreitunsgeschiwndigkeit aus der Dispersionsrelation @dispersion, die sich in @sec-dispersionsrelation, die aus dem Einsetzen einer ebenen Welle in die Wellengleichung folgte. Die gleiche Vorgehensweise läasst sich auch auf das Gitter anwenden @Taflove2005[S.108]: In die verlustfreien Update-Gleichungen wird die an Gitterpunkten abgetastete ebene Welle 

$ E_z bar.v_(i\,j)^n = E_0 thin
  e^(upright(i)(omega n Delta t - tilde(k)_x i Delta - tilde(k)_y j Delta)) $ <tf1-ansatz>

sowie ein entsprechender Ansatz mit den Amplituden $H_(x 0)$ und $H_(y 0)$ für
$H_x$ und $H_y$ eingesetzt. Nach Vereinfachung bleiben aus den drei Update-Gleichungen drei
algebraische Gleichungen für die Amplituden @Taflove2005[S. 109, mod. eq.(4.4a-c)]:

$ H_(x 0) = (Delta t thin E_(z 0))/(mu_0 Delta)
  dot sin(tilde(k)_y Delta\/2)/sin(omega Delta t\/2) , $ <tf1-hx0>

$ H_(y 0) = -(Delta t thin E_(z 0))/(mu_0 Delta)
  dot sin(tilde(k)_x Delta\/2)/sin(omega Delta t\/2) , $ <tf1-hy0>

$ E_(z 0) sin((omega Delta t)/2)
  = (Delta t)/(epsilon_0 Delta) [H_(x 0) sin((tilde(k)_y Delta)/2)
    - H_(y 0) sin((tilde(k)_x Delta)/2)] . $ <tf1-ez0>

 Setzt man @tf1-hx0 und @tf1-hy0 in @tf1-ez0 ein,
kürzt sich $E_(z 0)$ heraus, und mit $c = 1\/sqrt(mu epsilon)$ ergibt sich
für das quadratische Gitter die _numerische Dispersionsrelation_
@Taflove2005[S. 109, Gl. (4.5)]

$ [1/(c Delta t) sin((omega Delta t)/2)]^2
  = [1/Delta sin((tilde(k)_x Delta)/2)]^2
  + [1/Delta sin((tilde(k)_y Delta)/2)]^2 . $ <num-dispersion>


Für $Delta, Delta t -> 0$ folgt aus der Kleinwinkelnäherung $sin x approx x$ und @num-dispersion geht in $omega^2\/c^2 = k_x^2+k_y^2$ sprich in die Dispersionsrelation des Kontinuums über. Bei endlicher Schrittweite dagegen hängt $tilde(k)$ nichtlinear von $omega$ ab:
$tilde(v)_p = omega\/tilde(k)$ wird frequenzabhängig, und das Gitter verhält sich
wie ein dispersives Medium, entgegen der Vereinfachung die im Kontinuum festgelegt wurde.
Über
$tilde(k)_x = tilde(k) cos theta$ und $tilde(k)_y = tilde(k) sin theta$ geht
zudem die Richtung ein; das Gitter weißt folglich eine Richtungsabhängigkeit auf  @Taflove2005[S. 107].

Mit $S_c = S\/sqrt(2)$ aus @sec-cfl und $N_lambda$ lässt sich @num-dispersion
umformen zu @Taflove2005[S. 109, Gl. (4.6)]

$ 1/S_c^2 sin^2((pi S_c)/N_lambda)
  = sin^2((tilde(k) Delta cos theta)/2) + sin^2((tilde(k) Delta sin theta)/2) . $ <num-dispersion-n>

Diese Gleichung ist nicht geschlossen nach $tilde(k)$ auflösbar und wird daher
numerisch gelöst @Taflove2005[vgl. S. 114]. Die so berechneten Werte stimmen
mit den Messwerten aus @abb:tf1-dispersion und @abb:tf1-anisotropie näherungsweise überein.
Die
Abweichung von $c$ ist damit vollständig durch die numerische Dispersion
erklärt; das Programm setzt das Yee-Verfahren korrekt um.

Woher die drei Befunde stammen, zeigt eine Taylor-Entwicklung von
@num-dispersion-n für $N_lambda >> 1$ mit $sin x approx x - x^3\/6$:

$ tilde(v)_p/c approx 1 - pi^2/(6 N_lambda^2)
  (cos^4 theta + sin^4 theta - S_c^2) . $ <tf1-taylor>

Der Fehler fällt mit $N_lambda^(-2)$ und sinkt mit wachsendem $S_c$; der
Klammerausdruck ist längs der Achsen am größten und längs der Diagonalen am
kleinsten. Für $S = 1$, also $S_c = 1\/sqrt(2)$, verschwindet er dort sogar:
Längs der Gitterdiagonalen ist die Ausbreitung dann dispersionsfrei
@Taflove2005[S. 111]. Das Zerfließen des Pulses in @abb:tf1-puls ist die direkte
Folge: Seine hochfrequenten, schlechter aufgelösten Anteile laufen langsamer und
bleiben als Nachlauf zurück @Taflove2005[S. 107].

Die numerische Dispersion ist Folge des Abbruchsfehlers aus @sec-finite-differenzen, der zweiter Ordnung war --
die Auswirkungen dieses Fehlers zeigen sich also direkt in der Physik.
Während die Abweichung der Geschwindigkeit im niederschwelligen Prozentbereich als sehr gering erscheint, ergibt sich daraus ein Fehler der sich mit der Zeit linear aufsummiert. Nach längerer Simulationsdauer kann es somit zu einer völlig verfälschten Phase kommen, weshalb die Kentniss über die Größe des Fehlers von großer Relevanz für die sinngemäße Anwendung des FDTD-Verfahrens ist @Schneider2010[S.166].


== Testfall 2 --- Energieerhaltung und Langzeitstabilität <sec-tf2>
In diesem Testfall werden insbesondere die Langzeitstabilität und die Fähigkeit der Simulation Energie zu konservieren auf die Probe gestellt. Von besonderer Bedeutung wird dabei auch die in Abschnitt @sec-cfl definierte Stabilitätszahl $S$ sein. 

=== Aufbau und Messgröße


Ein quadratisches Gebiet von $101 times 101$ quadratischen Zellen mit $Delta = 1 "mm"$ wird angesetzt und an die Stelle der PML als Randschicht tritt ein ideal leitender Rand (@sec-pec). In der Mitte wird eine Stromquelle $J_Z$ platziert, die das Simulationsgebiet einmalig mit einem Ricker-Puls ansetzt. Danach wird dem System keinerlei Energie mehr zugeführt und das Feld wird an den Simulationsrändern vollständig reflektiert, sodass keine Energie das System verlässt. Gerechnet wird mit
$S = 0,5$, $0,99$ und $1$ über $20 thin 000$ Zellen Laufweg, also rund $200$
Durchquerungen des Hohlraums, sowie zum Vergleich mit $S = 1,01$ und $1,05$
knapp oberhalb der Stabilitätsgrenze @cfl.

Gemessen wird dabei die Energie, die im Gitter gespeichert ist. Für den diskreten Fall ist die Gesamtenergie in einem Gitter für den TM#sub[z]-Fall durch @Schuhmann2001[Gl.(27)]:

$ W^n = Delta^2/2 sum_(i,j) [epsilon_0 (E_z bar.v_(i\,j)^n)^2
  + mu_0 (H_x^(n-1/2) H_x^(n+1/2) + H_y^(n-1/2) H_y^(n+1/2))] . $ <tf2-energie>


Durch den in @sec-yee definierten Zeitversatz zwischen $E$ und $H$ um $1\/2$ wird $H$ zum Zeitpunkt $n$ durch Mittelung der beiden benachbarten Werte berechnet. 

Verglichen wird $W^n$ mit der Energie, die die Quelle bis dahin eingespeist
hat. Nach dem Satz von Poynting ändert sich die Feldenergie nur durch die Arbeit
$-integral bold(E) dot bold(J) thin d A$ der Ströme und durch Energiefluss über
den Rand @GriffithsDE[Gl. (8.9)]. An der ideal leitenden Wand verschwindet
$E_z$ und damit der Energiefluss. Dieselbe Bilanz gilt auch im Gitter
@Schuhmann2001[vgl. S. 307, Gl. (30)]; für die zeitdiskrete Form an der
Quellzelle ergibt sich (eigene Rechnung)

$ W_Q^n = - sum_(k=0)^(n-1) Delta^2 Delta t thin J_z^(k+1/2)
  (E_z^k + E_z^(k+1))/2 . $ <tf2-einspeisung>

Sollgröße ist daher

$ W^n = W_Q^n , $ <tf2-soll>

 also eine konstante Energie, sobald die Quelle abgeschaltet ist.

=== Ergebnis

#figure(
  image("../abbildungen/testfall_b_energie.pdf", width: 100%),
  caption: [Energie im PEC-Hohlraum, bezogen auf die eingespeiste Energie
    $W_Q$ bei $S = 0,99$. Links: Feldenergie $W$ für stabile und instabile
    Courant-Zahlen, Kreise: $W_Q$ bei $S = 0,99$. Rechts: größte bis zum
    Zeitpunkt $t$ aufgetretene relative Abweichung zwischen $W$ und $W_Q$ bei
    $S <= 1$ über die gesamte Rechendauer.],
) <abb:tf2-energie>


Für $S <= 1$ steigt die Energie anfangs aufgrund der Arbeit der Quelle an und bleibt nach dem Abschalten der Quelle weitesgehend konstant. Über die gesamte Simulationsdauer weichen $W$ und $W_Q$ dabei um maximal $10 ^(-15)$ voneinander ab, ohne merklichen zeitlichen Zuwachs der Abweichung. 

Oberhalb von $S=1$ ändert sich diese Stabilität fundamental. Mit dem Fortschreiten der Zeit wächst die Energie exponetiell an -- im logarithmischen Diagramm als Gerade -- und übertrifft die eingespeiste Energie nach nur kurzer Laufdauer um mehrere Größenordnungen. 

=== Deutung

=== Deutung

Dass $W_Q$ bei $S = 0,5$ knapp $1 %$ über dem Wert bei $S = 0,99$ liegt, ist
eine Folge der diskreten Quelle: Der Puls wird nur zu den Zeitpunkten
$n + 1\/2$ abgetastet und wirkt auf das Gitterfeld statt auf das
Kontinuumsfeld. Für die Bilanz @tf2-soll spielt das keine Rolle. Dass die
eingespeiste Energie danach erhalten bleibt, zeigen Schuhmann und Weiland: Unter
der Courant-Bedingung $abs(Delta t thin omega_i) <= 2$ haben alle Eigenwerte
des Zeitschritts den Betrag $1$, jede Gitterschwingung wird also weder gedämpft
noch verstärkt @Schuhmann2001[S. 310, Gl. (44)--(46)]; für das quadratische
Gitter entspricht das $S <= 1$ (eigene Umrechnung). Die verbleibenden
$10^(-15)$ sind Rundungsfehler der rund $16$-stelligen Gleitkommarechnung
(vgl. @sec-vorgehen). Für $S > 1$ wird für die kürzesten Gitterwellen ein
Eigenwert betragsmäßig größer als $1$ (eigene Folgerung aus Gl. (44)); diese
von Rundungsfehlern angeregten Wellen wachsen mit jedem Zeitschritt um einen
festen Faktor, also exponentiell @Taflove2005[S. 107]. Die Grenze @cfl ist damit
scharf, und der Fehler aus Testfall 1 betrifft für $S <= 1$ allein die Phase,
nicht die Energie.