// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Validierung an analytisch lösbaren Anordnungen

In Kapitel 4 wurden die zuvor hergeleiteten Update-Gleichungen in ein Python-Programm überführt. Offen bleibt jedoch weiterhin, inwiefern daraus gewonnene 
Lösungen tatsächlich aussagekräftig sind und vor allem, ob sie Ausbreitung, Reflexion und Transmission ebener Wellen so wiedergeben, wie es die analytischen Lösungen vorhersagen. Vor diesem Hintergrund werden im Folgenden verschiedene Testfälle mittels fdtd_core.py simuliert und anschließend mit analytischen Lösungen verglichen. 
Ziel ist es, zu beurteilen, ob die Simulationsergebnisse eine physikalische Aussagekraft besitzen und wovon mögliche Abweichungen abhängen.

== Fragestellung und Vorgehen <sec-vorgehen>

Die numerische Validierung gliedert sich in drei Testfälle, die jeweils eine andere Kenngröße erfassen und dem Kontinuum gegenüberstellen:
#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: (left, left, left),
    table.header([Testfall], [Geprüfter Bestandteil], [Analytische Sollgröße]),
    [1 --- Ausbreitung im Vakuum], [Update-Gleichungen, Zeitschritt],
      [Phasengeschwindigkeit],
    [2 --- Materialgrenzfläche], [Materialkoeffizienten],
      [Fresnel-Koeffizienten],
    [3 --- Absorbierender Rand], [Randschicht],
      [Restreflexion $abs(R(theta))$ (@pml-restreflexion-profil)],
  ),
  caption: [Die drei Testfälle und die jeweils geprüften Bestandteile.],
)

Von einem Laborversuch unterscheidet sich diese Prüfung in einem wesentlichen
Punkt: Es gibt keine Zufallsfehler. Dasselbe Programm liefert bei jedem Lauf
denselben Wert, und die Sollseite ist eine Formel und keine Messung. Eine
Fehlerrechnung im üblichen Sinn entfällt damit. Unsicherheiten gibt es dennoch
--- unter anderem durch Rundungsfehler infolge der begrenzten Genauigkeit von Gleitkommazahlen oder Frequenz- und Phasenfehler bei der diskreten Fourier-Transformation. Die Frequenzauflösung der Fourier-Transformation ist dabei durch die Zahl $N_T$ der aufgezeichneten Zeitschritte auf $Delta f = 1\/(N_T Delta t)$ begrenzt @lit:Schneider2010[S. 121, Gl. (5.31)].

Von besonderer Bedeutung für das Kapitel wird ohnehin nicht der absolute Fehler sein, sondern vielmehr, wie sich dieser verhält --- beispielsweise 
in Abhängigkeit von der Gitterbreite. Dies ermöglicht es, differenzierter die Zuverlässigkeit des Verfahrens zu beurteilen und die enstehenden Abweichungen zu klassifizieren. Abhängig ist das Fehlerverhalten vor allem von zwei
Größen. Nach @lit:Schneider2010[S. 113] bestimmen nicht die absoluten Schrittweiten
$Delta$ und $Delta t$ das Verhalten der Felder im Gitter, sondern allein die
Courant-Zahl $S$, die nur ihr Verhältnis enthält, und die _Auflösung_

$ N_lambda = lambda_0/Delta = c/(f Delta), $ <nlambda>

also die Zahl der Zellen je Wellenlänge. Ein in diesen beiden Größen angegebenes
Ergebnis gilt somit für jede Frequenz und jede Gittergröße gleichermaßen; das
Fehlerverhalten wird deshalb in allen Testfällen über $N_lambda$ und $S$
untersucht.

Um nicht für jedes $N_lambda$ eigene Testläufe mit jeweils einer anderen Frequenz durchführen zu müssen, tritt in den Testfällen an die Stelle der harmonischen Quelle ein Puls. Konkret wird im Folgenden der _Ricker-Puls_ verwendet, also die zweite Ableitung der Gaußkurve @lit:Schneider2010[S. 118, Gl. (5.20)]. Da nach @nlambda jede Frequenz bei gleicher Gitterbreite einer anderen Auflösung entspricht, liefert ein einziger Lauf nach der Fourier-Transformation die Ergebnisse für verschiedene $N_lambda$ zugleich. 


Der Code zu allen Testfällen ist über das im Anhang verlinkte _GitHub_-Archiv einsehbar.

== Testfall 1 --- Ebene Welle im Vakuum <sec-tf1>

Im ersten Testfall steht alleine die Ausbreitung im Vakuum im Zentrum. Dementsprechend werden 
lediglich die Update-Gleichungen @update-hx bis @update-ez mit  $C_a = 1$, $C_b = Delta t\/(epsilon_0 Delta)$ und die Wahl 
des Zeitschrittes geprüft; die PML begrenzt nur das Gebiet. 

=== Aufbau und Messgröße


Nach @sec-quellen wird eine Gitterlinie mit mehreren Stromquellen $J_z$ besetzt (@abb:tf1-aufbau), welche durch gleichphasige Einspeisung eine ebene Welle in $x$-Richtung erzeugen ($theta=0$). Für eine schräge Ausbreitung  um $theta$ wird die jeweilige Einspeisung der Punktquellen um 

$ t_d = (y sin theta)/c $ <tf1-verzoegerung>

verzögert, sodass sich die Elementarwellen zu
einer um $theta$ gekippten Front überlagern. Die Quelllinie durchspannt das gesamte Gebiet; bei schräger Ausbreitung werden die Quellen zu den Enden hin abgeschwächt.

#figure(
  image("../abbildungen/testfall_a_aufbau.pdf", width: 55%),
  caption: [Versuchsaufbau (schematisch). Die um $t_d$ verzögerten
    Stromquellen erzeugen eine unter $theta$ laufende ebene Welle; gemessen
    wird an $P_1$ und $P_2$ im Abstand $L$. (KI-generiert mit Claude Opus 5, Anthropic)],
) <abb:tf1-aufbau>



Parallel zur Ausbreitungsrichtung werden zwei Sonden $P_1$ und $P_2$ mit dem Abstand $L$ registriert; die _Fast Fourier Transformation_ liefert für jegliche Frequenzen jeweils Amplitude und Phase. Da eine Fläche konstanter Phase die Strecke $L$ in der Zeit $L\/tilde(v)_p$ zurücklegt, ergibt sich die Phasendifferenz $Delta phi(f) = 2 pi f L\/tilde(v)_p$, woraus für die Phasengeschwindigkeit folgt:

$ tilde(v)_p (f) = (2 pi f L)/(Delta phi(f)). $ <tf1-vp>

Die Tilde kennzeichnet Größen im Gitter.
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


Aus @abb:tf1-dispersion und @abb:tf1-anisotropie folgen drei systematische Erkenntnisse.
Erstens ist die Phasengeschwindigkeit im Gitter stets langsamer als im Kontinuum, wobei die Abweichung mit abnehmender Auflösung zunimmt.
Zweitens verläuft die Abweichung im doppelt-logarithmischen Diagramm annähernd parallel mit einer Geraden der Steigung $-2$, sodass sie mit $N_lambda^(-2)$ fällt. Drittens zeigt sich eine Anisotropie von $tilde(v)_p$: Längs der Gitterachsen ist die Welle am langsamsten, längs der Diagonalen schneller und erreicht dort für $S = 0,99$ beinahe $c$.
Visuell lässt sich das Phänomen der Frequenzabhängigkeit besonders bei der Betrachtung eines Pulses beobachten (@abb:tf1-puls): 
Bei grober Auflösung zerfällt der Puls in einen langen Nachlauf, bei feiner bleibt er erhalten. 

 Diese Ergebnisse widersprechen zunächst @phasengeschwindigkeit, nach der die Phasengeschwindigkeit im Kontinuum weder von der Frequenz noch von der Richtung abhängt.

#figure(
  image("../abbildungen/testfall_a_pulszerfliessen.pdf", width: 80%),
  caption: [$E_z$ längs der Ausbreitungsrichtung nach $280$ Zellen Laufweg,
    normiert auf die Anfangsamplitude, für grobe (links) und feine (rechts)
    Auflösung.],
) <abb:tf1-puls>

=== Deutung: numerische Dispersion


Ziel der folgenden Herleitung ist eine Formel für die Phasengeschwindigkeit
$tilde(v)_p$ im Gitter, die mit den Messwerten verglichen werden kann. Im Kontinuum folgt die Ausbreitungsgeschwindigkeit aus der Dispersionsrelation @dispersion, die in @sec-dispersionsrelation aus dem Einsetzen einer ebenen Welle in die Wellengleichung folgte. Die gleiche Vorgehensweise lässt sich auch auf das Gitter anwenden @lit:Taflove2005[S. 108]: In die verlustfreien Update-Gleichungen wird die an Gitterpunkten abgetastete ebene Welle 

$ E_z bar.v_(i\,j)^n = E_0 thin
  e^(upright(i)(omega n Delta t - tilde(k)_x i Delta - tilde(k)_y j Delta)) $ <tf1-ansatz>

sowie ein entsprechender Ansatz mit den Amplituden $H_(x 0)$ und $H_(y 0)$ für
$H_x$ und $H_y$ eingesetzt. Nach Vereinfachung bleiben aus den drei Update-Gleichungen drei
algebraische Gleichungen für die Amplituden @lit:Taflove2005[S. 109, mod. Gl. (4.4a--c)]:

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
@lit:Taflove2005[S. 109, Gl. (4.5)]

$ [1/(c Delta t) sin((omega Delta t)/2)]^2
  = [1/Delta sin((tilde(k)_x Delta)/2)]^2
  + [1/Delta sin((tilde(k)_y Delta)/2)]^2 . $ <num-dispersion>


Für $Delta, Delta t -> 0$ geht @num-dispersion mit der Kleinwinkelnäherung $sin x approx x$ in $omega^2\/c^2 = k_x^2+k_y^2$, sprich in die Dispersionsrelation des Kontinuums, über. Bei endlicher Schrittweite dagegen hängt $tilde(k)$ nichtlinear von $omega$ ab:
$tilde(v)_p$ wird frequenzabhängig, und das Gitter verhält sich
wie ein dispersives Medium, entgegen der Vereinfachung, die im Kontinuum festgelegt wurde.
Über
$tilde(k)_x = tilde(k) cos theta$ und $tilde(k)_y = tilde(k) sin theta$ geht
zudem die Richtung ein; das Gitter weist folglich eine Richtungsabhängigkeit auf  @lit:Taflove2005[S. 107].

Mit $S_c = S\/sqrt(2)$ aus @sec-cfl und $N_lambda$ lässt sich @num-dispersion
umformen zu @lit:Taflove2005[S. 109, Gl. (4.6)]

$ 1/S_c^2 sin^2((pi S_c)/N_lambda)
  = sin^2((tilde(k) Delta cos theta)/2) + sin^2((tilde(k) Delta sin theta)/2) . $ <num-dispersion-n>

Diese Gleichung ist nicht geschlossen nach $tilde(k)$ auflösbar und wird daher
numerisch gelöst @lit:Taflove2005[vgl. S. 114]. Die so berechneten Werte stimmen
mit den Messwerten aus @abb:tf1-dispersion und @abb:tf1-anisotropie näherungsweise überein.
Die
Abweichung von $c$ ist damit vollständig durch die numerische Dispersion
erklärt; das Programm setzt das Yee-Verfahren korrekt um.

Woher die drei Befunde stammen, zeigt eine eigene Taylor-Entwicklung von
@num-dispersion-n für $N_lambda >> 1$ mit $sin x approx x - x^3\/6$:

$ tilde(v)_p/c approx 1 - pi^2/(6 N_lambda^2)
  (cos^4 theta + sin^4 theta - S_c^2) . $ <tf1-taylor>

Der Fehler fällt mit $N_lambda^(-2)$ und sinkt mit wachsendem $S_c$; der
Klammerausdruck ist längs der Achsen am größten und längs der Diagonalen am
kleinsten. Für $S = 1$, also $S_c = 1\/sqrt(2)$, verschwindet er dort sogar:
Längs der Gitterdiagonalen ist die Ausbreitung dann dispersionsfrei
@lit:Taflove2005[S. 111]. Das Zerfließen des Pulses in @abb:tf1-puls ist die direkte
Folge: Seine hochfrequenten, schlechter aufgelösten Anteile laufen langsamer und
bleiben als Nachlauf zurück @lit:Taflove2005[S. 107].

Während die Abweichung der Geschwindigkeit im niederschwelligen Prozentbereich als sehr gering erscheint, ergibt sich daraus ein Fehler, der sich mit der Zeit linear aufsummiert. Nach längerer Simulationsdauer kann es somit zu einer völlig verfälschten Phase kommen, weshalb die Kenntnis über die Größe des Fehlers von großer Relevanz für die sachgemäße Anwendung des FDTD-Verfahrens ist @lit:Schneider2010[S. 166].


== Testfall 2 --- Reflexion und Transmission an einer dielektrischen Grenzfläche <sec-tf3>

In diesem Testfall trifft eine Welle auf ein Dielektrikum und wird teilweise reflektiert. Dabei wird insbesondere geprüft, ob die Reflexions- und Transmissionskoeffizienten mit dem Kontinuum übereinstimmen. 

=== Aufbau und Messgröße

Eine ebene Welle trifft aus dem Vakuum ($n_1 = 1$) senkrecht auf einen
verlustfreien Halbraum mit $epsilon_r in {2, 4, 9}$, also der Brechzahl $n_2 = sqrt(epsilon_r)$. Mit $mu = mu_0$ und
@impedanz ergeben sich die beiden Fresnel-Formeln für senkrechten Einfall
@lit:Schneider2010[S. 136, Gl. (5.94), (5.95)], die jeweils die anteilige _Reflexion_ bzw. _Transmission_ angeben:

$ Gamma = (n_1 - n_2)/(n_1 + n_2) ,
  quad
  T = (2 n_1)/(n_1 + n_2)  $ <tf3-fresnel>

In der Versuchsanordnung wird das Dielektrikum mit `set_material_box` ab der Spalte
$i_G$ gesetzt. Die Grenzfläche liegt damit auf dem $H_y$-Knoten bei $(i_G - 1\/2) Delta$. Angeregt wird wie in Testfall 1 über eine Linie von Stromquellen $J_z$ mit $theta = 0$. In beiden Materialien befindet sich jeweils eine Messsonde.
Ergänzend wird eine weitere Simulation durchgeführt, die sich nur darin unterscheidet, dass diesmal kein Dielektrikum vorhanden ist, die Impedanz also durchweg gleich bleibt. 

Nach @lit:Schneider2010[S. 136 f.] wird zweimal gerechnet: Aus Lauf A ohne Dielektrikum ergibt sich die einfallende Welle, während aus Lauf B mit Dielektrikum das Gesamtfeld folgt. 
Da $Gamma$ und $T$ im
Frequenzbereich definiert sind, werden die Spektren
frequenzweise geteilt @lit:Schneider2010[S. 137]. Aus den beiden Testläufen folgt so in Abhängigkeit von der Frequenz für die Fresnel-Koeffizienten @lit:Schneider2010[vgl. S. 144]:



$ Gamma_"Sim"= abs(hat(E)_B - hat(E)_A)/abs(hat(E)_A) bar.v_(P_r) ,
  quad
  T_"Sim" = abs(hat(E)_B)/abs(hat(E)_A) bar.v_(P_t) . $ <tf3-messung>


=== Ergebnis

#figure(
  image("../abbildungen/testfall_d_courant.pdf", width: 100%),
  caption: [Reflexions- und Transmissionsfaktor für $epsilon_r = 4$ bei
    verschiedenen Courant-Zahlen. Links: $abs(Gamma_"Sim")$ und $T_"Sim"$ mit
    den Fresnel-Werten @tf3-fresnel (punktiert). Rechts: relative Abweichung
    von @tf3-fresnel mit einer Vergleichsgeraden der Steigung $-2$.],
) <abb:tf3-courant>

Beide numerischen Koeffizienten liegen über den Fresnel-Werten; $abs(Gamma_"Sim")$ beträgt bei $N_lambda = 10$ etwa $0,42$ statt $1\/3$. Mit wachsender Auflösung konvergieren sie gegen die aus @tf3-fresnel berechneten Werte. Aus @abb:tf3-courant folgt zudem eine quadratische Abnahme dieser Abweichung mit $N_lambda^(-2)$, wohingegen die Courant-Zahl anders als in Testfall 1 kaum Einfluss nimmt. Mit der Permittivität wächst die Abweichung dagegen deutlich (@abb:tf3-material). 
#figure(
  image("../abbildungen/testfall_d_material.pdf", width: 55%),
  caption: [$abs(Gamma_"Sim")$ bezogen auf den Fresnel-Wert für
    verschiedene $epsilon_r$, punktiert $tilde(Gamma)\/Gamma$ nach
    @tf3-gitter; $S = 0,99$.],
) <abb:tf3-material>



=== Deutung


In dem Testfall hängt die ebene Welle mit $theta=0$ nicht von $y$ ab, sodass $H_x$ verschwindet und der Testfall so der Ausbreitung in einem eindimensionalen Gitter ähnelt.

Im Kontinuum braucht man für die beiden Unbekannten $Gamma$ und $T$ zwei
Gleichungen: $E_z$ ist an der Grenzfläche stetig, und $H_y$ ist an der
Grenzfläche stetig. Im Yee-Gitter gibt es an der Grenzfläche jedoch nur einen
$H_y$-Knoten. Die nächsten $E_z$-Knoten liegen eine halbe Zelle links und rechts
davon.

Die Stetigkeit von $H_y$ kann deshalb direkt übernommen werden. Links der
Grenzfläche überlagern sich dabei die einfallende Welle mit der Amplitude $1$ und
die reflektierte Welle mit der Amplitude $tilde(Gamma)$, rechts läuft allein die
transmittierte Welle mit der Amplitude $tilde(T)$. Da die Impedanz im Gitter
exakt ist @lit:Schneider2010[Gl. (7.46)], gilt für jede Welle wie im Kontinuum
$H_y prop n E_z$; bei der reflektierten Welle kehrt sich das Vorzeichen um, da
sie in $-x$-Richtung läuft. Die Stetigkeit von $H_y$ lautet damit
@lit:Schneider2010[Gl. (7.53)]

$ n_1 (1 - tilde(Gamma)) = n_2 tilde(T) . $ <tf3-hstetig>

Für $E_z$ geht das nicht, da es an der Grenzfläche keinen $E_z$-Wert gibt. Stattdessen wird die
Update-Gleichung @update-hy des $H_y$-Knotens auf der Grenzfläche verwendet. Sie
verbindet den $E_z$-Wert links der Grenzfläche mit dem $E_z$-Wert rechts davon
und übernimmt so die Rolle der Stetigkeit von $E_z$ @lit:Schneider2010[S. 169--171].

Der Unterschied zum Kontinuum besteht also darin, dass die Wellen nicht an der
Grenzfläche selbst verglichen werden, sondern jeweils eine halbe Zelle daneben.
Über diese halbe Zelle hat sich die Phase jeder Welle bereits um
$tilde(k)_ell Delta\/2$ weitergedreht, wobei $tilde(k)_ell$ die numerische
Wellenzahl im Medium $ell$ ist ($ell = 1$: Vakuum, $ell = 2$: Dielektrikum).
Die Update-Gleichung hat deshalb dieselbe Form wie die Stetigkeit von $E_z$ im
Kontinuum, $1 + Gamma = T$ (einfallende plus reflektierte Welle ergibt die
transmittierte Welle), nur trägt jede Welle zusätzlich den Faktor ihrer
Phasendrehung @lit:Schneider2010[vgl. Gl. (7.58), (7.64)]:

$ e^(upright(i) tilde(k)_1 Delta\/2) + tilde(Gamma) e^(-upright(i) tilde(k)_1 Delta\/2)
  = tilde(T) e^(upright(i) tilde(k)_2 Delta\/2) . $ <tf3-update>

Löst man die Stetigkeit von $H_y$ (@tf3-hstetig) und die Update-Gleichung
(@tf3-update) gemeinsam nach $tilde(Gamma)$ und $tilde(T)$ auf, bleibt
von dieser Phasendrehung jeweils nur der Kosinus übrig @lit:Schneider2010[S. 172, Gl. (7.71), (7.72)]:

$ tilde(Gamma) = (n_1 cos(tilde(k)_2 Delta\/2) - n_2 cos(tilde(k)_1 Delta\/2))
    /(n_1 cos(tilde(k)_2 Delta\/2) + n_2 cos(tilde(k)_1 Delta\/2)) ,
  quad
  tilde(T) = (2 n_1 cos(tilde(k)_1 Delta\/2))
    /(n_1 cos(tilde(k)_2 Delta\/2) + n_2 cos(tilde(k)_1 Delta\/2)) . $ <tf3-gitter>



Bei kleiner Zellbreite nähern sich die Kosinusfaktoren eins und @tf3-gitter geht in @tf3-fresnel über. 


Die quadratische Ordnung des aus den Kosinusfaktoren folgenden Fehlers lässt sich in eigener Rechnung mittels $cos x approx 1 - x^2 \/2$ und mittels Kleinwinkelnäherung von @num-dispersion-n,  aus welcher $tilde(k)_ell Delta\/2 approx pi n_ell\/N_lambda$ folgt, zeigen: 

$ abs(tilde(Gamma))/abs(Gamma) approx 1 + (n_1 n_2 pi^2)/N_lambda^2 . $ <tf3-naeherung>


Damit sind die Befunde aus @abb:tf3-courant und @abb:tf3-material erklärt, dass die Koeffizienten stets über der Sollgröße waren. Für große $N_lambda$ geht der zweite Summand gegen null und es folgt die analytische Sollgröße. Zudem zeigt sich auch die Zunahme des Fehlers bei großer Permittivität.


== Testfall 3 --- Absorbierender Rand <sec-tf4>

In diesem Testfall steht die PML aus @sec-pml im Vordergrund. Nach
@pml-restreflexion-profil hängt die Reflexion im Kontinuum nicht von der
Frequenz ab; sie beträgt bei senkrechtem Einfall genau die vorgegebene
Restreflexion $abs(R(0))$ und wächst bei schrägem Einfall auf
$abs(R(0))^(cos theta)$. Geprüft wird, ob dies im Gitter gilt und welche Rolle
die Dicke der Schicht spielt.

=== Aufbau und Messgröße

In der Mitte des Gebiets wird nach @sec-quellen
eine Gitterlinie mit Stromquellen $J_z$ besetzt, die wie in Testfall 1
gleichphasig einen Ricker-Puls einspeisen, sodass eine ebene Welle senkrecht
auf die rechte Schicht zuläuft. Wenige Zellen vor der Schicht liegt die Sonde $P$. Die Reflexion wird wie in
Testfall 2 (@sec-tf3) aus zwei Läufen bestimmt; an die Stelle des Laufs ohne
Dielektrikum tritt ein breiteres Gebiet, dessen rechter Rand
während der Messdauer nicht erreicht wird @lit:Taflove2005[S. 314, Abb. 7.2].
$abs(R_"Sim")$ folgt damit aus @tf3-messung. Da die Amplitude einer ebenen
Welle mit dem Laufweg nicht abnimmt, ist keine Korrektur nötig; die Höhe des
Gebiets ist so gewählt, dass Störungen vom oberen und unteren Rand und von den
Enden der Quelllinie die Sonde während der Messdauer nicht erreichen.

Variiert werden die Vorgabe $abs(R(0))$ zwischen $10^(-1)$ und $10^(-20)$ für
$5$, $10$ und $20$ Zellen bei $N_lambda = 10$. Für schrägen Einfall werden die
Quellen nach @tf1-verzoegerung verzögert, zu den Enden hin abgeschwächt und $theta$ bei $14$ Zellen und den
Vorgaben $10^(-3)$ und $10^(-6)$ zwischen $0°$ und $60°$ variiert; die Gebietshöhe wächst dabei mit
$theta$. Eine Kontrolle mit größerer Höhe bei $45°$ und $60°$ ändert die Werte nicht. Das Profil bleibt kubisch ($m = 3$), die Courant-Zahl beträgt
$S = 0,99$.

=== Ergebnis

#figure(
  image("../abbildungen/testfall_c_winkel.pdf", width: 100%),
  caption: [Gemessene Reflexion der PML bei $N_lambda = 10$. Links: über der
    Vorgabe $abs(R(0))$ für drei Dicken bei senkrechtem Einfall, punktiert
    $abs(R_"Sim") = abs(R(0))$, gestrichelt der Standardwert von
    `fdtd_core.py`. Rechts: über dem Einfallswinkel $theta$ für $14$ Zellen
    und zwei Vorgaben, gestrichelt $abs(R(0))^(cos theta)$ nach
    @pml-restreflexion-profil.],
) <abb:tf4-restreflexion>

Für schwache Vorgaben liegt die gemessene Reflexion unterhalb der Vorgabe, und
zwar umso deutlicher, je strenger die Vorgabe und je dünner die Schicht ist
(@abb:tf4-restreflexion, links); bei $10$ Zellen etwa beträgt sie für
$abs(R(0)) = 10^(-1)$ das $0,75$-Fache, für $10^(-3)$ nur noch das $0,11$-Fache
der Vorgabe. Für strengere Vorgaben steigt sie wieder an; eine kleinere Vorgabe führt also nicht
zwangsläufig zu weniger Reflexion: Bei $10$ Zellen liegt der Tiefstwert im Einbruch bei
$10^(-8)$ bei $7,0 dot 10^(-6)$, für $abs(R(0)) = 10^(-20)$ dagegen bei
$5,2 dot 10^(-4)$. Mit
wachsender Dicke sinkt das Minimum deutlich, von
$3,5 dot 10^(-3)$ bei $5$ auf $1,4 dot 10^(-6)$ bei $20$ Zellen.

Mit dem Einfallswinkel steigt die Reflexion für $abs(R(0)) = 10^(-3)$ von $1,8 dot 10^(-4)$ bei $0°$ auf
$1,4 dot 10^(-2)$ bei $60°$ und damit noch stärker als
$abs(R(0))^(cos theta)$, bleibt aber darunter (@abb:tf4-restreflexion, rechts). Für
den Standardwert $10^(-6)$ bleibt sie bis $45°$ unter $10^(-5)$ und steigt erst
bei $60°$ auf $3 dot 10^(-5)$.

=== Deutung

Die gemessene Reflexion setzt sich aus zwei Anteilen zusammen
@lit:Taflove2005[S. 293]. Der _Kontinuumsanteil_ ist die Welle, die die Schicht
zweimal durchläuft und dazwischen am ideal leitenden Rand zurückgeworfen wird;
er entspricht nach @pml-restreflexion-profil der Vorgabe $abs(R(0))$, bei
schrägem Einfall $abs(R(0))^(cos theta)$, und hängt nicht von der Frequenz ab. Der _Diskretisierungsanteil_ entsteht, weil die
Leitfähigkeit im Gitter nicht stetig ansteigt, sondern von Zelle zu Zelle
springt: Jeder Sprung wirft wie eine kleine Grenzfläche einen Teil der Welle
zurück @lit:Schneider2010[S. 309]. Der Potenzverlauf @pml-profil hält diese Sprünge
am Eingang klein und legt die großen in die Tiefe, wo die Felder bereits
abgeklungen sind @lit:Taflove2005[S. 293].

Wie groß der Diskretisierungsanteil ausfällt, hängt davon ab, wie gut das
stufige Profil den stetigen Verlauf des Kontinuums annähert. Mit wachsender
Auflösung $N_lambda$ entfallen mehr Zellen auf eine Wellenlänge, sodass das
Profil für die Welle glatter wird; zusätzliche Läufe für $N_lambda$ zwischen $7$
und $30$ (nicht abgebildet) bestätigen das jedoch nur bei Schichten bis etwa
$6$ Zellen; bei dickeren Schichten hängt die Reflexion nicht eindeutig von
$N_lambda$ ab. Dieselbe Wirkung hat eine dickere Schicht:
Je mehr Zellen sie umfasst, desto flacher steigt die Leitfähigkeit von null auf
$sigma_max$ an, desto kleiner ist jeder einzelne Sprung und desto geringer die
Reflexion (@abb:tf4-restreflexion, links).

Aus beiden Anteilen folgt auch das Minimum in @abb:tf4-restreflexion (links).
Da $sigma_max$ nach @pml-sigmamax mit $-ln abs(R(0))$ wächst, senkt eine
strengere Vorgabe zwar den Kontinuumsanteil, macht das Profil aber steiler und
vergrößert so die Sprünge. Das Optimum gleicht beide Anteile aus
@lit:Taflove2005[S. 293]; weil eine dickere Schicht die Sprünge verkleinert, sinkt
es mit der Dicke und verschiebt sich zu strengeren Vorgaben. Warum die
Reflexion bei schwachen Vorgaben sogar unter der Vorgabe liegt, wurde nicht
abschließend geklärt; eine mögliche Ursache ist die rekursive Auswertung nach
@pml-bc, die die kontinuierliche Dämpfung nur näherungsweise wiedergibt. Der
Anstieg mit $theta$ in @abb:tf4-restreflexion (rechts) folgt aus dem Kontinuumsanteil nach
@pml-restreflexion-profil; bei $10^(-6)$ ist dieser Anteil klein, sodass der
Anstieg erst bei $60°$ sichtbar wird.

Das von @lit:Taflove2005[S. 294] für $10$ Zellen als optimal genannte
$abs(R(0)) approx e^(-16) approx 10^(-7)$ liegt in @abb:tf4-restreflexion
(links) im flachen Bereich um das Minimum. Auch die Größenordnung von rund
$10^(-5)$ stimmt mit @lit:Taflove2005[S. 315, Abb. 7.3a] überein.

Für senkrechten Einfall ist damit zuerst die Dicke der Schicht maßgeblich; `r0` ist
anschließend nahe dem Minimum und nicht möglichst klein zu wählen. Die in den
übrigen Testfällen verwendeten $14$ Zellen erreichen mit dem Standardwert bei
$N_lambda = 10$ bis $60°$ höchstens etwa $3 dot 10^(-5)$; flacherer Einfall wurde nicht gemessen.

Insgesamt zeigen die drei Testfälle Abweichungen vom Kontinuum, die sich bis
auf einzelne Eigenschaften der PML auf die Diskretisierung zurückführen lassen;
ihre Einordnung folgt im Fazit.
