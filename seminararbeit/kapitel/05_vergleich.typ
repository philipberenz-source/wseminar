#import "../template.typ": *

= Validierung an analytisch lösbaren Anordnungen

In Kapitel 4 wurden die zuvor hergeleiteten Update-Gleichungen in ein Python-Programm überführt. Offen bleibt jedoch weiterhin, inwiefern daraus gewonnene 
Lösungen tatsächlich aussagekräftig sind und vor allem, ob sie Ausbreitung, Reflexion und Transmission ebener Wellen so wiedergeben, wie es die analytischen Lösungen vorhersagen. Vor diesem Hintergrund werden im Folgenden verschiedene Testfälle mit `fdtd_core.py` simuliert und anschließend mit analytischen Lösungen verglichen. 
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
in Abhängigkeit von der Gitterbreite. Dies ermöglicht es, differenzierter die Zuverlässigkeit des Verfahrens zu beurteilen und die entstehenden Abweichungen zu klassifizieren. Abhängig ist das Fehlerverhalten vor allem von zwei
Größen. Nach @lit:Schneider2010[S. 113] bestimmen nicht die absoluten Schrittweiten
$Delta$ und $Delta t$ das Verhalten der Felder im Gitter, sondern allein die
Courant-Zahl $S$, die nur ihr Verhältnis enthält, und die _Auflösung_

$ N_lambda = lambda_0/Delta = c/(f Delta), $ <nlambda>

also die Zahl der Zellen je Wellenlänge. Ein in diesen beiden Größen angegebenes
Ergebnis gilt somit für jede Frequenz und jede Gittergröße gleichermaßen; das
Fehlerverhalten wird deshalb in allen Testfällen über $N_lambda$ und $S$
untersucht.

Um nicht für jedes $N_lambda$ eigene Testläufe mit jeweils einer anderen Frequenz durchführen zu müssen, tritt in den Testfällen an die Stelle der harmonischen Quelle ein Puls. Konkret wird im Folgenden der _Ricker-Puls_ verwendet, also die zweite Ableitung der Gaußkurve @lit:Schneider2010[S. 118, Gl. (5.20)]. Da nach @nlambda jede Frequenz bei gleicher Gitterbreite einer anderen Auflösung entspricht, liefert ein einziger Lauf nach der Fourier-Transformation die Ergebnisse für verschiedene $N_lambda$ zugleich. 

Der Code zu allen Testfällen ist über das im Anhang verlinkte GitHub-Archiv einsehbar.

== Testfall 1 --- Ebene Welle im Vakuum <sec-tf1>

Im ersten Testfall steht alleine die Ausbreitung im Vakuum im Zentrum. Dementsprechend werden 
lediglich die Update-Gleichungen @update-hx bis @update-ez mit  $C_a = 1$, $C_b = Delta t\/(epsilon_0 Delta)$ und die Wahl 
des Zeitschrittes geprüft; die PML begrenzt nur das Gebiet. 

=== Aufbau und Messgröße

Nach @sec-quellen wird eine Gitterlinie mit mehreren Stromquellen $J_z$ besetzt (@abb:tf1-aufbau), welche durch gleichphasige Einspeisung eine ebene Welle in $x$-Richtung erzeugen ($theta=0$). Für eine schräge Ausbreitung  um $theta$ wird die jeweilige Einspeisung der Punktquellen um 

$ t_d = (y sin theta)/c $ <tf1-verzoegerung>

verzögert, sodass sich die Elementarwellen zu
einer um $theta$ gekippten Front überlagern @lit:Schneider2010[S. 312, mod. Gl. (11.24)]. Dabei bezeichnet $y$ den Abstand der jeweiligen Punktquelle vom oberen Ende der Quelllinie, welche das gesamte Gebiet durchspannt; bei schräger Ausbreitung werden die Quellen zu den Enden hin abgeschwächt.

#figure(
  image("../abbildungen/testfall_a_aufbau.pdf", width: 55%),
  caption: [Versuchsaufbau (schematisch). Die um $t_d$ verzögerten
    Stromquellen erzeugen eine unter $theta$ laufende ebene Welle; gemessen
    wird an $P_1$ und $P_2$ im Abstand $L$. (KI-generiert mit Claude Opus 5, Anthropic)],
) <abb:tf1-aufbau>

Parallel zur Ausbreitungsrichtung werden zwei Sonden $P_1$ und $P_2$ mit dem Abstand $L$ angeordnet; die _Fast Fourier Transformation_ liefert für jede Frequenz Amplitude und Phase. Da eine Fläche konstanter Phase die Strecke $L$ in der Zeit $L\/tilde(v)_p$ zurücklegt, ergibt sich die Phasendifferenz $Delta phi(f) = 2 pi f L\/tilde(v)_p$, woraus für die Phasengeschwindigkeit folgt:

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
Zweitens verläuft die Abweichung im doppelt-logarithmischen Diagramm annähernd parallel zu einer Geraden der Steigung $-2$, sodass sie mit $N_lambda^(-2)$ fällt. Drittens zeigt sich eine Anisotropie von $tilde(v)_p$: Längs der Gitterachsen ist die Welle am langsamsten, längs der Diagonalen schneller und erreicht dort für $S = 0,99$ beinahe $c$.
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
numerisch gelöst @lit:Taflove2005[S. 114]. Die so berechneten Werte stimmen
mit den Messwerten aus @abb:tf1-dispersion und @abb:tf1-anisotropie näherungsweise überein.
Die
Abweichung von $c$ lässt sich damit auf die numerische Dispersion
zurückführen; ein Fehler in der Umsetzung des Yee-Verfahrens ist nicht erkennbar.

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

Die Abweichung der Geschwindigkeit ist zwar klein, der daraus folgende Phasenfehler wächst jedoch linear mit der Laufstrecke. Bei langen Laufwegen kann die Phase daher erheblich verfälscht sein, weshalb die Kenntnis über die Größe des Fehlers von großer Relevanz für die sachgemäße Anwendung des FDTD-Verfahrens ist.

== Testfall 2 --- Reflexion und Transmission an einer dielektrischen Grenzfläche <sec-tf3>

In diesem Testfall trifft eine Welle auf ein Dielektrikum und wird teilweise reflektiert. Dabei wird insbesondere geprüft, ob die Reflexions- und Transmissionskoeffizienten mit dem Kontinuum übereinstimmen. 

=== Aufbau und Messgröße

Eine ebene Welle trifft aus dem Vakuum ($n_1 = 1$) senkrecht auf einen
verlustfreien Halbraum mit $epsilon_r in {2, 4, 9}$, also der Brechzahl $n_2 = sqrt(epsilon_r)$. Mit $mu = mu_0$ und
@impedanz ergeben sich die beiden Fresnel-Formeln für senkrechten Einfall
@lit:Schneider2010[S. 136, Gl. (5.94), (5.95)], die jeweils die anteilige Reflexion bzw. Transmission angeben:

$ Gamma = (n_1 - n_2)/(n_1 + n_2) ,
  quad
  T = (2 n_1)/(n_1 + n_2)  $ <tf3-fresnel>

In der Versuchsanordnung wird das Dielektrikum mit `set_material_box` ab der Spalte
$i_G$ gesetzt. Die Grenzfläche liegt damit auf dem $H_y$-Knoten bei $(i_G - 1\/2) Delta$. Angeregt wird wie in Testfall 1 über eine Linie von Stromquellen $J_z$ mit $theta = 0$. In beiden Materialien befindet sich jeweils eine Messsonde.
Ergänzend wird eine weitere Simulation durchgeführt, die sich nur darin unterscheidet, dass diesmal kein Dielektrikum vorhanden ist, die Impedanz also durchweg gleich bleibt. 

Nach @lit:Schneider2010[S. 136 f.] wird zweimal gerechnet: Aus Lauf A ohne Dielektrikum ergibt sich die einfallende Welle, während aus Lauf B mit Dielektrikum das Gesamtfeld folgt. 
Da $Gamma$ und $T$ im
Frequenzbereich definiert sind, werden die Spektren
frequenzweise geteilt @lit:Schneider2010[S. 137]. Aus den beiden Testläufen folgt so in Abhängigkeit von der Frequenz für die Fresnel-Koeffizienten @lit:Schneider2010[S. 144]:

$ Gamma_"Sim"= abs(hat(E)_B - hat(E)_A)/abs(hat(E)_A) bar.v_(P_r) ,
  quad
  T_"Sim" = abs(hat(E)_B)/abs(hat(E)_A) bar.v_(P_t) . $ <tf3-messung>

=== Ergebnis

#figure(
  image("../abbildungen/testfall_d_courant.pdf", width: 100%),
  caption: [Reflexions- und Transmissionsfaktor für $epsilon_r = 4$ bei
    verschiedenen Courant-Zahlen. Links: $abs(Gamma_"Sim")$ und $T_"Sim"$ mit
    den Fresnel-Werten @tf3-fresnel (punktiert). Rechts: relative Abweichung
    von @tf3-fresnel mit einer Vergleichsgeraden der Steigung $-2$. ],
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
Gleichungen: Sowohl $E_z$ als auch $H_y$ sind an der Grenzfläche stetig. Im Yee-Gitter gibt es an der Grenzfläche jedoch nur einen
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
Phasendrehung @lit:Schneider2010[Gl. (7.58), (7.64)]:

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

Die PML aus @sec-pml ist über @pml-sigmamax so ausgelegt, dass eine senkrecht
einfallende Welle im Kontinuum genau mit der vorgegebenen Restreflexion
$abs(R(0))$ zurückkehrt; bei schrägem Einfall wächst diese nach
@pml-restreflexion-profil auf $abs(R(0))^(cos theta)$. Geprüft wird, wie gut die
Schicht diesen Auslegungswert im Gitter erreicht. Variiert werden dazu die
Vorgabe $abs(R(0))$, die über @pml-sigmamax die Leitfähigkeit $sigma_max$
festlegt, sowie die Dicke der Schicht, die Auflösung $N_lambda$ und der
Einfallswinkel $theta$.

=== Aufbau und Messgröße

In der Mitte des Gebiets wird nach @sec-quellen eine Gitterlinie mit
Stromquellen $J_z$ besetzt, die wie in Testfall 1 einen Ricker-Puls einspeisen,
sodass eine ebene Welle auf die rechte Schicht zuläuft. Fünf Zellen vor der
Schicht liegt die Sonde $P$. Wie in Testfall 2 (@sec-tf3) wird die Reflexion aus
zwei Läufen bestimmt; an die Stelle des Laufs ohne Dielektrikum tritt ein
breiteres Gebiet, dessen rechter Rand während der Messdauer nicht erreicht wird
@lit:Taflove2005[S. 314, Abb. 7.2]. $abs(R_"Sim")$ folgt damit aus
@tf3-messung, und da der Puls breitbandig ist, liefert jeder Lauf die Reflexion
für alle $N_lambda$ zugleich. Die Höhe des Gebiets ist so gewählt, dass Störungen
vom oberen und unteren Rand und von den Enden der Quelllinie die Sonde während
der Messdauer nicht erreichen. Für schrägen Einfall werden die Quellen nach
@tf1-verzoegerung verzögert und zu den Enden hin abgeschwächt.

Variiert wird die Vorgabe zwischen $10^(-1)$ und $10^(-16)$ für Schichten aus
$5$ und $20$ Zellen, die Auflösung zwischen $N_lambda = 6$ und $30$ beim
Standardwert $abs(R(0)) = 10^(-6)$ von `fdtd_core.py` sowie der Winkel zwischen
$0°$ und $60°$ für $20$ Zellen. Wo nicht variiert, gelten $N_lambda = 20$, das
kubische Profil ($m = 3$) und $S = 0,99$.

=== Ergebnis

#figure(
  image("../abbildungen/testfall_pml.pdf", width: 100%),
  caption: [Gemessene Reflexion der PML. (a) Über der Vorgabe $abs(R(0))$ für
    zwei Dicken, punktiert $abs(R_"Sim") = abs(R(0))$. (b) Über der Auflösung
    $N_lambda$ beim Standardwert $abs(R(0)) = 10^(-6)$ (punktiert). (c) Über dem
    Einfallswinkel $theta$ für $20$ Zellen. ],
) <abb:tf4-restreflexion>

Mit strengerer Vorgabe sinkt die Reflexion zunächst und steigt dann langsam
wieder an, sodass jede Dicke ein Optimum besitzt (@abb:tf4-restreflexion a).
Solange die Vorgabe schwach ist, folgt die Messung ihrem Verlauf, liegt aber
stets darunter; bei $20$ Zellen etwa beträgt sie $3,5 dot 10^(-4)$ statt $10^(-3)$.
Mit wachsender Dicke sinkt das Optimum deutlich, von rund $10^(-3)$ bei $5$ auf rund $10^(-6)$ bei $20$ Zellen.

Über der Auflösung (b) fällt die Reflexion der dünnen Schicht zwischen
$N_lambda = 6$ und $20$ im Mittel etwa mit $N_lambda^(-2)$, danach flacher; die der
$20$-Zellen-Schicht hängt deutlich schwächer von $N_lambda$ ab und fällt nur von
$4,4 dot 10^(-6)$ auf $1,4 dot 10^(-6)$.

Mit dem Einfallswinkel (c) steigt die Reflexion für beide Vorgaben an; für
$10^(-6)$ bleibt sie bis
$45°$ nahe $10^(-6)$, obwohl nach @pml-restreflexion-profil dort bereits
$5,7 dot 10^(-5)$ zu erwarten wären, und steigt erst bei $60°$ deutlich an.
Für $10^(-3)$ liegt sie bei allen Winkeln unter $abs(R(0))^(cos theta)$, bei
$60°$ etwa mit $1,9 dot 10^(-2)$ statt $3,2 dot 10^(-2)$.

=== Deutung

Die gemessene Reflexion setzt sich aus zwei Anteilen zusammen
@lit:Taflove2005[S. 293]. Der _Kontinuumsanteil_ ist die Welle, die die Schicht
zweimal durchläuft und am ideal leitenden Rand zurückgeworfen wird; er
entspricht im Kontinuum @pml-restreflexion-profil. Der _Diskretisierungsanteil_
entsteht, weil die Leitfähigkeit im Gitter nicht stetig ansteigt, sondern von
Zelle zu Zelle springt. Im Kontinuum ist die Schicht an jeder Stelle angepasst, gleich wie
groß $sigma$ ist (@sec-pml). Im Gitter dagegen verwenden zwei benachbarte Zellen
unterschiedliche Koeffizienten $b_w$ und $C_w$ aus @pml-bc, und die Felder auf
beiden Seiten eines Sprungs werden wie an der Grenzfläche in Testfall 2
(@sec-tf3) nicht am selben Ort, sondern eine halbe Zelle versetzt miteinander
verknüpft. Ein solcher Sprung wirkt ähnlich wie eine kleine Grenzfläche und wirft
einen Teil der Welle zurück @lit:Schneider2010[S. 217]; wie an der
Grenzfläche in Testfall 2 wächst dieser Anteil mit der Größe des Sprungs. Ein allmählicher Anstieg verringert diese
Reflexion @lit:Schneider2010[S. 312], weshalb der Potenzverlauf @pml-profil die
Sprünge am Eingang der Schicht klein hält @lit:Taflove2005[S. 293].

Damit lassen sich alle vier Abhängigkeiten qualitativ zuordnen. Eine strengere Vorgabe
senkt den Kontinuumsanteil, erhöht aber $sigma_max$: Nach @pml-sigmamax ist
$sigma_max$ proportional zu $-ln abs(R(0))$, sodass etwa eine Vorgabe von
$10^(-16)$ statt $10^(-4)$ $sigma_max$ und damit jeden Sprung des Profils
vervierfacht. Das Optimum in (a) ist daher die Abwägung zwischen beiden Anteilen
@lit:Taflove2005[S. 293]. Da $sigma_max$ bei gleicher Vorgabe umgekehrt
proportional zur Dicke ist, verteilt eine dickere Schicht den Anstieg auf kleinere
Sprünge und senkt so das Optimum.

Die Auflösung wirkt auf den Diskretisierungsanteil der Reflexion, der durch
die gestufte Abbildung des Leitfähigkeitsprofils entsteht. Mit wachsendem
$N_lambda$ wird dieser Anteil kleiner, da das gestufte Profil gemessen an der
Wellenlänge feiner aufgelöst ist. Bei der dünnen Schicht mit den größten
Sprüngen sinkt die Reflexion in (b) daher deutlich mit wachsendem $N_lambda$.

Der Einfallswinkel (c) wirkt sich auf den Kontinuumsanteil aus. Nach
@pml-daempfung wirkt die Dämpfung nur über $k_x = (omega\/c) cos theta$, sodass
eine schräg einfallende Welle in der Schicht weniger abklingt. Da $cos theta$ im
Exponenten von @pml-restreflexion-profil steht, wächst dieser Anteil mit
abnehmendem $cos theta$ exponentiell und damit über Größenordnungen.

Offen bleibt, weshalb die Messung in (a) und (c) teils deutlich unter dem
Kontinuumswert nach @pml-restreflexion-profil liegt, bei $20$ Zellen und
$abs(R(0)) = 10^(-3)$ etwa um den Faktor $3$. Eine mögliche Erklärung, die
eine eigene, nicht durch Literatur belegte Überlegung ist, liefert die Näherung
in @pml-zerlegt, die Ableitung über jedes Zeitintervall als konstant
anzunehmen: Sie könnte bewirken, dass die Schicht im Gitter stärker dämpft als
ausgelegt, umso mehr, je größer der Zeitschritt ist. Dazu passt, dass ein
Kontrolllauf mit kleinerem Zeitschritt ($S = 0,5$) näher an der Vorgabe liegt.
Belegt ist die Erklärung damit jedoch nicht.

In der Praxis ist damit zwischen Dicke, Vorgabe und Einfallswinkel abzuwägen.
Eine dickere Schicht senkt die Reflexion, vergrößert aber das Rechengebiet; die
Vorgabe ist nahe dem Optimum der gewählten Dicke und nicht möglichst klein zu
wählen; und je schräger Wellen auf den Rand treffen können, desto mehr Reserve
brauchen Dicke und Vorgabe.

Insgesamt zeigen die drei Testfälle Abweichungen vom Kontinuum, die sich weitgehend aus der Diskretisierung ergeben;
ihre Einordnung folgt im Fazit.
