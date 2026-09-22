
== Testfall 2 --- Energieerhaltung und Langzeitstabilität <sec-tf2>

Testfall 1 hat gezeigt, dass das Gitter die _Phase_ einer Welle verfälscht.
Offen bleibt, ob es auch ihre _Energie_ verändert und ob die Rechnung über lange
Zeiten stabil bleibt. Beides wird in einem geschlossenen Hohlraum geprüft, aus
dem keine Energie entweichen kann.

=== Aufbau und Messgröße

Ein quadratisches Gebiet von $101 times 101$ Zellen mit $Delta = 1 "mm"$ wird
ringsum von einem ideal leitenden Rand (@sec-pec) begrenzt; eine PML gibt es
nicht. Im Zentrum regt eine Stromquelle $J_z$ nach @sec-vorgehen einen
Ricker-Puls mit $N_lambda = 20$ an der Spitzenfrequenz an, der nach etwa $40$
Zellen Laufweg abgeklungen ist. Danach wird das Feld an den Wänden fortwährend
reflektiert. Gerechnet wird mit $S = 0,5$, $0,99$ und $1$ über $20 thin 000$
Zellen Laufweg, also rund $200$ Durchquerungen des Hohlraums.

Gemessen wird die im Gitter gespeicherte Energie. Im Kontinuum beträgt ihre
Dichte $u = 1\/2 (epsilon_0 E^2 + B^2\/mu_0)$ @GriffithsDE[Gl. (8.13)]; summiert
man sie über alle Zellen, folgt für den TM#sub[z]-Fall die Energie je
Längeneinheit in $z$-Richtung @Schuhmann2001[vgl. Gl. (27)]

$ W^n = Delta^2/2 sum_(i,j) [epsilon_0 (E_z bar.v_(i\,j)^n)^2
  + mu_0 (H_x^(n-1/2) H_x^(n+1/2) + H_y^(n-1/2) H_y^(n+1/2))] . $ <tf2-energie>

Weil $H$ nach @sec-yee um einen halben Zeitschritt versetzt zu $E_z$ vorliegt,
wird $H^2$ zum Zeitpunkt $n$ durch das Produkt der beiden benachbarten Werte
ersetzt (eigene Wahl); beide Anteile beziehen sich so auf denselben Zeitpunkt.

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

insbesondere also eine konstante Energie, sobald die Quelle abgeschaltet ist.

=== Ergebnis

#figure(
  image("../abbildungen/testfall_b_energie.pdf", width: 100%),
  caption: [Energie im PEC-Hohlraum. Links: Feldenergie $W$ (Linien) und
    eingespeiste Energie $W_Q$ (Kreise) während der Anregung. Rechts: größte
    bis zum Zeitpunkt $t$ aufgetretene relative Abweichung zwischen $W$ und
    $W_Q$ über die gesamte Rechendauer.],
) <abb:tf2-energie>

Während der Anregung folgt die Feldenergie für alle drei Courant-Zahlen der
eingespeisten (@abb:tf2-energie, links); danach bleibt sie konstant. Über die
gesamte Rechendauer weichen $W$ und $W_Q$ um höchstens etwa $10^(-15)$
voneinander ab, ohne dass die Abweichung mit der Zeit merklich anwächst
(@abb:tf2-energie, rechts).

=== Deutung

Die Quelle speist je nach Zeitschritt geringfügig unterschiedlich viel Energie
ein: Bei $S = 0,5$ liegt $W_Q$ knapp $1 %$ über dem Wert bei $S = 0,99$. Im
Kontinuum hinge die eingespeiste Energie nicht vom Zeitschritt ab. Ursache ist,
dass der Puls nur zu den Zeitpunkten $n + 1\/2$ abgetastet wird und die Quelle
auf das Gitterfeld mit seiner numerischen Dispersion statt auf das
Kontinuumsfeld wirkt. Wie groß die Abweichung vom Kontinuumswert ist, wird hier
nicht bestimmt.

Was einmal eingespeist ist, bleibt dagegen erhalten. Schuhmann und Weiland
zeigen, dass das Leapfrog-Verfahren unter der Courant-Bedingung jede
Eigenschwingung des Gitters mit einem Faktor vom Betrag $1$ fortschreibt; die
diskrete Energie wird also weder gedämpft noch verstärkt
@Schuhmann2001[S. 310, Gl. (44)--(47)]. Ihre Courant-Bedingung
$abs(Delta t thin omega_i) <= 2$ entspricht für das quadratische
TM#sub[z]-Gitter gerade $S <= 1$ (eigene Umrechnung). Die verbleibende
Abweichung von etwa $10^(-15)$ stammt aus der Rundung: Python rechnet mit
Gleitkommazahlen doppelter Genauigkeit, also rund $16$ Dezimalstellen
(vgl. @sec-vorgehen). Dass sie über $200$ Durchquerungen nicht anwächst, zeigt
die Langzeitstabilität der Rechnung für $S <= 1$.

Der Fehler aus Testfall 1 betrifft damit allein die Phase einer Welle, nicht
ihre Energie. Die Energieerhaltung gilt allerdings nur im verlustfreien,
geschlossenen Gebiet; Leitfähigkeit und PML entziehen dem Feld gewollt Energie.
