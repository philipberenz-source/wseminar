
== Testfall 2 --- Energieerhaltung und Langzeitstabilität <sec-tf2>

Testfall 1 hat gezeigt, dass das Gitter die _Phase_ einer Welle verfälscht.
Offen bleibt, ob es auch ihre _Energie_ verändert und ob die Rechnung über lange
Zeiten stabil bleibt. Beides wird in einem geschlossenen Hohlraum geprüft, aus
dem keine Energie entweichen kann.

=== Aufbau und Messgröße

Ein quadratisches Gebiet von $101 times 101$ Zellen mit $Delta = 1 "mm"$ wird
ringsum von einem ideal leitenden Rand (@sec-pec) begrenzt; eine PML gibt es
nicht. Im Zentrum regt eine Stromquelle $J_z$ nach @sec-vorgehen einen
einzelnen Ricker-Puls mit $N_lambda = 20$ an der Spitzenfrequenz an. Nach der
doppelten Pulsverzögerung $2 t_0$, entsprechend etwa $50$ Zellen Laufweg, wird
die Quelle abgeschaltet ($J_z = 0$); danach wird keine Energie mehr zugeführt,
und das Feld wird an den Wänden fortwährend reflektiert. Gerechnet wird mit
$S = 0,5$, $0,99$ und $1$ über $20 thin 000$ Zellen Laufweg, also rund $200$
Durchquerungen des Hohlraums, sowie zum Vergleich mit $S = 1,01$ und $1,05$
knapp oberhalb der Stabilitätsgrenze @cfl.

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
  caption: [Energie im PEC-Hohlraum, bezogen auf die eingespeiste Energie
    $W_Q$ bei $S = 0,99$. Links: Feldenergie $W$ für stabile und instabile
    Courant-Zahlen, Kreise: $W_Q$ bei $S = 0,99$. Für $S > 1$ ist $H^2$ zum
    Zeitpunkt $n + 1\/2$ genommen, da das Produkt in @tf2-energie beim
    exponentiellen Anstieg das Vorzeichen wechselt. Rechts: größte bis zum
    Zeitpunkt $t$ aufgetretene relative Abweichung zwischen $W$ und $W_Q$ bei
    $S <= 1$ über die gesamte Rechendauer.],
) <abb:tf2-energie>

Für $S <= 1$ steigt die Feldenergie während des Pulses an und bleibt nach dem
Abschalten der Quelle konstant (@abb:tf2-energie, links). Sie stimmt dabei zu
jedem Zeitpunkt mit der eingespeisten Energie überein: Über die gesamte
Rechendauer weichen $W$ und $W_Q$ um höchstens etwa $10^(-15)$ voneinander ab,
ohne dass die Abweichung mit der Zeit merklich anwächst (@abb:tf2-energie,
rechts).

Oberhalb von $S = 1$ ändert sich das Bild grundlegend. Obwohl nach dem Puls
keine Energie mehr zugeführt wird, wächst $W$ exponentiell an --- im
logarithmischen Diagramm als Gerade --- und übersteigt die eingespeiste
Energie bei $S = 1,01$ nach rund $130$ und bei $S = 1,05$ nach rund $60$ Zellen
Laufweg um mehr als $30$ Größenordnungen. Je weiter $S$ die Grenze
überschreitet, desto steiler ist der Anstieg.

=== Deutung

Die Quelle speist je nach Zeitschritt geringfügig unterschiedlich viel Energie
ein: Bei $S = 0,5$ liegt $W_Q$ knapp $1 %$ über dem Wert bei $S = 0,99$. Im
Kontinuum hinge die eingespeiste Energie nicht vom Zeitschritt ab. Ursache ist,
dass der Puls nur zu den Zeitpunkten $n + 1\/2$ abgetastet wird und die Quelle
auf das Gitterfeld mit seiner numerischen Dispersion statt auf das
Kontinuumsfeld wirkt. Wie groß die Abweichung vom Kontinuumswert ist, wird hier
nicht bestimmt; für die Erhaltung spielt sie keine Rolle, denn geprüft wird
die Bilanz @tf2-soll.

Was einmal eingespeist ist, bleibt für $S <= 1$ erhalten. Schuhmann und Weiland
schreiben das Leapfrog-Verfahren als Iteration $bold(y)^(m+1) = bold(B) bold(y)^m$
des aus $E$ und $H$ zusammengesetzten Feldvektors $bold(y)$ und bestimmen die
Eigenwerte der Iterationsmatrix $bold(B)$ @Schuhmann2001[S. 309 f., Gl. (39)--(44)]:

$ lambda_B = (2 - (Delta t thin omega_i)^2)/2
  plus.minus sqrt(((2 - (Delta t thin omega_i)^2)/2)^2 - 1) , $ <tf2-eigenwerte>

wobei $omega_i$ die Eigenfrequenzen des räumlich diskretisierten Gitters sind.
Solange die Courant-Bedingung $abs(Delta t thin omega_i) <= 2$ erfüllt ist, ist
der Ausdruck unter der Wurzel negativ, beide Eigenwerte sind komplex und haben
den Betrag $abs(lambda_B) = 1$ @Schuhmann2001[S. 310, Gl. (45)--(46)]. Jede
Eigenschwingung des Gitters wird damit pro Zeitschritt nur in ihrer Phase
gedreht, nicht verstärkt oder gedämpft. Für das quadratische TM#sub[z]-Gitter
ist die größte Eigenfrequenz $omega_i = 2 sqrt(2) c\/Delta$ (kürzeste
Gitterwelle, $tilde(k)_x Delta = tilde(k)_y Delta = pi$), sodass die
Bedingung gerade $S <= 1$ entspricht (eigene Umrechnung). Die verbleibende
Abweichung von etwa $10^(-15)$ stammt aus der Rundung: Python rechnet mit
Gleitkommazahlen doppelter Genauigkeit, also rund $16$ Dezimalstellen
(vgl. @sec-vorgehen). Dass sie über $200$ Durchquerungen nicht anwächst, zeigt
die Langzeitstabilität der Rechnung.

Für $S > 1$ ist $abs(Delta t thin omega_i) > 2$ für die kürzesten Gitterwellen.
Dann wird der Ausdruck unter der Wurzel in @tf2-eigenwerte positiv, beide
Eigenwerte sind reell, und einer davon ist betragsmäßig größer als $1$ (eigene
Folgerung aus @tf2-eigenwerte). Die zugehörige Schwingung wächst mit jedem
Zeitschritt um diesen festen Faktor, also exponentiell --- wie in
@abb:tf2-energie beobachtet. Je größer $S$, desto größer der Faktor und desto
steiler der Anstieg. Eine physikalische Ursache braucht es dazu nicht: Der Puls
selbst enthält diese kurzen Wellen kaum, doch schon Rundungsfehler regen sie
an, und sie werden aufgeschaukelt, bis die Ergebnisse unbegrenzt anwachsen
@Taflove2005[S. 107]. Die in @sec-cfl angegebene Grenze ist damit scharf: Schon
$1 %$ darüber ist die Rechnung nach wenigen hundert Zellen Laufweg unbrauchbar.

Der Fehler aus Testfall 1 betrifft damit für $S <= 1$ allein die Phase einer
Welle, nicht ihre Energie. Die Energieerhaltung gilt allerdings nur im
verlustfreien, geschlossenen Gebiet; Leitfähigkeit und PML entziehen dem Feld
gewollt Energie.
