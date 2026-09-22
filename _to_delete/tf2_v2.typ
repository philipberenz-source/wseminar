
== Testfall 2 --- Energieerhaltung und Stabilität <sec-tf2>

Testfall 1 hat gezeigt, dass das Gitter die _Phase_ einer Welle verfälscht.
Offen bleibt, ob es auch ihre _Energie_ verändert und ob die Rechnung über
lange Zeiten stabil bleibt. Beides wird in einem geschlossenen Hohlraum geprüft,
aus dem keine Energie entweichen kann.

=== Aufbau und Messgröße

Ein quadratisches Gebiet von $101 times 101$ Zellen mit $Delta = 1 "mm"$ wird
ringsum von einem ideal leitenden Rand (@sec-pec) begrenzt; eine PML gibt es
nicht. Im Zentrum regt eine Stromquelle $J_z$ nach @sec-vorgehen einen
Ricker-Puls an, der nach etwa $50$ Zellen Laufweg abgeklungen ist. Danach
bleibt das Feld sich selbst überlassen und wird an den Wänden fortwährend
reflektiert.

Gemessen wird die im Hohlraum gespeicherte Energie. Im Kontinuum beträgt ihre
Dichte $u = 1\/2 (epsilon_0 E^2 + B^2\/mu_0)$ @GriffithsDE[Gl. (8.13)]; mit
$B = mu_0 H$ und dem Gitter als Integrationsraster folgt für den
TM#sub[z]-Fall die Energie je Längeneinheit in $z$-Richtung

$ W^n = Delta^2/2 sum_(i,j) [epsilon_0 (E_z bar.v_(i\,j)^n)^2
  + mu_0 (H_x^(n-1/2) H_x^(n+1/2) + H_y^(n-1/2) H_y^(n+1/2))] . $ <tf2-energie>

Weil $H$ nach @sec-yee um einen halben Zeitschritt versetzt zu $E_z$ vorliegt,
wird $H^2$ zum Zeitpunkt $n$ durch das Produkt der beiden benachbarten Werte
ersetzt (eigene Wahl); beide Anteile beziehen sich so auf denselben Zeitpunkt.
Zum Vergleich wird auch die „naive“ Summe aus $E_z^n$ und $H^(n+1/2)$
ausgewertet.

Die Sollgröße liefert der Satz von Poynting @GriffithsDE[Gl. (8.9)]: Die
Energie in einem Volumen ändert sich nur durch Arbeit an Ladungen und durch
Energiefluss $bold(S) = bold(E) times bold(H)$ über die Oberfläche. Nach dem
Abschalten der Quelle gibt es keine Ströme mehr, und an der ideal leitenden
Wand verschwindet $E_z$ und damit $bold(S) dot bold(n)$. Es muss also

$ W = "const." $ <tf2-soll>

gelten. Variiert wird die Courant-Zahl: $S = 0,5$, $0,99$ und $1$ liegen
innerhalb der Stabilitätsgrenze @cfl, $S = 1,001$, $1,01$ und $1,05$ knapp
darüber.

=== Ergebnis

#figure(
  image("../abbildungen/testfall_b_energie.pdf", width: 100%),
  caption: [Energie im PEC-Hohlraum, bezogen auf den Wert $W_0$ nach
    Abschalten der Quelle. Links: Verlauf für stabile und instabile
    Courant-Zahlen (logarithmisch). Rechts: größte bis zum Zeitpunkt $t$
    aufgetretene Abweichung von der Erhaltung bei $S <= 1$, zeitzentriert nach
    @tf2-energie und naiv.],
) <abb:tf2-energie>

Für $S <= 1$ bleibt die Energie nach dem Abschalten der Quelle konstant
(@abb:tf2-energie, links). Rechts zeigt sich, dass @tf2-energie über
$20 thin 000$ Zellen Laufweg --- rund $200$ Durchquerungen des Hohlraums ---
bis auf die Rechengenauigkeit von etwa $10^(-15)$ erhalten bleibt, ohne mit der
Zeit anzuwachsen. Die naive Summe schwankt dagegen um einige Prozent, wächst
aber ebenso wenig. Schon wenig oberhalb von $S = 1$ steigt die Energie dagegen
exponentiell an (Gerade im logarithmischen Diagramm), und zwar umso schneller,
je größer $S$ ist. Bereits bei $S = 1,001$ übersteigt sie ihren Anfangswert
nach wenigen hundert Zellen Laufweg um viele Größenordnungen.

=== Deutung

Die Schwankung der naiven Summe ist kein Energieverlust, sondern ein
Messartefakt: $E_z$ und $H$ gehören zu Zeitpunkten, die $Delta t\/2$
auseinanderliegen, und werden daher nicht zum selben Zeitpunkt summiert. Bezieht
man beide auf denselben Zeitpunkt, bleibt die Energie bis auf Rundungsfehler
erhalten. Das
Yee-Verfahren dämpft oder verstärkt eine Welle für $S <= 1$ also nicht; der
Fehler aus Testfall 1 betrifft allein die Phase.

Die Instabilität für $S > 1$ lässt sich aus der numerischen Dispersionsrelation
@num-dispersion ablesen (eigene Überlegung). Ihre rechte Seite wird für
$tilde(k)_x Delta = tilde(k)_y Delta = pi$, also für die kürzeste im Gitter
darstellbare Welle, maximal; dort verlangt die Gleichung
$sin(omega Delta t\/2) = sqrt(2) S_c = S$. Für $S <= 1$ existiert dazu ein
reelles $omega$, jede Gitterwelle schwingt also mit konstanter Amplitude. Für
$S > 1$ ist die Gleichung nur mit komplexem $omega$ erfüllbar, und eine
Lösung wächst mit jedem Zeitschritt um einen festen Faktor --- exponentiell,
wie in @abb:tf2-energie beobachtet. Eine physikalische Ursache braucht es dazu
nicht: Schon Rundungsfehler enthalten diese Anteile und werden aufgeschaukelt,
bis die Ergebnisse unbegrenzt anwachsen @Taflove2005[S. 107]. Die in
@sec-cfl angegebene Grenze @cfl ist damit scharf @Taflove2005[vgl. Kap. 4.7].
