
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

Die Sollgröße liefert der Satz von Poynting @GriffithsDE[Gl. (8.9)]: Die
Energie in einem Volumen ändert sich nur durch Arbeit an Ladungen und durch
Energiefluss $bold(S) = bold(E) times bold(H)$ über die Oberfläche. Nach dem
Abschalten der Quelle gibt es keine Ströme mehr, und an der ideal leitenden
Wand verschwindet $E_z$ und damit $bold(S) dot bold(n)$. Es muss also

$ W = "const." $ <tf2-soll>

gelten. Variiert werden nach @sec-vorgehen die Courant-Zahl --- $S = 0,5$,
$0,99$ und $1$ innerhalb der Stabilitätsgrenze @cfl, $S = 1,001$, $1,01$ und
$1,05$ knapp darüber --- sowie die Auflösung $N_lambda = 5$ bis $40$ an der
Spitzenfrequenz des Pulses.

=== Ergebnis

#figure(
  image("../abbildungen/testfall_b_energie.pdf", width: 100%),
  caption: [Energie im PEC-Hohlraum, bezogen auf den Wert $W_0$ nach
    Abschalten der Quelle. (a) Verschiedene Courant-Zahlen bei
    $N_lambda = 20$, punktiert der Anstieg nach @tf2-wachstum. (b) Verschiedene
    Auflösungen bei $S = 0,99$ und $S = 1,001$. (c) Größte bis zum Zeitpunkt $t$
    aufgetretene Abweichung von @tf2-soll bei $S <= 1$. In (a) und (b) ist
    $H^2$ zum Zeitpunkt $n + 1\/2$ genommen, da das Produkt in @tf2-energie
    beim exponentiellen Anstieg das Vorzeichen wechselt; auf der
    logarithmischen Skala ist der Unterschied nicht sichtbar.],
) <abb:tf2-energie>

Für $S <= 1$ bleibt die Energie nach dem Abschalten der Quelle konstant, und
zwar unabhängig von der Auflösung (@abb:tf2-energie a, b). Über
$20 thin 000$ Zellen Laufweg --- rund $200$ Durchquerungen des Hohlraums ---
weicht sie um höchstens etwa $10^(-15)$ von ihrem Anfangswert ab
(@abb:tf2-energie c). Schon wenig oberhalb von $S = 1$ steigt die Energie
dagegen exponentiell an (Gerade im logarithmischen Diagramm), umso schneller,
je größer $S$ ist. Die Steigung hängt nicht von $N_lambda$ ab; nur bei grober
Auflösung ($N_lambda = 5$) setzt der Anstieg früher ein.

=== Deutung

Das Yee-Verfahren dämpft oder verstärkt eine Welle für $S <= 1$ also nicht,
bei keiner Auflösung; der Fehler aus Testfall 1 betrifft allein die Phase. Die
verbleibende Abweichung von etwa $10^(-15)$ stammt aus der Rundung: Python
rechnet mit Gleitkommazahlen doppelter Genauigkeit, die nur rund $16$
Dezimalstellen speichern (vgl. @sec-vorgehen). Jede Rechenoperation rundet in
der letzten Stelle, also um etwa $10^(-16)$ relativ. Bei rund $10^4$
Gitterpunkten und einigen $10^4$ Zeitschritten heben sich diese Rundungen
größtenteils auf, sammeln sich aber langsam an; deshalb wächst die Abweichung
in @abb:tf2-energie c stufenweise und ist bei $S = 0,5$, das doppelt so viele
Zeitschritte benötigt, am größten. Ein systematischer Energieverlust ist nicht
erkennbar.

Die Instabilität für $S > 1$ lässt sich aus der numerischen Dispersionsrelation
@num-dispersion ablesen (eigene Überlegung). Ihre rechte Seite wird für
$tilde(k)_x Delta = tilde(k)_y Delta = pi$, also für die kürzeste im Gitter
darstellbare Welle, maximal; dort verlangt die Gleichung
$sin(omega Delta t\/2) = sqrt(2) S_c = S$. Für $S <= 1$ existiert dazu ein
reelles $omega$, jede Gitterwelle schwingt also mit konstanter Amplitude. Für
$S > 1$ ist die Gleichung nur mit komplexem $omega$ erfüllbar, und eine
Lösung wächst mit jedem Zeitschritt um einen festen Faktor. Mit
$omega Delta t\/2 = pi\/2 - upright(i) alpha$ und $cosh alpha = S$ wächst die
Feldamplitude je Zeitschritt um $e^(2 alpha)$, die Energie also um
$e^(4 alpha)$; umgerechnet auf den Laufweg $x = c t$ mit $sqrt(2)\/S$
Zeitschritten je Zelle folgt

$ W prop exp((4 sqrt(2) alpha)/S dot x/Delta), quad cosh alpha = S . $ <tf2-wachstum>

Dieser Anstieg stimmt mit der Steigung in @abb:tf2-energie a bis auf wenige
Prozent überein. Weil er allein von der kürzesten Gitterwelle bestimmt wird,
hängt er nicht von $N_lambda$ ab (@abb:tf2-energie b). Eine physikalische Ursache braucht es dazu
nicht: Schon Rundungsfehler enthalten diese Anteile und werden aufgeschaukelt,
bis die Ergebnisse unbegrenzt anwachsen @Taflove2005[S. 107]. Nur bei
$N_lambda = 5$ enthält bereits der Puls selbst merkliche Anteile nahe der
kürzesten Gitterwelle; der Anstieg startet dann von einem höheren Niveau und
setzt entsprechend früher ein. Die in
@sec-cfl angegebene Grenze @cfl ist damit scharf @Taflove2005[vgl. Kap. 4.7].
