// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Validierung an analytisch lösbaren Anordnungen

Kapitel 4 hat beschrieben, wie der Löser rechnet. Ob er richtig rechnet, ist
damit nicht gezeigt. Dieses Kapitel prüft ihn an vier Anordnungen, deren Lösung
sich in geschlossener Form angeben lässt. Sein Ergebnis sind keine neuen
physikalischen Erkenntnisse --- die Vergleichsgrößen sind seit Langem bekannt
---, sondern Fehlergrenzen, die Angabe ihrer Ursachen und Bedingungen, unter
denen die Rechnung verlässlich bleibt.

== Fragestellung und Vorgehen

Ein Programm kann an mehreren Stellen zugleich fehlerhaft sein; eine Anordnung,
die alle Bestandteile gleichzeitig beansprucht, könnte im Fehlerfall nicht sagen,
welcher versagt hat. Die vier Testfälle rücken deshalb jeweils einen anderen in
den Vordergrund: Testfall 1 ist der einzige, in dem sich das Feld über eine
größere Strecke und in alle Richtungen ausbreitet; Testfall 2 ist der einzige
geschlossene und prüft dadurch die Zeitintegration für sich; Testfall 3 ist der
einzige mit Materialkontrast; Testfall 4 prüft den Rand, der in allen
offenen Anordnungen mitwirkt.

#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: (left, left, left),
    table.header([Testfall], [Geprüfter Bestandteil], [Analytische Sollgröße]),
    [1 --- Ausbreitung im Vakuum], [Update-Gleichungen, Zeitschritt],
      [Phasengeschwindigkeit, Amplitudengesetz, Feldverlauf],
    [2 --- PEC-Hohlraumresonator], [Randbedingung, Langzeitverhalten],
      [Eigenfrequenzen, Energieerhaltung, Stabilitätsgrenze],
    [3 --- Materialgrenzfläche], [Materialkoeffizienten],
      [Fresnel-Koeffizienten],
    [4 --- Absorbierender Rand], [Randschicht],
      [Reflexionsfreiheit],
  ),
  caption: [Die vier Testfälle und die jeweils geprüften Bestandteile.],
)

Von einem Laborversuch unterscheidet sich diese Prüfung in zwei Punkten: Es gibt
keine Zufallsfehler, weil das Verfahren deterministisch ist, und die Sollseite
ist exakt. Die gesamte Abweichung geht damit auf das Verfahren oder auf das
Messverfahren zurück, und beide zu trennen ist die eigentliche Aufgabe. An die
Stelle der Fehlerrechnung treten drei andere Mittel.

*Konvergenzstudie.* Aussagekräftig ist nicht der einzelne Zahlenwert, sondern
_wie_ die Abweichung mit feinerem Gitter abnimmt. Kapitel 3 sagt aus den
zentralen Differenzen einen Fehler zweiter Ordnung voraus, wie ihn auch Taflove
und Hagness @Taflove2005 für das Yee-Verfahren angeben:

$ epsilon prop (Delta x)^p quad "mit" quad p = 2 $ <ordnung>

$p$ wird aus einer Ausgleichsgeraden in doppelt logarithmischer Auftragung
bestimmt; ein Wert nahe zwei bestätigt nicht nur das Ergebnis, sondern die
Herleitung.

*Kontrollrechnung.* Auch die Auswertung hat ihre eigene Ungenauigkeit. Wo das zu
erwarten ist, wird dasselbe Verfahren zusätzlich auf die _exakte_ Lösung
angewandt; was dabei herauskommt, ist der Beitrag der Auswertung allein.

*Nachweisgrenze.* Wo eine Anordnung mit Sollwert null zur Verfügung steht, wird
daran bestimmt, welche Abweichung das Messverfahren überhaupt noch auflösen kann.
Hinzu kommt in Testfall 2 ein Versuch, bei dem ein Versagen erwartet
wird --- eine Prüfung, die nur bestandene Fälle enthält, sagt wenig darüber aus,
ob sie einen Fehler überhaupt anzeigen könnte.

Als Maß für die Auflösung dient die bei Taflove und Hagness @Taflove2005
gebräuchliche Zahl der Gitterzellen je Wellenlänge,

$ N_lambda = lambda / (Delta x) $ <nlambda>

worin $lambda$ die Wellenlänge _im jeweiligen Medium_ ist; Standardwert sind 20
Zellen. Die Frequenz von 1 GHz ist ohne Belang, weil alle Ergebnisse über
$N_lambda$ angegeben und damit übertragbar sind. Zweiter Parameter ist die
Courant-Zahl $S$, hier so normiert, dass die Stabilitätsgrenze bei $S = 1$ liegt
--- genau der Parameter `courant` aus Abschnitt 4.2; sofern nichts anderes
angegeben ist, gilt $S$ = 0,99. Alle Zahlen stammen aus eigenen Läufen von
`fdtd_core.py`; Skripte und Rohdaten liegen im Anhang.

== Testfall 1 --- Ausbreitung im Vakuum

=== Sollgrößen

Im Vakuum breitet sich eine Welle mit der Lichtgeschwindigkeit aus, und zwar in
jede Richtung gleich schnell. Formal folgt das aus der Dispersionsrelation des
Kontinuums, die für eine ebene Welle mit den Komponenten $beta_x$, $beta_y$ der
Phasenkonstante lautet @Griffiths1999

$ mu epsilon omega^2 = beta_x^2 + beta_y^2 $ <dispkont>

Ein Punkt fester Phase erfüllt $omega t - beta x = "const"$; Ableiten nach der
Zeit liefert daraus die Phasengeschwindigkeit. Mit $beta = omega sqrt(mu
epsilon)$ aus @dispkont folgt nach Schneider @Schneider2010

$ c_p = (d x)/(d t) = omega/beta = 1/sqrt(mu epsilon) = c/sqrt(mu_r epsilon_r) $ <cpkont>

Entscheidend daran ist, dass rechts keine Frequenz mehr steht: Im Kontinuum
laufen alle Frequenzen gleich schnell. Für das Vakuum mit $mu_r = epsilon_r = 1$
lautet der Sollwert also $tilde(c)_p \/ c = 1$, und zwar unabhängig davon, wie
fein das Gitter ist und unter welchem Winkel die Welle läuft. Beides wird
geprüft.

Ein zweiter Sollwert betrifft die Amplitude. Eine Welle, die von einer Quelle
ausgeht, wird mit zunehmendem Abstand schwächer, weil sich dieselbe Leistung auf
einen immer größeren Umfang verteilt. In zwei Dimensionen ist dieser Umfang
$2 pi r$, sodass die Amplitude mit $r^(-1\/2)$ abfallen muss; in drei Dimensionen
wäre es eine Kugelfläche und damit das bekannte $1\/r$-Gesetz. Beide
unterscheiden sich so deutlich, dass die Messung zugleich prüft, ob die
Simulation wirklich das zweidimensionale Problem löst.

Die zugehörige exakte Lösung geben Taflove und Hagness @Taflove2005 an. Da eine
Kreiswelle in alle Richtungen gleich aussieht, hängt das Feld nur vom Abstand $r$
zur Quelle ab und nicht von $x$ und $y$ einzeln:

$ E_z (r, t) = "Re"{ A dot H_0^((1))(k r) dot e^(-i omega t) } $ <hankel>

Darin ist $A$ die Stärke der Quelle. Die _Hankel-Funktion_ $H_0^((1))$ übernimmt
dabei genau die Rolle, die im eindimensionalen Fall der Kosinus hat: Sie
beschreibt eine auslaufende Welle, nur dass ihre Amplitude zugleich nach außen
hin abnimmt. Für große Abstände geht sie in $r^(-1\/2) cos(k r - pi\/4)$ über und
liefert damit das oben genannte Gesetz; in Quellnähe weicht sie davon ab, weil
die Wellenfront dort noch stark gekrümmt ist. Zusammengesetzt ist sie aus zwei
Bessel-Funktionen, $H_0^((1)) = J_0 + i Y_0$, die Real- und Imaginärteil und
damit Amplitude und Phase liefern.

Diese beiden Bessel-Funktionen werden ohne Fremdbibliothek berechnet: $J_0$ über
die Trapezregel für $J_0 (z) = 1\/pi integral_0^pi cos(z sin t) dif t$, die bei
periodischem Integranden bis auf Maschinengenauigkeit konvergiert, und $Y_0$ über
die Polynomnäherung von Abramowitz und Stegun @AbramowitzStegun1964. Gegen die
exakte Reihe in Bruchrechnung geprüft liegen ihre Fehler bei $10^(-16)$
beziehungsweise $5 dot 10^(-9)$ --- weit unter allen gemessenen Abweichungen.

=== Aufbau und Messverfahren

Die beiden Sollgrößen verlangen verschiedene Geometrien. Die
Phasengeschwindigkeit wird an einer *ebenen Welle* gemessen: Eine Linienquelle
regt einen Kanal an, in dem das Feld nicht von $y$ abhängt, sodass kein
Krümmungseinfluss entsteht --- dafür ist nur eine Richtung zugänglich. Für
Richtungsabhängigkeit, Amplitudengesetz und Feldvergleich dient eine
*Kreiswelle* aus einer Punktquelle. Das Messverfahren ist in beiden Fällen
dasselbe: Nach dem Einschwingen läuft über acht Perioden eine
Fourier-Transformation bei der Anregungsfrequenz mit.

#figure(
  image("../abbildungen/abb_kreiswelle.pdf", width: 46%),
  caption: [Simuliertes Feld der Kreiswelle bei 20 Zellen je Wellenlänge.
    Eingezeichnet sind die beiden Richtungen, in denen der Fehler extremal wird.],
)

Aus Betrag und Phase folgen alle Messgrößen: die Phasenkonstante als Steigung
der entfalteten Phase über dem Abstand, das Amplitudengesetz aus dem Betrag und
das Differenzfeld aus dem Vergleich mit @hankel. Bei der Kreiswelle liegen die
Abtastpunkte auf Gitterpunkten --- für eine Richtung mit ganzzahligem Verhältnis
$(p, q)$ sind alle Punkte $(c + m p,\, c + m q)$ exakt ---, sodass keine
Interpolation nötig ist; zusätzlich wird über alle acht symmetrisch
gleichwertigen Strahlen gemittelt.

Die Kontrollrechnung ist hier unentbehrlich, denn auf einer gekrümmten
Wellenfront enthält die örtliche Phasensteigung einen Zusatz der Ordnung
$1\/(k r)$. Auf die exakte Lösung @hankel angewandt, ergibt dasselbe Verfahren
statt eins den Wert 0,99986: Die Auswertung allein täuscht eine um 140 ppm zu
kleine Geschwindigkeit vor. Da dieser Beitrag richtungsunabhängig ist, lässt er
sich herausrechnen; alle Werte der Kreiswelle sind entsprechend korrigiert.

=== Ergebnis: der Sollwert wird verfehlt

#figure(
  table(
    columns: (auto, auto, auto),
    align: (center, center, center),
    table.header([$N_lambda$], [$tilde(c)_p \/ c$ gemessen],
                 [Abweichung vom Sollwert]),
    [5], [0,9609396 ± 0,000033], [−3,91 %],
    [10], [0,9913194 ± 0,0000082], [−0,868 %],
    [20], [0,9978899 ± 0,000013], [−0,211 %],
    [40], [0,9994806 ± 0,000012], [−0,0519 %],
    [80], [0,9998755 ± 0,0000084], [−0,0125 %],
  ),
  caption: [Phasengeschwindigkeit der ebenen Welle längs der Gitterachse. Die
    Unsicherheit ist der Standardfehler der Steigung aus der Ausgleichsrechnung.],
) <tab:tf1>

Der Sollwert wird nicht getroffen. Die Welle läuft im Gitter durchgängig zu
langsam, bei fünf Zellen je Wellenlänge um 3,9 %, bei zwanzig noch um 0,21 %.
Das Vorzeichen ist in allen Läufen dasselbe und die Abweichung streut nicht ---
ein Ablesefehler scheidet damit aus.

Hinzu kommt eine zweite Verletzung des Sollwerts: Die gemessene Geschwindigkeit
hängt von der Ausbreitungsrichtung ab. Längs der Gitterachse beträgt der Fehler
bei 20 Zellen je Wellenlänge 0,211 %, in der Diagonalen dagegen nur 0,0041 % ---
ein Verhältnis von 51. Im Vakuum darf es beides nicht geben. Die Ursache muss
also im Verfahren liegen.

=== Erklärung: die Dispersionsrelation des Gitters

@dispkont gilt für das Kontinuum. Die entsprechende Beziehung für das Gitter
gewinnt Schneider @Schneider2010, indem er dieselbe ebene Welle in die diskreten
Update-Gleichungen einsetzt; die folgende Rechnung überträgt sein
eindimensionales Vorgehen auf die zwei Dimensionen dieser Arbeit.

Ausgangspunkt ist, dass eine zentrale Differenz auf eine ebene Welle wie eine
Multiplikation wirkt. Für die Zeitableitung gilt

$ tilde(partial)_t e^(i omega q Delta t) = (e^(i omega Delta t \/ 2) -
  e^(-i omega Delta t \/ 2))/(Delta t) e^(i omega q Delta t)
  = i 2/(Delta t) sin((omega Delta t)/2) e^(i omega q Delta t) $ <schiebung>

und für die beiden Ortsrichtungen dasselbe mit $Delta x$ beziehungsweise
$Delta y$. Mit den Abkürzungen

$ Omega = 2/(Delta t) sin((omega Delta t)/2), quad
  K_x = 2/(Delta x) sin((tilde(beta)_x Delta x)/2), quad
  K_y = 2/(Delta y) sin((tilde(beta)_y Delta y)/2) $ <ersetzung>

bedeutet eine Differenz nach der Zeit demnach eine Multiplikation mit
$i Omega$ und eine Differenz nach dem Ort eine mit $-i K$. Darin sind
$tilde(beta)_x$ und $tilde(beta)_y$ die Komponenten der Phasenkonstante im
Gitter. Jede Richtung bekommt ihren eigenen Ausdruck, und jeder enthält die
Zellweite _dieser_ Richtung --- genau daher rührt die Richtungsabhängigkeit.

Damit werden aus den drei Update-Gleichungen des Yee-Schemas drei algebraische
Beziehungen zwischen den Amplituden $hat(E)_0$, $hat(H)_(x,0)$ und
$hat(H)_(y,0)$. Das Durchflutungsgesetz liefert

$ epsilon Omega hat(E)_0 = K_y hat(H)_(x,0) - K_x hat(H)_(y,0) $ <ampere-disk>

und die beiden Teile des Induktionsgesetzes

$ mu Omega hat(H)_(x,0) = K_y hat(E)_0, quad
  mu Omega hat(H)_(y,0) = -K_x hat(E)_0 $ <faraday-disk>

Setzt man @faraday-disk in @ampere-disk ein, so kürzt sich die Feldamplitude $hat(E)_0$
heraus, und es bleibt die Dispersionsrelation des Gitters:

$ mu epsilon Omega^2 = K_x^2 + K_y^2 $ <disprel>

Der Vergleich mit @dispkont zeigt die Verwandtschaft: Es ist dieselbe Gleichung
mit denselben Rollen, und $mu epsilon$ steht an derselben Stelle --- nur sind die
Größen jetzt die des Gitters. Dass das Kontinuum darin als Grenzfall enthalten
ist, zeigt die Reihe $sin(xi) approx xi$ für kleine Argumente: Sie führt
$Omega$ in $omega$ und $K$ in $beta$ über, womit @disprel in @dispkont übergeht.

Läuft die Welle unter dem Winkel $theta$ zur $x$-Achse,
gilt $tilde(beta)_x = tilde(beta) cos theta$ und $tilde(beta)_y = tilde(beta)
sin theta$; mit gleichen Zellweiten $Delta x = Delta y = Delta$ und der Abkürzung
$u = tilde(beta) Delta \/ 2$ nimmt @disprel die von Schneider @Schneider2010
angegebene Form

$ sin^2( (pi S)/(sqrt(2) N_lambda) )
  = S^2/2 [ sin^2(u cos theta) + sin^2(u sin theta) ] $ <disp2d>

an. Daraus folgt genau die Größe, die gemessen wurde, nämlich das Verhältnis der
Geschwindigkeit im Gitter zur Lichtgeschwindigkeit:

$ tilde(c)_p / c = pi / (N_lambda u) $ <verhaeltnis>

Bemerkenswert ist, dass darin weder Frequenz noch Zellweite einzeln auftreten,
sondern nur ihr Verhältnis $N_lambda$. Für beliebige Winkel wird @disp2d
numerisch aufgelöst; für die beiden Sonderfälle gelingt die Auflösung nach
Schneider @Schneider2010 geschlossen:

$ u_"Achse" = arcsin[ sqrt(2)/S sin( (pi S)/(sqrt(2) N_lambda) ) ], quad
  u_"Diag" = sqrt(2) arcsin[ 1/S sin( (pi S)/(sqrt(2) N_lambda) ) ] $ <disp>

Für $theta = 0$ verschwindet der zweite Summand in @disp2d, für
$theta = 45 degree$ sind beide gleich groß.

Damit lässt sich die Messung prüfen --- nun aber gegen einen anderen Sollwert:
nicht mehr gegen $c$, sondern gegen das, was das Gitter selbst vorhersagt.
@disp liefert für die fünf Auflösungen aus @tab:tf1 die Werte 0,9609386,
0,9913168, 0,9978851, 0,9994746 und 0,9998689 --- die Messung trifft sie auf 1,0
bis 6,6 ppm. Die Abweichungen von bis zu vier Prozent aus @tab:tf1 und diese
wenigen ppm widersprechen einander also nicht, sondern messen Verschiedenes: jene
den Abstand zur Physik, diese den Abstand zur Vorhersage für das Gitter. Die
Restdifferenz wächst dabei nicht mit dem Fehler selbst, sondern bleibt bei allen
Auflösungen bei etwa 5 ppm; sie ist eine Eigenschaft der Auswertung und nicht des
Verfahrens.

Die Abweichung vom Sollwert ist damit erklärt: Sie ist kein Fehler der Umsetzung,
sondern eine vorhersagbare Eigenschaft des Verfahrens, die _numerische
Dispersion_.

#figure(
  image("../abbildungen/abb_dispersion.pdf", width: 88%),
  caption: [Verhältnis von numerischer zu wahrer Phasengeschwindigkeit nach
    @verhaeltnis. Links die Auflösungsabhängigkeit für beide Sonderfälle und
    drei Courant-Zahlen, rechts die Richtungsabhängigkeit bei 20 Zellen je
    Wellenlänge.],
) <abb:disp>

Damit ist auch die Konvergenz belegt. Verdoppelt man die Auflösung, so sinkt der
Fehler auf ein Viertel; gemessen wurden bei den vier Verdopplungsschritten die
Verhältnisse 4,50, 4,11, 4,06 und 4,17. Genau das ist mit _zweiter Ordnung_
gemeint: Der Fehler ist proportional zu $Delta x^2$, und halbiert man $Delta x$,
so viertelt sich dessen Quadrat. Die Ausgleichsgerade liefert entsprechend
$p$ = 2,06 ± 0,02.

Dass es gerade das Quadrat ist, folgt aus @ersetzung. Entwickelt man den Sinus
in seine Reihe $sin(u) = u - u^3\/6 plus.minus dots$, so wird daraus

$ K_x = beta_x (1 - (beta_x Delta x)^2/24 plus.minus dots) $ <reihe>

Das erste Glied ist die gesuchte Größe, das zweite der Fehler --- und dieses
enthält die Zellweite im Quadrat. Ein Glied _erster_ Ordnung fehlt, weil die
zentrale Differenz ihren Punkt symmetrisch einschließt, sodass sich die linearen
Anteile beider Nachbarn gegenseitig aufheben. Ebendieses Ausbleiben des linearen
Glieds macht das Verfahren von zweiter statt erster Ordnung.

Dass die Messung leicht über zwei liegt, ist ebenfalls vorhergesagt: Wertet man
@disp selbst auf denselben fünf Auflösungen aus, ergibt sich $p$ = 2,05. Die
Reihe bricht nach dem quadratischen Glied eben nicht ab, und auf groben Gittern
tragen die folgenden Glieder noch spürbar bei. Aus der Reihe fällt allein der
letzte Wert 4,17: Bei 80 Zellen je Wellenlänge beträgt der Fehler nur noch
125 ppm, sodass die oben genannte Restdifferenz der Auswertung von etwa 6 ppm ihn
bereits merklich verschiebt.

Ebenso ist die Richtungsabhängigkeit geklärt: Der rechte Teil von @abb:disp zeigt
die aus @disp2d berechnete Kurve über den gesamten Winkelbereich, und die
Messwerte folgen ihr; die größte Abweichung beträgt 111 ppm bei einer Streuung
von 30 ppm und ist damit rund zwanzigmal kleiner als der Effekt selbst.

Schließlich enthält @disp2d auch die Courant-Zahl, und deren Einfluss läuft der
Anschauung zuwider: Je größer der Zeitschritt, desto kleiner die Abweichung.
Längs der Achse steigt $tilde(c)_p \/ c$ bei 20 Zellen je Wellenlänge von
0,9961161 bei $S$ = 0,35 auf 0,9978899 bei $S$ = 0,99, und jeder dieser Werte
stimmt mit @disp auf 5 ppm überein. In der Diagonalen hebt bei $S = 1$ der
Arkussinus in @disp den Sinus auf, sodass $tilde(c)_p \/ c = 1$ folgt --- exakt
und für jede Auflösung. An der Stabilitätsgrenze läuft eine diagonal
fortschreitende Welle demnach genau so schnell wie in Wirklichkeit; ein solcher
Zeitschritt heißt bei Taflove und Hagness _magic time-step_ @Taflove2005. Längs
der Achse gilt das nicht, weil dort der Faktor $sqrt(2)$ _innerhalb_ des
Arkussinus steht. Da das Programm nur $S < 1$ zulässt, ist er nicht exakt
erreichbar; mit $S$ = 0,99 liegt die Diagonale mit 0,0041 % aber schon nahe
daran. Ein großer Zeitschritt ist damit doppelt erstrebenswert: Die Rechnung wird
schneller und zugleich genauer.

=== Weitere Messungen: Laufweg, Amplitude und Feld

Eine zu kurze Wellenlänge ist kein einmaliger Fehler, sondern einer, der sich
mit jeder Periode erneut aufaddiert. Gemessen wurde das in einem Kanal von 56
Wellenlängen Länge.

#figure(
  image("../abbildungen/abb_verzug.pdf", width: 52%),
  caption: [Verzug gegenüber der exakten Welle über der zurückgelegten
    Strecke. Die gestrichelten Geraden sind nicht angepasst, sondern folgen
    aus @disp.],
) <abb:verzug>

Wie @abb:verzug zeigt, wächst der Verzug streng linear. Gemessen wurden 0,002119 Wellenlängen je
durchlaufener Wellenlänge bei 20 Zellen und 0,008750 bei zehn; die aus @disp
berechneten Werte lauten 0,002119 und 0,008759. Eine halbe Wellenlänge Verzug
--- die simulierte Welle schwingt dann gegenphasig zur wirklichen --- ist
demnach nach 236 beziehungsweise 57 Wellenlängen erreicht. Zu jeder
Genauigkeitsangabe gehört deshalb die Größe des Rechengebiets.

Die zweite Messung betrifft nicht die Phase, sondern die Amplitude und damit den
zweiten Sollwert aus @hankel.

#figure(
  image("../abbildungen/abb_amplitude.pdf", width: 48%),
  caption: [Amplitude der Kreiswelle über dem Abstand von der Quelle,
    verglichen mit dem zwei- und dem dreidimensionalen Potenzgesetz.],
)

Die Ausgleichsrechnung über 1,5 bis 7,5 Wellenlängen liefert für den Exponenten
in $A prop r^(-p)$ den Wert 0,5002 ± 0,0001. Der Vergleichswert ist dabei nicht
exakt 0,5: Dieselbe Rechnung auf der exakten Lösung @hankel ergibt 0,4997, weil
diese noch Terme der Ordnung $1\/(k r)$ enthält. Das dreidimensionale Gesetz
liegt weit außerhalb; die Simulation löst also nachweislich das zweidimensionale
Problem.

#figure(
  image("../abbildungen/abb_differenzfeld.pdf", width: 88%),
  caption: [Betrag der Differenz zwischen simuliertem und exaktem Feld für drei
    Auflösungen, bezogen auf die größte im Gebiet auftretende Amplitude. Die
    Quellstärke ist auf dem Ring bei zwei Wellenlängen angepasst.],
) <abb:diff>

Zieht man @hankel punktweise vom simulierten Feld ab, erhält man ein Bild des
Fehlers über das ganze Gebiet. @abb:diff fasst die bisherigen Aussagen zusammen:
Der Fehler wächst nach außen, weil sich der Verzug mit jeder durchlaufenen
Wellenlänge aufaddiert; er ist längs der Gitterachsen am größten und in den
Diagonalen am kleinsten, woraus das vierzählige Muster entsteht; und er nimmt mit
feinerem Gitter quadratisch ab. Gemittelt über den Ring bei sechs Wellenlängen
beträgt er 15,1 %, 3,8 % und 1,2 % der dortigen Amplitude, was Verhältnissen von
3,9 und 3,2 entspricht.

Die Größenordnung lässt sich vorhersagen: Eine Phasenverschiebung $Delta phi$
erzeugt eine Felddifferenz vom Betrag $2 abs(sin(Delta phi \/ 2))$. Setzt man
dafür den über vier Wellenlängen aufsummierten Verzug aus @disp2d ein und
mittelt über alle Richtungen, ergeben sich 13,1 %, 3,2 % und 0,8 % --- 15 bis
35 % unter der Messung, obwohl nur der Phasen- und nicht der Amplitudenfehler
eingeht. Der Feldfehler ist damit keine eigenständige Fehlerquelle, sondern
dieselbe wie zuvor.

=== Diskussion

Zunächst passte das Differenzfeld nicht zum erwarteten Bild: Bei 40 Zellen sank
der Fehler nicht weiter, sondern blieb bei etwa 2,5 % stehen, und war überdies
nahezu richtungsunabhängig --- gerade nicht das, was numerische Dispersion
erzeugt. Die Aufspaltung in Betrag und Phase zeigte eine Welligkeit des Betrags
von 4,5 %, radial mit halber Wellenlänge periodisch: die Signatur einer stehenden
Welle und damit einer Reflexion am Gebietsrand. Eine Verdopplung der
absorbierenden Schicht von 0,75 auf 1,5 Wellenlängen senkte sie auf 0,7 % und
stellte das quadratische Verhalten wieder her. Testfall 1 hat damit
unbeabsichtigt einen Mangel aufgedeckt, der Gegenstand von Testfall 4 ist.

== Testfall 2 --- Der PEC-Hohlraumresonator

=== Sollgrößen

Da an einer ideal leitenden Wand nach Abschnitt 3.6 das elektrische Feld
verschwindet, können nur solche Wellen bestehen, die an allen vier Wänden einen
Knoten besitzen --- wie bei einer eingespannten Saite, nur in zwei Richtungen
zugleich. Daraus folgen die von Griffiths @Griffiths1999 hergeleiteten
Eigenfrequenzen

$ f_(m n) = c/2 sqrt((m/L_x)^2 + (n/L_y)^2), quad m, n >= 1 $ <eigenf>

worin $L_x$, $L_y$ die Kantenlängen sind und $m$, $n$ angeben, wie viele
Halbwellen zwischen die Wände passen. Die Einschränkung $m, n >= 1$ ist
wesentlich: Bei $m = 0$ gäbe es überall einen Knoten und damit gar kein Feld.

@eigenf gilt im Kontinuum. Auf dem Gitter sind die Wellenzahlen durch die Wände
ebenso exakt vorgegeben, die zugehörige Frequenz folgt aber aus @disprel statt
aus der Kontinuumsbeziehung. Setzt man $tilde(beta)_x = m pi \/ L_x$ und
$tilde(beta)_y = n pi \/ L_y$ dort ein und löst nach der Frequenz auf, erhält
man einen zweiten Sollwert $tilde(f)_(m n)$: @eigenf sagt voraus, was die
_Physik_ verlangt, $tilde(f)_(m n)$, was das _Verfahren_ liefern muss --- der
Vergleich mit beiden trennt Diskretisierungs- von Umsetzungsfehler. Als dritte
Sollgröße dient die Energie, denn ein geschlossener, verlustfreier Hohlraum kann
keine abgeben.

=== Aufbau und Messverfahren

Ein Rahmen aus ideal leitenden Zellen umschließt einen Bereich von 0,600 m ×
0,400 m. Der absorbierende Rand liegt außerhalb und bleibt damit wirkungslos,
sodass dieser Testfall als einziger vollständig geschlossen ist. Angeregt wird
mit einem kurzen Gauß-Puls an drei Punkten, aufgezeichnet an drei weiteren.

Quellen und Sonden liegen bewusst unsymmetrisch: Läge eine Quelle auf einer
Knotenlinie, würde die betreffende Mode gar nicht angeregt und fehlte im
Spektrum, ohne dass dies ein Fehler des Verfahrens wäre. Die Peaksuche läuft
unabhängig von der Vorhersage ab --- erst werden alle Linien gesucht, dann
zugeordnet ---, denn nur so kann sowohl eine fehlende als auch eine überzählige
Linie auffallen. Zwei Feinheiten sind nötig: Die Nebenmaxima des Hann-Fensters,
2,4 Frequenzstufen neben jeder Linie mit 2,7 % ihrer Höhe, werden verworfen, und
die Lage jeder Linie wird durch eine Parabel durch das Maximum und seine
Nachbarn verfeinert.

=== Ergebnis: Eigenfrequenzen

Im Band von 0,3 bis 1,2 GHz sagt @eigenf acht Frequenzen voraus. Das gemessene
Spektrum enthält genau acht Linien, von denen sich jede einer Mode zuordnen
lässt; zusätzliche treten nicht auf. Das ist ebenso aussagekräftig wie die
Genauigkeit der vorhandenen, denn ein falsch umgesetzter Rand würde
Schwingungen erzeugen, die es nach @eigenf gar nicht geben darf.

#figure(
  image("../abbildungen/abb_resonator.pdf", width: 88%),
  caption: [Links das gemessene Spektrum; die gepunkteten Linien sind die
    Eigenfrequenzen nach @eigenf, ihr sichtbarer Versatz zu den Messlinien ist
    der zu erklärende Fehler. Rechts der Betrag der relativen Abweichung jeder
    Mode gegenüber beiden Sollwerten, aufgetragen über der Zahl der Zellen, die
    auf ihre Wellenlänge entfallen.],
) <abb:res>

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (center, center, center, center, center),
    table.header([Mode], [$N_lambda$], [gemessen / MHz], [gegen @eigenf],
                 [gegen $tilde(f)_(m n)$]),
    [(1,1)], [20,0], [450,2257], [−0,035 %], [−0,03 ppm],
    [(2,1)], [14,4], [624,3224], [−0,039 %], [+0,04 ppm],
    [(1,2)], [11,4], [786,6743], [−0,424 %], [+0,00 ppm],
    [(3,1)], [10,7], [835,6438], [−0,275 %], [−0,03 ppm],
    [(2,2)], [10,0], [899,4992], [−0,140 %], [+0,00 ppm],
    [(3,2)], [8,5], [1059,6791], [−0,023 %], [+0,27 ppm],
    [(4,1)], [8,4], [1059,8416], [−0,695 %], [−0,15 ppm],
    [(1,3)], [7,8], [1138,2420], [−1,164 %], [−0,01 ppm],
  ),
  caption: [Eigenfrequenzen bei einer Zellweite von $L_x \/ 18$, entsprechend
    20 Zellen je Wellenlänge für die Grundmode. Aufzeichnungsdauer 30 µs.],
) <tab:res>

@tab:res enthält das zentrale Ergebnis, und der rechte Teil von @abb:res zeigt
es auf einen Blick. Gegenüber der Kontinuumsformel weichen die Frequenzen um bis
zu 1,16 % ab, durchgängig nach unten und umso stärker, je weniger Zellen auf die
Wellenlänge der Mode entfallen --- die roten Punkte steigen nach links an.
Gegenüber der Gitterlösung bleibt dagegen höchstens 0,27 ppm, im Mittel 0,07 ppm;
die blauen Punkte liegen ohne erkennbaren Trend am unteren Rand und damit vier
bis sechs Größenordnungen tiefer. Die gesamte Abweichung von der Physik wird also
durch die Dispersionsrelation erklärt, und es bleibt nichts übrig, was auf einen
Umsetzungsfehler hindeuten könnte.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (center, center, center, center),
    table.header([$N_lambda$ der Grundmode], [Gitter], [mittlere Abweichung],
                 [Verhältnis]),
    [10,0], [30 × 27], [2,00 %], [---],
    [20,0], [39 × 33], [0,183 %], [11,0],
    [39,9], [57 × 45], [0,0453 %], [4,03],
    [79,9], [93 × 69], [0,0113 %], [4,01],
  ),
  caption: [Auflösungsstudie: mittlere Abweichung der fünf Moden unterhalb
    1 GHz von der Kontinuumsformel @eigenf.],
) <tab:res_konv>

Ab 20 Zellen je Wellenlänge betragen die Verhältnisse bei Verdopplung 4,03 und
4,01, entsprechend einer Ordnung von 2,00. Der erste Schritt fällt mit 11,0 aus
der Reihe, weil bei zehn Zellen für die Grundmode auf die höheren Moden nur noch
vier bis fünf entfallen; dort ist die Reihenentwicklung, aus der @ordnung folgt,
nicht mehr gültig, und die Peaksuche findet nur vier statt fünf Linien. Die
Studie zeigt damit zugleich, wo das Verfahren aufhört, sich gutartig zu
verhalten.

=== Ergebnis: Energie und Stabilität

Zunächst ist zu klären, welche Größe überhaupt erhalten sein muss. Die
naheliegende Summe $W_"naiv"^n = epsilon\/2 abs(E_z^n)^2 + mu\/2
abs(bold(H)^(n+1\/2))^2$ ist es nicht, denn beide Felder sind um einen halben
Zeitschritt versetzt, und die Summe erfasst deshalb stets einen etwas falschen
Punkt des Pendelns zwischen ihnen. Erhalten ist die von Taflove und Hagness
@Taflove2005 angegebene Größe, in der beide Magnetfeldwerte symmetrisch um den
Zeitpunkt des elektrischen Feldes liegen:

$ W^n = epsilon/2 abs(E_z^n)^2
      + mu/2 bold(H)^(n-1\/2) dot bold(H)^(n+1\/2) $ <energie>

#figure(
  image("../abbildungen/abb_stabilitaet.pdf", width: 88%),
  caption: [Links die beiden Energiegrößen über der Zeit, rechts die Hüllkurve
    der größten Feldamplitude für fünf Courant-Zahlen.],
) <abb:stab>

Über 300 000 Zeitschritte, entsprechend 10 517 Perioden der Grundmode, schwankt
$W_"naiv"$ um 34,6 %, während $W$ nach @energie um $1 dot 10^(-12)$ % schwankt;
zehn über die Laufzeit verteilte Abschnittsmittel stimmen auf vierzehn
Nachkommastellen überein. Die Energie ist damit nicht näherungsweise, sondern bis
auf Rundungsfehler erhalten --- eine stärkere Aussage als eine kleine gemessene
Abweichung, denn sie bestätigt, dass das Schema die richtige Erhaltungsgröße
besitzt. Zugleich ist die scheinbare Schwankung von 34,6 % gar kein Fehler,
sondern ein Artefakt der gewählten Auswertungsgröße: Die größte in diesem Kapitel
beobachtete Abweichung stammt nicht aus der Rechnung.

Die Stabilität hängt daran, dass der Zeitschritt die Bedingung aus Abschnitt 3.5
einhält. Da die Konstruktion aus Abschnitt 4.2 eine Verletzung verhindert, musste
er nachträglich überschrieben werden; weil die Materialkoeffizienten und die
Randprofile aus Abschnitt 4.5 von $Delta t$ abhängen, sind sie danach neu zu
berechnen. Als Anfangsbedingung dient die höchste im Resonator mögliche Mode,
damit die Messung nicht darauf angewiesen ist, dass Rundungsfehler sie zufällig
anregen.

Die CFL-Bedingung gilt nach Taflove und Hagness @Taflove2005 für das unendliche
Gitter, auf dem Wellenzahlen bis zur Nyquist-Grenze vorkommen. In einem endlichen Resonator ist das
Spektrum diskret und erreicht diese Grenze nicht ganz; die tatsächliche Grenze
liegt deshalb etwas über $S = 1$. Berechnet man sie aus der höchsten wirklich
vorhandenen Mode, ergibt sich für die verwendeten 17 × 11 freien Zellen der Wert
1,00622.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (center, center, center),
    table.header([$S$], [Wachstum je Zeitschritt], [Verhalten]),
    [0,99000], [unter $10^(-7)$], [stabil],
    [1,00400], [unter $10^(-7)$], [stabil],
    [1,00620], [unter $10^(-7)$], [stabil],
    [1,00625], [$+1,6 dot 10^(-2)$], [Überlauf nach 10 918 Schritten],
    [1,00700], [$+7,9 dot 10^(-2)$], [Überlauf nach 2 300 Schritten],
    [1,01000], [$+1,7 dot 10^(-1)$], [Überlauf nach 1 051 Schritten],
  ),
  caption: [Verhalten an der Stabilitätsgrenze. Unterhalb des Umschlags schwankt
    die gemessene Rate um null und wechselt von Lauf zu Lauf das Vorzeichen;
    angegeben ist deshalb nur eine obere Schranke. Der Umschlag liegt zwischen
    1,00620 und 1,00625, die Vorhersage bei 1,00622.],
) <tab:stab>

Der gemessene Umschlag stimmt mit der Vorhersage auf besser als $5 dot 10^(-5)$
überein; oberhalb der Grenze wächst die Amplitude exponentiell und überschreitet
innerhalb weniger tausend Schritte jeden darstellbaren Wert.

=== Diskussion

Der Stabilitätsversuch ist der einzige, bei dem ein Versagen erwartet wird, und
gerade deshalb wichtig: Er zeigt, dass die Prüfanordnung einen Fehler überhaupt
sichtbar machen kann. Für den allgemeinen Gebrauch bleibt es dennoch bei
$S < 1$, denn ein offenes Gebiet enthält im Gegensatz zum Resonator Wellenzahlen
bis zur Nyquist-Grenze.

Die Genauigkeitsgrenze liegt hier nicht beim Gitter, sondern bei der Auswertung:
Eine Fourier-Transformation über 30 µs löst 33 kHz auf, also $7 dot 10^(-5)$ der
Grundfrequenz; dass die Frequenzen dennoch auf 0,3 ppm stimmen, verdankt sich der
Parabelverfeinerung. Zwei Moden liegen nur 0,16 MHz auseinander und verschmolzen
bei der zunächst gewählten Dauer von 3 µs zu einer Linie --- ein feineres Gitter
hätte hier nichts genützt, eine längere Simulation dagegen alles.

== Testfall 3 --- Reflexion an einer Materialgrenzfläche

=== Sollgröße

Trifft eine Welle senkrecht auf die Grenze zweier Medien, so folgt der
reflektierte Anteil aus der Fresnel-Formel, die Griffiths @Griffiths1999 daraus
herleitet, dass elektrisches und magnetisches Feld an einer Grenzfläche keinen
Sprung machen dürfen:

$ r = (n_1 - n_2)/(n_1 + n_2), quad "mit" quad n = sqrt(epsilon_r) $ <fresnel>

Darin ist $r$ das Verhältnis der reflektierten zur einfallenden Feldamplitude.
Ist das zweite Medium leitfähig, tritt an die Stelle der Permittivität die
komplexe Größe, mit der Sullivan @Sullivan2013 verlustbehaftete Medien
beschreibt:

$ underline(epsilon)_(r,2) = epsilon_(r,2) - i sigma/(omega epsilon_0) $ <epskomplex>

und @fresnel gilt unverändert weiter. Der Quotient $sigma \/ (omega epsilon_0
epsilon_r)$ heißt _Verlusttangens_ und misst, wie stark ein Material dämpft.

=== Aufbau und Messverfahren

Gerechnet wird in einem Kanal von drei Wellenlängen Höhe: Die ersten zwölf
Wellenlängen sind Vakuum, dahinter füllt das zu prüfende Material die restlichen
acht. Die Quelle belegt den gesamten linken Querschnitt, sodass eine ebene Welle
entsteht, die senkrecht auf die Grenzfläche trifft --- genau der Fall, für den
@fresnel gilt. Die Sonde liegt etwa auf halbem Weg zwischen Quelle und
Grenzfläche.

Der Reflexionskoeffizient lässt sich dort nicht unmittelbar ablesen, denn im
Gitter steht an jeder Stelle nur die Summe aus einfallender und reflektierter
Welle. Verwendet wird deshalb ein _Differenzverfahren_: Derselbe Aufbau wird
zweimal gerechnet, einmal mit und einmal ohne Materialsprung. Da beide Läufe bis
auf das Material identisch sind, ist die Differenz beider Aufzeichnungen an
derselben Sonde exakt das reflektierte Feld; aus dem Verhältnis seiner
Fourier-Transformierten zu der des Referenzlaufs folgt der Reflexionsgrad.

Angeregt wird mit einem Gauß-Puls über ein Band von 0,35 bis 2,4 GHz. Ein
einziger Lauf liefert damit den Reflexionsgrad für alle enthaltenen Frequenzen
--- und weil eine Frequenz zugleich eine Auflösung bedeutet, $N_(lambda,2) = c
\/ (f n_2 Delta x)$, ist die Auflösungsstudie in derselben Messung bereits
enthalten. Das Band entspricht bei $epsilon_r = 4$ Auflösungen von 8 bis 57
Zellen je Wellenlänge im Medium.

Als Bezugsfall dient durchgehend $epsilon_r = 4$: Die Zellweite ist so gewählt,
dass auf eine Wellenlänge in diesem Material genau 20 Zellen entfallen, und sie
bleibt in allen Läufen dieselbe. Für die übrigen Materialien folgt $N_(lambda,2)$
daraus zwangsläufig. Da die Wellenlänge im Medium um den Faktor $sqrt(epsilon_r)$
kürzer ist als im Vakuum, bei fester Zellweite also weniger Zellen auf sie
entfallen, sinkt $N_(lambda,2)$ mit wachsendem $epsilon_r$ --- von 32,7 bei
$epsilon_r$ = 1,5 auf 13,3 bei $epsilon_r$ = 9.

#figure(
  image("../abbildungen/abb_grenzflaeche.pdf", width: 88%),
  caption: [Links die Trennung von einfallendem und reflektiertem Feld durch
    Differenzbildung, rechts der gemessene Reflexionsgrad über der Auflösung im
    Medium; punktiert die Fresnel-Werte nach @fresnel.],
) <abb:grenz>

Die Nachweisgrenze ergibt sich aus dem Referenzlauf selbst: Nach dem Durchgang
des Pulses darf dort kein Signal mehr eintreffen. Was dennoch ankommt, stammt vom
absorbierenden Rand und beträgt 2,6 $dot 10^(-9)$ der Pulsamplitude --- sechs
Größenordnungen unter allen gemessenen Abweichungen.

=== Ergebnis

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (center, center, center, center, center),
    table.header([$epsilon_r$], [$N_(lambda,2)$], [$abs(r)$ gemessen],
                 [nach @fresnel], [Abweichung]),
    [1,50], [32,7], [0,10179], [0,10102], [+0,76 %],
    [2,25], [26,7], [0,20187], [0,20000], [+0,93 %],
    [4,00], [20,0], [0,33750], [0,33333], [+1,25 %],
    [6,00], [16,3], [0,42667], [0,42020], [+1,54 %],
    [9,00], [13,3], [0,50948], [0,50000], [+1,90 %],
  ),
  caption: [Reflexionsgrad bei 1 GHz für fünf verlustfreie Dielektrika, alle mit
    derselben Zellweite gerechnet. $N_(lambda,2)$ nimmt deshalb mit wachsendem
    $epsilon_r$ ab.],
) <tab:fresnel>

Der Reflexionsgrad wird auf ein bis zwei Prozent getroffen; die Abweichung ist in
allen fünf Fällen positiv und wächst von 0,76 auf 1,90 %, also systematisch und
nicht zufällig. Zwei Ursachen fallen darin allerdings zusammen, denn mit
$epsilon_r$ wächst nicht nur der Kontrast an der Grenzfläche, es sinkt zugleich
die Auflösung im Medium. Die folgende Messreihe trennt beides: Sie hält
$epsilon_r = 4$ fest und verändert allein die Auflösung, indem sie denselben Lauf
bei verschiedenen Frequenzen auswertet.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (center, center, center),
    table.header([$N_(lambda,2)$], [$abs(r)$ gemessen], [Abweichung]),
    [10], [0,35069], [+5,21 %],
    [15], [0,34082], [+2,25 %],
    [20], [0,33750], [+1,25 %],
    [30], [0,33517], [+0,55 %],
    [40], [0,33437], [+0,31 %],
  ),
  caption: [Auflösungsstudie für $epsilon_r = 4$, gewonnen aus demselben Lauf
    wie @tab:fresnel durch Auswertung bei verschiedenen Frequenzen. Die
    Ausgleichsrechnung liefert die Ordnung 2,03.],
) <tab:fresnel_konv>

Auch dieser Fehler ist damit ein Diskretisierungsfehler zweiter Ordnung und
verschwindet mit feinerem Gitter.

Bisher waren alle Materialien verlustfrei, sodass die Leitfähigkeit im Programm
gar nicht beansprucht wurde. Die dritte Messreihe holt das nach: Sie behält
$epsilon_r = 4$ bei und schaltet zusätzlich eine Leitfähigkeit $sigma$ hinzu.
Damit tritt an die Stelle der reellen Permittivität die komplexe Größe aus
@epskomplex, und geprüft werden nun auch die Koeffizienten `Ca` und `Cb` aus
Abschnitt 4.5.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (center, center, center, center, center),
    table.header([$sigma$ / (S/m)], [Verlusttangens], [$abs(r)$ gemessen],
                 [nach @fresnel], [Abweichung]),
    [0,01], [0,045], [0,33160], [0,33374], [−0,64 %],
    [0,05], [0,225], [0,34251], [0,34322], [−0,21 %],
    [0,20], [0,899], [0,44065], [0,43723], [+0,78 %],
    [1,00], [4,494], [0,71104], [0,69795], [+1,88 %],
  ),
  caption: [Reflexionsgrad bei 1 GHz für $epsilon_r = 4$ mit Leitfähigkeit,
    verglichen mit @fresnel bei komplexer Permittivität nach @epskomplex.],
)

Der Verlusttangens überstreicht dabei zwei Größenordnungen, vom nahezu
verlustfreien Dielektrikum bis zu einem Material, dessen Leitungsstrom den
Verschiebungsstrom um das Viereinhalbfache übertrifft. Die Abweichung bleibt
unter zwei Prozent und damit in derselben Größenordnung wie im verlustfreien
Fall; die Koeffizienten `Ca` und `Cb` sind damit ebenfalls geprüft.

=== Diskussion

Die Abweichung hat mindestens zwei Ursachen, die sich nicht vollständig trennen
lassen. Die erste ist die numerische Dispersion aus Testfall 1: Die Welle wird im
dichteren Medium stärker verzögert, weil dort weniger Zellen auf eine Wellenlänge
entfallen, sodass der wirksame Brechungsindex unterschiedlich stark zu groß
ausfällt. Setzt man die daraus folgenden wirksamen Indizes $n = c \/ tilde(c)_p$ nach
@verhaeltnis in @fresnel ein, ergibt sich
+0,25 % für $epsilon_r$ = 1,5 bis +0,63 % für $epsilon_r$ = 9 --- Vorzeichen und
Trend stimmen, doch erklärt das nur etwa ein Drittel des Betrags. Der Rest geht
auf die Darstellung der Grenzfläche zurück: Die Permittivität springt zwischen
zwei Zellen, während der zugehörige $E_z$-Knotenpunkt genau auf der Sprungstelle
liegt, sodass deren wirksame Lage um bis zu eine halbe Zelle unbestimmt ist.
Beide Beiträge sind zweiter Ordnung und erzeugen gemeinsam die gemessene
Konvergenz.

== Testfall 4 --- Der absorbierende Rand

=== Sollgröße und Messverfahren

Die absorbierende Schicht aus Abschnitt 4.5 soll den Rand des Rechengebiets
unsichtbar machen; der Sollwert ist also null. Wie gut ihr das gelingt, ist nicht
vorab bekannt, denn das Modul verbindet die Formulierung von Sullivan
@Sullivan2013 mit einem Dämpfungsparameter nach Taflove und Hagness
@Taflove2005.

Gemessen wird gegen eine _Referenzlösung_: Derselbe Puls läuft einmal in einem so
großen Gebiet, dass dessen Rand während der Beobachtungsdauer keine Rolle spielt,
und einmal im kleinen Gebiet mit dem zu prüfenden Rand; die Differenz beider
Felder ist genau der Fehler, den der Rand verursacht. Ausgewertet wird
frequenzaufgelöst, denn die Dicke der Schicht bemisst sich an der Wellenlänge:
Ein und dieselbe Schicht ist für hohe Frequenzen dick und für tiefe dünn.
Angeregt wird mit einem Ricker-Puls, dessen Spektrum bei der Frequenz null
verschwindet --- ein Gauß-Puls hinterließe ein statisches Restfeld, das kein Rand
absorbieren kann.

Geprüft wird in zwei Anordnungen, weil der Einfallswinkel eine Rolle spielt: eine
ebene Welle, die senkrecht auf die Schicht trifft, und eine Punktquelle, deren
Wellenfront den Rand unter allen Winkeln zugleich trifft --- einschließlich
streifendem Einfall und den vier Ecken.

=== Ergebnis

#figure(
  image("../abbildungen/abb_pml.pdf", width: 88%),
  caption: [Links die Energie im Gebiet über der Zeit für vier Schichtdicken,
    rechts die Restreflexion bei 1 GHz über der Dicke, für beide Anordnungen.],
) <abb:pml>

Links in @abb:pml fällt die Energie nach dem Aufbau des Pulses stufenweise ab,
weil die verbleibende Welle zwischen den Rändern hin und her läuft und bei jedem
Auftreffen erneut gedämpft wird. Jede Stufe entspricht also einem Auftreffen, und
ihre Höhe zeigt bereits qualitativ, was die folgende Messung beziffert: Je dicker
die Schicht, desto tiefer der Abfall.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (center, center, center, center),
    table.header([Dicke $d$], [Zellen], [senkrechter Einfall],
                 [alle Winkel (Punktquelle)]),
    [0,2 $lambda$], [4], [0,094], [0,62],
    [0,25 $lambda$], [5], [0,053], [0,44],
    [0,5 $lambda$], [10], [0,0077], [0,17],
    [0,75 $lambda$], [15], [0,00020], [0,090],
    [1,0 $lambda$], [20], [0,0010], [0,051],
    [1,5 $lambda$], [30], [0,00020], [0,016],
    [2,0 $lambda$], [40], [0,00018], [0,0039],
  ),
  caption: [Restreflexion $abs(R)$ des absorbierenden Randes bei 1 GHz, gemessen
    gegen die Referenzlösung bei 20 Zellen je Wellenlänge.],
) <tab:pml>

Der Einfallswinkel entscheidet: Trifft die Welle senkrecht auf, bleibt bei zehn
Zellen weniger als ein Prozent zurück; trifft sie unter allen Winkeln auf, sind
es 17 %. Der Unterschied beträgt 27 dB und wiegt damit
schwerer als eine Verdopplung der Schichtdicke, die je nach Ausgangsdicke 11 bis
22 dB bringt. Streifender Einfall und die Ecken, in denen
zwei Schichten überlappen, sind also der Schwachpunkt --- und weil in einer
offenen Anordnung stets beide Fälle vorkommen, ist der Wert für die Punktquelle
der maßgebliche.

Zu beachten ist außerdem die Einheit der ersten Spalte von @tab:pml --- die Dicke
ist in Wellenlängen angegeben und nicht in Zellen. Das ist plausibel, denn die Schicht
muss die Welle über eine hinreichende Strecke dämpfen, und „hinreichend“ bemisst
sich an der Wellenlänge. Die Folge daraus ist jedoch bemerkenswert: Bleibt die
Zellenzahl bei einer Verfeinerung des Gitters unverändert, schrumpft die Dicke in
Wellenlängen und der Rand wird _schlechter_.

Eine eigene Messreihe bestätigt das: Hält man die Schicht bei zehn Zellen fest
und verfeinert nur das Gitter, so steigt $abs(R)$ für die Punktquelle von 0,056
bei zehn Zellen je Wellenlänge über 0,17 bei zwanzig auf 0,21 bei vierzig. Damit ist auch die
Beobachtung aus Testfall 1 erklärt: Ein erster Durchlauf der dortigen
Auflösungsstudie zeigte keine quadratische Abnahme, sondern eine Sättigung, weil
das feinere Gitter den einen Fehler verbesserte und zugleich einen anderen
verschlechterte.

=== Diskussion

Der absorbierende Rand ist damit der schwächste Bestandteil des Programms ---
allerdings nur für schräg auftreffende Wellen. Mit den voreingestellten zehn
Zellen bleiben bei senkrechtem Einfall 0,8 % zurück, was neben der numerischen
Dispersion von zwei Promille nicht ins Gewicht fiele. Unter allen Winkeln sind es
dagegen 17 %, und damit ist der Rand um zwei Größenordnungen der größere Fehler.
Wer Genauigkeit im Promillebereich anstrebt, muss die Schicht auf mindestens zwei
Wellenlängen verdicken. Eine Variation des Dämpfungsparameters brachte gut
5 dB; die Grenze liegt also in der Formulierung der Schicht, nicht in der Wahl
ihrer Parameter.

== Zusammenfassung

#figure(
  table(
    columns: (auto, 1fr, auto, auto),
    align: (left, left, center, center),
    table.header([Testfall], [Prüfgröße], [Abweichung bei $N_lambda = 20$],
                 [Ordnung $p$]),
    [1], [Phasengeschwindigkeit, längs der Achse], [0,21 %], [2,06 ± 0,02],
    [1], [Phasengeschwindigkeit, in der Diagonalen], [0,0041 %], [---],
    [1], [Amplitudengesetz $A prop r^(-1\/2)$], [0,1 %], [---],
    [2], [Eigenfrequenzen gegen @eigenf (Mittel)], [0,18 %], [2,00],
    [2], [Eigenfrequenzen gegen die Gitterlösung], [0,000007 %], [---],
    [2], [Energieerhaltung über 10 517 Perioden], [$10^(-12)$ %], [---],
    [2], [Stabilitätsgrenze], [0,005 %], [---],
    [3], [Reflexionsgrad, verlustfrei], [1,3 %], [2,03],
    [3], [Reflexionsgrad, verlustbehaftet], [< 1,9 %], [---],
    [4], [Restreflexion des Randes, senkrecht (10 Zellen)], [0,8 %], [---],
    [4], [Restreflexion des Randes, alle Winkel], [17 %], [---],
  ),
  caption: [Übersicht der geprüften Größen. Die Ordnung $p$ stammt jeweils aus
    einer eigenen Auflösungsstudie über mindestens vier Auflösungen.],
) <tab:uebersicht>

Drei Beobachtungen ergeben sich aus @tab:uebersicht. Erstens liegen die
Abweichungen bei 20 Zellen je Wellenlänge zwischen zwei Promille und zwei Prozent
--- mit einer Ausnahme, dem absorbierenden Rand bei schrägem Einfall. Zweitens folgt jede der drei mit
einer Auflösungsstudie geprüften Größen einem Fehler zweiter Ordnung; die
unabhängig gemessenen Werte 2,06, 2,00 und 2,03 bestätigen die Herleitung aus
Kapitel 3 unmittelbar. Drittens sind die Abweichungen dort am kleinsten, wo sich
der Diskretisierungsfehler vollständig vorhersagen lässt: Phasengeschwindigkeit
und Eigenfrequenzen stimmen mit der _Gitterlösung_ um Größenordnungen besser
überein als mit der Kontinuumslösung.

#figure(
  image("../abbildungen/abb_konvergenz.pdf", width: 48%),
  caption: [Die drei Auflösungsstudien in einer Auftragung. Alle drei Größen
    fallen mit derselben Steigung, obwohl sie in drei verschiedenen Anordnungen
    an drei verschiedenen Sollwerten gemessen wurden.],
) <abb:konv>

Dass drei so unterschiedliche Größen --- eine Geschwindigkeit, eine Frequenz und
ein Amplitudenverhältnis --- in @abb:konv parallel verlaufen, ist die eigentliche
Bestätigung: Alle drei gehen auf dieselbe Ursache zurück, den Abbruchfehler der
zentralen Differenzen aus Abschnitt 3.1. Bei einem Umsetzungsfehler wäre ein so
gleichförmiges Verhalten in drei verschiedenen Anordnungen nicht zu erwarten.
Bemerkenswert ist außerdem, dass die größte beobachtete Abweichung gar kein
Fehler des Verfahrens war: Der Zeitversatz des Leapfrog-Schemas täuschte eine
Energieschwankung von 35 % vor. Bei der Beurteilung einer Simulation ist deshalb
zuerst zu prüfen, ob eine Abweichung überhaupt aus der Rechnung stammt und nicht
aus der Art, wie gemessen wird.

=== Vorgaben für den weiteren Gebrauch

+ *Mindestens 20 Zellen je Wellenlänge*, bezogen auf die Wellenlänge im
  dichtesten vorkommenden Material und bei breitbandiger Anregung auf den oberen
  Rand des benötigten Bandes.
+ *Absorbierende Schicht von mindestens einer Wellenlänge Dicke*, für
  Genauigkeit im Promillebereich von zwei --- und in Wellenlängen vorzugeben,
  nicht in Zellen.
+ *Zeitschritt dicht unterhalb der Stabilitätsgrenze*, weil das die Rechnung
  nicht nur beschleunigt, sondern auch genauer macht.
+ *Geometrie so ausrichten, dass die maßgebliche Ausbreitungsrichtung nicht auf
  einer Gitterachse liegt*, sofern es auf Laufzeiten oder Winkel ankommt.
+ *Laufstrecke zu jeder Genauigkeitsangabe nennen*, da sich der Phasenfehler
  linear aufsummiert.

=== Grenzen dieser Validierung

Die Grenzfläche in Testfall 3 verläuft gitterparallel. Über den Fehler an
gekrümmten Flächen, die durch Treppenstufen angenähert werden müssen, sagt das
Kapitel nichts aus; dort würde zusätzlich die Richtungsabhängigkeit aus
Testfall 1 wirksam. Ungeprüft bleiben ferner der schräge Einfall auf eine
Materialgrenzfläche, frequenzabhängige Materialien und der Übergang zu drei
Dimensionen. Alle Aussagen gelten
ausschließlich für den in Abschnitt 3.4 festgelegten Modellrahmen: zwei
Dimensionen, $"TM"_z$-Polarisation, lineare, isotrope und nichtdispersive Medien.

Schließlich sind alle Testfälle einfach genug, um analytisch lösbar zu sein, und
das ist zugleich ihr Zweck und ihre Beschränkung: Eine Validierung kann nur
solche Fehler ausschließen, die sie zu prüfen imstande ist. Ihre Aussagekraft
beruht deshalb weniger auf der Größe der Abweichungen als auf der Anlage der
Prüfung --- jeder Testfall beansprucht einen anderen Bestandteil, verglichen wird
nicht eine einzelne Zahl, sondern ein ganzer Verlauf, der Fehler wurde nicht nur
gemessen, sondern auch vorhergesagt, und ein Versuch ist so angelegt, dass er
scheitern muss. Hinzu kommt,
dass zwei Untersuchungen zunächst nicht das erwartete Ergebnis lieferten --- die
Auflösungsstudie in Testfall 1 und die Energiemessung in Testfall 2 --- und sich
beide Male die Ursache benennen ließ; gerade sie haben mehr über das Programm
gezeigt als die gelungenen.

Damit ist das Programm nicht bloß lauffähig, sondern in seinem Fehlerverhalten
bekannt: Man weiß, wie groß der Fehler bei gegebener Auflösung ist, wie er mit
ihr abnimmt, wovon er sonst abhängt und welche Abweichungen gar nicht im
Verfahren liegen. Eine selbst geschriebene FDTD-Simulation gibt die analytisch
bekannten Eigenschaften elektromagnetischer Wellen damit auf etwa ein Prozent
genau wieder --- allerdings nur unter Bedingungen, die man kennen muss.
