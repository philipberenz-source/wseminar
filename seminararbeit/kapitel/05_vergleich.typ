// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Validierung an analytisch lösbaren Anordnungen

Kapitel 4 hat beschrieben, wie der Löser rechnet. Ob er richtig rechnet, ist
damit nicht gezeigt. Dieses Kapitel prüft ihn an vier Anordnungen, deren Lösung
sich in geschlossener Form angeben lässt, und liefert Fehlergrenzen, ihre
Ursachen und die Bedingungen, unter denen die Rechnung verlässlich bleibt.

== Fragestellung und Vorgehen

Ein Programm kann an mehreren Stellen zugleich fehlerhaft sein; eine Anordnung,
die alle Bestandteile gleichzeitig beansprucht, könnte im Fehlerfall nicht sagen,
welcher versagt hat. Die vier Testfälle rücken deshalb jeweils einen anderen in
den Vordergrund: Nur in Testfall 1 breitet sich das Feld über eine größere
Strecke und in alle Richtungen aus, nur Testfall 2 ist geschlossen und prüft
dadurch die Zeitintegration für sich, nur Testfall 3 enthält einen
Materialkontrast, und Testfall 4 gilt dem Rand, der in allen offenen Anordnungen
mitwirkt.

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

Von einem Laborversuch unterscheidet sich diese Prüfung in einem wesentlichen
Punkt: Es gibt keine Zufallsfehler. Dasselbe Programm liefert bei jedem Lauf
denselben Wert, und die Sollseite ist eine Formel und keine Messung. Eine
Fehlerrechnung im üblichen Sinn entfällt damit. Unsicherheiten gibt es dennoch
--- Werte aus einer Ausgleichsgeraden tragen deren Standardfehler, und auch das
Auswerteverfahren rechnet nicht beliebig genau. Beides wird dort genannt, wo es
eine Aussage berührt, und sonst nicht wiederholt.

Die eigentliche Aufgabe liegt woanders: Jede Abweichung kann aus dem Verfahren
stammen oder aus der Art, wie gemessen wird. Drei Mittel halten beides
auseinander. Die *Konvergenzstudie* fragt nicht nach dem einzelnen Zahlenwert,
sondern danach, _wie_ die Abweichung kleiner wird, wenn man das Gitter
verfeinert. Aus den zentralen Differenzen folgt in Kapitel 3 ein Fehler zweiter
Ordnung, wie ihn auch Taflove und Hagness @Taflove2005 angeben:

$ epsilon prop (Delta x)^p quad "mit" quad p = 2 $ <ordnung>

Anschaulich: Halbiert man die Zellweite, so sollte der Fehler auf ein Viertel
sinken. Kommt für $p$ aus einer Ausgleichsgeraden zwei heraus, so ist nicht nur
das Ergebnis bestätigt, sondern die Herleitung dahinter. Die *Kontrollrechnung*
wendet dasselbe Messverfahren auf die _exakte_ Lösung an, deren Wert man kennt;
was es dort zu viel oder zu wenig anzeigt, stammt aus der Auswertung. Die
*Nachweisgrenze* schließlich zeigt an einer Anordnung mit Sollwert null, wie
klein eine Abweichung noch sein darf, um überhaupt erkennbar zu bleiben. Beide
werden dort ausgeführt, wo sie gebraucht werden. Hinzu kommt in Testfall 2 ein
Versuch, bei dem ein Versagen erwartet wird --- eine Prüfung, die nur bestandene
Fälle enthält, sagt wenig darüber aus, ob sie einen Fehler anzeigen könnte.

Als Maß für die Auflösung dient die Zahl der Gitterzellen je Wellenlänge,

$ N_lambda = lambda / (Delta x) $ <nlambda>

worin $lambda$ die Wellenlänge _im jeweiligen Medium_ ist; Standardwert sind 20
Zellen. Die Frequenz von 1 GHz ist ohne Belang, weil alle Ergebnisse über
$N_lambda$ angegeben und damit übertragbar sind. Zweiter Parameter ist die
Courant-Zahl $S$, so normiert, dass die Stabilitätsgrenze bei $S = 1$ liegt ---
der Parameter `courant` aus Abschnitt 4.2; sofern nichts anderes angegeben ist,
gilt $S$ = 0,99. Alle Zahlen stammen aus eigenen Läufen von `fdtd_core.py`;
Skripte und Rohdaten liegen im Anhang.

== Testfall 1 --- Ausbreitung im Vakuum

=== Anordnung und Sollgrößen

Im Vakuum breitet sich eine Welle mit Lichtgeschwindigkeit aus, und zwar in jede
Richtung gleich schnell. Formal folgt das aus der Dispersionsrelation des
Kontinuums, die für eine ebene Welle mit den Komponenten $beta_x$, $beta_y$ der
Phasenkonstante lautet @Griffiths1999

$ mu epsilon omega^2 = beta_x^2 + beta_y^2 $ <dispkont>

Daraus ergibt sich die Phasengeschwindigkeit zu $c_p = omega \/ beta = 1 \/
sqrt(mu epsilon)$. Entscheidend daran ist, dass rechts keine Frequenz mehr
steht: Im Kontinuum laufen alle Frequenzen gleich schnell. Für das Vakuum lautet
der Sollwert also $tilde(c)_p \/ c = 1$, unabhängig von Gitterweite und Winkel.
Beides wird geprüft.

Ein zweiter Sollwert betrifft die Amplitude. Eine Welle wird mit wachsendem
Abstand schwächer, weil sich dieselbe Leistung auf einen immer größeren Umfang
verteilt. In zwei Dimensionen ist dieser Umfang $2 pi r$, sodass die Amplitude
mit $r^(-1\/2)$ abfallen muss --- in drei Dimensionen wäre es eine Kugelfläche
und damit das bekannte $1\/r$-Gesetz. Beide unterscheiden sich so deutlich, dass
die Messung zugleich prüft, ob die Simulation wirklich das zweidimensionale
Problem löst. Die zugehörige exakte Lösung geben Taflove und Hagness
@Taflove2005 an; da eine Kreiswelle in alle Richtungen gleich aussieht, hängt
das Feld nur vom Abstand ab:

$ E_z (r, t) = "Re"{ A dot H_0^((1))(k r) dot e^(-i omega t) } $ <hankel>

Die _Hankel-Funktion_ $H_0^((1))$ übernimmt darin die Rolle des Kosinus im
eindimensionalen Fall, nur dass ihre Amplitude zugleich nach außen hin abnimmt:
Für große Abstände geht sie in $r^(-1\/2) cos(k r - pi\/4)$ über und liefert
damit das genannte Gesetz. Sie ist eigens programmiert, unter anderem nach einer
Polynomnäherung von Abramowitz und Stegun @AbramowitzStegun1964.

Die beiden Sollgrößen verlangen verschiedene Geometrien. Die
Phasengeschwindigkeit wird an einer *ebenen Welle* gemessen: Eine Linienquelle
regt einen Kanal an, in dem das Feld nicht von $y$ abhängt, sodass kein
Krümmungseinfluss entsteht --- dafür ist nur eine Richtung zugänglich. Für
Richtungsabhängigkeit, Amplitudengesetz und Feldvergleich dient eine
*Kreiswelle* aus einer Punktquelle, deren Messpunkte genau auf Gitterpunkte
fallen, sodass nirgends interpoliert werden muss. In beiden Fällen läuft nach dem
Einschwingen über acht Perioden eine Fourier-Transformation mit; aus Betrag und
Phase folgen dann alle Messgrößen.

#figure(
  placement: auto,
  image("../abbildungen/abb_kreiswelle.pdf", width: 26%),
  caption: [Simuliertes Feld der Kreiswelle bei 20 Zellen je Wellenlänge.
    Eingezeichnet sind die beiden Richtungen, in denen der Fehler extremal wird.],
)

Bei der Kreiswelle wird nun die Kontrollrechnung gebraucht, und zwar aus
folgendem Grund. Das Messverfahren liest die Geschwindigkeit aus der Steigung
der Phase über dem Abstand ab und unterstellt damit stillschweigend eine gerade
Wellenfront. Die einer Kreiswelle ist aber gekrümmt, und zwar umso stärker, je
näher man der Quelle kommt. Die abgelesene Steigung fällt dadurch etwas zu groß
aus, was einer etwas zu kleinen Geschwindigkeit entspricht --- ein Effekt, der
nichts mit der Simulation zu tun hat, sondern allein mit der Auswertung.

Wie groß er ausfällt, lässt sich unmittelbar bestimmen: Man wendet dasselbe
Verfahren auf die exakte Lösung @hankel an, deren Geschwindigkeit ja genau $c$
beträgt. Herauskommen müsste eins, heraus kommt 0,99986. Diese fehlenden 140 ppm
sind der Preis des Messverfahrens und kein Fehler der Rechnung; da sie in alle
Richtungen gleich hoch ausfallen, sind sie von allen Werten der Kreiswelle
abgezogen.

=== Ergebnis und Erklärung

Der Sollwert wird nicht getroffen. Gemessen wurde bei fünf Auflösungen von 5 bis
80 Zellen je Wellenlänge, und die Welle läuft in jedem Lauf zu langsam: um 3,9 %
bei fünf Zellen, um 0,21 % bei zwanzig und um 0,013 % bei achtzig. Das Vorzeichen
ist stets dasselbe, der Betrag sinkt gleichmäßig --- ein Ablesefehler scheidet
damit aus. Hinzu kommt eine zweite Verletzung: Die Geschwindigkeit hängt von der
Ausbreitungsrichtung ab. Längs der Gitterachse beträgt der Fehler bei 20 Zellen
0,21 %, in der Diagonalen dagegen nur 0,0041 % --- ein Verhältnis von 51. Im
Vakuum darf es weder das eine noch das andere geben. Die Ursache muss folglich
im Verfahren liegen.

@dispkont gilt für das Kontinuum; die entsprechende Beziehung für das Gitter
gewinnt Schneider @Schneider2010, indem er dieselbe ebene Welle in die diskreten
Update-Gleichungen einsetzt. Ausgangspunkt ist, dass eine zentrale Differenz auf
eine ebene Welle wie eine Multiplikation wirkt; für die Zeitableitung etwa gilt

$ tilde(partial)_t e^(i omega q Delta t) = i 2/(Delta t)
  sin((omega Delta t)/2) e^(i omega q Delta t) $ <schiebung>

und für die Ortsrichtungen dasselbe mit $Delta x$ beziehungsweise $Delta y$. Mit
den Abkürzungen

$ Omega = 2/(Delta t) sin((omega Delta t)/2), quad
  K_x = 2/(Delta x) sin((tilde(beta)_x Delta x)/2), quad
  K_y = 2/(Delta y) sin((tilde(beta)_y Delta y)/2) $ <ersetzung>

bedeutet eine Zeitableitung eine Multiplikation mit $i Omega$ und eine
Ortsableitung eine mit $-i K$. Jede Richtung bekommt dabei ihren eigenen
Ausdruck mit der Zellweite _dieser_ Richtung --- daher rührt die
Richtungsabhängigkeit. Mehr wird nicht gebraucht: Im Kontinuum wird aus jeder
Zeitableitung ein Faktor $i omega$ und aus jeder Ortsableitung ein Faktor
$-i beta$, und genau so entsteht @dispkont. Die diskreten Gleichungen haben
denselben Bau, nur mit $Omega$ und $K$ an diesen Stellen. Dieselbe Rechnung
liefert deshalb dieselbe Gleichung mit den ersetzten Größen:

$ mu epsilon Omega^2 = K_x^2 + K_y^2 $ <disprel>

Dass das Kontinuum darin als Grenzfall steckt, zeigt $sin(xi) approx xi$ für
kleine Argumente: Damit geht @disprel unmittelbar in @dispkont über. Der ganze
Unterschied zwischen Physik und Verfahren besteht also darin, dass dort, wo im
Kontinuum $omega$ und $beta$ stehen, im Gitter deren durch den Sinus verzerrte
Entsprechungen auftreten.

Läuft die Welle unter dem Winkel $theta$ zur $x$-Achse, so nimmt @disprel mit
$Delta x = Delta y = Delta$ und der Abkürzung $u = tilde(beta) Delta \/ 2$ die
von Schneider @Schneider2010 angegebene Form

$ sin^2( (pi S)/(sqrt(2) N_lambda) )
  = S^2/2 [ sin^2(u cos theta) + sin^2(u sin theta) ] $ <disp2d>

an, woraus genau die gemessene Größe folgt, nämlich $tilde(c)_p \/ c = pi \/
(N_lambda u)$. Bemerkenswert ist, dass darin weder Frequenz noch Zellweite
einzeln auftreten, sondern nur ihr Verhältnis. Für beliebige Winkel wird @disp2d
numerisch aufgelöst, für die beiden Sonderfälle geschlossen @Schneider2010:

$ u_"Achse" = arcsin[ sqrt(2)/S sin( (pi S)/(sqrt(2) N_lambda) ) ], quad
  u_"Diag" = sqrt(2) arcsin[ 1/S sin( (pi S)/(sqrt(2) N_lambda) ) ] $ <disp>

#figure(
  placement: auto,
  image("../abbildungen/abb_dispersion.pdf", width: 66%),
  caption: [Verhältnis von numerischer zu wahrer Phasengeschwindigkeit. Links
    die Auflösungsabhängigkeit für beide Sonderfälle und drei Courant-Zahlen,
    rechts die Richtungsabhängigkeit bei 20 Zellen je Wellenlänge.],
) <abb:disp>

Damit lässt sich die Messung ein zweites Mal prüfen, nun gegen einen anderen
Sollwert: Der erste war die Physik, also $c$; der zweite ist das, was das Gitter
nach @disp von sich aus vorhersagt. Über alle fünf Auflösungen trifft die Messung
diese zweite Vorhersage auf höchstens 6,6 ppm. Zu den bis zu vier Prozent von
eben steht das nicht im Widerspruch, denn beide Zahlen messen Verschiedenes ---
der Löser rechnet nicht falsch, sondern genau das aus, was das Gitter hergibt,
und das Gitter gibt eben nicht ganz die Physik her. Auf den einzelnen ppm-Wert
sollte man dabei nicht zu viel geben, da die Ausgleichsrechnung selbst in dieser
Größenordnung streut; aussagekräftig ist, dass alle fünf Auflösungen dasselbe
Bild ergeben und die Restdifferenz nicht mitwächst. Die Abweichung vom Sollwert
ist damit erklärt: Sie ist kein Umsetzungsfehler, sondern eine vorhersagbare
Eigenschaft des Verfahrens, die _numerische Dispersion_.

Auch die Konvergenz ist damit belegt. Verdoppelt man die Auflösung, sinkt der
Fehler auf ein Viertel; gemessen wurden Verhältnisse zwischen 4,06 und 4,50, und
die Ausgleichsgerade liefert $p$ = 2,06. Dass gerade das Quadrat auftritt, folgt
aus @ersetzung: Die Reihe $sin(u) = u - u^3\/6 plus.minus dots$ hat kein
lineares Glied, weil die zentrale Differenz ihren Punkt symmetrisch einschließt
--- ebendies begründet die zweite Ordnung aus Abschnitt 3.1. Dass die Messung
leicht darüber liegt, ist ebenfalls vorhergesagt: Wertet man @disp selbst auf
denselben Auflösungen aus, ergibt sich 2,05, denn die Reihe bricht nach dem
quadratischen Glied nicht ab.

Die Richtungsabhängigkeit klärt der rechte Teil von @abb:disp: Die Messwerte
liegen auf der aus @disp2d berechneten Kurve, die größte Abweichung beträgt rund
ein Zwanzigstel des Effekts selbst. Auch die Courant-Zahl steht in @disp2d, und
ihr Einfluss läuft der Anschauung zuwider: Je größer der Zeitschritt, desto
_kleiner_ die Abweichung --- längs der Achse sinkt der Fehler bei 20 Zellen von
0,39 % bei $S$ = 0,35 auf 0,21 % bei $S$ = 0,99.
Ein Sonderfall verdient Beachtung: Läuft die Welle diagonal und ist $S = 1$, so
heben sich in @disp Sinus und Arkussinus gerade auf, sodass $tilde(c)_p \/ c = 1$
übrig bleibt --- exakt und für jede Auflösung. An der Stabilitätsgrenze läuft
eine diagonale Welle also genau so schnell wie in Wirklichkeit; Taflove und
Hagness nennen das den _magic time-step_ @Taflove2005. Längs der Achse gilt es
nicht, weil dort der Faktor $sqrt(2)$ _innerhalb_ des Arkussinus steht. Ein
großer Zeitschritt ist damit doppelt erstrebenswert: Die Rechnung wird schneller
und zugleich genauer.

Eine zu kurze Wellenlänge ist kein einmaliger Fehler, sondern einer, der sich
mit jeder Periode erneut aufaddiert. Wie @abb:verzug zeigt, wächst der Verzug
streng linear, und die gemessene Steigung stimmt auf ein Promille mit der aus
@disp berechneten überein. Anschaulicher als die Steigung ist ihre Folge: Eine
halbe Wellenlänge Verzug --- die simulierte Welle schwingt dann gegenphasig zur
wirklichen --- ist bei 20 Zellen nach 236 durchlaufenen Wellenlängen erreicht,
bei zehn Zellen schon nach 57. Zu jeder Genauigkeitsangabe gehört deshalb die
Größe des Rechengebiets.

#figure(
  placement: auto,
  grid(
    columns: 2, gutter: 8pt,
    image("../abbildungen/abb_verzug.pdf", width: 100%),
    image("../abbildungen/abb_amplitude.pdf", width: 100%),
  ),
  caption: [Links der Verzug gegenüber der exakten Welle über der zurückgelegten
    Strecke; die gestrichelten Geraden sind nicht angepasst, sondern folgen aus
    @disp. Rechts die Amplitude der Kreiswelle über dem Abstand von der Quelle,
    verglichen mit dem zwei- und dem dreidimensionalen Potenzgesetz.],
) <abb:verzug>

Die Amplitude prüft den zweiten Sollwert. Die Ausgleichsrechnung über 1,5 bis
7,5 Wellenlängen liefert für den Exponenten in $A prop r^(-p)$ den Wert 0,5002.
Der Vergleichswert ist dabei nicht exakt 0,5: Dieselbe Rechnung auf @hankel
ergibt 0,4997, weil diese noch Terme der Ordnung $1\/(k r)$ enthält. Messung und
Kontrollrechnung liegen damit fünf Zehntausendstel auseinander, während das
dreidimensionale Gesetz mit $p = 1$ weit außerhalb liegt --- die Simulation löst
also nachweislich das zweidimensionale Problem.

#figure(
  placement: auto,
  image("../abbildungen/abb_differenzfeld.pdf", width: 66%),
  caption: [Betrag der Differenz zwischen simuliertem und exaktem Feld für drei
    Auflösungen, bezogen auf die größte im Gebiet auftretende Amplitude.],
) <abb:diff>

Zieht man @hankel punktweise vom simulierten Feld ab, erhält man ein Bild des
Fehlers über das ganze Gebiet. @abb:diff fasst alle bisherigen Aussagen zusammen:
Der Fehler wächst nach außen, weil sich der Verzug aufaddiert; er ist längs der
Gitterachsen am größten und in den Diagonalen am kleinsten, woraus das
vierzählige Muster entsteht; und er nimmt mit feinerem Gitter quadratisch ab ---
auf dem Ring bei sechs Wellenlängen von 15,1 % über 3,9 % auf 1,2 % der dortigen
Amplitude. Selbst diese Größenordnung lässt sich aus dem Verzug vorhersagen: Zwei
um $Delta phi$ verschobene Schwingungen unterscheiden sich um $2 abs(sin(Delta
phi \/ 2))$, was auf 13,1 %, 3,2 % und 0,8 % führt --- 13 bis 34 % unter der
Messung, obwohl nur der Phasen- und nicht der Amplitudenfehler eingeht. Der
Feldfehler ist damit keine eigene Fehlerquelle, sondern derselbe Fehler, nur
anders sichtbar gemacht.

== Testfall 2 --- Der PEC-Hohlraumresonator

=== Anordnung und Sollgrößen

Da an einer ideal leitenden Wand nach Abschnitt 3.6 das elektrische Feld
verschwindet, können nur solche Wellen bestehen, die an allen vier Wänden einen
Knoten besitzen --- wie bei einer eingespannten Saite, nur in zwei Richtungen
zugleich. Daraus folgen die von Griffiths @Griffiths1999 hergeleiteten
Eigenfrequenzen

$ f_(m n) = c/2 sqrt((m/L_x)^2 + (n/L_y)^2), quad m, n >= 1 $ <eigenf>

worin $L_x$, $L_y$ die Kantenlängen sind und $m$, $n$ angeben, wie viele
Halbwellen zwischen die Wände passen. Bei $m = 0$ gäbe es überall einen Knoten
und damit gar kein Feld.

@eigenf gilt allerdings im Kontinuum, und auf dem Gitter ändert sich davon nur
ein Teil. Wie viele Halbwellen zwischen die Wände passen, geben diese nach wie
vor exakt vor; welche Frequenz dazu gehört, folgt jedoch aus @disprel statt aus
der Kontinuumsbeziehung. Setzt man dort $tilde(beta)_x = m pi \/ L_x$ und
$tilde(beta)_y = n pi \/ L_y$ ein, so erhält man einen zweiten Sollwert
$tilde(f)_(m n)$. Für jede Mode stehen damit zwei Vergleichswerte bereit, und
darin liegt der eigentliche Nutzen dieses Testfalls: @eigenf sagt, was die
_Physik_ verlangt, $tilde(f)_(m n)$, was das _Verfahren_ liefern muss. Weicht
die Messung von der ersten Zahl ab, nicht aber von der zweiten, so liegt es an
der Diskretisierung und an nichts sonst. Als dritte Sollgröße dient die Energie,
denn ein geschlossener, verlustfreier Hohlraum kann keine abgeben.

Ein Rahmen aus ideal leitenden Zellen umschließt einen Bereich von 0,600 m ×
0,400 m; der absorbierende Rand liegt außerhalb und bleibt wirkungslos, sodass
dieser Testfall als einziger vollständig geschlossen ist. Angeregt wird mit
einem kurzen Gauß-Puls an drei Punkten, aufgezeichnet an drei weiteren. Quellen
und Sonden liegen bewusst unsymmetrisch: Läge eine Quelle auf einer Knotenlinie,
so würde die betreffende Mode gar nicht erst angeregt und fehlte im Spektrum,
ohne dass dies ein Fehler des Verfahrens wäre. Auch die Suche nach den Linien
läuft bewusst unabhängig von der Vorhersage ab --- erst werden alle Linien
gesucht, dann zugeordnet ---, denn nur so fällt auf, wenn eine fehlt oder eine
zu viel da ist. Zwei Zusätze sichern sie ab: Die Nebenmaxima des Fensters werden
an ihrem festen Abstand erkannt und verworfen, und die Lage jeder Linie wird
durch eine Parabel durch das Maximum und seine Nachbarn verfeinert.

=== Ergebnis

Im Band von 0,3 bis 1,2 GHz sagt @eigenf acht Frequenzen voraus. Das gemessene
Spektrum enthält genau acht Linien, jede einer Mode zuzuordnen; zusätzliche
treten nicht auf. Das ist ebenso aussagekräftig wie ihre Genauigkeit, denn ein
falsch umgesetzter Rand würde Schwingungen erzeugen, die es nach @eigenf gar
nicht geben darf.

#figure(
  placement: auto,
  image("../abbildungen/abb_resonator.pdf", width: 66%),
  caption: [Links das gemessene Spektrum; die gepunkteten Linien sind die
    Eigenfrequenzen nach @eigenf. Rechts der Betrag der relativen Abweichung
    jeder Mode gegenüber beiden Sollwerten über der Zahl der Zellen, die auf ihre
    Wellenlänge entfallen.],
) <abb:res>

Der rechte Teil von @abb:res enthält das zentrale Ergebnis. Gemessen wurde bei
einer Zellweite von $L_x \/ 18$, entsprechend 20 Zellen je Wellenlänge für die
Grundmode; auf die höchste erfasste Mode entfallen nur noch 7,8. Gegenüber der
Kontinuumsformel weichen die acht Frequenzen um bis zu 1,16 % ab, ausnahmslos
nach unten, und grob gilt: je weniger Zellen, desto größer die Abweichung.
Streng ist das nicht, denn nach Testfall 1 hängt sie auch von der Richtung ab.
Die Mode $(3,2)$ etwa hat wegen $3 pi \/ L_x = 2 pi \/ L_y$ beide
Wellenzahlkomponenten gleich groß, liegt also genau in der Diagonalen --- dort,
wo die numerische Dispersion am kleinsten ist --- und fällt deshalb trotz
geringer Auflösung als genaueste aller acht Moden aus.

Gegenüber der Gitterlösung bleibt dagegen höchstens 0,27 ppm, im Mittel
0,07 ppm; die blauen Punkte liegen ohne erkennbaren Trend drei bis sechs
Größenordnungen unter den roten. Die Abweichung von der Physik geht also
vollständig auf die Dispersionsrelation zurück, und es bleibt nichts übrig, was
auf einen Umsetzungsfehler hindeuten könnte. Eine Auflösungsstudie über vier
Gitter bestätigt auch hier die zweite Ordnung: Die mittlere Abweichung der fünf
Moden unter 1 GHz sinkt von 2,0 % bei zehn Zellen für die Grundmode auf 0,011 %
bei achtzig, mit Verhältnissen nahe vier ab 20 Zellen, entsprechend $p$ = 2,00.
Der erste Schritt fällt mit elf aus der Reihe, weil dort auf die höheren Moden
nur noch fünf bis sieben Zellen entfallen; die Reihenentwicklung hinter @ordnung
gilt dann nicht mehr. Die Studie zeigt damit zugleich, wo das Verfahren aufhört,
sich gutartig zu verhalten.

Bei der Energie ist zunächst zu klären, welche Größe überhaupt erhalten sein
muss. Elektrische und magnetische Energie einfach zu addieren, genügt nicht: Die
beiden Felder gehören zu Zeitpunkten, die einen halben Schritt auseinanderliegen,
und da die Energie zwischen ihnen pendelt, erwischt eine solche Summe stets einen
etwas falschen Punkt dieses Pendelns --- sie schwankt, obwohl nichts verloren
geht. Erhalten ist stattdessen die von Taflove und Hagness @Taflove2005
angegebene Größe, in der die beiden Magnetfeldwerte symmetrisch um den Zeitpunkt
des elektrischen Feldes liegen:

$ W^n = epsilon/2 abs(E_z^n)^2
      + mu/2 bold(H)^(n-1\/2) dot bold(H)^(n+1\/2) $ <energie>

#figure(
  placement: auto,
  image("../abbildungen/abb_stabilitaet.pdf", width: 66%),
  caption: [Links die beiden Energiegrößen über der Zeit, rechts die Hüllkurve
    der größten Feldamplitude für fünf Courant-Zahlen.],
) <abb:stab>

Der Unterschied ist drastisch. Über 300 000 Zeitschritte, entsprechend 10 517
Perioden der Grundmode, schwankt die naive Summe um 34,6 %, $W$ nach @energie
dagegen um $10^(-12)$ %. Die Energie ist damit nicht bloß näherungsweise,
sondern bis auf Rundungsfehler erhalten --- eine stärkere Aussage als eine klein
gemessene Abweichung, denn sie zeigt, dass das Schema die richtige
Erhaltungsgröße überhaupt besitzt. Bemerkenswert ist die Kehrseite: Die 34,6 %
sind kein Fehler der Rechnung, sondern ein Artefakt der falsch gewählten
Auswertungsgröße. Wer nur diese Zahl sähe, hielte den Löser für grob fehlerhaft.

Die Stabilität hängt daran, dass der Zeitschritt die Bedingung aus Abschnitt 3.5
einhält. Um sie überhaupt verletzen zu können, musste er nachträglich
überschrieben und die von ihm abhängigen Koeffizienten neu berechnet werden; als
Anfangszustand dient die höchste im Resonator mögliche Mode, damit die Messung
nicht darauf warten muss, dass Rundungsfehler sie zufällig anregen. Zu erwarten
ist dabei nicht genau $S = 1$: Die CFL-Bedingung gilt nach Taflove und Hagness
@Taflove2005 für das unendliche Gitter, auf dem Wellen bis hinunter zur
Nyquist-Grenze vorkommen. In einem endlichen Resonator ist das Spektrum diskret
und reicht nicht ganz heran, sodass die Grenze etwas darüber liegen muss ---
berechnet aus der höchsten vorhandenen Mode bei 1,00622.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (center, center, center),
    table.header([$S$], [Wachstum je Zeitschritt], [Verhalten]),
    [0,99000], [unter $10^(-7)$], [stabil],
    [1,00400], [unter $10^(-7)$], [stabil],
    [1,00620], [unter $10^(-7)$], [stabil],
    [1,00625], [$+1,7 dot 10^(-2)$], [Überlauf nach 10 918 Schritten],
    [1,00700], [$+7,9 dot 10^(-2)$], [Überlauf nach 2 300 Schritten],
    [1,01000], [$+1,7 dot 10^(-1)$], [Überlauf nach 1 051 Schritten],
  ),
  caption: [Verhalten an der Stabilitätsgrenze. Unterhalb des Umschlags schwankt
    die Rate um null und wechselt das Vorzeichen; angegeben ist deshalb nur eine
    obere Schranke. Der Umschlag liegt zwischen 1,00620 und 1,00625, die
    Vorhersage bei 1,00622.],
) <tab:stab>

Der gemessene Umschlag trifft die Vorhersage auf besser als $5 dot 10^(-5)$;
oberhalb der Grenze wächst die Amplitude exponentiell und überschreitet innerhalb
weniger tausend Schritte jeden darstellbaren Wert. Dieser Versuch ist der
einzige, bei dem ein Versagen erwartet wird, und gerade deshalb wichtig: Er
zeigt, dass die Prüfanordnung einen Fehler überhaupt sichtbar machen kann. Für
den allgemeinen Gebrauch bleibt es dennoch bei $S < 1$, denn ein offenes Gebiet
enthält Wellenzahlen bis zur Nyquist-Grenze. Die Genauigkeitsgrenze liegt hier
ohnehin bei der Auswertung: Zwei Moden liegen nur 0,16 MHz auseinander und
verschmolzen bei der zunächst gewählten Aufzeichnungsdauer zu einer Linie --- ein
feineres Gitter hätte nichts genützt, eine längere Simulation dagegen alles.

== Testfall 3 --- Reflexion an einer Materialgrenzfläche

=== Anordnung und Sollgröße

Trifft eine Welle senkrecht auf die Grenze zweier Medien, so folgt der
reflektierte Anteil aus der Fresnel-Formel, die Griffiths @Griffiths1999 daraus
herleitet, dass elektrisches und magnetisches Feld an einer Grenzfläche keinen
Sprung machen dürfen:

$ r = (n_1 - n_2)/(n_1 + n_2), quad "mit" quad n = sqrt(epsilon_r) $ <fresnel>

Darin ist $r$ das Verhältnis der reflektierten zur einfallenden Feldamplitude.
Ist das zweite Medium leitfähig, so tritt an die Stelle der Permittivität die
komplexe Größe $underline(epsilon)_(r,2) = epsilon_(r,2) - i sigma \/ (omega
epsilon_0)$, mit der Sullivan @Sullivan2013 verlustbehaftete Medien beschreibt,
und @fresnel gilt unverändert weiter. Der Quotient $sigma \/ (omega epsilon_0
epsilon_r)$ heißt _Verlusttangens_ und misst, wie stark ein Material dämpft.

Gerechnet wird in einem Kanal von drei Wellenlängen Höhe: Die ersten zwölf
Wellenlängen sind Vakuum, dahinter füllt das zu prüfende Material die restlichen
acht. Die Quelle belegt den gesamten linken Querschnitt, sodass eine ebene Welle
senkrecht auf die Grenzfläche trifft --- genau der Fall, für den @fresnel gilt.
Der Reflexionsgrad lässt sich an der Sonde nicht unmittelbar ablesen, denn dort
steht nur die Summe aus einfallender und reflektierter Welle. Verwendet wird
deshalb ein _Differenzverfahren_: Derselbe Aufbau wird zweimal gerechnet, mit und
ohne Materialsprung, sodass die Differenz beider Aufzeichnungen exakt das
reflektierte Feld ist.

Angeregt wird mit einem Gauß-Puls über ein Band von 0,35 bis 2,4 GHz. Ein
einziger Lauf liefert damit den Reflexionsgrad für alle enthaltenen Frequenzen
--- und weil eine Frequenz zugleich eine Auflösung bedeutet, steckt die
Auflösungsstudie bereits in derselben Messung. Als Bezugsfall dient durchgehend
$epsilon_r = 4$ mit 20 Zellen je Wellenlänge im Material; da die Zellweite in
allen Läufen dieselbe bleibt, sinkt $N_(lambda,2)$ mit wachsendem $epsilon_r$
zwangsläufig von 32,7 auf 13,3. Die Nachweisgrenze liefert der Referenzlauf
selbst: Was dort nach dem Durchgang des Pulses noch ankommt, stammt vom
absorbierenden Rand und beträgt $2{,}6 dot 10^(-9)$ der Pulsamplitude, sechs
Größenordnungen unter allen gemessenen Abweichungen.

=== Ergebnis

#figure(
  placement: auto,
  image("../abbildungen/abb_grenzflaeche.pdf", width: 66%),
  caption: [Links die Trennung von einfallendem und reflektiertem Feld durch
    Differenzbildung, rechts der gemessene Reflexionsgrad über der Auflösung im
    Medium; punktiert die Fresnel-Werte nach @fresnel.],
) <abb:grenz>

Der rechte Teil von @abb:grenz enthält das Ergebnis für alle fünf Dielektrika
zugleich, denn jede Kurve durchläuft ein ganzes Band von Auflösungen. Bei den
20 Zellen des Bezugsfalls wird der Reflexionsgrad auf ein bis zwei Prozent
getroffen; die Abweichung ist in allen fünf Fällen positiv und wächst von 0,76 %
bei $epsilon_r$ = 1,5 auf 1,90 % bei $epsilon_r$ = 9, ist also systematisch.

Zwei Ursachen fallen darin zusammen, denn mit $epsilon_r$ wächst nicht nur der
Kontrast an der Grenzfläche, es sinkt zugleich die Auflösung im Medium. Trennen
lassen sie sich an einer einzelnen Kurve: Verfolgt man den Bezugsfall
$epsilon_r = 4$ über das Band, so bleibt der Kontrast fest und nur die Auflösung
ändert sich. Die Abweichung sinkt dabei von 5,2 % bei zehn Zellen auf 0,31 % bei
vierzig, und die Ausgleichsrechnung liefert die Ordnung 2,03. Auch dieser Fehler
ist damit ein Diskretisierungsfehler zweiter Ordnung --- sichtbar daran, dass
alle fünf Kurven nach rechts hin auf ihre punktierte Fresnel-Linie zulaufen.

Bisher waren alle Materialien verlustfrei, sodass die Leitfähigkeit im Programm
gar nicht beansprucht wurde. Die dritte Messreihe holt das nach und prüft damit
zugleich die Koeffizienten `Ca` und `Cb` aus Abschnitt 4.5.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (center, center, center),
    table.header([$sigma$ / (S/m)], [Verlusttangens], [Abweichung von @fresnel]),
    [0,01], [0,045], [−0,64 %],
    [0,05], [0,225], [−0,21 %],
    [0,20], [0,899], [+0,78 %],
    [1,00], [4,494], [+1,88 %],
  ),
  caption: [Reflexionsgrad bei 1 GHz für $epsilon_r = 4$ mit Leitfähigkeit. Als
    einzige Messreihe des Kapitels besitzt diese keine eigene Abbildung.],
)

Der Verlusttangens überstreicht dabei zwei Größenordnungen, vom nahezu
verlustfreien Dielektrikum bis zu einem Material, dessen Leitungsstrom den
Verschiebungsstrom um das Viereinhalbfache übertrifft. Die Abweichung bleibt
unter zwei Prozent und damit in derselben Größenordnung wie im verlustfreien
Fall.

Sie hat mindestens zwei Ursachen, die sich nicht vollständig trennen lassen. Die
erste ist die numerische Dispersion aus Testfall 1: Im dichteren Medium
entfallen weniger Zellen auf eine Wellenlänge, sodass der wirksame
Brechungsindex unterschiedlich stark zu groß ausfällt. Setzt man die daraus
folgenden Indizes in @fresnel ein, so ergeben sich +0,25 % bei $epsilon_r$ = 1,5
bis +0,63 % bei $epsilon_r$ = 9 --- Vorzeichen und Trend stimmen, doch erklärt
das nur etwa ein Drittel des Betrags. Der Rest geht auf die Darstellung der
Grenzfläche zurück: Die Permittivität springt zwischen zwei Zellen, während der
zugehörige $E_z$-Knotenpunkt genau auf der Sprungstelle liegt, sodass deren
wirksame Lage um bis zu eine halbe Zelle unbestimmt ist. Beide Beiträge sind
zweiter Ordnung und erzeugen gemeinsam die gemessene Konvergenz.

== Testfall 4 --- Der absorbierende Rand

=== Anordnung und Sollgröße

Die absorbierende Schicht aus Abschnitt 4.5 soll den Rand des Rechengebiets
unsichtbar machen; der Sollwert ist also null. Wie gut ihr das gelingt, ist
nicht vorab bekannt, denn das Modul verbindet die Formulierung von Sullivan
@Sullivan2013 mit einem Dämpfungsparameter nach Taflove und Hagness
@Taflove2005.

Gemessen wird gegen eine _Referenzlösung_: Derselbe Puls läuft einmal in einem so
großen Gebiet, dass dessen Rand keine Rolle spielt, und einmal im kleinen Gebiet
mit dem zu prüfenden Rand; die Differenz beider Felder ist genau der Fehler, den
der Rand verursacht. Ausgewertet wird frequenzaufgelöst, denn ein und dieselbe
Schicht ist für hohe Frequenzen dick und für tiefe dünn. Der anregende
Ricker-Puls hat bei der Frequenz null kein Spektrum --- ein Gauß-Puls hinterließe
ein statisches Restfeld, das kein Rand absorbieren kann. Geprüft wird in zwei
Anordnungen: eine ebene Welle, die senkrecht auftrifft, und eine Punktquelle,
deren Wellenfront den Rand unter allen Winkeln zugleich trifft, einschließlich
streifendem Einfall und den vier Ecken.

=== Ergebnis

#figure(
  placement: auto,
  image("../abbildungen/abb_pml.pdf", width: 66%),
  caption: [Links die Energie im Gebiet über der Zeit für vier Schichtdicken,
    rechts die Restreflexion bei 1 GHz über der Dicke, für beide Anordnungen.],
) <abb:pml>

Links in @abb:pml fällt die Energie stufenweise ab, weil die verbleibende Welle
zwischen den Rändern hin und her läuft und bei jedem Auftreffen erneut gedämpft
wird; jede Stufe entspricht also einem Auftreffen. Entscheidend ist jedoch nicht
die Dicke allein, sondern der Einfallswinkel. Der rechte Teil zeigt es an zwei
Kurven, die bei keiner Dicke näher als um den Faktor sechs beieinanderliegen:
Trifft die Welle senkrecht auf, so bleibt bei zehn Zellen weniger als ein
Prozent zurück; trifft sie unter allen Winkeln auf, sind es 17 %. Das sind 27 dB
Unterschied --- mehr, als eine Verdopplung der Schichtdicke an irgendeiner
Stelle einbringt. Streifender Einfall und die Ecken, in denen zwei Schichten
überlappen, sind also der Schwachpunkt, und weil in einer offenen Anordnung
stets beide Fälle vorkommen, ist die Kurve der Punktquelle die maßgebliche.

Zu beachten ist außerdem die Einheit der Abszisse: Die Dicke ist in Wellenlängen
aufgetragen und nicht in Zellen. Das ist plausibel, denn die Schicht muss die
Welle über eine hinreichende Strecke dämpfen, und „hinreichend“ bemisst sich an
der Wellenlänge. Die Folge ist bemerkenswert: Bleibt die Zellenzahl bei einer
Verfeinerung des Gitters unverändert, schrumpft die Dicke in Wellenlängen und
der Rand wird _schlechter_. Eine eigene Messreihe bestätigt das --- hält man die
Schicht bei zehn Zellen fest und verfeinert nur das Gitter, so steigt $abs(R)$
für die Punktquelle von 0,056 bei zehn Zellen je Wellenlänge über 0,17 bei
zwanzig auf 0,21 bei vierzig. Damit ist auch die Beobachtung aus Testfall 1
erklärt: Die dortige Auflösungsstudie zeigte zunächst eine Sättigung statt einer
quadratischen Abnahme, weil das feinere Gitter den einen Fehler verbesserte und
zugleich einen anderen verschlechterte.

Der absorbierende Rand ist damit der schwächste Bestandteil des Programms ---
allerdings nur für schräg auftreffende Wellen, wo er die numerische Dispersion
um zwei Größenordnungen übertrifft. Wer Genauigkeit im Promillebereich anstrebt,
muss die Schicht auf mindestens zwei Wellenlängen verdicken; eine Variation des
Dämpfungsparameters brachte demgegenüber nur wenige Dezibel. Die Grenze liegt
also in der Formulierung der Schicht und nicht in der Wahl ihrer Parameter.

== Zusammenfassung

#figure(
  table(
    columns: (auto, 1fr, auto, auto),
    align: (left, left, center, center),
    table.header([Testfall], [Prüfgröße], [Abweichung bei $N_lambda = 20$],
                 [Ordnung $p$]),
    [1], [Phasengeschwindigkeit, längs der Achse], [0,21 %], [2,06],
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

#figure(
  placement: auto,
  image("../abbildungen/abb_konvergenz.pdf", width: 32%),
  caption: [Die drei Auflösungsstudien in einer Auftragung. Alle drei Größen
    fallen mit derselben Steigung, obwohl sie in drei verschiedenen Anordnungen
    an drei verschiedenen Sollwerten gemessen wurden.],
) <abb:konv>

Bei 20 Zellen je Wellenlänge bleibt keine der geprüften Größen über zwei Prozent
--- mit der einzigen Ausnahme des absorbierenden Randes bei schrägem Einfall.
Jede der drei mit einer Auflösungsstudie geprüften Größen folgt einem Fehler
zweiter Ordnung, und dass diese drei --- eine Geschwindigkeit, eine Frequenz und
ein Amplitudenverhältnis --- in @abb:konv parallel verlaufen, ist die eigentliche
Bestätigung: Alle gehen auf dieselbe Ursache zurück, den Abbruchfehler der
zentralen Differenzen aus Abschnitt 3.1. Bei einem Umsetzungsfehler wäre ein so
gleichförmiges Verhalten in drei verschiedenen Anordnungen nicht zu erwarten. Am
kleinsten sind die Abweichungen dort, wo sich der Diskretisierungsfehler
vollständig vorhersagen lässt: Phasengeschwindigkeit und Eigenfrequenzen stimmen
mit der _Gitterlösung_ um Größenordnungen besser überein als mit der
Kontinuumslösung.

Eine zweite Lehre ziehen die Fälle, in denen nicht das Verfahren danebenlag,
sondern die Messung: Der Zeitversatz des Leapfrog-Schemas täuschte eine
Energieschwankung von 35 % vor, die gekrümmte Wellenfront der Kreiswelle eine um
140 ppm zu kleine Geschwindigkeit. In beiden Fällen hätte man dem Löser einen
Fehler zugeschrieben, den er nicht gemacht hat. Bei der Beurteilung einer
Simulation ist deshalb zuerst zu prüfen, ob eine Abweichung überhaupt aus der
Rechnung stammt --- und nicht aus der Art, wie gemessen wird.

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

Die Grenzfläche in Testfall 3 verläuft gitterparallel; über den Fehler an
gekrümmten Flächen, die durch Treppenstufen angenähert werden müssen, sagt das
Kapitel nichts aus. Ungeprüft bleiben ferner der schräge Einfall auf eine
Materialgrenzfläche, frequenzabhängige Materialien und der Übergang zu drei
Dimensionen. Alle Aussagen gelten ausschließlich für den in Abschnitt 3.4
festgelegten Modellrahmen: zwei Dimensionen, $"TM"_z$-Polarisation, lineare,
isotrope und nichtdispersive Medien. Schließlich sind alle Testfälle einfach
genug, um analytisch lösbar zu sein --- was das Kapitel zeigt, ist deshalb nicht,
dass der Löser jedes Problem richtig rechnet, sondern dass er dort, wo sich sein
Ergebnis prüfen lässt, genau das liefert, was das Verfahren hergibt.