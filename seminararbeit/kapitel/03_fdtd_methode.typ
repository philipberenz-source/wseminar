// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Die FDTD-Methode
Ziel dieses Kapitel ist es die zuvor beleuchteten Maxwell-Gleichungen in ein System diskreter Gleichungen zu überführen. Dabei werden Schrittweise  drei Update-Gleichungen hergeleitet, die zu jedem Zeitpunkt an einem jeden Ort die Feldverteilung liefern. Diese bilden dann das Fundament für die tatsächliche maschinelle Implementierung in Kapitel 4.
== Grundidee: Finite-Differenzen-Approximation <sec-finite-differenzen>

Da ein Rechner den Grenzübergang $Delta x arrow.r 0$ nicht ausführen kann,
müssen die in @tmz-hx bis @tmz-ez auftretenden Ableitungen durch Quotienten
endlicher --- _finiter_ --- Differenzen ersetzt werden. Wie genau eine solche
Näherung ist, zeigt die Taylor-Entwicklung einer hinreichend oft stetig
differenzierbaren Funktion $u$ um die Stelle $x$
@Taflove2005[eq(2.10a) und eq(2.10b)]:

$ u(x plus.minus Delta x) = u(x) plus.minus Delta x (d u)/(d x)
  + (Delta x)^2/2 (d^2 u)/(d x^2)
  plus.minus (Delta x)^3/6 (d^3 u)/(d x^3) + ... $ <taylor>

Alle Ableitungen der rechten Seite sind an der Stelle $x$ zu nehmen, die Zeit
bleibt festgehalten. Da @tmz-hx bis @tmz-ez --- anders als die Wellengleichung
@welle --- nur erste Ableitungen enthalten, werden die beiden Reihen @taylor
voneinander subtrahiert: Dabei heben sich sämtliche Glieder gerader Ordnung
auf, und Division durch $2 Delta x$ liefert die _zentrale Differenz_
@Taflove2005[analog zu eq(2.12)]

$ (d u)/(d x)
  = (u(x + Delta x) - u(x - Delta x))/(2 Delta x)
  + O((Delta x)^2) $ <zentraldifferenz-kont>

Die Herleitung setzt dabei nicht voraus, dass die beiden Auswertungsstellen
gerade um $Delta x$ von $x$ entfernt liegen. Ersetzt man $Delta x$ durch
$Delta x\/2$, so folgt die für das FDTD-Verfahren maßgebliche Form

$ (d u)/(d x)
  = (u(x + Delta x\/2) - u(x - Delta x\/2))/(Delta x)
  + O((Delta x)^2) $ <zentraldifferenz-halb>

Der Differenzenquotient bleibt damit genau der Stelle $x$ zugeordnet, an der
die Ableitung gebraucht wird: Er ist der Änderungsquotient über ein
symmetrisch um $x$ gelegtes Intervall der Breite $Delta x$ und nicht, wie bei
einer einseitigen Differenz, über ein zu einer Seite verschobenes. Der
verbleibende _Abbruchfehler_ ist von zweiter Ordnung --- eine Halbierung von
$Delta x$ viertelt ihn ---, wobei erst in Kapitel 5 an einer analytisch
bekannten Lösung gemessen wird, ob das Verfahren als Ganzes diese Ordnung
tatsächlich erreicht.

== Das Yee-Gitter (2D) <sec-yee>

Die zentrale Differenz @zentraldifferenz-halb lässt sich auf @tmz-hx bis
@tmz-ez nicht unmittelbar anwenden, weil diese drei Gleichungen wechselseitig
gekoppelt sind: Die Zeitableitung von $H_x$ und $H_y$ hängt nach @tmz-hx und
@tmz-hy von der Ortsableitung von $E_z$ ab, die Zeitableitung von $E_z$ nach
@tmz-ez umgekehrt von den Ortsableitungen von $H_x$ und $H_y$. Keine der drei
Komponenten lässt sich damit für sich allein fortschreiben. Hinzu kommt eine
zweite, rein numerische Schwierigkeit: @zentraldifferenz-halb verlangt zwei
Funktionswerte je eine halbe Zelle beiderseits der Stelle, an der die
Ableitung gebraucht wird. Würden alle drei Komponenten an denselben
Gitterpunkten abgelegt, stünden dort gerade keine solchen Zwischenwerte zur
Verfügung; die Ableitung müsste entweder über eine ganze Zelle hinweg nach
@zentraldifferenz-kont --- mit dem vierfachen Fehleranteil --- oder aus
interpolierten Werten gebildet werden. Beides gäbe die in
@sec-finite-differenzen gewonnene Genauigkeit wieder preis.

Erst an dieser Stelle wird das Rechengebiet in ein regelmäßiges Gitter
zerlegt: Ort und Zeit werden an den äquidistanten Stützstellen
$x_i = i Delta x$, $y_j = j Delta y$ und $t_n = n Delta t$ abgetastet, wobei
$Delta x$ und $Delta y$ die Kantenlängen einer Gitterzelle und $Delta t$ die
Länge eines Zeitschritts bezeichnen. Yee löste 1966 beide oben genannten
Schwierigkeiten mit einer einzigen Festlegung darüber, _wo_ innerhalb dieses
Gitters die einzelnen Feldkomponenten abgelegt werden @Yee1966:
Sie liegen nicht am selben Ort, sondern räumlich um
eine halbe Zelle gegeneinander _versetzt_ (englisch _staggered_) --- und zwar
so, dass jede Komponente genau dort liegt, wo die Ableitung der jeweils anderen
Feldgröße benötigt wird. Für den in @sec-zweidimensional festgelegten
TM#sub[z]-Fall bedeutet das @Schneider2010[Abschn. 8.3, S. 185f.]:

/ $E_z$: liegt auf den ganzzahligen Gitterpunkten $(i Delta x, j Delta y)$;
/ $H_x$: liegt um eine halbe Zelle in $y$-Richtung versetzt bei $(i Delta x, (j+1\/2) Delta y)$;
/ $H_y$: liegt um eine halbe Zelle in $x$-Richtung versetzt bei $((i+1\/2) Delta x, j Delta y)$.

#figplaceholder(caption: "Yee-Gitter im TM_z-Fall: E_z auf den Gitterpunkten, H_x und H_y jeweils um eine halbe Zelle versetzt. Gestrichelt eine Zelle, deren drei Komponenten im Programm dieselben Indizes tragen.")

Dass diese Anordnung genau aufgeht, zeigt der Blick auf die einzelnen
Gleichungen. @tmz-hx wird am Ort von $H_x$ ausgewertet; die dort benötigte
Ableitung $partial E_z\/partial y$ stützt sich auf die beiden $E_z$-Werte bei
$j$ und $j+1$, die diesen Punkt gerade eine halbe Zelle nach oben und unten
einschließen. Für @tmz-hy gilt dasselbe in $x$-Richtung. Umgekehrt liegen die
für @tmz-ez am Ort von $E_z$ benötigten Werte $H_y$ bei $i plus.minus 1\/2$ und
$H_x$ bei $j plus.minus 1\/2$ ebenfalls symmetrisch um diesen Punkt herum. Jede
der vier auftretenden Ortsableitungen ist damit unmittelbar eine zentrale
Differenz nach @zentraldifferenz-halb, gebildet aus den zwei nächstgelegenen
Nachbarwerten und ohne jede Interpolation.

Die wechselseitige Abhängigkeit beider Felder ist damit allerdings noch nicht
aufgelöst. Yee versetzt deshalb zusätzlich die Zeitachse: Das elektrische Feld
wird zu ganzzahligen Vielfachen $n Delta t$ ausgewertet, das magnetische zu den
dazwischenliegenden halbzahligen Zeitpunkten $(n+1\/2) Delta t$
@Schneider2010[S. 186]. Beide Felder werden dann abwechselnd fortgeschrieben
--- ein Vorgehen, das nach dem Bocksprungspiel _leapfrog_ genannt wird: Aus dem
bekannten $E_z^n$ folgt über @tmz-hx und @tmz-hy das magnetische Feld zum
Zeitpunkt $n+1\/2$, aus diesem über @tmz-ez das elektrische Feld zum Zeitpunkt
$n+1$, und so fort. Die Kopplung ist damit kein Zirkelschluss mehr, sondern
eine Abfolge: Auf der rechten Seite jeder Aktualisierungsvorschrift stehen
ausschließlich bereits berechnete Werte. Das Verfahren ist _explizit_, es muss
zu keinem Zeitpunkt ein Gleichungssystem gelöst werden --- der wesentliche
Grund für den geringen Rechenaufwand je Zeitschritt. Zu beachten ist dabei, dass
nur $bold(E)$ gegen $bold(H)$ zeitlich versetzt sein muss; $H_x$ und $H_y$
untereinander liegen auf demselben Zeitpunkt @Schneider2010[S. 186].

Die halbzahligen Indizes sind eine Schreibweise für den Ort, keine zusätzlichen
Speicherplätze: Im Programm (Kapitel 4) werden $E_z$, $H_x$ und $H_y$ als drei
Felder mit denselben ganzzahligen Indizes $(i,j)$ geführt, und die Versetzung um
eine halbe Zelle ist dabei stillschweigend mitgedacht
@Schneider2010[Abb. 8.1, S. 187].
// LITERATURLUECKE: Ein weiteres, hier bewusst weggelassenes Argument fuer das
// Yee-Gitter ist, dass die quellfreien Gleichungen (Gauss, div B = 0)
// automatisch erhalten bleiben, wenn sie zu t=0 erfuellt sind. Beleg dafuer
// steht im Yee-Originalpapier; die Projektdatei "YEE 1966.pdf" enthaelt keine
// Textebene und war nicht lesbar. Vor Aufnahme am Original pruefen.

== Herleitung der Update-Gleichungen (2D, $"TM"_z$-Fall) <sec-update>

#todo[
  Ausgangspunkt: die drei relevanten Komponenten von Faraday- und
  Ampère-Gesetz für $E_z, H_x, H_y$ herleiten. Explizite Diskretisierung mit
  Leapfrog-Zeitschema (E auf ganzzahligen, H auf halbzahligen Zeitschritten)
  Schritt für Schritt zeigen.
]

$ H_x^(n+1/2) = H_x^(n-1/2) - (Delta t)/(mu Delta y) (E_z^n (i,j+1) - E_z^n (i,j)) $

$ H_y^(n+1/2) = H_y^(n-1/2) + (Delta t)/(mu Delta x) (E_z^n (i+1,j) - E_z^n (i,j)) $

$ E_z^(n+1) = E_z^n + (Delta t)/epsilon dot (
  (H_y^(n+1/2)(i,j) - H_y^(n+1/2)(i-1,j)) / (Delta x)
  - (H_x^(n+1/2)(i,j) - H_x^(n+1/2)(i,j-1)) / (Delta y)
) $

#note[
  Erweiterung um einen $sigma E_z$-Verlustterm ergänzen, falls in Kapitel 6
  eine Leitfähigkeit $sigma > 0$ (Wände, Möbel) verwendet wird — dann auch
  hier die Update-Gleichung entsprechend anpassen und im Text vermerken.
]

== Modellannahmen und Einschränkungen

#todo[
  Für jeden Punkt kurz begründen, warum die Annahme nötig ist und was sie
  ausschließt.
]

- *2D-Reduktion:* Translationsinvarianz in $z$-Richtung, Trennung in
  $"TM"_z$-/$"TE"_z$-Moden, keine $z$-Feldkomponenten-Kopplung.
- *Lineare Medien:* $epsilon, mu$ unabhängig von der Feldstärke (keine
  nichtlinearen Effekte wie der Kerr-Effekt).
- *Isotrope Medien:* $epsilon, mu$ skalar statt tensoriell — keine
  Richtungsabhängigkeit (schließt z. B. Kristalle mit anisotroper
  Permittivität aus).
- *Nicht-frequenzabhängige (nichtdispersive) Medien:* $epsilon_r, mu_r$
  konstant über das gesamte im Puls enthaltene Frequenzband — kein
  Debye- oder Lorentz-Modell für reale, frequenzabhängige Materialien
  (z. B. feuchte Baustoffe).

== Stabilität: die Courant-Friedrichs-Lewy-Bedingung (CFL)

#todo[
  Anschauliche Herleitung: Die numerische Ausbreitungsgeschwindigkeit im
  Gitter darf die physikalische nicht unterschreiten können, d. h. eine
  Information darf pro Zeitschritt höchstens eine Gitterzelle „überspringen“.
  Formale Herleitung optional über von-Neumann-Stabilitätsanalyse (Ansatz
  $E prop e^(i(bold(k) dot bold(r) - omega t))$ in die Update-Gleichungen
  einsetzen, Bedingung für $|"Verstärkungsfaktor"| <= 1$).
]

Ergebnis für das 2D-Gitter:

$ Delta t <= 1 / (c sqrt(1/(Delta x)^2 + 1/(Delta y)^2)) $

Für $Delta x = Delta y$ vereinfacht sich dies zu:

$ Delta t <= (Delta x) / (c sqrt(2)) $

Courant-Zahl: $S = (c Delta t) / (Delta x)$.

== Randbedingungen

#todo[
  Problem der endlichen Simulationsdomäne einführen: unphysikalische
  Reflexionen am Gebietsrand. Verwendetes Verfahren benennen (z. B.
  Mur-Randbedingung 1. Ordnung oder PEC-Rand, je nach Testfall aus
  Kapitel 5). Verweis auf die genauere PML (Perfectly Matched Layer) als
  Standardverfahren setzen, siehe Ausblick (Kapitel 7).
]
