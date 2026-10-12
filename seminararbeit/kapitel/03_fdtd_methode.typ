// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Die FDTD-Methode
Ziel dieses Kapitels ist es, die zuvor beleuchteten Maxwell-Gleichungen in ein
System diskreter Gleichungen zu überführen. Dabei werden schrittweise drei
Update-Gleichungen hergeleitet, die zu jedem Zeitpunkt an einem jeden Ort die
Feldverteilung liefern. Diese bilden dann das Fundament für die tatsächliche
 Implementierung in Kapitel 4.

== Grundidee: Finite-Differenzen-Approximation <sec-finite-differenzen>

Da ein Rechner den Grenzübergang $Delta x arrow.r 0$ nicht berechnen kann,
müssen die in @tmz-hx bis @tmz-ez auftretenden Ableitungen durch Quotienten
endlicher --- _finiter_ --- Differenzen umschrieben werden. Die Genauigkeit
dieser Methode zeigt sich in der Taylor-Entwicklung einer hinreichend oft stetig
differenzierbaren Funktion $u$ um die Stelle $x$
@Taflove2005[Gl. (2.10a) und Gl. (2.10b)]:

$ u(x plus.minus Delta x) = u(x) plus.minus Delta x (d u)/(d x)
  + (Delta x)^2/2 (d^2 u)/(d x^2)
  plus.minus (Delta x)^3/6 (d^3 u)/(d x^3) + ... $ <taylor>

Anschließend werden die beiden in @taylor zusammengefassten Reihen
voneinander subtrahiert: Dabei heben sich sämtliche Glieder gerader Ordnung
auf, und Division durch $2 Delta x$ liefert die _zentrale Differenz_
@Taflove2005[analog zu Gl. (2.12)]

$ (d u)/(d x)
  = (u(x + Delta x) - u(x - Delta x))/(2 Delta x)
  + O((Delta x)^2) $ <zentraldifferenz-kont>

Ersetzt man nun $Delta x$ durch $Delta x\/2$, so folgt die für das
FDTD-Verfahren vorausgesetzte Form

$ (d u)/(d x)
  = (u(x + Delta x\/2) - u(x - Delta x\/2))/(Delta x)
  + O((Delta x)^2) $ <zentraldifferenz-halb>

Der Differenzenquotient bleibt damit genau der Stelle $x$ zugeordnet, an der
die Ableitung gebraucht wird: Er beschreibt die Änderungsrate über ein
symmetrisch um $x$ gelegtes Intervall der Breite $Delta x$. Der
verbleibende _Abbruchfehler_ ist von zweiter Ordnung --- eine Halbierung von
$Delta x$ geht somit mit einer Viertelung des Fehlers einher.

== Das Yee-Gitter (2D) <sec-yee>


Die zentrale Differenz @zentraldifferenz-halb lässt sich auf @tmz-hx bis
@tmz-ez zwar anwenden, ändert aber nichts an dem Problem der gegenseitigen Kopplung: Die Zeitableitung von $H_x$ und $H_y$ hängt nach @tmz-hx und
@tmz-hy von der Ortsableitung von $E_z$ ab, die Zeitableitung von $E_z$ nach
@tmz-ez umgekehrt von den Ortsableitungen von $H_x$ und $H_y$. Keine der drei
Komponenten lässt sich damit als alleinstehende Größe lösen.

An dieser Stelle wird das Rechengebiet in ein regelmäßiges Gitter
zerlegt: Ort und Zeit werden an den diskreten Stützstellen
$x_i = i Delta x$, $y_j = j Delta y$ und $t_n = n Delta t$ abgetastet, wobei
$Delta x$ und $Delta y$ die Kantenlängen einer Gitterzelle und $Delta t$ einen
Zeitschritt beschreiben. Ausgehend von dieser Diskretisierung löste Yee die
oben beschriebene Schwierigkeit mit einer einzigen Festlegung darüber, _wo_
innerhalb dieses Gitters die einzelnen Feldkomponenten gespeichert werden
@Yee1966: Sie liegen nicht am selben Ort, sondern räumlich um
eine halbe Zelle gegeneinander _versetzt_ (englisch _staggered_) --- und zwar
so, dass jede Komponente genau dort liegt, wo die Ableitung der jeweils anderen
Feldgröße benötigt wird. Für den in @sec-zweidimensional festgelegten
TM#sub[z]-Fall bedeutet das @Yee1966[S. 303]:

/ $E_z$: liegt auf den ganzzahligen Gitterpunkten $(i Delta x, j Delta y)$;
/ $H_x$: liegt um eine halbe Zelle in $y$-Richtung versetzt bei $(i Delta x, (j+1\/2) Delta y)$;
/ $H_y$: liegt um eine halbe Zelle in $x$-Richtung versetzt bei $((i+1\/2) Delta x, j Delta y)$.

#figure(
  image("../abbildungen/abb_yee.png", width: 55%),
  caption: [Anordnung der Feldkomponenten im Yee-Gitter für den
    TM#sub[z]-Fall: $E_z$ liegt auf den ganzzahligen Gitterpunkten, $H_x$ um
    eine halbe Zelle in $y$-Richtung versetzt, $H_y$ um eine halbe Zelle in
    $x$-Richtung. Abbildung entnommen aus @Liu2022[Abb. 4].],
) <abb:yee>

Der Sinn der in @abb:yee dargestellten Anordnung zeigt sich beim Blick auf die
einzelnen Gleichungen.
@tmz-hx wird am Ort von $H_x$ ausgewertet; die dort benötigte
Ableitung $partial E_z\/partial y$ stützt sich auf die beiden $E_z$-Werte bei
$j$ und $j+1$, die diesen Punkt gerade eine halbe Zelle nach oben und unten
einschließen. Für @tmz-hy gilt dasselbe in $x$-Richtung. Umgekehrt liegen die
für @tmz-ez am Ort von $E_z$ benötigten Werte $H_y$ bei $i plus.minus 1\/2$ und
$H_x$ bei $j plus.minus 1\/2$ ebenfalls symmetrisch um diesen Punkt herum. Jede
der vier auftretenden Ortsableitungen ist damit unmittelbar eine zentrale
Differenz nach @zentraldifferenz-halb, gebildet aus den zwei nächstgelegenen
Nachbarwerten.

Die wechselseitige zeitliche Abhängigkeit beider Felder ist damit jedoch noch
nicht aufgelöst. Yee versetzt deshalb zusätzlich die Anordnung auf der
Zeitachse: Das elektrische Feld wird zu ganzzahligen Vielfachen $n Delta t$
ausgewertet, das magnetische zu den dazwischenliegenden halbzahligen Zeitpunkten
$(n+1\/2) Delta t$ @Taflove2005[S. 60 f.]. Beide Felder werden dann abwechselnd
fortgeschrieben --- ein Vorgehen, das nach dem Bocksprungspiel _leapfrog_
benannt ist: Aus dem bekannten $E_z^n$ folgt über @tmz-hx und @tmz-hy das
magnetische Feld zum Zeitpunkt $n+1\/2$, aus diesem über @tmz-ez das elektrische
Feld zum Zeitpunkt $n+1$, und so fort. Auf der rechten Seite jeder
Aktualisierungsvorschrift stehen damit ausschließlich bereits berechnete Werte;
das Verfahren ist _explizit_, es muss zu keinem Zeitpunkt ein Gleichungssystem
gelöst werden.

Durch den zeitlichen Versatz ist die
Differenz zweier um $Delta t$ auseinanderliegender Werte ist der halbzahligen
Zwischenzeit zugeordnet, an der die jeweils andere Feldgröße gerade vorliegt.
Damit ist auch die Zeitableitung eine zentrale Differenz nach
@zentraldifferenz-halb und das Verfahren ist in
Ort _und_ Zeit von zweiter Ordnung, also $O((Delta x)^2, (Delta y)^2,
(Delta t)^2)$.

== Herleitung der Update-Gleichungen (2D, $"TM"_z$-Fall) <sec-update>

Unter Verwendung des Yee-Gitters aus @sec-yee und Bildung zentraler Differenzen
nach @zentraldifferenz-halb ergeben sich aus @tmz-hx bis @tmz-ez unmittelbar
drei explizite Aktualisierungsvorschriften, im Folgenden _Update-Gleichungen_
genannt: Für ein _quadratisches_ Gitter mit
$Delta x = Delta y = Delta$ ergeben sich folgende Gleichungen @Taflove2005[mod. Gl. (3.41)]:

$ H_x bar.v_(i\,j+1/2)^(n+1/2) = H_x bar.v_(i\,j+1/2)^(n-1/2)
  - (Delta t)/(mu_0 Delta)
    (E_z bar.v_(i\,j+1)^n - E_z bar.v_(i\,j)^n) $ <update-hx>

$ H_y bar.v_(i+1/2\,j)^(n+1/2) = H_y bar.v_(i+1/2\,j)^(n-1/2)
  + (Delta t)/(mu_0 Delta)
    (E_z bar.v_(i+1\,j)^n - E_z bar.v_(i\,j)^n) $ <update-hy>

$ E_z bar.v_(i\,j)^(n+1) &= C_a bar.v_(i\,j) dot E_z bar.v_(i\,j)^n \
  & quad + C_b bar.v_(i\,j) dot (H_y bar.v_(i+1/2\,j)^(n+1/2)
    - H_y bar.v_(i-1/2\,j)^(n+1/2)
    + H_x bar.v_(i\,j-1/2)^(n+1/2)
    - H_x bar.v_(i\,j+1/2)^(n+1/2)
    - J_z bar.v_(i\,j)^(n+1/2) Delta) $ <update-ez>

mit den allein vom Material abhängigen und deshalb vor Beginn der Zeitschleife
einmalig berechneten Koeffizienten @Taflove2005[Gl. (3.31a) und Gl. (3.31b)]

$ C_a bar.v_(i\,j) = (1 - (sigma Delta t)/(2 epsilon))/(1 + (sigma Delta t)/(2 epsilon)),
  quad
  C_b bar.v_(i\,j) = ((Delta t)/(epsilon Delta))/(1 + (sigma Delta t)/(2 epsilon)) $ <koeffizienten>

wobei $epsilon = epsilon_0 epsilon_r bar.v_(i\,j)$ und $sigma = sigma bar.v_(i\,j)$
die am Ort des betreffenden $E_z$-Wertes gültigen Materialgrößen nach @material
bezeichnen.

Ziel der folgenden Umformung ist eine Gleichung, die den neuen Wert
$E_z bar.v_(i\,j)^(n+1)$ allein aus bereits bekannten Größen berechnet. Die Form von @koeffizienten geht auf zwei Festlegungen zurück, die in @tmz-ez
noch offen sind. So wird die dort auftretende Stromdichte aufgeteilt in
einen von außen eingeprägten Quellstrom $J_z$ und den vom Feld selbst
getriebenen Leitungsstrom nach dem Ohmschen Gesetz @material:

$ J_z^"ges" = J_z + sigma E_z $ <stromaufteilung>

Damit tritt $E_z$ in @tmz-ez zweimal auf --- einmal unter der Zeitableitung,
einmal im Verlustterm. Die Zeitableitung ist nach @zentraldifferenz-halb dem
Zeitpunkt $n+1\/2$ zugeordnet; da eine Gleichung nur dann erfüllt sein kann,
wenn all ihre Glieder für _denselben_ Zeitpunkt gelten, muss auch der
Verlustterm dort ausgewertet werden. Genau dort liegt jedoch kein gespeicherter
Wert von $E_z$ vor, weil das elektrische Feld nach @sec-yee nur zu ganzzahligen
Zeitpunkten geführt wird. Der Verlustterm wird deshalb durch das arithmetische
Mittel der beiden umgebenden Werte ersetzt (_semi-implizites_ Vorgehen) @Taflove2005[S. 63, mod. Gl. (3.26)]:

$ E_z bar.v_(i\,j)^(n+1/2) approx
  (E_z bar.v_(i\,j)^(n+1) + E_z bar.v_(i\,j)^n)/2 $ <semiimplizit>

Setzt man @stromaufteilung und @semiimplizit in @tmz-ez ein, zerfällt der
Verlustterm in zwei Teile: Der eine enthält den gesuchten Wert
$E_z bar.v_(i\,j)^(n+1)$, der andere den bereits bekannten Wert
$E_z bar.v_(i\,j)^n$. Bringt man den ersten auf die linke Seite und lässt den
zweiten rechts stehen @Taflove2005[mod. Gl. (3.28)], so folgt 

$ (1 + (sigma Delta t)/(2 epsilon)) E_z bar.v_(i\,j)^(n+1)
  &= (1 - (sigma Delta t)/(2 epsilon)) E_z bar.v_(i\,j)^n \
  & quad + (Delta t)/(epsilon Delta) (H_y bar.v_(i+1/2\,j)^(n+1/2)
    - H_y bar.v_(i-1/2\,j)^(n+1/2)
    + H_x bar.v_(i\,j-1/2)^(n+1/2)
    - H_x bar.v_(i\,j+1/2)^(n+1/2)
    - J_z bar.v_(i\,j)^(n+1/2) Delta) $ <sortiert>

Division durch die Klammer der linken Seite isoliert $E_z bar.v_(i\,j)^(n+1)$
und liefert @update-ez: Die Vorfaktoren vereinfachen sich dabei zu 
$C_a$ und $C_b$.

Für verlustfreie Gebiete ($sigma = 0$) entfallen beide Hälften des
Verlustterms, $C_a = 1$ und $C_b = Delta t\/(epsilon Delta)$, und @update-ez
geht in die ungedämpfte Form über.

Mit @update-hx bis @update-ez ist die Feldentwicklung vollständig festgelegt,
sobald ein Anfangszustand vorliegt. Dieser ist in allen folgenden Rechnungen das
feldfreie Gebiet, $E_z = H_x = H_y = 0$ zum Zeitpunkt $n = 0$.

== Anregung des Feldes <sec-quellen>

Ein solches Gitter bleibt unter @update-hx bis @update-ez für alle Zeiten
feldfrei; Energie muss deshalb während der Zeitschleife eingebracht werden.
Verwendet werden dafür durchgehend _additive_ Quellen @Schneider2010[S. 46f.],
die eine Quellfunktion $f$ an einem Gitterpunkt aufaddieren --- entweder als
_weiche Quelle_ unmittelbar auf $E_z$ oder als eingeprägte Stromdichte $J_z$,
die nach @stromaufteilung über @update-ez wirkt:

$ "weich:" quad E_z bar.v_(i\,j)^(n+1) arrow.l
  E_z bar.v_(i\,j)^(n+1) + f^(n+1)
  quad quad quad
  "Strom:" quad J_z bar.v_(i\,j)^(n+1/2) = f^(n+1/2) $ <quellen>

Nur die zweite Form ist physikalisch deutbar, da $J_z$ als Stromdichte unmittelbar in @tmz-ez auftritt; beiden gemeinsam ist, dass der
Quellpunkt regulär weiter aktualisiert wird und eine zurücklaufende Welle ihn
daher ungehindert passiert.

Beide wirken an einem einzigen Gitterpunkt und erzeugen damit eine
Zylinderwelle. Für eine ebene Welle wird eine das Gebiet durchspannende
Gitterlinie senkrecht zur Ausbreitungsrichtung mit identischen Quellen besetzt:
Nach dem _Huygensschen Prinzip_ überlagern sich deren Elementarwellen zu einer
ebenen Front, sodass das Feld längs dieser Linie konstant ist. Beide strahlen nach zwei
Seiten ab; ausgewertet wird in Kapitel 5 jeweils eine davon.

== Stabilität: die Courant-Friedrichs-Lewy-Bedingung (CFL) <sec-cfl>

Damit die Simulation numerisch stabil bleibt, ist die Wahl des Zeitschrittes nicht frei, sondern muss einem Kriterium folgen.
Der Grund liegt in der Nächste-Nachbar-Kopplung des Yee-Gitters: Ein Wert an
einem Gitterpunkt beeinflusst je Zeitschritt nur seine unmittelbaren Nachbarn,
sodass sich eine Information im Gitter höchstens um eine Zelle je Schritt
ausbreitet. Physikalisch läuft sie in derselben Zeit die Strecke $c Delta t$.
Stabil kann das Verfahren daher nur sein, wenn die numerische Ausbreitung der
physikalischen mindestens folgen kann. Durch eine von-Neumann-Analyse der Update-Gleichungen ergibt sich für den _2D_-Fall die _Courant-Friedrichs-Lewy-Bedingung (CFL)_, die die Wahl des Zeitschrittes in Abhängigkeit von der Gitterbreite einschränkt @Taflove2005:

$ Delta t <= 1/(c sqrt(1\/(Delta x)^2 + 1\/(Delta y)^2))
  quad arrow.r.double quad
  Delta t <= Delta/(c sqrt(2)) $ <cfl>

wobei rechts wieder das quadratische Gitter aus @sec-update vorausgesetzt ist.
Dass darin die Vakuumlichtgeschwindigkeit $c$ steht, ist keine Einschränkung:
In jedem Medium mit $epsilon_r >= 1$ ist die Ausbreitungsgeschwindigkeit
kleiner, sodass @cfl stets die schärfere Forderung stellt.

Wie nahe die Rechnung an dieser Grenze geführt wird, gibt die _Courant-Zahl_
$S = sqrt(2) c Delta t \/ Delta$ an, die hier so normiert ist, dass die
Stabilitätsgrenze @cfl gerade bei $S = 1$ liegt. Die in der Literatur übliche
Courant-Zahl $S_c = c Delta t\/Delta$ mit der Stabilitätsgrenze
$S_c <= 1\/sqrt(2)$ hängt damit über $S_c = S\/sqrt(2)$ zusammen. In den folgenden Kapiteln wird der Zeitschritt stets über $S$ festgelegt. Welchen Einfluss die _Courant-Zahl_ auf die Simulationsergebnisse hat, wird in Kapitel 5 näher untersucht. 

== Randbedingungen <sec-rand>

=== Ideal leitender Rand (PEC) <sec-pec>

Für den ideal leitenden Rand gilt nach @sec-randbedingungen
$bold(E)_parallel = 0$. Im TM#sub[z]-Fall ist $E_z$ die einzige elektrische
Komponente und liegt für jede senkrecht zur $x y$-Ebene verlaufende Wand stets
_parallel_ zu ihr; die Bedingung verkürzt sich damit auf $E_z = 0$. Bezeichnet
$M$ die Menge der als leitend markierten Gitterpunkte, so werden deren
$E_z$-Werte nach jeder Auswertung von @update-ez wieder auf null gesetzt
@Schneider2010[S. 171f.]:

$ E_z bar.v_(i\,j)^(n+1) = 0 quad "für alle" (i,j) in M $ <pec>

Die magnetischen Komponenten werden unverändert nach @update-hx und @update-hy
fortgeschrieben. 

=== Absorbierender Rand (_Perfectly Matched Layer_, PML) <sec-pml>

An den Rändern des Simulationsgebiets fehlen benachbarte $H_y$- bzw. $H_x$-Felder, um zentrale Differenzen zu bilden, weshalb $E_z$ an diesen Stellen auf null gesetzt wird. Damit wird das Feld nach @sec-pec vollständig zurückgeworfen. In vielen Fällen ist dieser Effekt jedoch unerwünscht, insbesondere wenn die Ausbreitung im unendlichen Raum simuliert wird.
Naheliegend wäre es deshalb, den eigentlichen Rändern eine Schicht vorzulagern, welche die auslaufende Welle dämpft, sodass keine Energie zurück in das eigentliche Simulationsgebiet geworfen wird (@abb:pml-schema).

#figure(
  image("../abbildungen/abb_pml_schema.pdf", width: 48%),
  caption: [Anordnung der absorbierenden Schicht: Sie umgibt das eigentliche
    Rechengebiet auf allen vier Seiten und ist außen weiterhin durch einen
    PEC-Rand abgeschlossen; in den Ecken überlappen die Schichten beider
    Richtungen.],
) <abb:pml-schema>

Eine gewöhnliche verlustbehaftete Schicht leistet das jedoch nicht: Sie dämpft in alle
Richtungen gleich, während eine reflexionsfreie Anpassung verlangt, dass
ausschließlich die Ausbreitung _senkrecht_ zur Grenzfläche gedämpft wird
@Schneider2010[Abschn. 11.2 bis 11.4]. Im Folgenden wird deshalb zunächst die
Koordinate komplex gestreckt, daraus die Restreflexion der Schicht bestimmt und
das Ergebnis anschließend in den Zeitbereich übertragen und in die
Update-Gleichungen eingebaut.

*Komplexe Koordinatenstreckung.* Der Ausweg besteht darin, nicht das Material
zu verändern, sondern die Koordinate. Die hier verwendete Formulierung
heißt _komplexe Koordinatenstreckung_,
weil die dabei eingeführten Faktoren die einzelnen Koordinatenrichtungen
skalieren @Schneider2010[Abschn. 11.1 und 11.5]. Zu jeder Raumrichtung $w in {x, y}$ wird dabei
ein komplexer Streckfaktor eingeführt @Schneider2010[Gl. (11.86)]:

$ S_w = 1 + sigma_w/(upright(i) omega epsilon_0) $ <pml-sw>



Der Streckfaktor wird nun jeder Ortsableitung vorangestellt, und zwar derart, dass
$1\/S_w$ stets an die Ableitung in _derselben_ Richtung $w$ gebunden ist. Dafür
wird der Nabla-Operator erweitert, hier auf zwei Dimensionen
eingeschränkt @Schneider2010[mod. Gl. (11.93)]

$ tilde(nabla) = hat(x) 1/S_x partial/(partial x)
  + hat(y) 1/S_y partial/(partial y) $ <pml-nabla>

Ausgeschrieben für den TM#sub[z]-Fall aus @sec-zweidimensional bedeutet das
@Schneider2010[mod. Gl. (11.90) bis Gl. (11.92)]

$ upright(i) omega epsilon E_z = 1/S_x (partial H_y)/(partial x)
    - 1/S_y (partial H_x)/(partial y),
  quad
  upright(i) omega mu_0 H_y = 1/S_x (partial E_z)/(partial x),
  quad
  upright(i) omega mu_0 H_x = -1/S_y (partial E_z)/(partial y) $ <pml-tmz>

Wichtig ist, dass sich diese Streckung richtungsweise vorschreiben lässt.
So wird die Richtungsselektivität gewährleistet: In den
Schichten links und rechts ist $sigma_x > 0$, während $sigma_y = 0$ und damit
$S_y = 1$ bleibt; in den
Schichten oben und unten ist es umgekehrt (@abb:pml-schema). Der
Reflexionsfaktor verschwindet dadurch im Kontinuum für _jeden_ Einfallswinkel und _jede_
Frequenz. Setzt man @pml-sw in eine unter $theta$ einlaufende ebene Welle nach
@ansatz ein, so tritt an die Stelle von $k_x$ das Produkt $S_x k_x$, und der
Ortsanteil der Welle zerfällt in zwei Faktoren
@Schneider2010[Gl. (11.82) und Gl. (11.83)]:

$ e^(-upright(i) S_x k_x x - upright(i) k_y y)
  = underbrace(e^(-(sigma_x k_x)/(omega epsilon_0) x), "reell, klingt ab")
  dot underbrace(e^(-upright(i) (k_x x + k_y y)), "unveränderte Phase") $ <pml-daempfung>

*Dicke und Verlauf der Schicht.* Hinter der dämpfenden Schicht bleibt der
PEC-Rand bestehen. Eine Welle durchläuft die Schicht der Dicke $d$, wird an
diesem Rand gemäß @pec zurückgeworfen und durchläuft sie ein zweites Mal;
maßgeblich für die verbleibende Reflexion ist deshalb die auf dem gesamten Hin-
und Rückweg angesammelte Dämpfung. Für den Verlauf innerhalb der Schicht ist
ein Potenzverlauf üblich @Taflove2005[S. 292, Gl. (7.60a)]

$ sigma_x (rho) = sigma_max rho^m,
  quad rho = a/d in [0, 1] $ <pml-profil>

worin $a$ den Abstand vom inneren Schichtrand bezeichnet. Mit $k_x = (omega\/c) cos theta$ und
$1\/(c epsilon_0) = eta_0$, dem Freiraumwellenwiderstand, folgt aus dem
Dämpfungsfaktor in @pml-daempfung für den Betrag des zurückkehrenden Anteils

$ abs(R(theta)) = exp(-(2 eta_0 sigma_max d cos theta)/(m+1)) $ <pml-restreflexion-profil>

@Taflove2005[S. 292, Gl. (7.61)]. Darin tritt keine Frequenz auf. 
Logarithmiert man @pml-restreflexion-profil für senkrechten Einfall
($theta = 0$) und löst nach der Leitfähigkeit auf, so folgt
@Taflove2005[S. 293, Gl. (7.62)]

$ sigma_max = -((m+1) ln abs(R(0)))/(2 eta_0 d) $ <pml-sigmamax>

*Rückkehr in den Zeitbereich.* Ziel der folgenden zwei Schritte ist es,
@pml-tmz in eine Form zu bringen, die sich wie @update-hx bis @update-ez
Zeitschritt für Zeitschritt auswerten lässt und dabei je Gitterpunkt nur einen
zusätzlichen gespeicherten Wert benötigt. @pml-tmz steht im Frequenzbereich, das
FDTD-Verfahren rechnet jedoch im Zeitbereich, und einem Produkt im
Frequenzbereich entspricht dort eine Faltung. Benötigt wird also die
Rücktransformierte $overline(S)_w$ von $1\/S_w$. Mit @pml-sw ist

$ 1/S_w = (upright(i) omega epsilon_0)/(sigma_w + upright(i) omega epsilon_0)
  = 1 - 1/(1 + upright(i) omega epsilon_0\/sigma_w) $ <pml-inverses>

Der erste Summand liefert die Delta-Distribution $delta(t)$, die in der Faltung
die ursprüngliche Ableitung unverändert reproduziert; für den zweiten genügt das
Korrespondenzpaar $(1 + upright(i) omega tau)^(-1) arrow.l.r tau^(-1)
e^(-t\/tau) u(t)$ @Schneider2010[Gl. (11.104)]. Mit $tau = epsilon_0\/sigma_w$ ergibt sich zusammen
@Schneider2010[mod. Gl. (11.107) und Gl. (11.108)]

$ overline(S)_w (t) = delta(t) + zeta_w (t),
  quad
  zeta_w (t) = -sigma_w/epsilon_0 e^(-(sigma_w\/epsilon_0) t) u(t) $ <pml-zeta>

wobei $u(t)$ die Sprungfunktion ist. Im Zeitbereich lautet damit etwa die
Gleichung für $H_x$ @Schneider2010[mod. Gl. (11.99)]

$ mu_0 (partial H_x)/(partial t) = -(partial E_z)/(partial y)
  - zeta_y (t) star (partial E_z)/(partial y) $ <pml-zeit>

für die übrigen Komponenten gilt dasselbe mit dem jeweils zugehörigen $zeta_w$.

*Rekursive Auswertung der Faltung.* Eine Faltung über die gesamte Vorgeschichte würde den in @sec-yee gewonnenen Vorteil des expliziten Verfahrens zunichtemachen. Der Kern in @pml-zeta ist jedoch eine reine
Exponentialfunktion, und für einen solchen Kern lässt sich die Faltung rekursiv
auswerten. Dazu wird die Ableitung über je einen Zeitschritt als konstant
angesehen, sodass jedes Teilintegral derselbe Faktor ist, multipliziert
mit der $i$-ten Potenz einer Konstanten @Schneider2010[mod. Gl. (11.125) bis
Gl. (11.127)]:

$ integral_(i Delta t)^((i+1) Delta t) zeta_w (tau) d tau = C_w (b_w)^i,
  quad
  b_w = e^(-(sigma_w Delta t)\/epsilon_0),
  quad
  C_w = b_w - 1 $ <pml-bc>

Für jede
der vier in @pml-tmz auftretenden Ortsableitungen wird deshalb eine Hilfsgröße
$Psi$ mitgeführt, die die gesamte Vorgeschichte in einer einzigen Zahl je
Gitterpunkt speichert.

*Modifizierte Update-Gleichungen.* An jeder Stelle, an der in @update-hx bis
@update-ez ein
Differenzenquotient steht, tritt nun dieser Quotient _zuzüglich_ der zugehörigen
Hilfsgröße auf. Für $H_x$ tritt an die Stelle
von @update-hx das Paar

$ Psi_(H_x) bar.v_(i\,j+1/2)^n = b_y bar.v_(j+1/2) thin
    Psi_(H_x) bar.v_(i\,j+1/2)^(n-1)
  + C_y bar.v_(j+1/2) thin
    (E_z bar.v_(i\,j+1)^n - E_z bar.v_(i\,j)^n)/Delta $ <pml-psi-hx>

$ H_x bar.v_(i\,j+1/2)^(n+1/2) = H_x bar.v_(i\,j+1/2)^(n-1/2)
  - (Delta t)/mu_0 [(E_z bar.v_(i\,j+1)^n - E_z bar.v_(i\,j)^n)/Delta
    + Psi_(H_x) bar.v_(i\,j+1/2)^n] $ <pml-update-hx>

Für $H_y$ gilt dasselbe mit $b_x$ und $C_x$. Für die elektrische Komponente sind zwei Hilfsgrößen nötig, weil in @pml-tmz
zwei Ableitungen mit unterschiedlichen Streckfaktoren auftreten:

$ Psi_(E\,x) bar.v_(i\,j)^(n+1/2) &= b_x bar.v_i thin
    Psi_(E\,x) bar.v_(i\,j)^(n-1/2)
  + C_x bar.v_i thin
    (H_y bar.v_(i+1/2\,j)^(n+1/2) - H_y bar.v_(i-1/2\,j)^(n+1/2))/Delta \
  Psi_(E\,y) bar.v_(i\,j)^(n+1/2) &= b_y bar.v_j thin
    Psi_(E\,y) bar.v_(i\,j)^(n-1/2)
  + C_y bar.v_j thin
    (H_x bar.v_(i\,j+1/2)^(n+1/2) - H_x bar.v_(i\,j-1/2)^(n+1/2))/Delta $ <pml-psi-ez>

Diese treten in @update-ez als zusätzlicher Summand
$Delta (Psi_(E\,x) - Psi_(E\,y))$ innerhalb der Klammer von $C_b$ auf.

Im Inneren des Rechengebiets ist
$sigma_w = 0$ und damit $b_w = 1$ und $C_w = 0$, sodass die Gleichungen dort exakt in @update-hx bis @update-ez übergehen. Die Bezeichnung _perfectly matched_ darf allerdings nicht darüber
hinwegtäuschen, dass die Anpassung nur im Kontinuum exakt ist. Auf dem
diskreten Gitter kommt es gegenüber dem Kontinuum zu systematischen Abweichungen, die im Rahmen der Validierung in Kapitel 5 
erfasst werden.