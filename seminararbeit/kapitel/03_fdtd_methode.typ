// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Die FDTD-Methode
Ziel dieses Kapitels ist es, die zuvor beleuchteten Maxwell-Gleichungen in ein
System diskreter Gleichungen zu überführen. Dabei werden schrittweise drei
Update-Gleichungen hergeleitet, die zu jedem Zeitpunkt an einem jeden Ort die
Feldverteilung liefern. Diese bilden dann das Fundament für die tatsächliche
maschinelle Implementierung in Kapitel 4.

== Grundidee: Finite-Differenzen-Approximation <sec-finite-differenzen>

Da ein Rechner den Grenzübergang $Delta x arrow.r 0$ nicht berechnen kann,
müssen die in @tmz-hx bis @tmz-ez auftretenden Ableitungen durch Quotienten
endlicher --- _finiter_ --- Differenzen umschrieben werden. Die Genauigkeit
dieser Methode zeigt sich in der Taylor-Entwicklung einer hinreichend oft stetig
differenzierbaren Funktion $u$ um die Stelle $x$
@Taflove2005[eq(2.10a) und eq(2.10b)]:

$ u(x plus.minus Delta x) = u(x) plus.minus Delta x (d u)/(d x)
  + (Delta x)^2/2 (d^2 u)/(d x^2)
  plus.minus (Delta x)^3/6 (d^3 u)/(d x^3) + ... $ <taylor>

Anschließend werden die beiden in @taylor zusammengefassten Reihen
voneinander subtrahiert: Dabei heben sich sämtliche Glieder gerader Ordnung
auf, und Division durch $2 Delta x$ liefert die _zentrale Differenz_
@Taflove2005[analog zu eq(2.12)]

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
symmetrisch um $x$ gelegtes Intervall der Breite $Delta x$ und nicht, wie bei
einer einseitigen Differenz, über ein zu einer Seite verschobenes. Der
verbleibende _Abbruchfehler_ ist von zweiter Ordnung --- eine Halbierung von
$Delta x$ geht somit mit einer Viertelung des Fehlers einher.

== Das Yee-Gitter (2D) <sec-yee>

Die zentrale Differenz @zentraldifferenz-halb lässt sich auf @tmz-hx bis
@tmz-ez nicht unmittelbar anwenden, weil diese drei Gleichungen miteinander
gekoppelt sind: Die Zeitableitung von $H_x$ und $H_y$ hängt nach @tmz-hx und
@tmz-hy von der Ortsableitung von $E_z$ ab, die Zeitableitung von $E_z$ nach
@tmz-ez umgekehrt von den Ortsableitungen von $H_x$ und $H_y$. Keine der drei
Komponenten lässt sich damit als alleinstehende Größe lösen.

An dieser Stelle wird das Rechengebiet in ein regelmäßiges Gitter
zerlegt: Ort und Zeit werden an den äquidistanten Stützstellen
$x_i = i Delta x$, $y_j = j Delta y$ und $t_n = n Delta t$ abgetastet, wobei
$Delta x$ und $Delta y$ die Kantenlängen einer Gitterzelle und $Delta t$ einen
Zeitschritt beschreiben. Ausgehend von dieser Diskretisierung löste Yee die
oben beschriebene Schwierigkeit mit einer einzigen Festlegung darüber, _wo_
innerhalb dieses Gitters die einzelnen Feldkomponenten abgelegt werden
@Yee1966: Sie liegen nicht am selben Ort, sondern räumlich um
eine halbe Zelle gegeneinander _versetzt_ (englisch _staggered_) --- und zwar
so, dass jede Komponente genau dort liegt, wo die Ableitung der jeweils anderen
Feldgröße benötigt wird. Für den in @sec-zweidimensional festgelegten
TM#sub[z]-Fall bedeutet das @Yee1966[S.303]:

/ $E_z$: liegt auf den ganzzahligen Gitterpunkten $(i Delta x, j Delta y)$;
/ $H_x$: liegt um eine halbe Zelle in $y$-Richtung versetzt bei $(i Delta x, (j+1\/2) Delta y)$;
/ $H_y$: liegt um eine halbe Zelle in $x$-Richtung versetzt bei $((i+1\/2) Delta x, j Delta y)$.

#figplaceholder(caption: "Yee-Gitter im TM_z-Fall: E_z auf den Gitterpunkten, H_x und H_y jeweils um eine halbe Zelle versetzt. Gestrichelt eine Zelle, deren drei Komponenten im Programm dieselben Indizes tragen.")

Ein Blick auf die einzelnen Gleichungen zeigt den Sinn hinter dieser Anordnung.
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
$(n+1\/2) Delta t$ @Taflove2005[S.60f.]. Beide Felder werden dann abwechselnd
fortgeschrieben --- ein Vorgehen, das nach dem Bocksprungspiel _leapfrog_
benannt ist: Aus dem bekannten $E_z^n$ folgt über @tmz-hx und @tmz-hy das
magnetische Feld zum Zeitpunkt $n+1\/2$, aus diesem über @tmz-ez das elektrische
Feld zum Zeitpunkt $n+1$, und so fort. Auf der rechten Seite jeder
Aktualisierungsvorschrift stehen damit ausschließlich bereits berechnete Werte;
das Verfahren ist _explizit_, es muss zu keinem Zeitpunkt ein Gleichungssystem
gelöst werden --- der wesentliche Grund für den geringen Rechenaufwand je
Zeitschritt.

Entscheidend ist dabei, dass auch die Zeitachse versetzt abgetastet wird: Die
Differenz zweier um $Delta t$ auseinanderliegender Werte ist der halbzahligen
Zwischenzeit zugeordnet, an der die jeweils andere Feldgröße gerade vorliegt.
Damit ist auch die Zeitableitung eine zentrale Differenz nach
@zentraldifferenz-halb und nicht etwa eine einseitige, und das Verfahren ist in
Ort _und_ Zeit von zweiter Ordnung, also $O((Delta x)^2, (Delta y)^2,
(Delta t)^2)$.

== Herleitung der Update-Gleichungen (2D, $"TM"_z$-Fall) <sec-update>

Unter Verwendung des Yee-Gitters aus @sec-yee und Bildung zentraler Differenzen
nach @zentraldifferenz-halb ergeben sich aus @tmz-hx bis @tmz-ez unmittelbar
drei explizite Aktualisierungsvorschriften, im Folgenden _Update-Gleichungen_
genannt: Für ein _quadratisches_ Gitter mit
$Delta x = Delta y = Delta$ ergeben sich folgende Gleichungen @Taflove2005[mod. eq(3.41)]:

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
einmalig berechneten Koeffizienten @Taflove2005[eq(3.31a) und eq(3.31b)]

$ C_a bar.v_(i\,j) = (1 - (sigma Delta t)/(2 epsilon))/(1 + (sigma Delta t)/(2 epsilon)),
  quad
  C_b bar.v_(i\,j) = ((Delta t)/(epsilon Delta))/(1 + (sigma Delta t)/(2 epsilon)) $ <koeffizienten>

wobei $epsilon = epsilon_0 epsilon_r bar.v_(i\,j)$ und $sigma = sigma bar.v_(i\,j)$
die am Ort des betreffenden $E_z$-Wertes gültigen Materialgrößen nach @material
bezeichnen.

Die Form von @koeffizienten geht auf zwei Festlegungen zurück, die in @tmz-ez
noch offen sind. So wird die dort auftretende Stromdichte aufgeteilt in
einen von außen eingeprägten Quellstrom $J_z$ und den vom Feld selbst
getriebenen Leitungsstrom nach dem Ohmschen Gesetz @material:

$ J = J_z + sigma E_z $ <stromaufteilung>

Damit tritt $E_z$ in @tmz-ez zweimal auf --- einmal unter der Zeitableitung,
einmal im Verlustterm. Die Zeitableitung ist nach @zentraldifferenz-halb dem
Zeitpunkt $n+1\/2$ zugeordnet; da eine Gleichung nur dann erfüllt sein kann,
wenn all ihre Glieder für _denselben_ Zeitpunkt gelten, muss auch der
Verlustterm dort ausgewertet werden. Genau dort liegt jedoch kein gespeicherter
Wert von $E_z$ vor, weil das elektrische Feld nach @sec-yee nur zu ganzzahligen
Zeitpunkten geführt wird. Der Verlustterm wird deshalb durch das arithmetische
Mittel der beiden umgebenden Werte ersetzt (_semi-implizites_ Vorgehen):

$ E_z bar.v_(i\,j)^(n+1/2) approx
  (E_z bar.v_(i\,j)^(n+1) + E_z bar.v_(i\,j)^n)/2 $ <semiimplizit>

Setzt man @stromaufteilung und @semiimplizit in @tmz-ez ein, zerfällt der
Verlustterm in zwei Teile: Der eine enthält den gesuchten Wert
$E_z bar.v_(i\,j)^(n+1)$, der andere den bereits bekannten Wert
$E_z bar.v_(i\,j)^n$. Bringt man den ersten auf die linke Seite und lässt den
zweiten rechts stehen, so folgt

$ (1 + (sigma Delta t)/(2 epsilon)) E_z bar.v_(i\,j)^(n+1)
  &= (1 - (sigma Delta t)/(2 epsilon)) E_z bar.v_(i\,j)^n \
  & quad + (Delta t)/(epsilon Delta) (H_y bar.v_(i+1/2\,j)^(n+1/2)
    - H_y bar.v_(i-1/2\,j)^(n+1/2)
    + H_x bar.v_(i\,j-1/2)^(n+1/2)
    - H_x bar.v_(i\,j+1/2)^(n+1/2)
    - J_z bar.v_(i\,j)^(n+1/2) Delta) $ <sortiert>

Division durch die Klammer der linken Seite isoliert $E_z bar.v_(i\,j)^(n+1)$
und liefert @update-ez: Aus dem Verhältnis der beiden runden Klammern wird
$C_a$, aus dem Vorfaktor $Delta t\/(epsilon Delta)$ wird $C_b$. Der Faktor $2$
in beiden Koeffizienten stammt damit allein aus der Mittelung @semiimplizit.

Für verlustfreie Gebiete ($sigma = 0$) entfallen beide Hälften des
Verlustterms, $C_a = 1$ und $C_b = Delta t\/(epsilon Delta)$, und @update-ez
geht in die ungedämpfte Form über.

Mit @update-hx bis @update-ez ist die Feldentwicklung vollständig festgelegt,
sobald ein Anfangszustand vorliegt. Dieser ist in allen folgenden Rechnungen das
feldfreie Gebiet, $E_z = H_x = H_y = 0$ zum Zeitpunkt $n = 0$; die eingeprägte
Stromdichte $J_z$ in @update-ez ist damit die einzige Größe, über die Energie in
das Rechengebiet gelangt. Ihre konkrete zeitliche und räumliche Gestalt hängt
vom jeweiligen Versuch ab und wird deshalb erst in Kapitel 4 festgelegt.

== Stabilität: die Courant-Friedrichs-Lewy-Bedingung (CFL) <sec-cfl>

Damit die Simulation numerisch stabil bleibt, ist die Wahl des Zeitschrittes nicht frei, sondern muss einem Kriterium folgen.
Der Grund liegt in der Nächste-Nachbar-Kopplung des Yee-Gitters: Ein Wert an
einem Gitterpunkt beeinflusst je Zeitschritt nur seine unmittelbaren Nachbarn,
sodass sich eine Information im Gitter höchstens um eine Zelle je Schritt
ausbreitet. Physikalisch läuft sie in derselben Zeit die Strecke $c Delta t$.
Stabil kann das Verfahren daher nur sein, wenn die numerische Ausbreitung der
physikalischen mindestens folgen kann. Durch eine von-Neumann-Analyse der Update-Gleichungen ergibt sich für den _2D_-Fall die _Courant-Lewy-Bedingung (CFL)_, die die Wahl des Zeitschrittes in Abhängigkeit von der Gitterbreite einschränkt @Taflove2005:

$ Delta t <= 1/(c sqrt(1\/(Delta x)^2 + 1\/(Delta y)^2))
  quad arrow.r.double quad
  Delta t <= Delta/(c sqrt(2)) $ <cfl>

wobei rechts wieder das quadratische Gitter aus @sec-update vorausgesetzt ist.
Dass darin die Vakuumlichtgeschwindigkeit $c$ steht, ist keine Einschränkung:
In jedem Medium mit $epsilon_r >= 1$ ist die Ausbreitungsgeschwindigkeit
kleiner, sodass @cfl stets die schärfere Forderung stellt.

Wie nahe die Rechnung an dieser Grenze geführt wird, gibt die _Courant-Zahl_
$S = sqrt(2) c Delta t \/ Delta$ an, die hier so normiert ist, dass die
Stabilitätsgrenze @cfl gerade bei $S = 1$ liegt. In den folgenden Kapiteln wird der Zeitschritt stets durch diese bestimmt. Welchen Einfluss die _Courant-Zahl_ auf die Simulationsergebnisse hat wird in Kapitel 5 näher untersucht. 

== Randbedingungen <sec-rand>

=== Ideal leitender Rand (PEC) <sec-pec>

Für den ideal leitenden Rand gilt nach @sec-randbedingungen
$bold(E)_parallel = 0$. Im TM#sub[z]-Fall ist $E_z$ die einzige elektrische
Komponente und liegt für jede in der $x y$-Ebene verlaufende Wand stets
_parallel_ zu ihr; die Bedingung verkürzt sich damit auf $E_z = 0$. Sie ist
keine Fortschreibungsvorschrift, sondern wird dem Gitter aufgeprägt: Bezeichnet
$M$ die Menge der als leitend markierten Gitterpunkte, so werden deren
$E_z$-Werte nach jeder Auswertung von @update-ez wieder auf null gesetzt
@Schneider2010[S. 171f.]:

$ E_z bar.v_(i\,j)^(n+1) = 0 quad "für alle" (i,j) in M $ <pec>

Die magnetischen Komponenten werden unverändert nach @update-hx und @update-hy
fortgeschrieben. 

=== Absorbierender Rand (PML) <sec-pml>

Während @sec-pec beschreibt, wie man bewusst eine Totalreflexion herbeiruft, findet eine solche an den Rändern des Simulationsgebietes ebenfalls statt. In vielen Fällen ist dieser Effekt jedoch unerwünscht, da so die Energie vollständig in den Simulationsraum zurückgeworfen wird, statt sich weiter im unendlichen Raum auszubreiten.
Die nahestehendeste Lösung wäre es deshalb den eigentlichen Rändern eine Schicht vorzulagern, welche die auslaufende Welle dämpft, sodass keine Energie zurück in das eigentliche Simulationsgebiet geworfen wird.

#figplaceholder(caption: [Anordnung der absorbierenden Schicht: Sie umgibt das eigentliche Rechengebiet auf allen vier Seiten und ist außen weiterhin durch einen PEC-Rand abgeschlossen. Rechts der Verlauf der Dämpfungsgröße $sigma_x$ quer zur Schicht; in den Ecken überlappen die Schichten beider Richtungen.])

Eine gewöhnliche verlustbehaftete Schicht leistet das jedoch nicht: Sie dämpft in alle
Richtungen gleich, während eine reflexionsfreie Anpassung verlangt, dass
ausschließlich die Ausbreitung _senkrecht_ zur Grenzfläche gedämpft wird ---
weshalb sie nur bei senkrechtem Einfall reflexionsfrei bleibt und bei jedem
schrägen Einfall zurückwirft @Schneider2010[Abschn. 11.2 bis 11.4].

*Komplexe Koordinatenstreckung.* Der Ausweg besteht darin, nicht das Material
zu verändern, sondern die Koordinate. Bérenger erreichte dies 1994 zunächst,
indem er jede Feldkomponente in zwei künstliche Anteile aufspaltete
@Berenger1994; die hier verwendete, ohne Aufspaltung auskommende Formulierung
geht auf Chew und Weedon zurück und heißt _komplexe Koordinatenstreckung_,
weil die dabei eingeführten Faktoren die einzelnen Koordinatenrichtungen
skalieren @Schneider2010[Abschn. 11.1 und 11.5]. Zu jeder Raumrichtung $w in {x, y}$ wird dabei
ein komplexer Streckfaktor eingeführt @Schneider2010[eq(11.84) und eq(11.86)]:

$ S_w = 1 + sigma_w/(upright(i) omega epsilon_0) $ <pml-sw>



Der Streckfaktor wird nun jeder Ortsableitung vorangestellt, und zwar derart, dass
$1\/S_w$ stets an die Ableitung in _derselben_ Richtung $w$ gebunden ist. Dafür
wird der Nabla-Operator erweitert, hier auf zwei Dimensionen
eingeschränkt @Schneider2010[mod. eq(11.93)]

$ tilde(nabla) = hat(x) 1/S_x partial/(partial x)
  + hat(y) 1/S_y partial/(partial y) $ <pml-nabla>

unter dessen Anwendung die beiden Rotationsgleichungen in der Frequenzdomäne die Gestalt $upright(i) omega epsilon
bold(E) = tilde(nabla) times bold(H)$ und $-upright(i) omega mu bold(H) =
tilde(nabla) times bold(E)$ annehmen @Schneider2010[eq(11.94) und eq(11.95)].
Ausgeschrieben für den TM#sub[z]-Fall aus @sec-zweidimensional bedeutet das
@Schneider2010[mod. eq(11.90) bis eq(11.92)]

$ upright(i) omega epsilon E_z = 1/S_x (partial H_y)/(partial x)
    - 1/S_y (partial H_x)/(partial y),
  quad
  upright(i) omega mu_0 H_y = 1/S_x (partial E_z)/(partial x),
  quad
  upright(i) omega mu_0 H_x = -1/S_y (partial E_z)/(partial y) $ <pml-tmz>

Wichtig ist, dass sich diese Streckung richtungsweise vorschreiben lässt.
So wird die Richtungsselektivität gewährleistet. Die Ausbreitung längs der Grenzfläche bleibt damit unverändert, gedämpft
wird allein die dazu senkrechte @Schneider2010[eq(11.72) und eq(11.73)]. Der
Reflexionsfaktor verschwindet dadurch für _jeden_ Einfallswinkel und _jede_
Frequenz. Setzt man @pml-sw in eine unter $theta$ einlaufende ebene Welle nach
@ansatz ein, so tritt an die Stelle von $k_x$ das Produkt $S_x k_x$, und der
Ortsanteil der Welle zerfällt in zwei Faktoren
@Schneider2010[eq(11.82) und eq(11.83)]:

$ e^(-upright(i) S_x k_x x - upright(i) k_y y)
  = underbrace(e^(-(sigma_x k_x)/(omega epsilon_0) x), "reell, klingt ab")
  dot underbrace(e^(-upright(i) (k_x x + k_y y)), "unveränderte Phase") $ <pml-daempfung>

Der erste Faktor ist reell und ändert deshalb allein die Amplitude; der zweite
ist Zeichen für Zeichen derselbe wie im angrenzenden verlustfreien Gebiet. Die
Welle wird also gedämpft, läuft in Richtung und Betrag des Wellenvektors aber
genau so weiter.

*Dicke und Verlauf der Schicht.* Hinter der dämpfenden Schicht bleibt weiterhin der PEC-Rand bestehen. Eine Welle durchläuft die Schicht der Dicke $d$, wird dort jedoch gemäß @pec
 zurückgeworfen und durchläuft sie nochmals; maßgeblich für
die verbleibende Reflexion ist deshalb die auf diesem Hin- und Rückweg
angesammelte Dämpfung. Mit $k_x = (omega\/c) cos theta$ und
$1\/(c epsilon_0) = eta_0$, dem Freiraumwellenwiderstand, folgt aus dem
Dämpfungsfaktor in @pml-daempfung für den Betrag des insgesamt zurückkehrenden
Anteils

$ abs(R(theta)) = exp(-2 eta_0 cos theta integral_0^d sigma_x (x) d x) $ <pml-restreflexion>

Daraus lassen sich drei Aussagen ablesen. Erstens enthält @pml-restreflexion
keine Frequenz --- im Kontinuum ist die Schicht für alle Frequenzen gleich
wirksam. Zweitens geht nur das Integral über $sigma_x$ ein, nicht aber, wie die
Dämpfung über die Dicke verteilt ist. Drittens steht $cos theta$ im Exponenten,
sodass die Dämpfung bei streifendem Einfall ($theta arrow.r 90 degree$)
verschwindet. Der letzte Punkt ist die grundsätzliche Schwäche jeder PML und
wird in Testfall 4 in Kapitel 5 quantitativ geprüft.

Auf dem diskreten Gitter gilt die zweite Aussage nicht mehr. Dort ist die
Schicht in Zellen zerlegt, und ein sprunghafter Anstieg von $sigma_x$ von einer
Zelle zur nächsten wirkt wie eine Materialgrenze und erzeugt eine rein
numerische Reflexion @Schneider2010[Abschn. 11.1 und 11.2]. $sigma_x$ wird
deshalb vom inneren Schichtrand zum Gitterrand hin allmählich gesteigert,
hier nach einem Potenzverlauf

$ sigma_x (rho) = sigma_max rho^m,
  quad rho = a/d in [0, 1] $ <pml-profil>

worin $a$ den Abstand vom inneren Schichtrand bezeichnet. Die beiden freien
Größen wirken gegenläufig: Ein großes $sigma_max$ senkt nach @pml-restreflexion die
physikalische Restreflexion, verschärft aber zugleich den Sprung je Zelle und damit die
numerische. Ihre konkrete Wahl wird deshalb erst in Kapitel 4 getroffen und in
Testfall 4 überprüft.

*Rückkehr in den Zeitbereich.* @pml-tmz steht im Frequenzbereich, das
FDTD-Verfahren rechnet jedoch im Zeitbereich, und einem Produkt im
Frequenzbereich entspricht dort eine Faltung. Benötigt wird also die
Rücktransformierte $overline(S)_w$ von $1\/S_w$; der Querstrich kennzeichnet
dabei die zugehörige Zeitfunktion, und zwar die zu $1\/S_w$ und nicht die zu
$S_w$ selbst, weil in @pml-tmz gerade der Kehrwert vor den Ableitungen steht.
Mit @pml-sw ist

$ 1/S_w = (upright(i) omega epsilon_0)/(sigma_w + upright(i) omega epsilon_0)
  = 1 - 1/(1 + upright(i) omega epsilon_0\/sigma_w) $ <pml-inverses>

Der erste Summand liefert die Delta-Distribution $delta(t)$, die in der Faltung
die ursprüngliche Ableitung unverändert reproduziert; für den zweiten genügt ein
_Korrespondenzpaar_. Darunter versteht man einen Ausdruck im Frequenzbereich
zusammen mit derjenigen Zeitfunktion, aus der die Fourier-Transformation ihn
erzeugt; weil beide eindeutig zusammengehören, lässt sich die eine aus der
anderen ablesen, statt das Transformationsintegral auszurechnen. Das hier
passende Paar lautet $(1 + upright(i) omega tau)^(-1) arrow.l.r tau^(-1)
e^(-t\/tau) u(t)$ @Schneider2010[eq(11.104)] --- einem einfachen Bruch im
Frequenzbereich entspricht also eine abklingende Exponentialfunktion im
Zeitbereich. Mit $tau = epsilon_0\/sigma_w$ ergibt sich zusammen
@Schneider2010[mod. eq(11.107) und eq(11.108)]

$ overline(S)_w (t) = delta(t) + zeta_w (t),
  quad
  zeta_w (t) = -sigma_w/epsilon_0 e^(-(sigma_w\/epsilon_0) t) u(t) $ <pml-zeta>

wobei $u(t)$ die Sprungfunktion ist. Im Zeitbereich lautet damit etwa die
Gleichung für $H_x$ @Schneider2010[mod. eq(11.99)]

$ mu_0 (partial H_x)/(partial t) = -(partial E_z)/(partial y)
  - zeta_y (t) star (partial E_z)/(partial y) $ <pml-zeit>

also die gewohnte Gleichung, ergänzt um genau einen Faltungsterm; für die
übrigen Komponenten gilt dasselbe mit dem jeweils zugehörigen $zeta_w$.

*Rekursive Auswertung der Faltung.* Eine Faltung über die gesamte
Vorgeschichte des Feldes wäre in der Zeitschleife nicht zu bezahlen --- sie
würde je Gitterpunkt und Zeitschritt einen mit $n$ wachsenden Aufwand
verursachen und den in @sec-yee gewonnenen Vorteil des expliziten Verfahrens
zunichtemachen. Der Kern in @pml-zeta ist jedoch eine reine
Exponentialfunktion, und für einen solchen Kern lässt sich die Faltung rekursiv
auswerten. Dazu wird die Ableitung über je einen Zeitschritt als konstant
angesehen --- dieselbe Treppenannäherung, die dem Verfahren ohnehin zugrunde
liegt ---, sodass das Faltungsintegral in eine Summe über Teilintegrale
zerfällt. Jedes dieser Teilintegrale ist wegen
$integral e^(-alpha t) d t prop e^(-alpha t)$ derselbe Faktor, multipliziert
mit der $i$-ten Potenz einer Konstanten @Schneider2010[mod. eq(11.125) bis
eq(11.127)]:

$ integral_(i Delta t)^((i+1) Delta t) zeta_w (tau) d tau = C_w (b_w)^i,
  quad
  b_w = e^(-(sigma_w Delta t)\/epsilon_0),
  quad
  C_w = b_w - 1 $ <pml-bc>

Zieht man in der so entstandenen Summe den Term $i = 0$ heraus, so ist der Rest
bis auf den Faktor $b_w$ genau dieselbe Summe zum vorangegangenen Zeitschritt
@Schneider2010[eq(11.128) bis eq(11.132)]. Der Wert der Faltung hängt damit nur
noch von seinem eigenen Vorgängerwert und der aktuellen Ableitung ab. Für jede
der vier in @pml-tmz auftretenden Ortsableitungen wird deshalb eine Hilfsgröße
$Psi$ mitgeführt, die die gesamte Vorgeschichte in einer einzigen Zahl je
Gitterpunkt speichert.

*Modifizierte Update-Gleichungen.* Damit lassen sich @update-hx bis
@update-ez unmittelbar erweitern: An jeder Stelle, an der dort ein
Differenzenquotient steht, tritt nun dieser Quotient _zuzüglich_ der zugehörigen
Hilfsgröße auf, und die Hilfsgröße wird unmittelbar davor nach @pml-bc
fortgeschrieben. Für die beiden magnetischen Komponenten treten an die Stelle
von @update-hx und @update-hy die Paare

$ Psi_(H_x) bar.v_(i\,j+1/2)^n = b_y bar.v_(j+1/2) thin
    Psi_(H_x) bar.v_(i\,j+1/2)^(n-1)
  + C_y bar.v_(j+1/2) thin
    (E_z bar.v_(i\,j+1)^n - E_z bar.v_(i\,j)^n)/Delta $ <pml-psi-hx>

$ H_x bar.v_(i\,j+1/2)^(n+1/2) = H_x bar.v_(i\,j+1/2)^(n-1/2)
  - (Delta t)/mu_0 [(E_z bar.v_(i\,j+1)^n - E_z bar.v_(i\,j)^n)/Delta
    + Psi_(H_x) bar.v_(i\,j+1/2)^n] $ <pml-update-hx>

$ Psi_(H_y) bar.v_(i+1/2\,j)^n = b_x bar.v_(i+1/2) thin
    Psi_(H_y) bar.v_(i+1/2\,j)^(n-1)
  + C_x bar.v_(i+1/2) thin
    (E_z bar.v_(i+1\,j)^n - E_z bar.v_(i\,j)^n)/Delta $ <pml-psi-hy>

$ H_y bar.v_(i+1/2\,j)^(n+1/2) = H_y bar.v_(i+1/2\,j)^(n-1/2)
  + (Delta t)/mu_0 [(E_z bar.v_(i+1\,j)^n - E_z bar.v_(i\,j)^n)/Delta
    + Psi_(H_y) bar.v_(i+1/2\,j)^n] $ <pml-update-hy>

Für die elektrische Komponente sind zwei Hilfsgrößen nötig, weil in @pml-tmz
zwei Ableitungen mit unterschiedlichen Streckfaktoren auftreten:

$ Psi_(E\,x) bar.v_(i\,j)^(n+1/2) &= b_x bar.v_i thin
    Psi_(E\,x) bar.v_(i\,j)^(n-1/2)
  + C_x bar.v_i thin
    (H_y bar.v_(i+1/2\,j)^(n+1/2) - H_y bar.v_(i-1/2\,j)^(n+1/2))/Delta \
  Psi_(E\,y) bar.v_(i\,j)^(n+1/2) &= b_y bar.v_j thin
    Psi_(E\,y) bar.v_(i\,j)^(n-1/2)
  + C_y bar.v_j thin
    (H_x bar.v_(i\,j+1/2)^(n+1/2) - H_x bar.v_(i\,j-1/2)^(n+1/2))/Delta $ <pml-psi-ez>

$ E_z bar.v_(i\,j)^(n+1) &= C_a bar.v_(i\,j) dot E_z bar.v_(i\,j)^n \
  & quad + C_b bar.v_(i\,j) dot (H_y bar.v_(i+1/2\,j)^(n+1/2)
    - H_y bar.v_(i-1/2\,j)^(n+1/2)
    + H_x bar.v_(i\,j-1/2)^(n+1/2)
    - H_x bar.v_(i\,j+1/2)^(n+1/2) \
  & quad quad quad quad quad
    + Delta (Psi_(E\,x) bar.v_(i\,j)^(n+1/2)
    - Psi_(E\,y) bar.v_(i\,j)^(n+1/2))
    - J_z bar.v_(i\,j)^(n+1/2) Delta) $ <pml-update-ez>

Die Koeffizienten $C_a$ und $C_b$ sind unverändert die aus @koeffizienten; der
Faktor $Delta$ vor der Klammer mit den Hilfsgrößen steht dort, weil $C_b$ nach
@koeffizienten bereits $1\/Delta$ enthält, die Klammer also bare Differenzen
und keine Differenzenquotienten aufnimmt. Die Reihenfolge innerhalb eines
Zeitschritts ist damit festgelegt: zuerst @pml-psi-hx und @pml-psi-hy, dann
@pml-update-hx und @pml-update-hy, anschließend @pml-psi-ez und zuletzt
@pml-update-ez.

Die Größen $b_w$ und $C_w$ hängen nach @pml-bc allein von $sigma_w$ am Ort der
jeweiligen Komponente ab. Sie sind daher --- wie $C_a$ und $C_b$ --- vor Beginn
der Zeitschleife einmalig zu berechnen, und zwar als eindimensionale Profile,
da $sigma_x$ nur von $x$ und $sigma_y$ nur von $y$ abhängt; für die um eine
halbe Zelle versetzten magnetischen Komponenten sind sie bei $i + 1\/2$
beziehungsweise $j + 1\/2$ auszuwerten. Im Inneren des Rechengebiets ist
$sigma_w = 0$ und damit $b_w = 1$ und $C_w = 0$: Die Hilfsgrößen bleiben, vom
feldfreien Anfangszustand aus, dauerhaft null, und @pml-psi-hx bis
@pml-update-ez gehen exakt in @update-hx bis @update-ez über. Die Schicht wirkt
also ausschließlich dort, wo $sigma_w > 0$ ist, und der Übergang ins Innere ist
ohne Fallunterscheidung stetig. In den vier Ecken des Gebiets überlappen die
Schichten beider Richtungen; dort sind $sigma_x$ und $sigma_y$ gleichzeitig von
null verschieden. Auch das erfordert keine Sonderbehandlung, da @pml-update-ez
ohnehin beide Hilfsgrößen enthält.

Der Preis der Schicht ist gering: je Ortsableitung --- also viermal --- eine
zusätzliche Zahl je Gitterpunkt sowie je Zeitschritt zwei zusätzliche
Multiplikationen und zwei Additionen. Die Bezeichnung _perfectly matched_ darf allerdings nicht darüber
hinwegtäuschen, dass die Anpassung nur im Kontinuum exakt ist. Auf dem
diskreten Gitter bleiben zwei Fehlerquellen: die stufenweise Änderung von
$sigma_w$ nach @pml-profil und die numerische Dispersion, derentwegen die
diskrete Welle der Dispersionsbeziehung des Kontinuums, aus der die Anpassung
folgt, nicht exakt gehorcht. Beide sind analytisch nicht zugänglich; wie groß
die verbleibende Reflexion tatsächlich ausfällt, muss deshalb gemessen werden
--- Testfall 4 in Kapitel 5 tut das.

#note[Vier Punkte für die Überarbeitung.

(i) *Was nicht aus @Schneider2010 stammt.* Eigen sind @pml-restreflexion und
die drei daraus abgelesenen Aussagen: Die Gleichung folgt aus dem Exponenten
von eq(11.83) mit $k_x = (omega\/c) cos theta$ und dem PEC-Rand hinter der
Schicht, steht so aber nicht bei Schneider. Der
Potenzverlauf @pml-profil steht ebenfalls nicht in dessen Kapitel 11 --- dort
wird nur gefordert, $sigma$ „allmählich“ zu steigern. Für @pml-profil und für
$sigma_max$ ist deshalb noch eine belegbare Quelle nachzutragen (Bérenger 1994
oder Taflove Kap. 7); bis dahin ist beides bewusst ohne Beleg formuliert.

Vor allem aber stehen die sechs Update-Gleichungen @pml-psi-hx bis
@pml-update-ez so *nicht* bei Schneider. Sein Kapitel 11 endet mit der
Rekursion eq(11.132) und einer in Worten beschriebenen Reihenfolge
(Abschn. 11.6); indizierte Update-Gleichungen für den TM#sub[z]-Fall gibt er
nicht an. Sie sind hier aus eq(11.132), @pml-bc und den eigenen @update-hx bis
@update-ez zusammengesetzt. Wörtlich übernommen ist allein die Rekursion. Dass
die Hilfsgröße *innerhalb* der $C_b$-Klammer steht, folgt aus eq(11.110), wo
$Psi$ neben der Rotation steht und damit denselben Faktor $Delta t\/epsilon$
trägt; Schneiders Formulierung „add or subtract $Psi$ to the electric field“
(Abschn. 11.6, Schritt 3) lässt das offen. Ein numerischer Vergleich mit
`fdtd_core.py` bestätigt die hier gewählte Einordnung: Beide Rechnungen stimmen
auf Maschinengenauigkeit überein.

(ii) *Warum „mod.“ bei @pml-zeta und @pml-bc.* Schneider führt in eq(11.87) den
allgemeineren Streckfaktor $S_w = kappa_w + sigma_w\/(a_w + upright(i) omega
epsilon_0)$ und schreibt eq(11.107), eq(11.126) und eq(11.127) damit an.
@pml-sw ist dessen Sonderfall $kappa_w = 1$, $a_w = 0$ und zugleich wörtlich
seine eq(11.86); die daraus folgenden @pml-zeta und @pml-bc sind deshalb als
„mod.“ zitiert.

(iii) *Verhältnis zum Programm.* Kapitel 4 ersetzt $b_w$ durch $(1-x)\/(1+x)$
mit $x = sigma_w Delta t\/(2 epsilon_0)$ und setzt $C_w = b_w - 1 =
-2x\/(1+x)$ (`g3` und `-g2`) --- strukturell also exakt @pml-bc mit der
bilinearen Näherung von $e^(-2x)$. Achtung: $x$ ist *nicht* klein. Bei
$Delta = 2$ mm, $S = "0,99"$ und zehn Zellen läuft $x$ von $0$ am inneren bis
$"1,12"$ am äußeren Schichtrand; dort liefert $(1-x)\/(1+x) = -"0,06"$ gegen
$e^(-2x) = "0,11"$, also sogar das falsche Vorzeichen. Bis $x approx "0,5"$
(inneres Drittel der Schicht) stimmen beide auf vier Nachkommastellen überein.
Ein numerischer Vergleich beider Varianten ergab eine Feldabweichung von rund
1 % der Maximalamplitude, bei praktisch gleicher Restreflexion des Randes
(Faktor 1,01) --- die Näherung ist also unschädlich, weil das Feld die
äußersten Zellen kaum noch erreicht. Der Satz dazu gehört in Abschnitt 4.5 und
darf *nicht* mit „Übereinstimmung bis zweiter Ordnung“ begründet werden.

(iv) *Kapitel 5.* Der Satz in Testfall 4, das Modul verbinde „die Formulierung
von Sullivan mit einem Dämpfungsparameter nach Taflove und Hagness“, trifft
nach dieser Herleitung nicht mehr zu und ist auf @Schneider2010 umzustellen.]
