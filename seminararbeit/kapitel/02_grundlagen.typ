// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Theoretische Grundlagen

Dieses Kapitel stellt knapp die Gleichungen bereit, aus denen in Kapitel 3 die
FDTD-Methode abgeleitet und mit denen in Kapitel 5 ihr Ergebnis überprüft wird:
die Maxwell-Gleichungen, die Materialgleichungen, die Reduktion auf zwei
Dimensionen, die Randbedingungen an Grenzflächen sowie die Wellengleichung und
ihre einfachste Lösung, die ebene Welle. Da es sich um Grundlagen und nicht um
den eigentlichen Gegenstand der Arbeit handelt, bleiben die Herleitungen kurz.

== Die Maxwell-Gleichungen <sec-maxwell>

Maxwell zeigte 1865, dass sich alle elektrischen und magnetischen
Erscheinungen durch ein Gleichungssystem beschreiben lassen @Maxwell1865, das
zusammen mit der Lorentz-Kraft den gesamten theoretischen Inhalt der
klassischen Elektrodynamik zusammenfasst @GriffithsDE[S. 417f.]. Auf die freien
Ladungen und Ströme bezogen @GriffithsDE[S. 420] lauten sie in differentieller
Form @GriffithsDE[eq(7.56)]

$ nabla dot bold(D) = rho_"frei" $ <gauss-e>
$ nabla dot bold(B) = 0 $ <gauss-b>
$ nabla times bold(E) = - (partial bold(B))/(partial t) $ <faraday>
$ nabla times bold(H) = bold(J)_"frei" + (partial bold(D))/(partial t) $ <ampere>

mit der elektrischen Feldstärke $bold(E)$, der elektrischen Flussdichte
$bold(D)$, der magnetischen Feldstärke $bold(H)$, der magnetischen Flussdichte
$bold(B)$ sowie der Dichte $rho_"frei"$ und $bold(J)_"frei"$ der freien
Ladungen und Ströme. @gauss-e und @gauss-b sind das Gaußsche Gesetz und sein
magnetisches Gegenstück; da nach allem, was bekannt ist, keine magnetischen
Ladungen existieren @GriffithsDE[S. 419], gilt @gauss-b uneingeschränkt. Beide
enthalten keine Zeitableitung und treten in der weiteren Rechnung nicht mehr
auf. Gegenstand der Arbeit ist allein die Kopplung von @faraday (Faradaysches
Induktionsgesetz: ein sich änderndes $bold(B)$ erzeugt ein umlaufendes
$bold(E)$) und @ampere (Ampère-Maxwell-Gesetz: ein umlaufendes $bold(H)$
entsteht durch einen Strom _oder_ ein sich änderndes $bold(D)$). Der Term
$partial bold(D)\/partial t$ ist Maxwells eigentliche Ergänzung des zuvor
bekannten Ampèreschen Gesetzes @GriffithsDE[S. 417] und heißt
_Verschiebungsstrom_ @GriffithsDE[S. 422]; erst er schließt den Kreis, und bildet die Grundlage für die Aubreitung _elektromagnetischer_ Wellen

== Materialgleichungen <sec-material>

Die vier Gleichungen verknüpfen sechs Feldgrößen und sind unterbestimmt; es
fehlt die Angabe, wie ein Material auf ein Feld reagiert. Diese
_Materialgleichung_ folgt nicht aus den Maxwell-Gleichungen und hängt vom
Material ab; in ihrer einfachsten, für _lineare_ Medien gültigen Form
@GriffithsDE[S. 422] wird hier durchgehend angesetzt

$ bold(D) = epsilon bold(E), quad bold(B) = mu bold(H), quad
  bold(J) = sigma bold(E) $ <material>

mit $epsilon = epsilon_0 epsilon_r$ und $mu = mu_0 mu_r$
sowie der Leitfähigkeit $sigma$ ; der letzte Ausdruck ist das Ohmsche Gesetz in
lokaler Form @GriffithsDE[S. 368].


Dass @material in dieser Form gilt, ist eine Annahme, die vier Vereinfachungen @Taflove2005[S. 52f.] vorrausetzt:

/ Linearität: $epsilon,mu,sigma$ unabhängig von der Feldstärke.
/ Isotropie: $epsilon,mu$ Zahlen statt Tensoren --- $bold(D) parallel bold(E)$.
/ Nichtdispersivität: $epsilon,mu,sigma$ unabhängig von der Frequenz.
/ Zeitinvarianz: $epsilon,mu,sigma$ ändern sich während der Simulation nicht.

Hinzu kommt die Festlegung auf _nichtmagnetische_ Materialien ($mu_r=1$,
$mu=mu_0$ überall); das Programm in Kapitel 4 kennt entsprechend nur eine
ortsabhängige Permittivität $epsilon_r(x,y)$ und Leitfähigkeit $sigma(x,y)$,
innerhalb einer Gitterzelle jedoch konstant (_stückweise homogen_).

== Reduktion auf zwei Dimensionen: der TM#sub[z]-Fall <sec-zweidimensional>

Die Arbeit rechnet nicht im dreidimensionalen Raum, sondern in einer Ebene, da
der Gesamtaufwand (Zellen mal Zeitschritte, Kapitel 3) bei $N$ Zellen je
Richtung in 2D bei $N^3$ liegt statt bei $N^4$ in 3D. Diese Reduktion ist
keine Näherung, sondern _exakt_, sofern weder Felder noch Materialverteilung
von der dritten Koordinate abhängen:

$ partial/(partial z) equiv 0 $ <zinvarianz>

Unter @zinvarianz  zerfallen @faraday und @ampere in zwei vollständig
entkoppelte, unabhängig lösbare Gruppen: den _TM#sub[z]-Fall_ (transversal
magnetisch, Komponenten $E_z,H_x,H_y$) und den _TE#sub[z]-Fall_ (transversal
elektrisch, $E_x,E_y,H_z$) @GriffithsDE[S.516].
// LITERATURLUECKE: Die Bezeichnungen TM_z/TE_z und die Entkopplung kommen bei
// Griffiths nicht vor. Als Beleg der Nomenklatur waeren Schneider2010
// (Abschn. 8.3, S. 185, und 8.7, S. 220) sowie Sullivan2013 (Kap. 3, S. 54 f.)
// geeignet -- Seitenangaben vor Zitat am Buch/Volltext pruefen.
Beide Fälle sind strukturell gleich gebaut; die Arbeit verwendet durchgehend
den TM#sub[z]-Fall, weil er mit nur einer elektrischen Feldkomponente Anregung
und Randbedingung auf eine skalare Größe zurückführt. Mit $mu=mu_0$ und unter Verwendung der in 2.2 angesprochenen Vereinfachungen vereinfachen sich die Gleichungen @Taflove2005[eq(3.13)] zu

$ mu_0 (partial H_x)/(partial t) = - (partial E_z)/(partial y) $ <tmz-hx>
$ mu_0 (partial H_y)/(partial t) = (partial E_z)/(partial x) $ <tmz-hy>
$ epsilon (partial E_z)/(partial t) + J
  = (partial H_y)/(partial x) - (partial H_x)/(partial y) $ <tmz-ez>


== Randbedingungen an Grenzflächen <sec-randbedingungen>

An der Grenze zweier Medien ändern sich $epsilon,mu,sigma$ sprunghaft; aus der
Integralform der Maxwell-Gleichungen folgen dafür _Randbedingungen_
@GriffithsDE[S. 424]. Maßgeblich ist die Stetigkeit der zur Grenzfläche
_parallelen_ Feldanteile: Eine sehr dünne, die Grenzfläche durchstoßende
Ampèresche Schleife liefert für $bold(E)$, da @faraday keinen Quellterm
enthält, uneingeschränkt

$ bold(E)_(1,parallel) = bold(E)_(2,parallel) $ <stetig>

während dieselbe Konstruktion, angewandt auf @ampere, für $bold(H)$ nur dann
Stetigkeit liefert, wenn kein freier Flächenstrom fließt @GriffithsDE[S. 424]
--- was für die in Abschnitt 5.3 betrachteten Dielektrika zutrifft.


Ein wichtiger Sonderfall ist der eines _ideal leitenden_ Randes (_perfect
electric conductor_, PEC): Mit $sigma arrow.r infinity$ würde im Leiter ein
unendlich großer Strom fließen, weshalb im Inneren $bold(E)$ verschwinden muss
@GriffithsDE[S.368]; mit @stetig folgt daraus für die Grenzfläche
$bold(E)_parallel = 0$.
Im TM#sub[z]-Fall verkürzt sich das auf $E_z=0$.

== Die elektromagnetische Wellengleichung <sec-wellengleichung>

Aus @faraday und @ampere folgt für ein quellfreies Gebiet
($rho_"frei"=bold(J)_"frei"=0$) mit konstantem $epsilon,mu$ eine einzelne
Gleichung für $bold(E)$ @GriffithsDE[S.481f.]: Bildet man die Rotation von
@faraday, nutzt die Vektoridentität @GriffithsDE[S.482]

$ nabla times (nabla times bold(E)) = nabla (nabla dot bold(E)) - nabla^2 bold(E) $

und setzt @gauss-e sowie @ampere ein, ergibt sich die _Wellengleichung_
@GriffithsDE[eq.(9.41)]

$ nabla^2 bold(E) - mu epsilon (partial^2 bold(E))/(partial t^2) = 0 $ <welle>

Dieselbe Gleichung folgt für $bold(H)$, wenn man stattdessen bei @ampere
ansetzt. In @welle kommen keine Ladungen oder Ströme mehr vor: Eine einmal
erzeugte Störung breitet sich mit der durch $mu epsilon$ festgelegten
Geschwindigkeit aus, im Vakuum @GriffithsDE[S.482]

$ c = 1/sqrt(mu_0 epsilon_0) approx 2,998 dot 10^8 " m/s" $

also der Lichtgeschwindigkeit --- woraus Maxwell schloss, dass Licht selbst
eine elektromagnetische Welle ist.

== Die ebene Welle als analytische Referenzlösung <sec-ebenewelle>

Die einfachste Lösung von @welle, die in Kapitel 5 als Vergleichsgröße dient,
ist die ebene Welle

$ bold(E)(bold(r), t) = E_0 thin e^(i(omega t - bold(k) dot bold(r))) hat(z) $ <ansatz>

Darin ist $E_0$ die (skalare) Amplitude, $hat(z)$ der Einheitsvektor in
$z$-Richtung --- der Vektor $bold(E)$ hat im TM#sub[z]-Fall ja nur die
Komponente $E_z$ ---, $bold(r) = (x,y,z)$ der (volle, dreidimensionale)
Ortsvektor des betrachteten Punktes, $t$ die Zeit, $omega=2pi f$ die
Kreisfrequenz und $bold(k)$ der Wellenvektor, dessen Richtung die
Ausbreitungsrichtung angibt @GriffithsDE[mod. eq.(1.43)]. Wegen @zinvarianz
muss $bold(k)$ dabei stets in der $x y$-Ebene liegen ($k_z=0$); an $bold(r)$
selbst ändert das nichts, es bleibt der volle Ortsvektor, aber sein
$z$-Anteil fällt im Skalarprodukt $bold(k) dot bold(r) = k_x x + k_y y$
heraus, sodass das Feld --- wie von @zinvarianz gefordert --- tatsächlich
nicht von $z$ abhängt. Da das Feld hier durch Verwendung der Eueleridentität komplex angesetzt ist,
beschreibt der Realteil das physikalische Feld. Eine ebene Welle heißt so,
weil zu jedem festen Zeitpunkt $t$ alle Orte $bold(r)$ gleicher Phase
$omega t - bold(k) dot bold(r) = "const."$ eine zu $bold(k)$ senkrechte Ebene
bilden, auf der Amplitude und Phase überall gleich sind @GriffithsDE[S.483].



== Grenzen analytischer Lösungsverfahren
Die in @sec-ebenewelle gezeigte ebene Welle stellt eine der wenigen exakten
Lösungen der Maxwell-Gleichungen dar. Analytische Lösungsverfahren stoßen in der
Praxis schnell an ihre Grenzen: Geschlossene mathematische Ausdrücke existieren
lediglich für stark symmetrische Spezialfälle --- etwa den unendlichen,
homogenen Raum, ebene Grenzflächen (Fresnelsche Formeln), ideale Zylinder oder
Kugeln (Mie-Streuung) sowie einfache Hohlleitergeometrien.

Sobald ein Problem reale Randbedingungen aufweist --- etwa komplexe,
unregelmäßige Objektgeometrien, mehrschichtige oder ortsabhängige Materialien
($epsilon_r(x,y)$, $sigma(x,y)$) --- ist eine geschlossene analytische
Berechnung nicht mehr möglich. Um die Feldverteilung für solche realistischen
Szenarien dennoch präzise zu bestimmen, muss das kontinuierliche System der
Maxwell-Gleichungen diskretisiert und numerisch gelöst werden. Dies motiviert den
Übergang zur _Finite-Difference Time-Domain_-Methode (FDTD), die im folgenden
Kapitel hergeleitet wird.
