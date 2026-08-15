// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Theoretische Grundlagen

Dieses Kapitel stellt die Gleichungen bereit, aus denen die späteren Kapitel
hervorgehen. Es beginnt mit den Maxwell-Gleichungen, ergänzt sie um die
Materialgleichungen und die Bedingungen an Materialgrenzen, leitet daraus die
Wellengleichung her und gibt schließlich die ebene Welle als deren einfachste
Lösung an. Sämtliche Aussagen gelten dabei zunächst im dreidimensionalen Raum;
die Beschränkung auf zwei Dimensionen, mit der die vorliegende Arbeit rechnet,
wird erst in Abschnitt 3.4 eingeführt und begründet.

== Die Maxwell-Gleichungen

James Clerk Maxwell hat 1865 gezeigt, dass sich elektrische und magnetische
Erscheinungen durch ein einziges Gleichungssystem beschreiben lassen
@Maxwell1865. In differentieller Form lauten die vier Gleichungen

$ nabla dot bold(D) = rho_"frei" $ <gauss-e>

$ nabla dot bold(B) = 0 $ <gauss-b>

$ nabla times bold(E) = - (partial bold(B))/(partial t) $ <faraday>

$ nabla times bold(H) = bold(J)_"frei" + (partial bold(D))/(partial t) $ <ampere>

Darin bezeichnen $bold(E)$ die elektrische Feldstärke, $bold(D)$ die elektrische
Flussdichte, $bold(H)$ die magnetische Feldstärke, $bold(B)$ die magnetische
Flussdichte, $rho_"frei"$ die Dichte der freien Ladungen und $bold(J)_"frei"$
die der freien Ströme. Die beiden Rechenoperationen sind die _Divergenz_
$nabla dot$, die angibt, wie stark ein Feld an einer Stelle entspringt oder
mündet, und die _Rotation_ $nabla times$, die angibt, wie stark es dort
umläuft. Jede der vier Gleichungen besitzt eine anschauliche Bedeutung
@Fleisch2008:

@gauss-e ist das _Gaußsche Gesetz_. Es besagt, dass elektrische Ladungen
Quellen des elektrischen Flusses sind: Wo eine positive Ladung sitzt, tritt
Fluss aus, wo eine negative sitzt, tritt er ein. Feldlinien beginnen und enden
also an Ladungen.

@gauss-b ist das entsprechende Gesetz für den Magnetismus, und dass rechts null
steht, ist die eigentliche Aussage: Es gibt keine magnetischen Ladungen.
Magnetische Feldlinien haben weder Anfang noch Ende, sondern sind stets
geschlossen.

@faraday ist das _Faradaysche Induktionsgesetz_. Ein sich zeitlich änderndes
Magnetfeld erzeugt ein umlaufendes elektrisches Feld --- der Vorgang, auf dem
jeder Generator beruht. Das Minuszeichen drückt aus, dass das entstehende Feld
der Änderung entgegenwirkt.

@ampere ist das _Ampère-Maxwell-Gesetz_. Ein umlaufendes Magnetfeld entsteht
nicht nur um einen Strom, sondern ebenso um ein sich änderndes elektrisches
Feld. Dieser zweite Anteil $partial bold(D) \/ partial t$ heißt
_Verschiebungsstrom_ und ist Maxwells eigentliche Ergänzung. Erst er schließt
den Kreis: Ein veränderliches $bold(E)$ erzeugt ein $bold(H)$, dieses wiederum
nach @faraday ein $bold(E)$, und so kann sich eine Störung ohne Ladungen und
ohne Ströme durch den leeren Raum fortpflanzen. Genau diese Kopplung von
@faraday und @ampere ist der Gegenstand der gesamten Arbeit; die beiden
Divergenzgleichungen treten in der Rechnung nicht mehr auf.

== Materialgleichungen

Die vier Gleichungen verknüpfen sechs Feldgrößen und sind damit unterbestimmt.
Was fehlt, ist die Angabe, wie ein Material auf ein Feld antwortet. Für die hier
betrachteten Medien gilt

$ bold(D) = epsilon bold(E), quad bold(B) = mu bold(H), quad
  bold(J) = sigma bold(E) $ <material>

mit der Permittivität $epsilon = epsilon_0 epsilon_r$, der Permeabilität
$mu = mu_0 mu_r$ und der elektrischen Leitfähigkeit $sigma$. Die
Naturkonstanten $epsilon_0$ und $mu_0$ beschreiben das Vakuum, die
dimensionslosen Faktoren $epsilon_r$ und $mu_r$ geben an, um welchen Faktor ein
Material die jeweilige Flussdichte gegenüber dem Vakuum verändert. Der letzte
Ausdruck ist das Ohmsche Gesetz in lokaler Form und beschreibt, dass in einem
leitfähigen Material ein Feld einen Strom antreibt.

Dass @material in dieser einfachen Gestalt gilt, ist eine Annahme und keine
Selbstverständlichkeit. Sie setzt voraus, dass $epsilon$, $mu$ und $sigma$ nicht
von der Feldstärke abhängen, dass sie Zahlen und nicht richtungsabhängige
Größen sind und dass sie nicht von der Frequenz abhängen. Die drei
Voraussetzungen werden in Abschnitt 3.4 als Modellannahmen ausdrücklich
festgehalten.

== Randbedingungen an Grenzflächen

An der Grenze zweier Medien ändern sich $epsilon$, $mu$ und $sigma$ sprunghaft.
Die Maxwell-Gleichungen enthalten Ableitungen und sind an einer solchen Stelle
nicht unmittelbar anwendbar; stattdessen folgen aus ihnen Bedingungen, die die
Felder auf beiden Seiten verknüpfen. Für die vorliegende Arbeit sind zwei davon
maßgeblich.

Erstens sind die zur Grenzfläche _parallelen_ Anteile von $bold(E)$ und
$bold(H)$ stetig. Bezeichnet man sie mit dem Index $parallel$, so gilt beim
Übergang von Medium 1 nach Medium 2

$ bold(E)_(1,parallel) = bold(E)_(2,parallel), quad
  bold(H)_(1,parallel) = bold(H)_(2,parallel) $ <stetig>

sofern an der Grenze kein Flächenstrom fließt. Anschaulich folgt das aus
@faraday und @ampere: Ein Umlaufintegral um einen flachen Weg, der die
Grenzfläche umschließt, umfasst im Grenzfall verschwindender Dicke keine Fläche
mehr und damit auch keinen Fluss, sodass die parallelen Anteile auf beiden
Seiten übereinstimmen müssen. Auf @stetig beruht die Berechnung der Reflexion an
einer Grenzfläche in Abschnitt 5.3.

Zweitens folgt daraus der Sonderfall eines _ideal leitenden_ Materials, in der
englischen Bezeichnung _perfect electric conductor_ und abgekürzt PEC. In einem
Leiter mit $sigma arrow.r infinity$ würde jedes Feld augenblicklich einen
unendlich großen Strom antreiben; im Inneren muss $bold(E)$ daher verschwinden.
Mit @stetig ergibt sich für die Grenzfläche

$ bold(E)_parallel = 0 $

Ein ideal leitender Rand erzwingt also, dass das elektrische Feld tangential zu
ihm null ist. Diese Bedingung bestimmt die Eigenfrequenzen des Hohlraums in
Abschnitt 5.2 und wird im Programm dadurch umgesetzt, dass die betreffenden
Feldwerte in jedem Zeitschritt auf null gesetzt werden.

== Die elektromagnetische Wellengleichung

Aus den gekoppelten Gleichungen @faraday und @ampere folgt eine einzelne
Gleichung für $bold(E)$ allein. Betrachtet wird dazu ein Gebiet ohne freie
Ladungen und ohne freie Ströme, also $rho_"frei" = 0$ und
$bold(J)_"frei" = 0$, mit räumlich konstantem $epsilon$ und $mu$.

Der erste Schritt besteht darin, auf beiden Seiten von @faraday die Rotation zu
bilden:

$ nabla times (nabla times bold(E)) = - partial/(partial t) (nabla times bold(B)) $

Auf der linken Seite lässt sich die Vektoridentität

$ nabla times (nabla times bold(E)) = nabla (nabla dot bold(E)) - nabla^2 bold(E) $

anwenden. Wegen $rho_"frei" = 0$ und konstantem $epsilon$ folgt aus @gauss-e
zunächst $nabla dot bold(E) = 0$, sodass der erste Term entfällt. Auf der
rechten Seite wird $bold(B) = mu bold(H)$ eingesetzt und anschließend @ampere
verwendet, das im quellfreien Fall zu
$nabla times bold(H) = epsilon thin partial bold(E) \/ partial t$ wird. Damit
bleibt

$ - nabla^2 bold(E) = - mu epsilon (partial^2 bold(E))/(partial t^2) $

und nach Umstellen die _Wellengleichung_

$ nabla^2 bold(E) - mu epsilon (partial^2 bold(E))/(partial t^2) = 0 $ <welle>

Eine völlig gleichartige Rechnung, die bei @ampere statt bei @faraday ansetzt,
liefert dieselbe Gleichung für $bold(H)$. Beide Felder gehorchen also derselben
Wellengleichung, was ihre gemeinsame Ausbreitung beschreibt.

@welle ist bemerkenswert, weil in ihr keine Ladungen und keine Ströme mehr
vorkommen. Eine einmal erzeugte Störung des Feldes bewegt sich fort, ohne dass
Materie daran beteiligt wäre --- und der Vorfaktor
$mu epsilon$ legt fest, wie schnell. Für das Vakuum ergibt sich

$ c = 1/sqrt(mu_0 epsilon_0) approx 2,998 dot 10^8 " m/s" $

also die Lichtgeschwindigkeit. Aus diesem Zusammentreffen schloss Maxwell, dass
Licht selbst eine elektromagnetische Welle ist.

== Die ebene Welle als analytische Referenzlösung

Die einfachste Lösung von @welle ist die ebene Welle. Sie dient in Kapitel 5
als Vergleichsgröße und wird deshalb hier vollständig angegeben. Der Ansatz
lautet

$ bold(E)(bold(r), t) = bold(E)_0 thin e^(i(omega t - bold(k) dot bold(r))) $ <ansatz>

mit der Kreisfrequenz $omega = 2 pi f$ und dem _Wellenvektor_ $bold(k)$, dessen
Richtung die Ausbreitungsrichtung angibt und dessen Betrag
$beta = abs(bold(k))$ als _Phasenkonstante_ bezeichnet wird. Die Schreibweise
mit der komplexen Exponentialfunktion ist eine Rechenerleichterung; gemeint ist
jeweils der Realteil.

Setzt man @ansatz in @welle ein, so ergibt jede Ableitung nach dem Ort einen
Faktor $-i bold(k)$ und jede nach der Zeit einen Faktor $i omega$. Es folgt

$ -beta^2 + mu epsilon thin omega^2 = 0
  quad ==> quad
  omega = beta / sqrt(mu epsilon) = c beta $

Diese Beziehung zwischen Frequenz und Phasenkonstante heißt
_Dispersionsrelation_. Aus ihr folgen unmittelbar die drei Größen, mit denen
Kapitel 5 arbeitet. Die _Phasengeschwindigkeit_, also die Geschwindigkeit, mit
der ein Punkt fester Phase fortschreitet, beträgt

$ c_p = omega/beta = c = 1/sqrt(mu_r epsilon_r) dot c_0 $

und ist damit unabhängig von der Frequenz: Im Kontinuum laufen alle Frequenzen
gleich schnell. Die _Wellenlänge_ ergibt sich zu

$ lambda = (2 pi)/beta = c/f $

Schließlich stehen die Beträge der beiden Felder in einem festen Verhältnis.
Setzt man @ansatz in @faraday ein, so folgt
$bold(k) times bold(E)_0 = omega mu bold(H)_0$ und daraus

$ abs(bold(E)_0) / abs(bold(H)_0) = (omega mu)/beta = sqrt(mu/epsilon) = Z $

Diese Größe heißt _Wellenwiderstand_; für das Vakuum beträgt sie
$Z_0 = sqrt(mu_0 \/ epsilon_0) approx 377 thin Omega$. Sie wird in
Abschnitt 5.3 benötigt, um aus den Stetigkeitsbedingungen @stetig die
Reflexions- und Transmissionskoeffizienten einer Grenzfläche zu bestimmen.

Damit sind alle analytischen Vergleichsgrößen der Arbeit bereitgestellt: die
Ausbreitungsgeschwindigkeit, die Wellenlänge, das Verhältnis der Feldstärken
sowie die Bedingungen an Grenzflächen und an ideal leitenden Rändern.

== Grenzen analytischer Lösungsverfahren

#todo[
  Argumentieren: Geschlossene Lösungen existieren nur bei hoher
  geometrischer Symmetrie (unendlicher homogener Raum, ebene Grenzfläche,
  Zylinder, Kugel, Hohlleiter mit konstantem Querschnitt). Sobald mehrere
  Objekte mit unterschiedlichen Materialeigenschaften und unregelmäßiger
  Form auftreten, ist keine geschlossene Lösung mehr angebbar, woraus die
  Notwendigkeit eines numerischen Verfahrens folgt --- Überleitung zu
  Kapitel 3.
]
