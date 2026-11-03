// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Theoretische Grundlagen

Dieses Kapitel stellt knapp die Gleichungen bereit, aus denen in Kapitel 3 die
FDTD-Methode abgeleitet und mit denen in Kapitel 5 ihr Ergebnis überprüft wird.
Abschließend werden die Grenzen analytischer Lösungsverfahren aufgezeigt. Da es
sich um Grundlagen und nicht um den eigentlichen Gegenstand der Arbeit
handelt, bleiben die Herleitungen kurz.

== Die Maxwell-Gleichungen <sec-maxwell>

Maxwell zeigte 1865, dass sich alle elektrischen und magnetischen
Erscheinungen durch ein Gleichungssystem beschreiben lassen, das
zusammen mit der Lorentz-Kraft den gesamten theoretischen Inhalt der
klassischen Elektrodynamik zusammenfasst @lit:GriffithsDE[S. 417 f.]. Auf die freien
Ladungen und Ströme bezogen @lit:GriffithsDE[S. 420] lauten sie in differentieller
Form @lit:GriffithsDE[Gl. (7.56)]

$ nabla dot bold(D) = rho_"frei" $ <gauss-e>
$ nabla dot bold(B) = 0 $ <gauss-b>
$ nabla times bold(E) = - (partial bold(B))/(partial t) $ <faraday>
$ nabla times bold(H) = bold(J)_"frei" + (partial bold(D))/(partial t) $ <ampere>

mit der elektrischen Feldstärke $bold(E)$, der elektrischen Flussdichte
$bold(D)$, der magnetischen Feldstärke $bold(H)$, der magnetischen Flussdichte
$bold(B)$ sowie der Dichte $rho_"frei"$ und $bold(J)_"frei"$ der freien
Ladungen und Ströme. @gauss-e und @gauss-b enthalten keine Zeitableitung und
werden im Zeitschrittverfahren in Kapitel 3 nicht ausgewertet; @gauss-e wird
lediglich einmal bei der Herleitung der Wellengleichung benötigt. Gegenstand der
Arbeit ist dementsprechend vorrangig die Kopplung von @faraday (Faradaysches
Induktionsgesetz: ein sich änderndes $bold(B)$ erzeugt ein umlaufendes
$bold(E)$) und @ampere (Ampère-Maxwell-Gesetz: ein umlaufendes $bold(H)$
entsteht durch einen Strom _oder_ ein sich änderndes $bold(D)$). Der Term
$partial bold(D)\/partial t$ heißt _Verschiebungsstrom_ @lit:GriffithsDE[S. 422];
erst er schließt den Kreis und bildet damit die Grundlage für die Ausbreitung
_elektromagnetischer_ Wellen.

== Materialgleichungen <sec-material>

Die vier Gleichungen verknüpfen sechs Feldgrößen und sind unterbestimmt; es
fehlt die Angabe, wie ein Material auf ein Feld reagiert. Diese
_Materialgleichung_ folgt nicht aus den Maxwell-Gleichungen und hängt vom
Material ab; in ihrer einfachsten, für _lineare_ Medien gültigen Form
@lit:GriffithsDE[S. 422] wird hier durchgehend angesetzt

$ bold(D) = epsilon bold(E), quad bold(B) = mu bold(H), quad
  bold(J) = sigma bold(E) $ <material>

mit $epsilon = epsilon_0 epsilon_r$ und $mu = mu_0 mu_r$
sowie der Leitfähigkeit $sigma$; der letzte Ausdruck ist das Ohmsche Gesetz in
lokaler Form @lit:GriffithsDE[S. 368].


Dass @material in dieser Form gilt, ist eine Annahme, die vier Vereinfachungen @lit:Taflove2005[S. 52 f.] voraussetzt:

/ Linearität: $epsilon,mu,sigma$ unabhängig von der Feldstärke.
/ Isotropie: $epsilon,mu,sigma$ Zahlen statt Tensoren --- $bold(D) parallel bold(E), bold(J) parallel bold(E)$.
/ Nichtdispersivität: $epsilon,mu,sigma$ unabhängig von der Frequenz.
/ Zeitinvarianz: $epsilon,mu,sigma$ ändern sich während der Simulation nicht.

Hinzu kommt die Festlegung auf _nichtmagnetische_ Materialien ($mu_r=1$,
$mu=mu_0$ überall); das Programm in Kapitel 4 kennt entsprechend nur eine
ortsabhängige Permittivität $epsilon_r(x,y)$ und Leitfähigkeit $sigma(x,y)$,
innerhalb einer Gitterzelle jedoch konstant (_stückweise homogen_).

== Reduktion auf zwei Dimensionen: der TM#sub[z]-Fall <sec-zweidimensional>

Eine weitere Vereinfachung ist die Beschränkung auf zwei Dimensionen: Felder und
Materialien werden nur in einer Ebene betrachtet.
Die Reduktion ist dabei keine Näherung, sondern _exakt_, sofern weder Felder
noch Materialverteilung von der dritten Koordinate abhängen:

$ partial/(partial z) equiv 0 $ <zinvarianz>

Unter @zinvarianz  zerfallen @faraday und @ampere in zwei vollständig
entkoppelte, unabhängig lösbare Gruppen: den _TM#sub[z]-Fall_ (transversal
magnetisch, Komponenten $E_z,H_x,H_y$) und den _TE#sub[z]-Fall_ (transversal
elektrisch, $E_x,E_y,H_z$) @lit:GriffithsDE[S. 516].

Beide Fälle sind strukturell gleich gebaut; die Arbeit verwendet durchgehend
den TM#sub[z]-Fall, weil er mit nur einer elektrischen Feldkomponente Anregung
und Randbedingung auf eine skalare Größe zurückführt. Mit $mu=mu_0$ und unter Verwendung der in @sec-material angesprochenen Vereinfachungen lauten die Gleichungen @lit:Taflove2005[Gl. (3.13)]

$ mu_0 (partial H_x)/(partial t) = - (partial E_z)/(partial y) $ <tmz-hx>
$ mu_0 (partial H_y)/(partial t) = (partial E_z)/(partial x) $ <tmz-hy>
$ epsilon (partial E_z)/(partial t) + J_z^"ges"
  = (partial H_y)/(partial x) - (partial H_x)/(partial y) $ <tmz-ez>

wobei $J_z^"ges"$ die gesamte Stromdichte in $z$-Richtung bezeichnet; ihre
Aufteilung folgt in @sec-update.


== Randbedingungen an Grenzflächen <sec-randbedingungen>

An der Grenze zweier Medien ändern sich $epsilon, mu, sigma$ sprunghaft. Aus der
Integralform der Maxwell-Gleichungen folgt, dass die zur Grenzfläche
_parallelen_ Feldanteile dabei stetig bleiben @lit:GriffithsDE[S. 424, Gl. (7.62),
Gl. (7.63)]:

$ bold(E)_(1,parallel) = bold(E)_(2,parallel), quad
  bold(H)_(1,parallel) = bold(H)_(2,parallel) $ <stetig>

Für $bold(H)$ gilt dies nur, solange auf der Grenzfläche kein Strom fließt; bei
den verlustfreien Dielektrika in @sec-tf3 ist das der Fall.

Ein wichtiger Sonderfall ist der _ideal leitende_ Rand (PEC): Für
$sigma arrow.r infinity$ bleibt die Stromdichte $bold(J) = sigma bold(E)$ nur
endlich, wenn $bold(E)$ im Inneren verschwindet @lit:GriffithsDE[S. 368]; mit
@stetig folgt an der Grenzfläche $bold(E)_parallel = 0$, im TM#sub[z]-Fall
also $E_z = 0$.

== Die elektromagnetische Wellengleichung <sec-wellengleichung>

Aus @faraday und @ampere folgt für ein quellen- und verlustfreies Gebiet
($rho_"frei"=0$ und $bold(J)_"frei"=sigma bold(E)=0$, also $sigma=0$) mit konstantem $epsilon,mu$ eine einzelne
Gleichung für $bold(E)$ @lit:GriffithsDE[S. 481 f.]: Bildet man die Rotation von
@faraday, nutzt die Vektoridentität @lit:GriffithsDE[S. 482]

$ nabla times (nabla times bold(E)) = nabla (nabla dot bold(E)) - nabla^2 bold(E) $

und setzt @gauss-e sowie @ampere ein, ergibt sich die _Wellengleichung_
@lit:GriffithsDE[mod. Gl. (9.41)]

$ nabla^2 bold(E) - mu epsilon (partial^2 bold(E))/(partial t^2) = 0 $ <welle>

Dieselbe Gleichung folgt für $bold(H)$, wenn man stattdessen bei @ampere
ansetzt. In @welle kommen keine Ladungen oder Ströme mehr vor: Eine einmal
erzeugte Störung breitet sich mit der durch $mu epsilon$ festgelegten
Geschwindigkeit aus, im Vakuum @lit:GriffithsDE[S. 482]

$ c = 1/sqrt(mu_0 epsilon_0) approx 2,998 dot 10^8 " m/s" $

also der Lichtgeschwindigkeit --- woraus Maxwell schloss, dass Licht selbst
eine elektromagnetische Welle ist.

== Die ebene Welle als analytische Referenzlösung <sec-ebenewelle>

Die einfachste Lösung von @welle, die in Kapitel 5 als Vergleichsgröße dient,
ist die ebene Welle

$ bold(E)(bold(r), t) = E_0 thin e^(upright(i)(omega t - bold(k) dot bold(r))) hat(z) $ <ansatz>

Darin ist $E_0$ die Amplitude, $hat(z)$ der Einheitsvektor in $z$-Richtung ---
im TM#sub[z]-Fall hat $bold(E)$ nur die Komponente $E_z$ ---, $bold(r)$ der
Ortsvektor, $omega=2pi f$ die Kreisfrequenz und $bold(k)$ der Wellenvektor,
dessen Richtung die Ausbreitungsrichtung angibt @lit:GriffithsDE[mod. Gl. (9.43)].
Wegen @zinvarianz liegt $bold(k)$ stets in der $x y$-Ebene, sodass
$bold(k) dot bold(r) = k_x x + k_y y$ gilt. Physikalisch ist der Realteil des
komplexen Ansatzes; $upright(i)$ bezeichnet hier und im Folgenden die imaginäre
Einheit. Ihren Namen trägt die ebene Welle, weil alle Orte gleicher Phase zu
jedem Zeitpunkt eine zu $bold(k)$ senkrechte Ebene bilden @lit:GriffithsDE[S. 483],
die sich mit der Zeit in Richtung $bold(k)$ verschiebt (@abb:ebene-welle).

#figure(
  image("../abbildungen/abb_ebene_welle.png", width: 58%),
  caption: [Ebene Welle: Der Feldvektor zeigt in $z$-Richtung, die Ausbreitung
    erfolgt längs $bold(k) parallel hat(x)$. Die eingezeichneten Wellenflächen
    stehen senkrecht auf $bold(k)$; auf jeder von ihnen hat das Feld zu einem
    festen Zeitpunkt überall denselben Wert. Abbildung entnommen aus
    @lit:Orosz2022.],
) <abb:ebene-welle>

== Die Dispersionsrelation <sec-dispersionsrelation>

Mit $omega$ und $bold(k)$ enthält @ansatz zwei zunächst frei wählbare Größen,
die die Wellengleichung jedoch aneinander bindet. Setzt man den Ansatz in
@welle ein, so reproduziert jede Ableitung den Exponentialausdruck bis auf
einen Faktor: Die zweifache Zeitableitung liefert $-omega^2$, der
Laplace-Operator $-k^2$ mit $k^2 = k_x^2 + k_y^2$, da $k_z$ nach @zinvarianz
entfällt. Der gemeinsame Faktor kürzt sich heraus, und aus der
Differentialgleichung wird eine algebraische Bedingung, die
_Dispersionsrelation_ @lit:GriffithsDE[vgl. S. 483]

$ k^2 = mu epsilon omega^2 $ <dispersion>

Eine ebene Welle löst @welle also genau dann, wenn $k$ und $omega$ in diesem
allein vom Medium bestimmten Verhältnis stehen. Das Verhältnis $omega\/k$ ist zugleich die
Geschwindigkeit, mit der sich eine Fläche fester Phase verschiebt --- die
_Phasengeschwindigkeit_

$ v_p = omega/k = 1/sqrt(mu epsilon) $ <phasengeschwindigkeit>

Entscheidend ist, was in @phasengeschwindigkeit _nicht_ auftritt: weder $omega$
noch die Richtung von $bold(k)$. Im Kontinuum laufen also alle Frequenzen und
alle Richtungen gleich schnell --- Ersteres folgt aus dem in @sec-material
frequenzunabhängig angesetzten $epsilon$, Letzteres aus der dort geforderten
Isotropie. 

== Die Wellenimpedanz <sec-impedanz>

In einer ebenen Welle sind $bold(E)$ und $bold(H)$ nicht unabhängig voneinander:
Setzt man @ansatz in @faraday ein, so folgt, dass beide senkrecht aufeinander
stehen und ihr Amplitudenverhältnis allein vom Medium abhängt. Dieses
Verhältnis heißt _Wellenimpedanz_ @lit:Schneider2010[Gl. (7.46), S. 169]

$ eta = E_0/H_0 = sqrt(mu/epsilon)
  quad "bzw. im Vakuum" quad
  eta_0 = sqrt(mu_0/epsilon_0) approx 376,7 thin Omega $ <impedanz>

Trifft eine Welle auf eine Grenzfläche, an der sich $eta$ ändert, so müsste sie
auf der anderen Seite ein anderes Verhältnis von $E$ zu $H$ haben, während
@stetig verlangt, dass beide Felder an der Grenzfläche übereinstimmen. Beides
zugleich ist nur möglich, wenn ein Teil der Welle reflektiert wird; darauf
beruht die Reflexion am Dielektrikum in @sec-tf3.

== Grenzen analytischer Lösungsverfahren
Die in diesem Kapitel behandelten Lösungen sind stark idealisiert: Die ebene
Welle setzt einen unendlichen, homogenen Raum voraus, die Reflexion eine ebene,
unendlich ausgedehnte Grenzfläche. Geschlossene analytische Lösungen existieren
nur für solche Sonderfälle.
Reale Anordnungen wie das _iPhone_ aus der Einleitung erfüllen keine dieser
Voraussetzungen: Die Geometrie ist gekrümmt und endlich, das Material
ortsabhängig ($epsilon_r (x,y)$, $sigma(x,y)$). Für sie muss das
kontinuierliche Gleichungssystem diskretisiert und numerisch gelöst werden ---
im folgenden Kapitel mit der _Finite-Difference Time-Domain_-Methode (FDTD).