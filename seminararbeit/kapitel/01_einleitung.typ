// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Einleitung

== Relevanz und Motivation

Als der Airbus A380 im Dezember 2006 die gemeinsame Musterzulassung von EASA
und FAA erhielt und zehn Monate später den Liniendienst aufnahm @EASA_A380, lag
die Vergabe des Auftrags für sein Wetterradar bereits mehr als vier Jahre
zurück @Honeywell_RDR4000. Wesentlicher Bestandteil des Radars ist die Antenne, die in der Bugspitze hinter dem
Radom sitzt, welches sie vor Fahrtwind und Niederschlag schützt. Über den gesamten
Zeitraum der Entwicklung hinweg existierte diese Antenne in ihrer endgültigen
Einbausituation nicht, und dennoch musste sie konzepiert und für den Airbus optimiert werden. Möglich war das nur mittels Computern.

Anspruchsvoll wird eine solche Simulation durch das Bauteil, das der Antenne
vorgelagert ist. Das Wetterradar arbeitet im X-Band bei etwa 9,4 GHz, und das
Radom davor ist keine einfache Trennwand, sondern eine gekrümmte,
mehrschichtige und verlustbehaftete Kunststoffschale. Deren Aufbau schwächt die
durchtretende Welle und verschiebt ihre Phase über die Fläche hinweg
unterschiedlich stark, sodass die Phasenfront verkippt und die Antennenkeule
gegenüber der mechanischen Achse abgelenkt wird. Diesen Winkelfehler bezeichnet
man in der Radartechnik als _boresight error_; er bildet zusammen mit der
Transmissionsdämpbfung und dem Nebenkeulenpegel die zentrale Entwurfsgröße eines
Radoms @Radome2015. Eine Gewitterzelle erscheint dann unter einer leicht
falschen Richtung. Da der Fehler erst aus dem Zusammenspiel von Antenne,
Radomgeometrie und Materialaufbau entsteht, lässt er sich an der Antenne allein
nicht bestimmen.

Dieser Fall ist jedoch keinesfalls lediglich Sonderfall der Luftfahrt, sondern steht stellvertretend
für eine ganze Klasse technischer Problemstellungen. Zu ihr gehören der Entwurf
von Antennenarrays hinter dielektrischen Abdeckungen, die Auslegung von Filtern
und Wellenleitern der integrierten Photonik sowie die Konstruktion von Sende-
und Empfangsspulen in Magnetresonanztomographen. Ihnen allen ist gemeinsam,
dass mehrere Materialien unterschiedlicher Permittivität und Leitfähigkeit auf
unregelmäßigen, häufig gekrümmten oder geschichteten Grenzflächen
zusammentreffen und dass das Feld gerade an diesen Übergängen bestimmt werden
muss.

An Konfigurationen dieser Art scheitert die analytische Behandlung. Die
klassische Elektrodynamik gilt zwar als abgeschlossene Theorie, denn die
Maxwell-Gleichungen legen seit 1865 das Feld bei gegebenen Quellen,
Materialeigenschaften und Randbedingungen eindeutig fest @Maxwell1865.
Geschlossene Lösungen sind jedoch nur für wenige idealisierte Konfigurationen
bekannt, etwa für die ebene Welle im homogenen unbegrenzten Raum, für die
Reflexion an einer unendlich ausgedehnten ebenen Grenzfläche, für den
rechteckigen Hohlraumresonator oder für die Streuung an einer Kugel. Allen
diesen Fällen ist gemeinsam, dass das Material homogen und die Geometrie
hochsymmetrisch und idealsiert ist. Sobald diese Voraussetzungen entfallen, sobald also
inhomogene, geschichtete oder verlustbehaftete Materialien mit beliebig
geformten Grenzflächen zusammentreffen, ist das entstehende System partieller
Differentialgleichungen geschlossen nicht mehr lösbar. Die entscheidende Grenze
verläuft damit nicht zwischen bekannter und unbekannter Physik, sondern
zwischen lösbaren und unlösbaren Randwertproblemen(gemeint sind Differentialgleichungen, die durch Randbedingungen näher beschrieben werden).

Sämtliche genannten Anwendungen liegen jenseits dieser Grenze, sodass das Feld
dort auf anderem Weg bestimmt werden muss. Naheliegend wäre, es zu messen statt
es zu berechnen. Dieser Weg scheitert jedoch an drei Punkten. Erstens setzt
jede Messung ein gefertigtes Exemplar voraus, das in den frühen
Entwurfsphasen gerade nicht zur Verfügung steht. Zweitens erfordert jede
Entwurfsänderung einen neuen Prototyp, gleich ob die Wandstärke, der
Schichtaufbau oder die Krümmung verändert wird, sodass sich mehrere Varianten
weder in vertretbarer Zeit noch zu vernünftigen Kosten systematisch vergleichen
lassen. Drittens liefert eine Messung das Feld nur dort, wo sich eine Sonde
anbringen lässt, die es zudem selbst stört, während der Verlauf innerhalb des
Materials unzugänglich bleibt. Gerade dort aber entsteht der gesuchte Effekt.
Die numerische Lösung ist deshalb kein bequemer Ersatz für das Experiment,
sondern in vielen Fällen der einzige gangbare Weg. In weiten Teilen der
Hochfrequenztechnik ist die Simulation an die Stelle des Prototyps getreten,
und Entwurfe werden konzepiert, bevor ein Bauteil überhaupt gefertigt wird.

Rechnerisch zugänglich werden diese Probleme durch Diskretisierung. Dabei wird
das Kontinuum durch ein endliches Gitter ersetzt und die
Differentialquotienten werden durch Differenzenquotienten approximiert. Das
verbreitetste Verfahren dieser Art ist die _Finite-Difference
Time-Domain_-Methode (FDTD), die die zeitliche Entwicklung des Feldes
schrittweise auf einem diskreten Raumgitter berechnet.

Durchgeführt werden solche Rechnungen in der Praxis fast ausschließlich mit
kommerzieller Software. Ein verbreitetes Paket ist CST Studio Suite (Dassault
Systèmes/SIMULIA) @CSTStudio, dessen Quellcode nicht öffentlich zugänglich ist.
Aus Sicht der Anwendung liegt damit ein Blackbox-System vor, denn Geometrie,
Materialparameter und Anregung werden vorgegeben und eine Feldverteilung wird
ausgegeben, ohne dass der Weg dazwischen von außen nachvollziehbar wäre.
Methodisch ist das bedeutsam, weil ein Simulationsergebnis die exakte Lösung
eines Systems von Differenzengleichungen darstellt und nicht die Lösung der
Maxwell-Gleichungen. Beide stimmen erst im Grenzfall verschwindender
Gitterweite überein. Auf einem endlichen Gitter treten deshalb systematische
Abweichungen auf, die sich im Ergebnis nicht ohne Weiteres von physikalischen
Effekten unterscheiden lassen.

Umso auffälliger ist, wie selten die naheliegendste Rückfrage gestellt wird,
nämlich inwiefern eine solche numerische Lösung die zugehörige analytische
überhaupt reproduziert. Das hat einen strukturellen Grund. Für genau jene
Konfigurationen, aufgrund derer simuliert wird, existiert definitionsgemäß keine
analytische Vergleichslösung, sodass sich die Übereinstimmung nur dort
überprüfen ließe, wo man gar nicht simulieren müsste. Die wenigen exakt
lösbaren Fälle gelten in der Anwendung als uninteressant und bleiben daher
unbeachtet. Gerade weil die Anzeige von Gewitterzellen sicherheitsrelevant und
ein gefertigtes Radom nur noch begrenzt korrigierbar ist, liegt es nahe, diese
Rückfrage einmal systematisch zu stellen. Die vorliegende Arbeit kehrt die
übliche Herangehensweise daher um und behandelt die analytisch lösbaren Fälle
nicht als Trivialfälle, sondern als die einzigen Anordnungen, an denen sich das
Verfahren überhaupt prüfen lässt. Dazu wird ein eigener, vollständig
dokumentierter FDTD-Löser in Python implementiert und gegen exakt berechenbare
Konfigurationen gehalten --- nicht um bekannte Ergebnisse zu bestätigen,
sondern um zu bestimmen, wie genau das Verfahren rechnet und wovon diese
Genauigkeit abhängt.

== Historischer Kontext

Die Entwicklung von der Feldtheorie zum numerischen Standardverfahren lässt
sich an drei Stationen nachzeichnen.

James Clerk Maxwell führte 1865 in _A Dynamical Theory of the Electromagnetic
Field_ @Maxwell1865 Elektrizität, Magnetismus und Optik in einem einheitlichen
Gleichungssystem zusammen und sagte die Existenz sich mit Lichtgeschwindigkeit
ausbreitender elektromagnetischer Wellen voraus. Der experimentelle Nachweis
gelang Heinrich Hertz rund zwanzig Jahre später. Der so gesetzte theoretische
Rahmen ist bis heute unverändert gültig und bildet die Grundlage der
vorliegenden Arbeit.

Kane S. Yee legte 1966 mit einer knappen Veröffentlichung @Yee1966 die
Grundlage der numerischen Umsetzung. Sein Ansatz besteht darin, die
elektrischen und magnetischen Feldkomponenten nicht am selben Ort und nicht zum
selben Zeitpunkt zu speichern, sondern räumlich und zeitlich gegeneinander
versetzt auf einem sogenannten versetzten Gitter, dem Yee-Gitter, anzuordnen.
Auf diese Weise lassen sich die gekoppelten Rotationsgleichungen durchgängig
mit zentralen Differenzenquotienten approximieren, und beide Feldgrößen können
alternierend fortgeschrieben werden (_leapfrog_-Schema). Die zum damaligen
Zeitpunkt verfügbare Rechenleistung beschränkte die Anwendung allerdings auf
sehr kleine Modellprobleme.

Erst das Wachstum der verfügbaren Rechenkapazitäten machte das Verfahren
praktisch nutzbar. Allen Taflove prägte ab Mitte der 1970er Jahre die
Bezeichnung FDTD und baute die Methode zu einem Standardwerkzeug der
_Computational Electromagnetics_ aus, insbesondere durch die Entwicklung
absorbierender Randbedingungen und die Erweiterung auf verlustbehaftete sowie
dispersive Materialien @Taflove2005. Gegenwärtig ist FDTD in nahezu allen
kommerziellen Programmpaketen zur Feldsimulation vertreten --- unter anderem in
den eingangs erwähnten, deren Innenleben dem Anwender verborgen bleibt.

== Persönliche Motivation

Für Simulationen habe ich mich erstmals in der zehnten Jahrgangsstufe
interessiert, als ich mich mit der numerischen Beschreibung mechanischer Wellen
beschäftigt habe. In der elften Jahrgangsstufe kam die Programmierseite hinzu:
Im Rahmen des Astro-Pi-Wettbewerbs der ESA#footnote[Die European Astro Pi
Challenge ist ein Bildungsprogramm der ESA in Zusammenarbeit
mit der Raspberry Pi Foundation. Ziel ist die Erstellung von
Python-Programmen, die auf zwei mit Sensoren ausgestatteten
Raspberry-Pi-Einheiten an Bord der Internationalen Raumstation ausgeführt
werden @AstroPi.] habe ich eigenen Code geschrieben und dabei erfahren, wie
viel Sorgfalt eine Aufgabenstellung verlangt, deren Ergebnis nicht im eigenen
Ermessen liegt. Die vorliegende Arbeit setzt beides
fort, denn sie führt von der mechanischen zur elektromagnetischen Welle und
verbindet die physikalische Beschreibung mit der eigenen Implementierung.

Hinzu kommt mein Interesse an technischen Anwendungen. In der Schule und in der
theoretischen Physik werden Szenarien weitgehend idealisiert: Man rechnet mit
homogenen Medien, unendlich ausgedehnten Grenzflächen, punktförmigen Quellen
und reibungsfreien Vorgängen. Didaktisch ist das notwendig, denn nur so werden
die zugrunde liegenden Gesetzmäßigkeiten überhaupt sichtbar. Zugleich entsteht
dabei der Eindruck, physikalische Probleme seien stets in geschlossener Form
lösbar, obwohl die Idealisierung in realen Anwendungen gerade dort entfällt, wo
es interessant wird. Ein Verfahren wie FDTD setzt an genau dieser Stelle an,
weil es beliebige Geometrien und Materialverteilungen behandeln kann. Diese
Arbeit bot mir daher die Gelegenheit, einmal ein Problem zu bearbeiten, das
sich nicht durch Vereinfachung auflösen lässt.

== Zentrale Fragestellung und Aufbau der Arbeit

Aus den dargestellten Überlegungen ergibt sich die folgende Leitfrage:

#align(center)[
  #block(width: 92%)[
    _Wie zuverlässig reproduziert eine selbst implementierte zweidimensionale
    FDTD-Simulation die analytisch bekannten Eigenschaften elektromagnetischer
    Wellen, und wovon hängt diese Zuverlässigkeit ab?_
  ]
]

Um sie zu beantworten, wird die FDTD-Methode zunächst aus den
Maxwell-Gleichungen hergeleitet und in ein Python-Programm überführt. Dieses
Programm wird anschließend an drei Konfigurationen geprüft, deren Lösung sich
in geschlossener Form angeben lässt: an der ebenen Welle im freien Raum, an der
Reflexion an einer ebenen Grenzfläche und am rechteckigen Hohlraumresonator.
Die drei Anordnungen sind so gewählt, dass jede von ihnen einen anderen der
drei Diskretisierungsschritte in den Vordergrund rückt --- die des Raums, die
des Materials und die der Zeit ---, sodass sich im Fehlerfall benennen lässt,
welcher der drei versagt hat.

Der Gegenstand der Arbeit ist damit ausdrücklich das Verfahren und nicht ein
bestimmtes Bauteil. Ihr Ergebnis sind keine physikalischen Erkenntnisse ---
die Vergleichsgrößen sind seit Langem bekannt --- sondern quantitative
Fehlergrenzen in Abhängigkeit der Diskretisierung, zusammen mit der Angabe,
woher die verbleibenden Abweichungen jeweils stammen. Beantwortet wird damit
genau jene Rückfrage, die in der Praxis gewöhnlich offenbleibt.

Aus der Leitfrage ergeben sich zwei Teilfragen:

+ Mit welcher relativen Abweichung gibt die Simulation die analytisch bekannten
  Größen wieder, namentlich die Ausbreitungsgeschwindigkeit und die Wellenlänge
  im freien Raum, das Reflexionsverhalten an einer Grenzfläche und die
  Resonanzfrequenzen eines Hohlraums?
+ Wie hängen diese Abweichungen von den Diskretisierungsparametern ab,
  insbesondere von der Anzahl der Gitterzellen pro Wellenlänge, und welche
  prinzipiellen Grenzen des Verfahrens folgen daraus?

Im Einzelnen stellt Kapitel 2 die theoretischen Voraussetzungen bereit: die
Maxwell-Gleichungen in differentieller Form, die Herleitung der
Wellengleichung, die ebene Welle als analytische Referenzlösung sowie die
Ursachen für das Versagen analytischer Verfahren bei realistischen Geometrien.
Aus diesen Grundlagen entwickelt Kapitel 3 die FDTD-Methode, vom
Differenzenquotienten über das Yee-Gitter bis zu den Update-Gleichungen des
zweidimensionalen $"TM"_z$-Falls, und ergänzt sie um die
Courant-Friedrichs-Lewy-Bedingung sowie die verwendeten Randbedingungen. Wie
sich diese Gleichungen in ein lauffähiges Programm übersetzen lassen,
dokumentiert Kapitel 4. Kapitel 5 unterzieht dieses Programm den drei
Testfällen und bestimmt die gesuchten Fehlergrenzen, bevor Kapitel 6 die
Ergebnisse zusammenführt und einen Ausblick darauf gibt, wie sich das
charakterisierte Verfahren auf Problemstellungen anwenden ließe, für die keine
analytische Lösung mehr existiert --- etwa auf den eingangs beschriebenen
Radomdurchgang.
