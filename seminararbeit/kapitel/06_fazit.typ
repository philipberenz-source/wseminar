// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Fazit und Ausblick

== Fazit

Ziel dieser Arbeit war es, die Maxwell-Gleichungen in eine zweidimensionale
FDTD-Simulation zu überführen, diese in Python umzusetzen und die Genauigkeit
des entstandenen Solvers zu untersuchen.

Den Kern bildet dabei die Herleitung. Ausgehend von den Rotationsgleichungen,
reduziert auf den $"TM"_z$-Fall, wurden die Ableitungen über die
Taylor-Entwicklung durch zentrale Differenzenquotienten zweiter Ordnung ersetzt.
Die räumlich und zeitlich versetzte Anordnung der Feldkomponenten im Yee-Gitter
sorgt dafür, dass jede dieser Differenzen genau dort ausgewertet wird, wo die
Ableitung benötigt wird; so entstehen explizite Update-Gleichungen, die $E$ und
$H$ abwechselnd fortschreiben. Die Stabilitätsbedingung folgt aus demselben
Gleichungssystem, und selbst die PML lässt sich über die komplexe
Koordinatenstreckung und die rekursive Auswertung der entstehenden Faltung in
diese explizite Form überführen. Die Umsetzung in Python benötigt daher nur
wenige Zeilen je Feldkomponente.

Die Validierung in Kapitel 5 bestätigt diese Herleitung in doppelter Hinsicht:
Der Solver reproduziert die analytisch bekannten Eigenschaften
elektromagnetischer Wellen, und auch seine Abweichungen lassen sich aus den
hergeleiteten diskreten Gleichungen vorhersagen (@sec-zwischenfazit). Dieselben
Gleichungen, aus denen das Verfahren aufgebaut ist, beschreiben somit zugleich
seine Grenzen. Gerade darin liegt der Wert einer vollständigen Herleitung: Sie
macht nicht nur nachvollziehbar, wie das Verfahren arbeitet, sondern auch, wann
seinen Ergebnissen zu trauen ist.

Einschränkend ist festzuhalten, dass sich die Untersuchung auf zwei
Dimensionen, frequenzunabhängige Materialien und Anordnungen mit bekannter
analytischer Lösung beschränkt --- also auf Fälle, in denen eine Simulation am
wenigsten erforderlich ist. Die Arbeit versteht sich deshalb als Grundlage: Sie
liefert das Verständnis des Verfahrens, auf dem ein begründeter Einsatz bei
komplexeren Problemen aufbauen kann.

== Ausblick

Die gesamte Arbeit beschränkt sich auf den zweidimensionalen $"TM"_z$-Fall mit
frequenzunabhängigen Materialien. Daraus ergeben sich zwei naheliegende
Erweiterungen.

*Dritte Dimension.* Das Yee-Gitter lässt sich unmittelbar auf drei Dimensionen
mit allen sechs Feldkomponenten erweitern. Damit
entfiele die Annahme einer in $z$-Richtung unendlich ausgedehnten Anordnung,
und $"TM"_z$- und $"TE"_z$-Anteile würden gemeinsam erfasst. Der Preis ist der
Rechenaufwand: Die Zahl der Zellen wächst dann mit $N_lambda^3$ statt mit
$N_lambda^2$.

*Dispersive Medien.* Die in @sec-material getroffene Vereinfachung einer
frequenzunabhängigen Permittivität trifft auf viele reale Stoffe, etwa
wasserhaltige, nicht zu. Ein Debye-Modell beschreibt solche Medien und lässt
sich im Zeitbereich über eine zusätzliche Differentialgleichung oder eine
rekursive Faltung umsetzen @Schneider2010[Abschn. 10.3, 10.7]. Letztere arbeitet
mit einem Dämpfungsfaktor @Schneider2010[S. 308] und folgt damit demselben Prinzip wie die Hilfsgrößen der PML in
@pml-bc.