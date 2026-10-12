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

Die Validierung in Kapitel 5 zeigt, dass der Solver die analytisch bekannten
Eigenschaften elektromagnetischer Wellen reproduziert, wenn auch nicht exakt.
Die verbleibenden Abweichungen sind jedoch keine Mängel der Implementierung,
sondern systematische Folgen der Diskretisierung, die sich aus den hergeleiteten
diskreten Gleichungen erklären und vorhersagen lassen (@sec-zwischenfazit).

Daraus ergibt sich eine zentrale Erkenntnis dieser Arbeit: Ein numerischer Fehler
ist nicht grundsätzlich problematisch, sofern seine Ursache und sein Verhalten
bekannt sind. Ist bekannt, wovon eine Abweichung abhängt und wie sie sich mit
Auflösung, Zeitschritt oder Randschicht verändert, so lassen sich die Parameter
einer Simulation gezielt so wählen, dass der Fehler unterhalb einer geforderten
Grenze bleibt. Die Zuverlässigkeit eines Ergebnisses hängt somit weniger davon
ab, ob ein Fehler auftritt, als davon, ob er verstanden und abgeschätzt werden
kann.

Nachgewiesen wurde dies an Anordnungen mit bekannter analytischer Lösung. Die
gewonnenen Zusammenhänge bilden jedoch die Grundlage dafür, das Verfahren auch
auf Probleme anzuwenden, für die keine analytische Lösung existiert, und die
Verlässlichkeit der Ergebnisse dort begründet einzuschätzen.

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
rekursive Faltung umsetzen. Letztere arbeitet
mit einem Dämpfungsfaktor @Schneider2010[S. 308] und folgt damit demselben Prinzip wie die Hilfsgrößen der PML in
@pml-bc.