// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Fazit und Ausblick

== Fazit

Die Herleitung in Kapitel 3 zeigt, dass sich das Verfahren aus wenigen Bausteinen ergibt:
Zentrale Differenzen ersetzen die Ableitungen, und erst die räumlich und
zeitlich versetzte Anordnung der Felder im Yee-Gitter macht daraus explizite
Update-Gleichungen, die elektrisches und magnetisches Feld abwechselnd
fortschreiben. Den Zeitschritt begrenzt die Courant-Friedrichs-Lewy-Bedingung,
und auch die absorbierende Randschicht fügt sich über die komplexe
Koordinatenstreckung und die rekursive Auswertung der entstehenden Faltung in
dieses Verfahren ein. Mit NumPy (Kapitel 4) lassen sich die Update-Gleichungen vektorisiert auf
ganze Felder anwenden, wodurch sich Simulationen effizient ausführen lassen.

Dieselbe Diskretisierung, die das Verfahren so einfach macht, bestimmt auch
seine Grenzen. Die Validierung bestätigt zwar, dass der Solver Ausbreitung, Reflexion
und Transmission ebener Wellen wiedergibt, aber nicht exakt: Der Abbruchfehler
der zentralen Differenz verlangsamt die Welle im Gitter, die Verknüpfung der
Felder eine halbe Zelle neben der Grenzfläche erhöht Reflexions- und
Transmissionskoeffizient, und das stufige Leitfähigkeitsprofil trägt zur
Restreflexion der PML bei. Die Dispersion hängt von Auflösung, Courant-Zahl und
Ausbreitungsrichtung ab, der Grenzflächenfehler von Auflösung und
Permittivität; beide lassen sich mit den diskreten Gleichungen
vorhersagen. Bei der Randschicht sind ihre Dicke, die vorgegebene
Restreflexion und der Einfallswinkel maßgeblich.

Ein Fehler ist damit nicht von vornherein ein Mangel der Simulation. Lässt sich
einordnen, woraus er sich ergibt und wovon er abhängt, so können Auflösung,
Zeitschritt und Randschicht vor einer Rechnung so gewählt werden, dass er unter
einer geforderten Grenze bleibt. Die Zuverlässigkeit eines Ergebnisses hängt
somit weniger davon ab, ob ein Fehler auftritt, als davon, ob er verstanden
ist.

Im Rahmen der Arbeit ergeben sich jedoch auch mehrere Grenzen. Das Modell beschränkt sich auf zwei
Dimensionen und frequenzunabhängige Materialien, und auch innerhalb dieses
Modells blieben Fehlerquellen und Testgrößen unbeleuchtet: Weder die
elektrische Leitfähigkeit und die daraus folgende Absorption noch die
Energieerhaltung in einem geschlossenen System oder Eigenfrequenzen wurden
untersucht. Ungeprüft blieben zudem Fehler durch die Geometrie,
etwa der _Staircasing_-Fehler, der entsteht, wenn das Gitter schräge oder
gekrümmte Grenzflächen nur stufenförmig auflöst. 
Reale Geräte wie das
Smartphone aus der Einleitung gehen jedoch über die getroffenen Vereinfachungen hinaus, und die
gewonnenen Zusammenhänge gelten zunächst nur innerhalb des untersuchten
Rahmens. Sie zeigen aber, welche Größen auch dort für die Verlässlichkeit einer
FDTD-Simulation zu beachten sind, wo keine analytische Lösung zum Vergleich
vorliegt.

== Ausblick

Die Einschränkungen aus dem Fazit geben die nächsten Schritte vor. Am
naheliegendsten ist es, die bisher ungeprüften Fehlerquellen mit dem
vorhandenen Solver zu untersuchen. Die Leitfähigkeit ist durch die Koeffizienten
$C_a$ und $C_b$ bereits implementiert, sodass sich die Reflexion an einem
verlustbehafteten Halbraum mit demselben Zwei-Lauf-Verfahren wie in Testfall 2
messen und mit der analytischen Lösung für leitfähige Medien vergleichen ließe.
Energieerhaltung und Eigenfrequenzen ließen sich in einem ideal leitenden
Hohlraum prüfen, und an einer schräg im Gitter liegenden Grenzfläche könnte abgeschätzt werden, wie stark der Staircasing-Fehler ins Gewicht fällt.

Zwei weitere Erweiterungen betreffen das Modell selbst. Die erste ist der
Übergang in die dritte Dimension: Das Yee-Gitter lässt sich unmittelbar auf drei
Dimensionen mit allen sechs Feldkomponenten übertragen. Damit entfiele die
Annahme einer in $z$-Richtung unendlich ausgedehnten Anordnung, und
$"TM"_z$- und $"TE"_z$-Anteile würden gemeinsam erfasst. Einhergehen würde dies mit erhöhtem
Rechenaufwand, da die Zahl der Zellen dann mit $N_lambda^3$ statt mit
$N_lambda^2$ wächst. Die zweite betrifft die Materialien: Die in @sec-material
getroffene Vereinfachung einer frequenzunabhängigen Permittivität trifft auf
viele reale Stoffe, etwa wasserhaltige, nicht zu. Polare Stoffe dieser Art
beschreibt ein Debye-Modell, das sich im Zeitbereich über eine zusätzliche
Differentialgleichung oder eine rekursive Faltung umsetzen lässt. Letztere
arbeitet mit einem Dämpfungsfaktor @lit:Schneider2010[S. 308] und folgt damit
demselben Prinzip wie die Hilfsgrößen der PML in @pml-bc.

Erst nach diesen Schritten wäre eine praxisnahe Anordnung ohne analytische
Lösung der nächste sinnvolle Gegenstand, an dem sich die hier gewonnenen
Zusammenhänge bewähren müssten.
