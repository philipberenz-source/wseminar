// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Fazit und Ausblick

== Fazit

Die vier Testfälle aus Kapitel 5 zeigen, dass alle Abweichungen des Solvers vom
Kontinuum systematisch sind und sich aus den diskreten Gleichungen selbst
vorhersagen lassen. Die Phasengeschwindigkeit liegt unter $c$ und hängt von
Frequenz und Richtung ab (@sec-tf1), die Reflexion an einer Grenzfläche fällt
größer aus als nach Fresnel (@sec-tf3). Beide Abweichungen gehen auf den
Abbruchfehler zweiter Ordnung der zentralen Differenz zurück und fallen daher mit
$N_lambda^(-2)$. Da die Messwerte mit der numerischen Dispersionsrelation
@num-dispersion-n und den Gitterkoeffizienten @tf3-gitter übereinstimmen, setzt
der Solver das Yee-Verfahren korrekt um. Die Courant-Zahl bestimmt dagegen vor
allem die Stabilität: Für $S <= 1$ bleibt die Energie bis auf Rundungsfehler
erhalten, darüber wächst sie exponentiell (@sec-tf2). Am absorbierenden Rand
entscheidet schließlich in erster Linie die Dicke der Schicht über die
verbleibende Reflexion (@sec-tf4).

Die numerische Genauigkeit hängt somit vor allem von der Auflösung $N_lambda$
ab, wobei die Abweichung nach @tf3-naeherung zusätzlich mit der Brechzahl der
beteiligten Medien wächst. Hinzu kommen die Ausbreitungsrichtung, die
Courant-Zahl und am Rand die Dicke der PML. Da sich der Phasenfehler mit der
Laufstrecke aufsummiert, genügt eine einzelne Fehlerangabe nicht; die Auflösung
ist vielmehr nach Größe und Laufzeit des jeweiligen Problems zu wählen.

Diese Fehler sind damit nicht grundsätzlich problematisch, entscheidend ist
jedoch, sie zu kennen: Da ihre Größe vorhersagbar ist, lässt sich vorab
abschätzen, wie verlässlich ein Ergebnis ist. Verringern lassen sie sich durch
eine ausreichende Auflösung im optisch dichtesten Medium, eine Courant-Zahl nahe
der Stabilitätsgrenze und, sofern wählbar, eine Ausbreitung entlang der
Gitterdiagonalen, die für $S = 1$ sogar dispersionsfrei ist (@tf1-taylor). Am
Rand ist zuerst die Dicke der PML zu erhöhen.

== Ausblick

Die gesamte Arbeit beschränkt sich auf den zweidimensionalen $"TM"_z$-Fall mit
frequenzunabhängigen Materialien. Daraus ergeben sich drei naheliegende
Erweiterungen.

*Dritte Dimension.* Das Yee-Gitter lässt sich unmittelbar auf drei Dimensionen
mit allen sechs Feldkomponenten erweitern @Schneider2010[S. 246--253]. Damit
entfiele die Annahme einer in $z$-Richtung unendlich ausgedehnten Anordnung,
und $"TM"_z$- und $"TE"_z$-Anteile würden gemeinsam erfasst. Der Preis ist der
Rechenaufwand: Die Zahl der Zellen wächst dann mit $N_lambda^3$ statt mit
$N_lambda^2$.

*Dispersive Medien.* Die in @sec-material getroffene Vereinfachung einer
frequenzunabhängigen Permittivität trifft auf viele reale Stoffe, etwa
wasserhaltige, nicht zu. Ein Debye-Modell beschreibt solche Medien und lässt
sich im Zeitbereich über eine zusätzliche Differentialgleichung oder eine
rekursive Faltung umsetzen @Schneider2010[Abschn. 10.3, 10.7]. Letztere arbeitet
mit einem Dämpfungsfaktor $e^(-Delta t\/tau_d)$ @Schneider2010[S. 308,
Gl. (10.120)] und folgt damit demselben Prinzip wie die Hilfsgrößen der PML in
@pml-bc, sodass sie sich in die bestehende Struktur von `fdtd_core.py` einfügen
ließe (eigene Einschätzung).

*Beschleunigung auf Grafikprozessoren.* Jede Update-Gleichung benötigt nur die
Werte benachbarter Gitterpunkte, sodass alle Zellen eines Zeitschritts
gleichzeitig berechnet werden können. Parallelisierung über Threads und mehrere
Prozesse behandelt @Schneider2010[Kap. 13]; auf einer GPU ließe sich dasselbe
Prinzip auf tausende Rechenkerne übertragen. Für diese Arbeit ist das unmittelbar
relevant: Halbiert man die Zellweite, so viertelt sich nach Kapitel 5 der
Fehler, im zweidimensionalen Gitter steigt der Aufwand wegen der vierfachen
Zellzahl und der nach @cfl doppelten Schrittzahl jedoch um den Faktor acht.
