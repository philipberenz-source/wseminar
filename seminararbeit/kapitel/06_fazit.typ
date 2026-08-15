// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Fazit und Ausblick

== Zusammenfassung

#todo[
  Die Leitfrage aus Abschnitt 1.4 anhand der Ergebnisse aus Kapitel 5
  beantworten. Kernpunkte, die hier zusammenlaufen sollten:
  - Bei zwanzig Zellen pro Wellenlänge liegen alle geprüften Größen innerhalb
    eines Prozents vom analytischen Wert.
  - Der Fehler fällt mit dem Quadrat der Zellweite, und diese Abhängigkeit
    ließ sich aus dem Abbruch der Taylorreihe herleiten, nicht nur messen.
  - Er wächst linear mit der Laufstrecke und hängt von der Ausbreitungsrichtung
    ab --- eine einzelne Fehlerangabe ist deshalb unvollständig.
  - Die beiden größten beobachteten Abweichungen stammten gar nicht aus dem
    Verfahren, sondern aus der Auswertung.
  Wichtig: Hier nicht die Ergebnisse aus Kapitel 5 wiederholen, sondern sie auf
  die Leitfrage beziehen und bewerten.
]

== Grenzen der Arbeit

#todo[
  Ehrlich benennen, was die Validierung nicht abdeckt --- die Aufstellung aus
  Abschnitt 5.5 lässt sich hier straffen:
  - nur gitterparallele Grenzflächen, keine gekrümmten,
  - keine verlustbehafteten und keine frequenzabhängigen Materialien,
  - kein schräger Einfall, keine dritte Dimension,
  - und grundsätzlich: Alle Testfälle sind einfach genug, um analytisch lösbar
    zu sein. Das ist ihr Zweck, aber auch ihre Beschränkung.
]

== Ausblick

#todo[
  Folgende Punkte ausformulieren --- der erste ist der wichtigste, weil er die
  Klammer zur Einleitung schließt:
]

- *Anwendung auf ein analytisch unzugängliches Problem.* Der eingangs
  beschriebene Radomdurchgang wäre der naheliegende nächste Schritt: ein
  zweidimensionaler Schnitt durch Antennenapertur und gekrümmte Radomschale,
  aus dem sich die Ablenkung der Hauptkeule bestimmen ließe. Erst die in
  Kapitel 5 bestimmten Fehlergrenzen machen aus dem dabei berechneten Winkel
  eine belastbare Angabe. Zu beachten wäre insbesondere, dass an einer
  gekrümmten Fläche die Treppenstufennäherung und die numerische Anisotropie
  zusammenwirken und beide einen scheinbaren Winkelfehler erzeugen.
- *Verlustbehaftete Materialien.* Ein vierter Testfall zur Eindringtiefe in
  einem leitfähigen Medium würde die in Abschnitt 4.1 eingeführten
  Materialkoeffizienten auch außerhalb ihres verlustfreien Grenzfalls prüfen.
- *Erweiterung auf drei Dimensionen*, mit dem vollen Feldsatz und deutlich
  höherem Rechenaufwand.
- *Frequenzabhängige Medien*, etwa über ein Debye-Modell, womit sich auch
  breitbandige Anregungen in realen Materialien behandeln ließen.
- *Effizienzsteigerung* durch Parallelisierung oder GPU-Beschleunigung, um
  feiner aufgelöste Gitter in vertretbarer Zeit rechnen zu können --- was nach
  den Ergebnissen aus Kapitel 5 unmittelbar die Genauigkeit erhöht.
