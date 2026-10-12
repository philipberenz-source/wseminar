// ---------------------------------------------------------------
// main.typ — Haupt-/Einstiegsdatei der Seminararbeit
// ---------------------------------------------------------------
// Zum Kompilieren immer DIESE Datei übersetzen (nicht die Kapitel
// einzeln). Metadaten unten anpassen, Kapiteltexte liegen unter
// kapitel/*.typ und werden per #include eingebunden — so lässt
// sich an einem Kapitel arbeiten, ohne den Rest zu verlieren.

#import "template.typ": *

#show: doc => project(
  title: "Die FDTD-Methode zur numerischen Lösung der Maxwell-Gleichungen",
  subtitle: "Theoretische Herleitung, Implementierung und analytische Validierung  ",
  author: "Philip Berenz",
  school: "Oskar-von-Miller-Gymnasium",
  course: "Wissenschaftspropädeutisches Seminar im Leitfach Physik",
  seminar-topic: "Wellen und Teilchen",
  supervisor: "Roman Störmer",
  date: datetime(year: 2026, month: 8, day: 1),
  doc,
)

#include "kapitel/01_einleitung.typ"
#include "kapitel/02_grundlagen.typ"
#include "kapitel/03_fdtd_methode.typ"
#include "kapitel/04_implementierung.typ"
#include "kapitel/05_vergleich.typ"
#include "kapitel/06_fazit.typ"

#pagebreak()
// Zitierstil: "ieee" — numerische Verweise in eckigen Klammern, [1, S. 303].
// Begründung: das ist die Konvention der Physik und der Elektrotechnik; die
// hier zitierten Arbeiten (Yee, Taflove/Hagness, Bérenger, Sullivan) sind
// selbst so gesetzt. Eckige Klammern kollidieren zudem nicht mit der
// Gleichungsnummerierung (34) und den Klammern im Fließtext, wie es die
// runden Klammern von "iso-690-numeric" taten.
// Umstellen auf Fußnoten (deutsche Zitierweise) erfordert nur den Austausch
// des style-Arguments gegen "chicago-notes"; die Zitate im Text
// (@Key bzw. @Key[S. 51]) bleiben dabei unverändert.
#bibliography("literatur.bib", title: "Literaturverzeichnis", style: "ieee")

// ---------------- Anhang ----------------
#counter(heading).update(0)
#set heading(numbering: "A.1", supplement: [Anhang])
#include "kapitel/anhang.typ"

// ---------------- Offene Punkte (vor der Abgabe löschen) ----------------
#pagebreak()
#todo[
  *Offene Punkte vor der Abgabe* (diesen Block danach aus `main.typ` löschen)

  + Literatur Liu 2022 (Yee-Abbildung, Kap. 3): Band, Seiten, DOI und ggf.
    weitere Autoren nachtragen.
  + Literatur Griffiths: Auflage und Erscheinungsjahr der verwendeten Ausgabe
    prüfen.
  + Kap. 3, CFL-Bedingung: Seitenzahl der Quelle bei Taflove ergänzen.
  + Kap. 5.1: drei bis vier Sätze zu Messfehlern ergänzen (warum Rundungs- und
    FFT-Fehler klein sind; Einfluss von Sondenabstand, $sqrt(r)$-Korrektur und
    Referenzgebiet).
  + Kap. 4: einen Satz am Anfang, dass die NumPy-Umsetzung ausdrücklich Teil der
    Aufgabenstellung ist; „Löser“ zu „Solver“ vereinheitlichen.
  + Ton vereinheitlichen: Kap. 2.4 und 2.8, Kap. 3.6.2 (PML), Kap. 5.3.3
    (Deutung Testfall 2), Kap. 5.5 (Testfall 4), Kap. 5.6 (Zwischenfazit) und
    Kap. 6.1 (Fazit) in eigene Sprache bringen.
  + Erklärung: Ort, Datum und Unterschrift eintragen; KI-Formulierung mit dem
    Kursleiter abstimmen.
]
