// ---------------------------------------------------------------
// main.typ — Haupt-/Einstiegsdatei der Seminararbeit
// ---------------------------------------------------------------
// Zum Kompilieren immer DIESE Datei übersetzen (nicht die Kapitel
// einzeln). Metadaten unten anpassen, Kapiteltexte liegen unter
// kapitel/*.typ und werden per #include eingebunden — so lässt
// sich an einem Kapitel arbeiten, ohne den Rest zu verlieren.

#import "template.typ": *
#import "@preview/alexandria:0.2.2": *

// Alexandria erlaubt ein in Gruppen gegliedertes Literaturverzeichnis
// (das eingebaute bibliography() kann nur eine Liste). Zitate im Text
// tragen deshalb das Präfix "lit:", z. B. @lit:Yee1966[S. 303].
#show: alexandria(prefix: "lit:", read: path => read(path))

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
// Die Nummern richten sich nach der ersten Zitierung im Text; innerhalb jeder
// Gruppe stehen die Einträge aufsteigend sortiert.
#load-bibliography("literatur.bib", style: "ieee")

// Zuordnung der Quellen zu den Gruppen. Neue Quellen hier eintragen —
// ein zitierter Schlüssel, der in keiner Gruppe steht, erscheint unter
// "Sonstige Quellen", damit kein Verweis ins Leere zeigt.
#let gruppen = (
  ([Primärliteratur (Originalarbeiten)],
    ("Maxwell1865", "Yee1966", "Berenger1994", "Schuhmann2001", "ShinFan2012", "Liu2022")),
  ([Sekundärliteratur (Lehrbücher und Lehrmaterialien)],
    ("Taflove2005", "Sullivan2013", "GriffithsDE", "Fleisch2008", "Schneider2010", "Orosz2022")),
  ([Software-Dokumentation und Online-Quellen],
    ("NumPyFFT", "CSTStudio", "AstroPi", "AppleSpecs")),
)

#heading(numbering: none)[Literaturverzeichnis]
#context {
  let bib = get-bibliography("lit:")
  let alle = gruppen.map(g => g.at(1)).flatten()
  let rest = bib.references.filter(r => r.key not in alle)
  for (titel, keys) in gruppen + (([Sonstige Quellen], rest.map(r => r.key)),) {
    let refs = bib.references.filter(r => r.key in keys)
    if refs.len() > 0 {
      heading(level: 2, numbering: none, outlined: false, titel)
      render-bibliography((..bib, references: refs), title: none)
    }
  }
}

// ---------------- Anhang ----------------
#counter(heading).update(0)
#set heading(numbering: "A.1", supplement: [Anhang])
#include "kapitel/anhang.typ"

// ---------------- Offene Punkte (vor der Abgabe löschen) ----------------