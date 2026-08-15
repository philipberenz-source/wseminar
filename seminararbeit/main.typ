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
  subtitle: "Theoretische Herleitung, analytische Validierung und numerische Fallstudie",
  author: "Philip Berenz",
  school: "Oskar-von-Miller-Gymnasium",
  course: "Wissenschaftspropädeutisches Seminar im Leitfach Physik",
  seminar-topic: "Wellen und Teilchen",
  supervisor: "Roman Störmer",
  date: datetime(year: 2026, month: 8, day: 1),
  doc,
)

// Notationskonventionen (einmal festlegen, dann durchgängig verwenden):
// E: elektrische Feldstärke, D: elektrische Flussdichte, H: magnetische
// Feldstärke, B: magnetische Flussdichte, J: Stromdichte, σ: Leitfähigkeit,
// ε = ε₀εᵣ: Permittivität, μ = μ₀μᵣ: Permeabilität, ρ_frei: freie Ladungsdichte.
// Werden Symbole später ergänzt oder umbenannt (z. B. zusätzliche
// Verlustterme), bitte hier vermerken und konsistent in allen Kapiteln anpassen.

#include "kapitel/01_einleitung.typ"
#include "kapitel/02_grundlagen.typ"
#include "kapitel/03_fdtd_methode.typ"
#include "kapitel/04_implementierung.typ"
#include "kapitel/05_vergleich.typ"
#include "kapitel/06_fazit.typ"

#pagebreak()
#bibliography("literatur.bib", title: "Literaturverzeichnis", style: "iso-690-numeric")
