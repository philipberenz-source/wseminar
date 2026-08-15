// ---------------------------------------------------------------
// template.typ — Formatvorlage für die W-Seminararbeit
// ---------------------------------------------------------------
// Hier werden Seitenlayout, Titelseite, Nummerierung und ein paar
// Hilfsfunktionen (todo, note) definiert. Inhaltliche Änderungen
// bitte in main.typ bzw. in den Dateien unter kapitel/ vornehmen,
// nicht hier.

#let project(
  title: "Titel der Arbeit",
  subtitle: none,
  author: "Vorname Nachname",
  school: "Schule",
  course: "W-Seminar",
  seminar-topic: "Rahmenthema des W-Seminars",
  supervisor: "Kursleitung",
  date: datetime.today(),
  doc,
) = {
  set document(title: title, author: author)

  set page(
    paper: "a4",
    margin: (left: 3cm, right: 2cm, top: 2.5cm, bottom: 2.5cm),
    numbering: "1",
    number-align: center,
  )

  // Schriftart bewusst nicht erzwungen — Typst nimmt seine mitgelieferte
  // Standardschrift. Falls eine bestimmte Schrift gefordert ist, hier
  // ergänzen, z. B.: font: "Libertinus Serif"
  set text(lang: "de", size: 11pt)
  set par(justify: true, leading: 0.65em, first-line-indent: 0em)
  set heading(numbering: "1.1")
  set math.equation(numbering: "(1)")

  // Überschriften 1. Ordnung beginnen auf neuer Seite
  show heading.where(level: 1): it => {
    pagebreak(weak: true)
    v(0.5em)
    text(size: 18pt, weight: "bold", it)
    v(0.3em)
    line(length: 100%, stroke: 0.6pt + luma(170))
    v(1em)
  }
  show heading.where(level: 2): it => {
    v(1.1em)
    text(size: 13pt, weight: "bold", it)
    v(0.45em)
  }
  show heading.where(level: 3): it => {
    v(0.85em)
    text(size: 11.5pt, weight: "bold", it)
    v(0.3em)
  }

  // ---------------- Tabellen: schlichter Buchsatz ----------------
  // Keine Gitterlinien, nur eine Linie unter der Kopfzeile; Kopfzeile hellgrau
  // hinterlegt und fett. Gilt automatisch für alle #table(...) im Dokument.
  set table(
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt + luma(100)) },
    fill: (x, y) => if y == 0 { luma(243) },
    inset: (x: 7pt, y: 5.5pt),
  )
  show table.cell.where(y: 0): set text(weight: "bold", size: 9.5pt)
  show table: set text(size: 9.5pt)

  // ---------------- Abbildungen und Bildunterschriften ----------------
  show figure: set block(above: 1.4em, below: 1.4em)
  show figure.caption: set text(size: 9.5pt, fill: luma(60))

  // ---------------- Codeblöcke ----------------
  show raw.where(block: true): it => block(
    fill: luma(247),
    stroke: 0.5pt + luma(215),
    radius: 3pt,
    inset: (x: 9pt, y: 8pt),
    width: 100%,
    text(size: 8.5pt, it),
  )

  // ---------------- Titelseite ----------------
  align(center)[
    #v(1.5cm)
    #text(size: 13pt)[#school]
    #v(0.3cm)
    #text(size: 12pt)[#course — #seminar-topic]
    #v(3cm)
    #text(size: 22pt, weight: "bold")[#title]
    #if subtitle != none [
      #v(0.4cm)
      #text(size: 14pt, style: "italic")[#subtitle]
    ]
    #v(3cm)
    #text(size: 13pt)[Seminararbeit von #author]
    #v(0.5cm)
    #text(size: 12pt)[Kursleitung: #supervisor]
    #v(1fr)
    #text(size: 11pt)[#date.display("[day].[month].[year]")]
    #v(1cm)
  ]
  pagebreak()

  // ---------------- Inhaltsverzeichnis ----------------
  outline(title: "Inhaltsverzeichnis", indent: auto)
  pagebreak()

  doc
}

// -----------------------------------------------------------
// Hilfsfunktionen für effizientes Arbeiten am Entwurf
// -----------------------------------------------------------

// Platzhalter für noch zu schreibenden Fließtext.
// Erscheint als gelber Kasten — vor der Abgabe alle #todo[...] auflösen.
// Tipp: im Editor nach "TODO" suchen, um offene Stellen zu finden.
#let todo(body) = block(
  fill: yellow.lighten(70%),
  inset: 8pt,
  radius: 3pt,
  stroke: 0.5pt + yellow.darken(30%),
  width: 100%,
)[#text(weight: "bold", fill: yellow.darken(40%))[TODO: ] #body]

// Kurze Randnotiz / Erinnerung (kursiv, rot), nicht Teil des finalen Texts
#let note(body) = text(fill: red.darken(10%), style: "italic")[Notiz: #body]

// Platzhalter für eine noch einzufügende Abbildung
#let figplaceholder(caption: "Abbildung folgt") = figure(
  block(
    width: 100%,
    height: 4cm,
    fill: luma(240),
    stroke: 0.5pt + luma(180),
    radius: 3pt,
    align(center + horizon)[#text(fill: luma(120))[Platzhalter für Abbildung]],
  ),
  caption: caption,
)
