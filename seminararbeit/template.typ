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
    numbering: none,
    number-align: center,
  )

  set text(lang: "de", size: 11pt)
  set par(justify: true, leading: 0.65em, first-line-indent: 0em)
  set heading(numbering: "1.1")
  set math.equation(numbering: n => numbering("(1.1)", counter(heading).get().first(), n))

  show heading.where(level: 1): it => {
    counter(math.equation).update(0)
    pagebreak(weak: true)
    v(0.5em)
    text(size: 18pt, weight: "bold", it)
    v(0.3em)
    line(length: 100%, stroke: 0.6pt + luma(170))
    v(1em)
  }
  let needspace(min) = {
    block(breakable: false, height: min, above: 0pt, below: 0pt)
    v(-min)
  }
  show heading.where(level: 2): it => {
    needspace(6.5cm)
    v(1.1em)
    block(sticky: true, text(size: 13pt, weight: "bold", it))
    v(0.45em)
  }
  show heading.where(level: 3): it => {
    needspace(4.5cm)
    v(0.85em)
    block(sticky: true, text(size: 11.5pt, weight: "bold", it))
    v(0.3em)
  }

  set text(costs: (orphan: 100%, widow: 100%))

  set table(
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt + luma(100)) },
    fill: (x, y) => if y == 0 { luma(243) },
    inset: (x: 7pt, y: 5.5pt),
  )
  show table.cell.where(y: 0): set text(weight: "bold", size: 9.5pt)
  show table: set text(size: 9.5pt)

  show figure: set block(above: 1.4em, below: 1.4em)
  show figure.caption: set text(size: 9.5pt, fill: luma(60))

  show raw.where(block: true): it => block(
    fill: luma(247),
    stroke: 0.5pt + luma(215),
    radius: 3pt,
    inset: (x: 6pt, y: 5pt),
    width: 100%,
    text(size: 7.4pt, it),
  )

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

  outline(title: "Inhaltsverzeichnis", indent: auto, depth: 2)
  pagebreak()

  set page(numbering: "1")
  counter(page).update(1)

  doc
}
