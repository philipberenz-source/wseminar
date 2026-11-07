// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Quellcode und Website <anhang-code>

Sämtlicher Quellcode dieser Arbeit ist im folgenden GitHub-Repository
abgelegt:

#align(center)[#link("https://github.com/philipberenz-source/wseminar")[`https://github.com/philipberenz-source/wseminar`]]

Es enthält den in Kapitel 4 beschriebenen Solver `fdtd_core.py` sowie die
Skripte zu allen drei Testfällen aus Kapitel 5, mit denen sämtliche Messwerte
und Diagramme erzeugt wurden (Ordner `code snippets`).

Ergänzend zur Arbeit ist eine Website zum W-Seminar entstanden, auf der der
Solver direkt im Browser visuell genutzt werden kann:

#align(center)[#link("https://glowing-crisp-3f1e8f.netlify.app")[`https://glowing-crisp-3f1e8f.netlify.app`]]

#heading(numbering: none, outlined: true)[Erklärung]

Ich erkläre hiermit, dass ich die Seminararbeit selbstständig angefertigt
und nur die im Literaturverzeichnis angeführten Quellen und Hilfsmittel benutzt
habe.

Als KI-basiertes Hilfsmittel wurde das Sprachmodell Opus 5 des Unternehmens Anthropic in folgendem Umfang verwendet:

- zur sprachlichen Überarbeitung des Textes und zur Umsetzung einzelner,
  von mir geprüfter Korrekturvorschläge im Text,
- zum Erstellen eines Layouts in Typst,
- zur Erstellung der begleitenden W-Seminar-Website, auf der der Solver visuell
  genutzt werden kann,
- zur Erstellung der als KI-generiert gekennzeichneten Abbildungen in Kapitel 3 und 5
Die durch das Hilfsmittel erzeugten Inhalte habe ich geprüft; für die
inhaltliche Richtigkeit der gesamten Arbeit trage ich die Verantwortung.

#v(2.5cm)

#grid(
  columns: (1fr, 1fr),
  column-gutter: 1.5cm,
  [#line(length: 100%, stroke: 0.5pt) #text(size: 9.5pt)[Ort, Datum]],
  [#line(length: 100%, stroke: 0.5pt) #text(size: 9.5pt)[Unterschrift]],
)
