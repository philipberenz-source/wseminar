# Arbeiten mit dieser Typst-Vorlage in VS Code (Tinymist)

## Wichtigster Punkt: Hauptdatei anpinnen

Tinymist kompiliert standardmäßig **die gerade geöffnete Datei**. Wenn du
`kapitel/03_fdtd_methode.typ` offen hast, wird nur diese Datei allein
gebaut — ohne Titelseite, ohne Literaturverzeichnis. Deshalb schlagen
Zitate wie `@Maxwell1865` dort fehl ("failed to resolve reference"), und
der Compiler bricht ab.

**Lösung, einmal pro VS-Code-Sitzung:**

1. `main.typ` öffnen
2. Befehlspalette: `Cmd+Shift+P`
3. `Typst Pin Main` eingeben und ausführen (Befehl: `tinymist.pinMainToCurrent`)

Ab dann wird immer `main.typ` kompiliert, egal welche Kapitel-Datei du
gerade bearbeitest. Die Vorschau (`Cmd+K V`) zeigt dann das ganze Dokument.

Hinweis: Tinymist merkt sich das Pinning nicht über einen Neustart hinweg
— nach dem Neuöffnen von VS Code also erneut ausführen. Wenn dich das
stört, kannst du in den VS-Code-Einstellungen `tinymist.typstExtraArgs`
auf `["main.typ"]` setzen.

## Fehler eingrenzen

Falls trotzdem etwas nicht rendert: `test.typ` öffnen und kompilieren.
Die Datei testet isoliert nur drei Dinge — die `#todo`-Funktion, das
Zitieren mit `@Yee1966`, und Mathematik. So lässt sich sagen, ob das
Problem am Template oder an einer Kapitel-Datei liegt.

## Dateistruktur

| Datei | Zweck |
|---|---|
| `main.typ` | Einstiegsdatei — **diese** kompilieren/anpinnen |
| `template.typ` | Layout, Titelseite, Hilfsfunktionen |
| `kapitel/*.typ` | je ein Kapitel; jede Datei importiert `template.typ` selbst |
| `literatur.bib` | Quellen, zitieren im Text mit `@Yee1966` |
| `test.typ` | Diagnose-Datei, gehört nicht zur Arbeit |

## Wichtig zu wissen

`#include` teilt in Typst **keine** Variablen. Jede Kapitel-Datei braucht
deshalb ihre eigene Zeile `#import "../template.typ": *` ganz oben —
sonst kennt sie `#todo` und `#note` nicht. Wenn du eine neue Kapitel-Datei
anlegst, diese Zeile nicht vergessen.

## Offene Stellen finden

Im Editor nach `TODO` suchen (`Cmd+Shift+F`). Jeder gelbe Kasten im PDF
ist eine noch zu schreibende Stelle. Vor der Abgabe müssen alle `#todo[]`
und `#note[]` verschwunden sein.
