# Arbeiten mit dieser Typst-Vorlage in VS Code (Tinymist)

## Wichtigster Punkt: Hauptdatei anpinnen

Tinymist kompiliert standardmäßig **die gerade geöffnete Datei**. Wenn du
`kapitel/03_fdtd_methode.typ` offen hast, wird nur diese Datei allein
gebaut — ohne Titelseite, ohne Literaturverzeichnis. Deshalb schlagen
Zitate wie `@Maxwell1865` dort fehl ("failed to resolve reference"), und
jedes `@Zitat` wird im Editor rot unterkringelt.

**Das ist bereits fest eingestellt:** `.vscode/settings.json` im
Workspace-Root setzt `tinymist.typstExtraArgs` auf
`["seminararbeit/main.typ"]`, sodass Tinymist immer `main.typ` kompiliert
— unabhängig davon, welche Kapitel-Datei gerade offen ist. Damit
verschwinden die roten Fehler unter den Zitaten dauerhaft, ganz ohne
manuelles Pinnen.

Falls die rote Unterkringelung trotzdem auftaucht (z. B. nach einem
Extension-Update): VS Code neu laden (`Cmd/Ctrl+Shift+P` →
`Developer: Reload Window`). Alternativ manuell pinnen: `main.typ`
öffnen, Befehlspalette, `Typst Pin Main` (`tinymist.pinMainToCurrent`) —
gilt dann aber nur bis zum nächsten Neustart von VS Code.

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