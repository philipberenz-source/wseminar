// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

= Die FDTD-Methode

== Grundidee: Finite-Differenzen-Approximation

#todo[
  Zentrale Differenzenquotienten 2. Ordnung als Ersatz für partielle
  Ableitungen einführen und den Diskretisierungsfehler benennen.
]

$ ((partial f)/(partial x))_i approx (f_(i+1/2) - f_(i-1/2)) / (Delta x) $

Diskretisierungsfehler: $O(Delta x^2)$.

== Das Yee-Gitter (2D)

#todo[
  Räumlich um eine halbe Zelle versetzte Anordnung von $E_z$, $H_x$, $H_y$
  ($"TM"_z$-Fall) bzw. $H_z$, $E_x$, $E_y$ ($"TE"_z$-Fall) erläutern. Begründen:
  zentrierte Differenzen für die Rotation werden ohne Interpolation exakt,
  die Divergenzbedingung bleibt automatisch erhalten, wenn sie zu $t=0$
  erfüllt ist.
]

#figplaceholder(caption: "Yee-Gitter im 2D-Fall (TM_z): Anordnung von E_z, H_x, H_y")

== Herleitung der Update-Gleichungen (2D, $"TM"_z$-Fall)

#todo[
  Ausgangspunkt: die drei relevanten Komponenten von Faraday- und
  Ampère-Gesetz für $E_z, H_x, H_y$ herleiten. Explizite Diskretisierung mit
  Leapfrog-Zeitschema (E auf ganzzahligen, H auf halbzahligen Zeitschritten)
  Schritt für Schritt zeigen.
]

$ H_x^(n+1/2) = H_x^(n-1/2) - (Delta t)/(mu Delta y) (E_z^n (i,j+1) - E_z^n (i,j)) $

$ H_y^(n+1/2) = H_y^(n-1/2) + (Delta t)/(mu Delta x) (E_z^n (i+1,j) - E_z^n (i,j)) $

$ E_z^(n+1) = E_z^n + (Delta t)/epsilon dot (
  (H_y^(n+1/2)(i,j) - H_y^(n+1/2)(i-1,j)) / (Delta x)
  - (H_x^(n+1/2)(i,j) - H_x^(n+1/2)(i,j-1)) / (Delta y)
) $

#note[
  Erweiterung um einen $sigma E_z$-Verlustterm ergänzen, falls in Kapitel 6
  eine Leitfähigkeit $sigma > 0$ (Wände, Möbel) verwendet wird — dann auch
  hier die Update-Gleichung entsprechend anpassen und im Text vermerken.
]

== Modellannahmen und Einschränkungen

#todo[
  Für jeden Punkt kurz begründen, warum die Annahme nötig ist und was sie
  ausschließt.
]

- *2D-Reduktion:* Translationsinvarianz in $z$-Richtung, Trennung in
  $"TM"_z$-/$"TE"_z$-Moden, keine $z$-Feldkomponenten-Kopplung.
- *Lineare Medien:* $epsilon, mu$ unabhängig von der Feldstärke (keine
  nichtlinearen Effekte wie der Kerr-Effekt).
- *Isotrope Medien:* $epsilon, mu$ skalar statt tensoriell — keine
  Richtungsabhängigkeit (schließt z. B. Kristalle mit anisotroper
  Permittivität aus).
- *Nicht-frequenzabhängige (nichtdispersive) Medien:* $epsilon_r, mu_r$
  konstant über das gesamte im Puls enthaltene Frequenzband — kein
  Debye- oder Lorentz-Modell für reale, frequenzabhängige Materialien
  (z. B. feuchte Baustoffe).

== Stabilität: die Courant-Friedrichs-Lewy-Bedingung (CFL)

#todo[
  Anschauliche Herleitung: Die numerische Ausbreitungsgeschwindigkeit im
  Gitter darf die physikalische nicht unterschreiten können, d. h. eine
  Information darf pro Zeitschritt höchstens eine Gitterzelle „überspringen“.
  Formale Herleitung optional über von-Neumann-Stabilitätsanalyse (Ansatz
  $E prop e^(i(bold(k) dot bold(r) - omega t))$ in die Update-Gleichungen
  einsetzen, Bedingung für $|"Verstärkungsfaktor"| <= 1$).
]

Ergebnis für das 2D-Gitter:

$ Delta t <= 1 / (c sqrt(1/(Delta x)^2 + 1/(Delta y)^2)) $

Für $Delta x = Delta y$ vereinfacht sich dies zu:

$ Delta t <= (Delta x) / (c sqrt(2)) $

Courant-Zahl: $S = (c Delta t) / (Delta x)$.

== Randbedingungen

#todo[
  Problem der endlichen Simulationsdomäne einführen: unphysikalische
  Reflexionen am Gebietsrand. Verwendetes Verfahren benennen (z. B.
  Mur-Randbedingung 1. Ordnung oder PEC-Rand, je nach Testfall aus
  Kapitel 5). Verweis auf die genauere PML (Perfectly Matched Layer) als
  Standardverfahren setzen, siehe Ausblick (Kapitel 7).
]
