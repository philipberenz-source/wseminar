#set page(paper: "a4", margin: 2.2cm)
#set text(lang: "de", size: 11pt)
#set par(justify: true)
#set math.equation(numbering: "(1)")

#align(center, text(14pt, weight: "bold")[Herleitung der Näherung für $abs(tilde(Gamma)) \/ abs(Gamma)$ (Gl. 69)])
#v(4pt)

Der Index $ell$ steht für das Medium: $ell = 1$ Vakuum ($n_1 = 1$), $ell = 2$
Dielektrikum ($n_2 = sqrt(epsilon_r)$).

*1. Argument der Kosinusfaktoren.* Ausgangspunkt ist die Dispersionsrelation
(Gl. 59 mit $theta = 0$ und $c -> c\/n_ell$):

$ sin((tilde(k)_ell Delta)/2) = n_ell/S_c sin((pi S_c)/N_lambda) . $ <disp>

Für feine Auflösung sind beide Argumente klein, sodass $sin x approx x$ auf beiden
Seiten gilt:

$ (tilde(k)_ell Delta)/2 approx n_ell/S_c dot (pi S_c)/N_lambda = (pi n_ell)/N_lambda . $ <arg>

Die Courant-Zahl $S_c$ kürzt sich dabei heraus.

*2. Kosinus nähern.* Mit $cos x approx 1 - x^2\/2$ und der Abkürzung
$alpha = pi^2 \/ (2 N_lambda^2)$ folgt aus @arg

$ cos((tilde(k)_ell Delta)/2) approx 1 - n_ell^2 alpha . $ <cos>

*3. Einsetzen in $tilde(Gamma)$.* Mit @cos werden Zähler und Nenner von Gl. 68
ausmultipliziert und zusammengefasst:

$ "Zähler:" quad
  n_1 (1 - n_2^2 alpha) - n_2 (1 - n_1^2 alpha)
  &= (n_1 - n_2) + n_1 n_2 alpha (n_1 - n_2) \
  &= (n_1 - n_2)(1 + n_1 n_2 alpha) , $

$ "Nenner:" quad
  n_1 (1 - n_2^2 alpha) + n_2 (1 - n_1^2 alpha)
  &= (n_1 + n_2) - n_1 n_2 alpha (n_1 + n_2) \
  &= (n_1 + n_2)(1 - n_1 n_2 alpha) . $

*4. Bruch bilden.* Der Vorfaktor $(n_1 - n_2)\/(n_1 + n_2)$ ist genau der
Fresnel-Wert $Gamma$. Mit $y = n_1 n_2 alpha$ gilt

$ tilde(Gamma) = Gamma dot (1 + y)/(1 - y) . $ <bruch>

*5. Kleine Terme vernachlässigen.* Für kleines $y$ gilt $1\/(1 - y) approx 1 + y$;
das Glied $y^2$ entfällt:

$ (1 + y)/(1 - y) approx (1 + y)(1 + y) approx 1 + 2 y
  = 1 + (n_1 n_2 pi^2)/N_lambda^2 . $

*6. Ergebnis.* Da $tilde(Gamma)$ und $Gamma$ dasselbe Vorzeichen haben, folgt Gl. 69:

#block(stroke: 0.6pt, inset: 8pt, width: 100%)[
$ abs(tilde(Gamma)) / abs(Gamma) approx 1 + (n_1 n_2 pi^2)/N_lambda^2 . $
]

Die Näherung gilt für feine Auflösung. Bei $epsilon_r = 4$ liefert sie für
$N_lambda = 30$ eine Abweichung von 2,2 % (exakt nach Gl. 68: 2,2 %), für
$N_lambda = 10$ dagegen 19,7 % (exakt: 25,1 %).
