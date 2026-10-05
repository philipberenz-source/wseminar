#set page(paper: "a4", margin: 2.2cm)
#set text(lang: "de", size: 11pt)
#set par(justify: true)
#set math.equation(numbering: "(1)")

#align(center, text(14pt, weight: "bold")[Herleitung der Näherung für $tilde(v)_p \/ c$ (Gl. 60)])
#v(4pt)

*1. Phasengeschwindigkeit im Gitter.* Im Kontinuum gilt $c = omega \/ k$. Im Gitter
ist die Frequenz $omega$ dieselbe, da sie von der Quelle vorgegeben wird; die
Wellenzahl ist jedoch $tilde(k)$ statt $k$. Damit folgt

$ tilde(v)_p / c = (omega \/ tilde(k)) / (omega \/ k) = k / tilde(k) . $ <vp>

*2. Ansatz.* Für feine Auflösung weicht $tilde(k)$ nur wenig von $k$ ab. Mit einer
kleinen Größe $delta$ wird gesetzt

$ tilde(k) = k (1 + delta) quad => quad
  tilde(v)_p / c = 1 / (1 + delta) approx 1 - delta . $ <ansatz>

Mit $k = 2 pi \/ lambda_0$ und $Delta = lambda_0 \/ N_lambda$ gilt außerdem

$ (k Delta) / 2 = pi / N_lambda =: u
  quad => quad
  (tilde(k) Delta) / 2 = u (1 + delta) . $ <u>

*3. Dispersionsrelation.* Ausgangspunkt ist die numerische Dispersionsrelation
(Gl. 59):

$ 1/S_c^2 sin^2(S_c u)
  = sin^2(u (1 + delta) cos theta) + sin^2(u (1 + delta) sin theta) . $ <disp>

*4. Taylor-Entwicklung.* Aus $sin x approx x - x^3\/6$ folgt
$sin^2 x approx x^2 - x^4\/3$. Für die linke Seite ergibt sich

$ 1/S_c^2 sin^2(S_c u) approx u^2 - (S_c^2 u^4)/3 , $ <links>

für die rechte Seite mit $cos^2 theta + sin^2 theta = 1$ und unter Vernachlässigung
kleiner Terme höherer Ordnung ($delta^2$, $delta u^4$)

$ u^2 (1 + 2 delta) - u^4/3 (cos^4 theta + sin^4 theta) . $ <rechts>

*5. Gleichsetzen.* Setzt man @links und @rechts gleich, hebt sich $u^2$ auf beiden
Seiten weg:

$ - (S_c^2 u^4)/3 = 2 delta u^2 - u^4/3 (cos^4 theta + sin^4 theta) . $

Division durch $2 u^2$ liefert

$ delta = u^2/6 (cos^4 theta + sin^4 theta - S_c^2)
        = pi^2/(6 N_lambda^2) (cos^4 theta + sin^4 theta - S_c^2) . $ <delta>

*6. Ergebnis.* Eingesetzt in @ansatz folgt Gl. 60:

#block(stroke: 0.6pt, inset: 8pt, width: 100%)[
$ tilde(v)_p / c approx 1 - pi^2/(6 N_lambda^2)
    (cos^4 theta + sin^4 theta - S_c^2) . $
]
