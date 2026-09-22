
Aussagekräftig ist das Fehlerverhalten vor allem in Abhängigkeit von zwei
Größen. Nach @Schneider2010[S. 113] bestimmen nicht die absoluten Schrittweiten
$Delta$ und $Delta t$ das Verhalten der Felder im Gitter, sondern allein die
Courant-Zahl $S$, die nur ihr Verhältnis enthält, und die _Auflösung_

$ N_lambda = lambda_0/Delta = c/(f Delta), $ <nlambda>

also die Zahl der Zellen je Wellenlänge. Ein in diesen beiden Größen angegebenes
Ergebnis gilt somit für jede Frequenz und jede Gittergröße gleichermaßen; das
Fehlerverhalten wird deshalb in allen Testfällen über $N_lambda$ und $S$
untersucht.

Damit dafür nicht für jedes $N_lambda$ ein eigener Lauf nötig ist, wird in allen
Testfällen mit einem Puls statt mit einer harmonischen Quelle angeregt
@Schneider2010[S. 118]. Verwendet wird der _Ricker-Puls_, die zweite Ableitung
einer Gaußkurve @Schneider2010[S. 118, Gl. (5.20)]. Er enthält ein breites
Frequenzspektrum, im Gegensatz zum Gaußpuls aber keinen Gleichanteil, der im
Gitter unphysikalische Artefakte erzeugen kann @Schneider2010[S. 118]. Da nach
@nlambda jeder Frequenz eine Auflösung entspricht, liefert ein einziger Lauf
nach der Fourier-Transformation aus @sec-auswertung die Ergebnisse für einen
ganzen Bereich von $N_lambda$ zugleich.
