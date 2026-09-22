sowie ein entsprechender Ansatz mit den Amplituden $H_(x 0)$ und $H_(y 0)$ für
$H_x$ und $H_y$ eingesetzt. Jede zentrale Differenz verknüpft dann zwei Werte
des Ansatzes, die eine halbe Schrittweite vor und nach dem Bezugspunkt liegen
und sich nur um den Faktor $e^(plus.minus upright(i) omega Delta t\/2)$
unterscheiden. Für die Zeitableitung in @update-hy folgt etwa

$ (H_y bar.v_(i+1/2\,j)^(n+1/2) - H_y bar.v_(i+1/2\,j)^(n-1/2))/(Delta t)
  = H_y bar.v_(i+1/2\,j)^n
    (e^(upright(i) omega Delta t\/2) - e^(-upright(i) omega Delta t\/2))/(Delta t)
  = (2 upright(i))/(Delta t) sin((omega Delta t)/2) dot H_y bar.v_(i+1/2\,j)^n , $

wobei $H_y bar.v^n$ den Wert des Ansatzes zum Zeitpunkt $n$ bezeichnet. An die
Stelle des Faktors $upright(i) omega$, den die Zeitableitung im Kontinuum
liefert, tritt also $(2 upright(i)\/Delta t) sin(omega Delta t\/2)$; im Ort wird
entsprechend aus $-upright(i) tilde(k)_x$ der Faktor
$-(2 upright(i)\/Delta) sin(tilde(k)_x Delta\/2)$. Nach dem Kürzen des
gemeinsamen Exponentialfaktors bleiben aus den drei Update-Gleichungen drei
algebraische Gleichungen für die Amplituden @Taflove2005[S. 109, Gl. (4.4a--c),
mit $Delta x = Delta y = Delta$]:

$ H_(x 0) = (Delta t thin E_(z 0))/(mu_0 Delta)
  dot sin(tilde(k)_y Delta\/2)/sin(omega Delta t\/2) , $ <tf1-hx0>

$ H_(y 0) = -(Delta t thin E_(z 0))/(mu_0 Delta)
  dot sin(tilde(k)_x Delta\/2)/sin(omega Delta t\/2) , $ <tf1-hy0>

$ E_(z 0) sin((omega Delta t)/2)
  = (Delta t)/(epsilon_0 Delta) [H_(x 0) sin((tilde(k)_y Delta)/2)
    - H_(y 0) sin((tilde(k)_x Delta)/2)] . $ <tf1-ez0>

@tf1-hx0 und @tf1-hy0 folgen aus den Update-Gleichungen für $H_x$ und $H_y$,
@tf1-ez0 aus der für $E_z$. Setzt man @tf1-hx0 und @tf1-hy0 in @tf1-ez0 ein,
kürzt sich $E_(z 0)$ heraus, und mit $c = 1\/sqrt(mu_0 epsilon_0)$ ergibt sich
für das quadratische Gitter die _numerische Dispersionsrelation_
@Taflove2005[S. 109, Gl. (4.5)]

$ [1/(c Delta t) sin((omega Delta t)/2)]^2
  = [1/Delta sin((tilde(k)_x Delta)/2)]^2
  + [1/Delta sin((tilde(k)_y Delta)/2)]^2 . $ <num-dispersion>

