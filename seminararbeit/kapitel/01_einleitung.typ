// Eigenes Import nötig: #include teilt keine Variablen mit main.typ.
#import "../template.typ": *

// ---------------------------------------------------------------
// Zitierkonvention dieser Arbeit (durchgängig anwenden):
//   1. Verweis auf ein Werk als Ganzes -> kein Zusatz, nur @Key.
//      Beispiel: der sechsseitige Aufsatz @Yee1966.
//   2. Verweis auf eine konkrete Aussage in einem umfangreichen Werk
//      -> Fundstelle im Zusatz: @Key[S. 51] bzw. @Key[S. 186 f.].
//      Enthält der Zusatz Kommata oder Klammern, die ausgeschriebene
//      Form nehmen: #cite(<Key>, supplement: [S. 51, Gl. (2.7)]).
//   3. Quellen mit eigener Abschnittszählung statt Paginierung
//      (Maxwell 1865) -> [Abschn. 97].
//   4. Internetquellen -> keine Seitenangabe; das Abrufdatum steht
//      als urldate im Eintrag in literatur.bib.
//   5. Ein blanker @Key bezeichnet eine sinngemäße Übernahme;
//      wörtliche Übernahmen stehen in Anführungszeichen. Ein
//      vorangestelltes "vgl." wird nicht verwendet, da es zur
//      deutschen Fußnotenzitierweise gehört und nicht zum hier
//      gesetzten numerischen IEEE-Stil.
// ---------------------------------------------------------------

= Einleitung

== Relevanz und Motivation

Wer telefoniert, hält einen Sender wenige Zentimeter neben den Kopf. Ein Teil
der abgestrahlten Leistung dringt dabei in das Gewebe ein und wird dort in
Wärme umgesetzt. Wie groß dieser Anteil ist, beschreibt die spezifische
Absorptionsrate (SAR), also die je Masse absorbierte Leistung, gemittelt über
eine festgelegte Gewebemasse; am Kopf gilt in Deutschland ein Grenzwert von
2 W/kg @BfS_SAR. Belastbar ist ein solcher Grenzwert allerdings nur so weit wie
das Verfahren, mit dem die zugehörige Größe bestimmt wird.
#note[Mittelungsmasse (nach ICNIRP 10 g) an der Richtlinie prüfen und
anstelle von "eine festgelegte Gewebemasse" einsetzen.]

Ein Kopf ist nämlich kein homogener Körper. Haut, Schädelknochen, Liquor und
Hirngewebe unterscheiden sich in Permittivität und Leitfähigkeit deutlich
voneinander, ihre Grenzflächen verlaufen gekrümmt und unregelmäßig, und die
Antenne steht so nah, dass sich das einfallende Feld nicht als ebene Welle
behandeln lässt. Gesucht ist zudem die Feldstärke im Gewebe selbst, denn erst
aus ihr folgt die absorbierte Leistung.

Behandelt wird der Fall deshalb numerisch. Taflove und Hagness führen ihn als
Standardbeispiel der _Computational Electromagnetics_ an: Ein Mobiltelefon wird
bis auf Zellweiten von 0,1 mm aufgelöst und neben ein aus
Magnetresonanzaufnahmen gewonnenes Kopfmodell aus 15 Gewebetypen gesetzt, in
dem sich anschließend die Spitzenwerte der SAR bestimmen lassen
@Taflove2005[S. 11]. Ein Sonderfall ist das nicht: Ob Antennen hinter
dielektrischen Abdeckungen, Wellenleiter der integrierten Photonik oder Spulen
in Magnetresonanztomographen --- stets treffen Materialien unterschiedlicher
Permittivität und Leitfähigkeit auf unregelmäßigen Grenzflächen zusammen, und
stets muss das Feld gerade dort bestimmt werden.

An Konfigurationen dieser Art scheitert die analytische Behandlung --- nicht
aus Unkenntnis der Physik. Die klassische Elektrodynamik gilt als abgeschlossene
und in ihrem Geltungsbereich vollständige Theorie @GriffithsDE[S. 18]: Die
Maxwell-Gleichungen legen das Feld bei gegebenen Quellen, Materialeigenschaften
und Randbedingungen eindeutig fest @Maxwell1865[Abschn. 70]. Geschlossene
Lösungen sind jedoch nur für wenige Konfigurationen bekannt, in denen das
Material homogen und die Geometrie hochsymmetrisch ist. Die entscheidende
Grenze verläuft damit nicht zwischen bekannter und unbekannter Physik, sondern
zwischen lösbaren und unlösbaren Randwertproblemen (also
Differentialgleichungen mit vorgegebenen Randbedingungen).

Messen statt rechnen scheitert hier ebenfalls: Eine Sonde stört das Feld und
lässt sich in lebendem Gewebe ohnehin nicht anbringen, und jede
Entwurfsänderung erforderte ein neues gefertigtes Exemplar, sodass sich
Varianten weder in vertretbarer Zeit noch zu vernünftigen Kosten vergleichen
ließen. Rechnerisch zugänglich wird das Problem erst durch Diskretisierung,
also dadurch, dass das Kontinuum durch ein endliches Gitter ersetzt wird. Das
verbreitetste Verfahren dieser Art ist die _Finite-Difference
Time-Domain_-Methode (FDTD).

Durchgeführt werden solche Rechnungen fast ausschließlich mit kommerzieller
Software wie CST Studio Suite @CSTStudio, deren Quellcode nicht offenliegt; aus
Sicht der Anwendung liegt damit ein Blackbox-System vor. Methodisch bedeutsam
ist dabei, dass ein Simulationsergebnis die exakte Lösung eines Systems von
Differenzengleichungen darstellt und nicht die der Maxwell-Gleichungen; beide
stimmen erst im Grenzfall verschwindender Gitterweite überein. Auf einem
endlichen Gitter treten deshalb systematische Abweichungen auf, die sich nicht
ohne Weiteres von physikalischen Effekten unterscheiden lassen.

Wie groß diese Abweichungen ausfallen, bleibt in der Anwendung meist offen. Wie
schwer das wiegt, zeigt der eingangs beschriebene Fall: An den ausgegebenen
Zahlen hängt dort eine Aussage über die Belastung eines Menschen, und ob eine
berechnete SAR knapp unter oder über dem Grenzwert liegt, ist ohne eine Angabe
zur numerischen Genauigkeit nicht zu entscheiden. Die Lücke hat einen
strukturellen Grund: Für jene Konfigurationen, aufgrund derer simuliert wird,
existiert definitionsgemäß keine Vergleichslösung. Die vorliegende Arbeit kehrt
die übliche Herangehensweise daher um und behandelt die analytisch lösbaren
Fälle nicht als Trivialfälle, sondern als die einzigen Anordnungen, in denen
das wahre Feld bekannt ist und sich das Verfahren deshalb überhaupt prüfen
lässt.

== Historischer Kontext

James Clerk Maxwell führte 1865 in _A Dynamical Theory of the Electromagnetic
Field_ Elektrizität, Magnetismus und Optik in einem einheitlichen
Gleichungssystem zusammen und sagte die Existenz sich mit Lichtgeschwindigkeit
ausbreitender elektromagnetischer Wellen voraus @Maxwell1865[Abschn. 97]; der
experimentelle Nachweis gelang Heinrich Hertz rund zwanzig Jahre später. Damit
wurde deutlich, "dass das sichtbare Licht nur ein winziges Fenster im breiten
Spektrum der elektromagnetischen Strahlung darstellt, das sich vom Radiobereich
über Mikrowellen, Infrarot und Ultraviolett bis zum Röntgen- und Gammabereich
erstreckt" @GriffithsDE[S. 19] --- die Strahlung eines Mobiltelefons unterliegt
also denselben Gleichungen wie das Licht einer Lampe. Dieser Rahmen gilt
unverändert und bildet die Grundlage der vorliegenden Arbeit.

// Bewusst ohne Seitenzusatz (Regel 1 oben): der Verweis gilt dem sechs
// Seiten kurzen Aufsatz als Ganzem, dessen Umfang S. 302--307 bereits im
// Literaturverzeichnis steht.
Kane S. Yee legte 1966 mit einer knappen Veröffentlichung @Yee1966 die
Grundlage der numerischen Umsetzung. Er ordnet die elektrischen und
magnetischen Feldkomponenten nicht gemeinsam an, sondern räumlich und zeitlich
gegeneinander versetzt auf einem sogenannten Yee-Gitter. So lassen sich die
gekoppelten Rotationsgleichungen durchgängig mit zentralen
Differenzenquotienten approximieren, und beide Feldgrößen können alternierend
fortgeschrieben werden (_leapfrog_-Schema). Die damals verfügbare
Rechenleistung beschränkte die Anwendung allerdings auf sehr kleine
Modellprobleme.

Praktisch nutzbar wurde das Verfahren erst mit wachsenden Rechenkapazitäten.
Allen Taflove führte die Bezeichnung FDTD ein und baute es zu einem
Standardwerkzeug der _Computational Electromagnetics_ aus, insbesondere durch
absorbierende Randbedingungen und die Erweiterung auf verlustbehaftete
Materialien. Getragen wurde diese Entwicklung zunächst von militärischen
Anforderungen; erst seit etwa 1990 verschiebt sich das Feld zu kommerziellen
Anwendungen in Kommunikation und Rechentechnik @Taflove2005[S. 1] --- zu denen
das eingangs beschriebene Mobiltelefon zählt.
#note[Jahr der Namensprägung und Seitenzahl des historischen Überblicks in
Kapitel 1 von Taflove und Hagness am Buch prüfen.]

== Persönliche Motivation

Für Simulationen habe ich mich erstmals in der zehnten Jahrgangsstufe
interessiert, als ich mich mit der numerischen Beschreibung mechanischer Wellen
beschäftigt habe. In der elften kam die Programmierseite hinzu: Im Rahmen des
Astro-Pi-Wettbewerbs der ESA#footnote[Die European Astro Pi Challenge ist ein
Bildungsprogramm der ESA gemeinsam mit der Raspberry Pi Foundation; die
eingereichten Python-Programme laufen auf zwei Raspberry-Pi-Einheiten an Bord
der Internationalen Raumstation @AstroPi.] habe ich eigenen Code geschrieben
und dabei erfahren, wie viel Sorgfalt eine Aufgabenstellung verlangt, deren
Ergebnis nicht im eigenen Ermessen liegt. Die vorliegende Arbeit führt beides
fort: von der mechanischen zur elektromagnetischen Welle und von der
physikalischen Beschreibung zur eigenen Implementierung.

Hinzu kommt mein Interesse an technischen Anwendungen. In der Schule und in der
theoretischen Physik werden Szenarien weitgehend idealisiert: homogene Medien,
unendlich ausgedehnte Grenzflächen, punktförmige Quellen. Didaktisch ist das
notwendig, denn nur so werden die Gesetzmäßigkeiten sichtbar; zugleich entsteht
der Eindruck, physikalische Probleme seien stets geschlossen lösbar, obwohl die
Idealisierung gerade dort entfällt, wo es interessant wird. Ein Verfahren wie
FDTD setzt genau dort an, und diese Arbeit bot mir die Gelegenheit, einmal ein
Problem zu bearbeiten, das sich nicht durch Vereinfachung auflösen lässt.

== Zentrale Fragestellung und Aufbau der Arbeit

Aus den dargestellten Überlegungen ergibt sich die folgende Leitfrage:

#align(center)[
  #block(width: 92%)[
    _Wie zuverlässig reproduziert eine selbst implementierte zweidimensionale
    FDTD-Simulation die analytisch bekannten Eigenschaften elektromagnetischer
    Wellen, und wovon hängt diese Zuverlässigkeit ab?_
  ]
]

Um sie zu beantworten, wird ein eigener, vollständig dokumentierter FDTD-Löser
in Python implementiert und an drei Konfigurationen geprüft, deren Lösung sich
in geschlossener Form angeben lässt: an der ebenen Welle im freien Raum, an der
Reflexion an einer ebenen Grenzfläche und am rechteckigen Hohlraumresonator.
Sie sind so gewählt, dass jede einen anderen der drei Diskretisierungsschritte
in den Vordergrund rückt --- die des Raums, die des Materials und die der
Zeit ---, sodass sich im Fehlerfall benennen lässt, welcher versagt hat.
Gegenstand ist damit das Verfahren und nicht ein bestimmtes Bauteil: Das
Ergebnis sind keine physikalischen Erkenntnisse, sondern Fehlergrenzen in
Abhängigkeit der Diskretisierung samt der Herkunft der verbleibenden
Abweichungen.

Aus der Leitfrage ergeben sich zwei Teilfragen:

+ Mit welcher relativen Abweichung gibt die Simulation die analytisch bekannten
  Größen wieder --- Ausbreitungsgeschwindigkeit und Wellenlänge im freien Raum,
  Reflexionsverhalten an einer Grenzfläche, Resonanzfrequenzen eines Hohlraums?
+ Wie hängen diese Abweichungen von den Diskretisierungsparametern ab,
  insbesondere von der Anzahl der Gitterzellen pro Wellenlänge, und welche
  prinzipiellen Grenzen des Verfahrens folgen daraus?

Kapitel 2 stellt dazu die theoretischen Voraussetzungen bereit, Kapitel 3
entwickelt daraus die FDTD-Methode bis zu den Update-Gleichungen des
zweidimensionalen $"TM"_z$-Falls, Kapitel 4 dokumentiert die Implementierung
und Kapitel 5 die drei Testfälle. Kapitel 6 führt die Ergebnisse zusammen und
gibt einen Ausblick auf Problemstellungen ohne analytische Lösung --- etwa auf
die eingangs beschriebene Absorption im Kopfgewebe.
