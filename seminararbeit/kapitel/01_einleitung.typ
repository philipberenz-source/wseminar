#import "../template.typ": *

= Einleitung

#place(hide[
  @lit:Maxwell1865 @lit:Hertz1889 @lit:Yee1966 @lit:Berenger1994
  @lit:ChewWeedon1994 @lit:Liu2022
  @lit:Taflove2005 @lit:GriffithsDE @lit:Orosz2022 @lit:Schneider2010 @lit:Teixeira2023
  @lit:Oskooi2010 @lit:Warren2016 @lit:AppleSpecs @lit:AnsysLumerical
])

== Relevanz und Motivation

Als Apple im September 2023 das iPhone 15 Pro vorstellte, dessen Antennen in den umlaufenden Titanrahmen des Gehäuses integriert sind, musste das Antennendesign auf Gehäuse und Innenaufbau abgestimmt sein. Die Schwierigkeit liegt dabei in den Multiband-Anforderungen: Auf wenigen Zentimetern Gehäuserand müssen Mobilfunkbänder von etwa 600 MHz bis 6 GHz, WLAN bei 2,4, 5 und 6 GHz, Bluetooth, GNSS und NFC bei 13,56 MHz gleichzeitig untergebracht werden @lit:AppleSpecs. Ein Teil der abgestrahlten Leistung wird dabei nicht nach außen abgegeben, sondern an den inneren Grenzflächen von Display, Platine, Akku und Rückglas reflektiert, im Material absorbiert und in Wärme umgewandelt oder in ungewollte Richtungen gestreut. Da diese Reflexions- und Absorptionseigenschaften stark von der Materialzusammensetzung, den Schichtdicken und der jeweiligen Frequenz abhängen, entsteht ein hochgradig gekoppeltes System. Eine tatsächliche Anfertigung verschiedener Prototypen mit jeweils unterschiedlichen Geometrien und Materialkonfigurationen ist somit aufgrund des finanziellen und zeitlichen Aufwands nur eingeschränkt möglich. Stattdessen erfolgt die Abstimmung über computergestützte Simulationen, die die elektromagnetischen Eigenschaften eines Geräts im Voraus berechnen, sodass bereits der erste Prototyp ein optimiertes Abstrahlverhalten besitzen kann.

Der Fall steht dabei stellvertretend für eine ganze Klasse technischer Problemstellungen --- von der Radarsignatur eines Kampfjets bis zur Spulenauslegung in Magnetresonanztomographen. Ihnen allen ist gemeinsam, dass Materialien unterschiedlicher Permittivität $epsilon_r$ und Leitfähigkeit $sigma$ auf unregelmäßigen, häufig gekrümmten oder geschichteten Grenzflächen zusammentreffen. Zwar beschreiben die Maxwell-Gleichungen das Feld bei gegebenen Quellen und Materialeigenschaften vollständig @lit:Maxwell1865, geschlossene Lösungen sind aber nur für wenige hochsymmetrische Idealfälle bekannt: die ebene Welle im homogenen Raum, die Reflexion an einer unendlich ausgedehnten ebenen Grenzfläche oder die Streuung an einer Kugel. Entfällt diese Symmetrie, wie in den genannten Beispielen, so ist für das entstehende System partieller Differentialgleichungen im Allgemeinen keine geschlossene analytische Lösung verfügbar.

Wo eine analytische Lösung nicht verfügbar ist, treten deshalb numerische Verfahren an ihre Stelle: Das Kontinuum wird in ein endliches Gitter überführt und Differentialquotienten durch Differenzenquotienten angenähert.

Eines der verbreitetsten Verfahren in diesem Zusammenhang ist die _Finite-Difference Time-Domain_-Methode (FDTD), welche in der Praxis meist über kommerzielle Software angewandt wird @lit:Teixeira2023. Weit verbreitet ist dabei etwa die Software Ansys Lumerical FDTD, deren Quellcode jedoch nicht öffentlich zugänglich ist @lit:AnsysLumerical. Mit MEEP @lit:Oskooi2010 oder gprMax @lit:Warren2016 existieren jedoch auch _Open-Source_-Implementierungen; die eigentliche Lücke liegt daher weniger in der Verfügbarkeit des Codes als im Nachvollzug der Herleitung: Gegenüber dem physikalischen Kontinuum treten auf einem diskreten Gitter systematische Abweichungen auf, deren Ursachen bei der alleinigen Verwendung kommerzieller Blackbox-Programme oder fremder Quelltexte nur eingeschränkt nachvollzogen werden können.

Dementsprechend verfolgt diese Arbeit den Ansatz, die Maxwell-Gleichungen in ein System diskreter Gleichungen zu überführen und daraus einen eigenen FDTD-Solver in Python zu implementieren, welcher anschließend gegen analytisch beschriebene Testfälle gehalten wird. Dadurch wird sowohl die numerische Genauigkeit des Verfahrens bestimmt als auch untersucht, welche Fehlerquellen sich durch die Diskretisierung ergeben. 

== Historischer Kontext

Maxwell führte 1865 Elektrizität, Magnetismus und Optik in einem einheitlichen Gleichungssystem zusammen und sagte die Existenz sich mit Lichtgeschwindigkeit ausbreitender elektromagnetischer Wellen voraus @lit:Maxwell1865; der experimentelle Nachweis gelang Hertz rund zwanzig Jahre später @lit:Hertz1889. Damit waren die grundlegenden Gesetzmäßigkeiten der klassischen Elektrodynamik in einem einheitlichen theoretischen Rahmen formuliert. Die Beschreibung komplexer elektromagnetischer Konfigurationen stellte jedoch weiterhin insbesondere aufgrund der zu lösenden Randwertprobleme eine Herausforderung dar.

Die numerische Umsetzung begründete Kane S. Yee 1966 @lit:Yee1966, indem er die elektrischen und magnetischen Feldkomponenten räumlich und zeitlich gegeneinander versetzt auf einem sogenannten Yee-Gitter anordnete. Dadurch lassen sich die gekoppelten Rotationsgleichungen durchgängig mit zentralen Differenzenquotienten approximieren und beide Feldgrößen alternierend fortschreiben (_leapfrog_-Schema).

Praktisch nutzbar wurde das Verfahren allerdings erst mit wachsender Rechenkapazität. Allen Taflove trug Ende des 20. Jahrhunderts maßgeblich zur Etablierung und Verbreitung der FDTD-Methode bei und baute die Methode zu einem Standardwerkzeug der _Computational Electrodynamics_ aus @lit:Taflove2005.

In der heutigen Zeit findet das Verfahren in vielen Bereichen Verwendung, etwa in der Medizintechnik, beim Entwurf von Antennenarrays, bei der Modellierung von Halbleiterbauelementen sowie in der Photonik und Optik @lit:Teixeira2023.

== Persönliche Motivation

Schon vor meiner Beschäftigung mit FDTD hat mich die Verbindung von Physik, Mathematik und Informatik besonders interessiert. Während meiner Teilnahme am Astro-Pi-Projekt der ESA konnte ich selbst erleben, wie sich physikalische Fragestellungen mit Programmierung und mathematischen Methoden bearbeiten lassen. Gerade diese Art des Arbeitens hat mir gezeigt, wie viel spannender ein Thema werden kann, wenn man nicht nur theoretische Zusammenhänge betrachtet, sondern sie auch selbst mit einem Programm untersuchen und sichtbar machen kann. Deshalb finde ich FDTD für meine Seminararbeit besonders interessant: Die Methode verbindet genau diese drei Bereiche und ermöglicht es mir, physikalische Vorgänge mathematisch herzuleiten und anschließend selbst numerisch zu simulieren. Gleichzeitig bietet mir die Arbeit die Möglichkeit, die in der Schule meist getrennt behandelten Fächer einmal in einer konkreten wissenschaftlichen Anwendung zusammenzuführen.

== Zentrale Fragestellung und Aufbau der Arbeit

Aus den dargestellten Überlegungen ergibt sich die folgende Leitfrage:

#align(center)[
  #block(width: 92%)[
    _Wie lassen sich die Maxwell-Gleichungen in eine zweidimensionale
    FDTD-Simulation überführen und in Python implementieren, und wovon hängt die
    numerische Genauigkeit des daraus entwickelten FDTD-Solvers beim Vergleich mit
    analytischen Lösungen für die Ausbreitung, Reflexion und Transmission ebener Wellen ab?_
  ]
]

Kapitel 2 stellt dazu die theoretischen Voraussetzungen bereit: die Maxwell-Gleichungen in differentieller Form, die Herleitung der Wellengleichung, die ebene Welle als analytische Referenzlösung sowie die Grenzen analytischer Verfahren. Aus diesen Grundlagen entwickelt Kapitel 3 die FDTD-Methode, vom Differenzenquotienten über das Yee-Gitter bis zu den Update-Gleichungen des zweidimensionalen $"TM"_z$-Falls, und ergänzt sie um die Courant-Friedrichs-Lewy-Bedingung sowie die verwendeten Randbedingungen. Die Übersetzung in ein lauffähiges Programm dokumentiert Kapitel 4. Kapitel 5 unterzieht den entwickelten Solver drei Testfällen, deren Lösungen sich jeweils in geschlossener Form angeben lassen, und leitet aus dem Vergleich die Fehlerquellen und Grenzen des Verfahrens ab. Kapitel 6 fasst die Erkenntnisse der Arbeit abschließend zusammen und gibt einen Ausblick.
