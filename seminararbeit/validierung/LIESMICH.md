# Validierungsrechnungen zu Kapitel 5

Alle Zahlen und Abbildungen des Kapitels 5 entstehen aus diesen Skripten.
Benötigt werden nur `numpy` und `matplotlib`.

| Datei | Inhalt |
|---|---|
| `fdtd_core.py` | Kopie des Lösers aus Kapitel 4 (unverändert) |
| `analytik.py` | analytische Vergleichsgrößen: Bessel-/Hankel-Funktionen ohne scipy, Dispersionsrelation des Yee-Gitters, Fresnel-Koeffizienten. Direkt ausführbar — dann läuft der Selbsttest gegen die exakte Reihe in Bruchrechnung. |
| `tf1_kanal.py` | Testfall 1a: ebene Welle — Auflösungsstudie, Courant-Studie, Verzug über 50 λ, Pulsverformung |
| `tf1_kreis.py` | Testfall 1b: Kreiswelle — Richtungsabhängigkeit, Amplitudengesetz, Differenzfeld |
| `tf2_resonator.py` | Geometrie und Quellpositionen des Resonators (wird von den beiden folgenden importiert) |
| `tf2a_spektrum.py` | Testfall 2: Eigenfrequenzen, Auflösungsstudie, Referenzlauf über 30 µs |
| `tf2b_energie.py` | Testfall 2: Energieerhaltung über 300 000 Schritte, Stabilitätsgrenze |
| `tf3_grenzflaeche.py` | Testfall 3: Reflexionsgrad, Verluste, Vergleich zweier Messverfahren |
| `tf4_pml.py` | Testfall 4: Restreflexion gegen Referenzlösung, Dämpfungsprofil |
| `abbildungen.py` | erzeugt alle PDF-Abbildungen aus den Daten in `daten/` |
| `ZAHLEN.txt` | Zusammenstellung aller im Kapitel zitierten Messwerte |
| `daten/` | Ergebnisse der Läufe (JSON und NPZ). Die Skripte schreiben in einen Unterordner `data/`; der Ordner `daten/` enthält die archivierte Fassung. |

Reihenfolge: erst die vier `tf*`-Skripte, dann `abbildungen.py`.
Der Pfad, in den `abbildungen.py` schreibt, steht dort oben in der Variablen `AUS`.
