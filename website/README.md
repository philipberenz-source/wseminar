# Website + FDTD Studio (Browser)

Die Website stellt die FDTD-Methode kurz vor (`index.html`). Dazu gehört das Web-Studio
(`studio.html`), in dem der Python-Kern **`fdtd_core.py`** über
[Pyodide](https://pyodide.org) (Python + NumPy als WebAssembly) direkt im Browser läuft,
lokal und ohne Installation.

## Starten

```bash
cd website
python -m http.server 8000
```

Danach <http://localhost:8000> öffnen. Beim ersten Start des Studios lädt der Browser
Pyodide einmalig (ca. 10 MB) von cdn.jsdelivr.net, danach kommt es aus dem Cache.

> Das Studio braucht einen Webserver (kein `file://`), weil Python in einem Web Worker läuft.
> Die Startseite funktioniert auch per Doppelklick.

## Aufbau

| Datei | Inhalt |
| --- | --- |
| `index.html` | Startseite: Idee, Yee-Gitter, Leapfrog, CFL, Dispersion, PML, Validierung |
| `studio.html` | Web-Studio (Oberfläche) |
| `py/fdtd_core.py` | **Kopie** des Kerns aus `code snippets/`. Nach Änderungen am Kern neu kopieren |
| `py/bridge.py` | Szenenaufbau (aus `fdtd_studio.py` übernommen, ohne Tk/matplotlib) |
| `assets/js/py-worker.js` | lädt Pyodide + NumPy, führt Kern und Bridge im Hintergrund-Thread aus |
| `assets/js/studio.js` | Bedienung, Feld- und Diagrammdarstellung (rechnet selbst nicht) |
| `assets/js/fdtd.js` | Szenenmodell/Presets, Farbskalen, dazu ein JS-Port des Kerns für die kleinen Animationen der Startseite |
| `assets/js/main.js` | Animationen und Demos der Startseite |

Kern nach Änderungen synchronisieren:

```bash
cp "code snippets/fdtd_core.py" website/py/
```

Szenen-JSON ist mit dem Desktop-Studio (`fdtd_studio.py`) kompatibel.

## Flackerschutz

Dauerstrich-Quellen lassen jeden Bildpunkt mit der Quellfrequenz schwingen. Die Startseite
rechnet deshalb nur einen Zeitschritt pro Bild und hat einen Button „Animation aus“; sie
beachtet außerdem `prefers-reduced-motion`. Im Studio begrenzt der Flackerschutz (standardmäßig
an) die Schritte pro Bild so, dass jeder Bildpunkt höchstens etwa zweimal pro Sekunde wechselt.
