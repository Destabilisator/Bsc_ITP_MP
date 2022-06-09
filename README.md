# Bsc_ITP_MP
Numerische Analyse von Spinsystemen

Zur Betrachtung von Spinsystemen können unterschiedliche numerische Verfahren angewendet werden.
Gegenwärtig werden hier einige Ansätze der exakten Diagionalisierung verfolgt.

Der Sourcecode muss als "Bsc_ITP_MP" im Ordner "cmake-build-release" kompiliert werden.
Anschließend können über die beiliegenden -sh-Dateien uinterschiedliche Plopts erzeigt werden.

Bis N = 30 liegen bereits alle notwendigen Ordner vor. Das Limit deses Programms liegt be N = 32.
Allerdings steigen der RAM-Bedarf und die Rechenzeit sehr schnell an.

Alle entstehenden Plots sind im Ordner "results/" beziehungsweise im Unterordner "results/3DData/" zu finden.
Plots unterschiedlicher Systemgrößen (N) überschreiben sich nicht gegenseitig.
Es können die Suszeptibilität, die speziifische Wärmekapazität, die Anregungsenergie und die Energiedispersion
für unterschiedliche Temperaturen und Kopplungskonstanten dargestellt werden.


# when compiling via the command line:
mkdir build <br />
cd build <br />
cmake -DCMAKE_BUILD_TYPE=Release .. <br />
make
