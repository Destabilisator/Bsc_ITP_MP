# Bsc_ITP_MP
Numerische Analyse von Spinsystemen

Zur Betrachtung von Spinsystemen können unterschiedliche numerische Verfahren angewendet werden.
Gegenwärtig werden hier einige Ansätze der exakten Diagionalisierung verfolgt.

Der Sourcecode muss als "Bsc_ITP_MP" im Ordner "cmake-build-release" kompiliert werden.
Anschließend können über die beiliegenden -sh-Dateien uinterschiedliche Plopts erzeigt werden.

Bis N = 32 liegen bereits alle notwendigen Ordner vor. Das Limit deses Programms liegt be N = 32 (int 32 bit). Um größere Systeme zu betrachten müssen sehr viele "int" in "long" geänder werden, viel Spaß...

Alle entstehenden Plots sowie die zu Grunde liegenden Daten werden im Ordner "results/" abgespeichert.
Es können die Suszeptibilität, die speziifische Wärmekapazität, die Anregungsenergie, die Energiedispersion und Spinlückenenergie
für unterschiedliche Temperaturen und Kopplungskonstanten dargestellt werden.


# when compiling via the command line:

manual: <br />

mkdir build <br />
cd build <br />
cmake -DCMAKE_BUILD_TYPE=Release .. <br />
make

or with the script: <br />

cd build <br />
bash compile.sh
