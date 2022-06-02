#!/bin/bash

start=0.01
end=5
h=1.0

noX=-X
#noX=noX

start_time=$SECONDS
#build=build
build=cmake-build-release

if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

./$build/$prgm 6 $start $end 10000 $h -1 $noX silent && echo "" # 10000
./$build/$prgm 8 $start $end 5000 $h -1 $noX silent && echo "" # 10000
./$build/$prgm 10 $start $end 500 $h -1 $noX silent && echo "" # 5000
./$build/$prgm 12 $start $end 500 $h -1 $noX silent && echo "" # 1000
# ./$build/$prgm 14 $start $end 50 $h -1 $noX silent && echo "" # 50
# ./$build/$prgm 16 $start $end 50 $h -1 $noX silent && echo ""
# ./$build/$prgm 18 $start $end 20 $h -1 $noX silent && echo ""
# ./$build/$prgm 20 $start $end 20 $h -1 $noX silent && echo ""

start_time_plots=$SECONDS
$pth ./plotting/plotDeltaE.py
$pth ./plotting/plotSpecificHeatTConst.py
$pth ./plotting/plotSpecificHeatJConst.py
$pth ./plotting/plotSusceptibilityTConst.py
$pth ./plotting/plotSusceptibilityJConst.py
$pth ./plotting/plotDispersion.py
$pth ./plotting/plotSpinGap.py

elapsed=$(( SECONDS - start_time ))
elapsed_plots=$(( SECONDS - start_time_plots ))
echo ""
echo "plotting done, this took: $elapsed_plots seconds"
echo "all done, total elapsed time: $elapsed seconds"
