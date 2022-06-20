#!/bin/bash

# args: low or high
if [[ "$#" -ge 1 ]]; then
	regime=$1
else
	regime=low
fi

# data
start=0.01
end=2.5
h=0.0

#plot
plotstart=25
plotend=$end

noX=-X
#noX=noX

start_time=$SECONDS
build=build
#build=cmake-build-release

if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

# #$pth plotting/deleteData.py 6 SG
./$build/$prgm 6 $start $end 50 $h -1 $noX silent && echo ""
# #$pth plotting/deleteData.py 8 SG
./$build/$prgm 8 $start $end 50 $h -1 $noX silent && echo ""
# #$pth plotting/deleteData.py 10 SG
./$build/$prgm 10 $start $end 50 $h -1 $noX silent && echo ""
# #$pth plotting/deleteData.py 12 SG
./$build/$prgm 12 $start $end 50 $h 25 $noX silent && echo ""
./$build/$prgm 14 $start $end 50 $h 10 $noX silent && echo "" # 25
./$build/$prgm 16 $start $end 50 $h 10 $noX silent && echo "" # 25
./$build/$prgm 18 $start $end 50 $h 10 $noX silent && echo "" # 25

# $pth plotting/plotStatisticsC.py 0 $plotend low && echo "" && echo ""
# $pth plotting/plotStatisticsX.py 0 $plotend low
# $pth plotting/plotSpinGapQT.py low && echo "" && echo ""
# $pth plotting/plotDeltaEQT.py low && echo "" && echo ""

./$build/$prgm 20 $start $end 50 $h 10 $noX silent && echo ""
./$build/$prgm 22 $start $end 50 $h 5 $noX silent && echo ""
./$build/$prgm 24 $start $end 50 $h 5 $noX silent && echo ""
./$build/$prgm 26 $start $end 50 $h 1 $noX silent && echo ""
./$build/$prgm 28 $start $end 50 $h 1 $noX silent && echo ""
./$build/$prgm 30 $start $end 50 $h 1 $noX silent && echo ""
./$build/$prgm 32 $start $end 50 $h 1 $noX silent && echo ""

#$pth plotting/plotSpinGapQT.py high && echo "" && echo ""

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
