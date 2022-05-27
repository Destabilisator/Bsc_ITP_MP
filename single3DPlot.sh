#!/bin/bash

start_time=$SECONDS

N=$1
J_START=$2
J_END=$3
J_COUNT=$4

T_START=$5
T_END=$6
T_COUNT=$7

if [[ "$#" -ge 8 ]]; then
	CORES=$8
else
	CORES=-1
fi

if [[ "$#" -ge 9 ]]; then
	noX=${9}
else
	noX=-
fi

if [[ "$#" -ge 10 ]]; then
	show=${10}
else
	show=-
fi

if [[ "$#" -ge 10 ]]; then
	SILENT=${11}
else
	SILENT=-
fi


if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

$pth ./plotting/deleteData.py $N
./cmake-build-release/$prgm 3D $N $J_START $J_END $J_COUNT $T_START $T_END $T_COUNT $CORES $noX $SILENT
echo ""
echo "plotting for N = $N:"
$pth ./plotting/plot3D.py $N $show

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
