#!/bin/bash

start_time=$SECONDS

if [[ "$#" -ge 4 ]]; then
	N=$1
	J_START=$2
	J_END=$3
	J_COUNT=$4
else
	N=4
	J_START=1
	J_END=2
	J_COUNT=30
fi

if [[ "$#" -ge 5 ]]; then
	CORES=$5
else
	CORES=-1
fi

if [[ "$#" -ge 6 ]]; then
	noX=$6
else
	noX=-X
fi

if [[ "$#" -ge 7 ]]; then
	show=$7
else
	show=-
fi

if [[ "$#" -ge 8 ]]; then
	SILENT=$8
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

./cmake-build-release/$prgm $N $J_START $J_END $J_COUNT $CORES $noX $SILENT
echo ""
echo "plotting for N = $N:"
$pth ./plotting/plot.py $N $show $noX

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
