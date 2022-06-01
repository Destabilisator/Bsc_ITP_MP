#!/bin/bash

start_time=$SECONDS

if [[ "$#" -ge 4 ]]; then
	N=$1
	START=$2
	END=$3
	COUNT=$4
else
	N=4
	START=1
	END=2
	COUNT=30
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

build=build
#build=cmake-build-release

./$build/$prgm $N $START $END $COUNT $CORES $noX $SILENT
echo ""
echo "plotting for N = $N:"
$pth ./plotting/plot.py $N $show $noX

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
