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
	h=$5
else
	h=0.0
fi

if [[ "$#" -ge 6 ]]; then
	CORES=$6
else
	CORES=-1
fi

if [[ "$#" -ge 7 ]]; then
	noX=$7
else
	noX=-X
fi

if [[ "$#" -ge 8 ]]; then
	show=$8
else
	show=-
fi

if [[ "$#" -ge 9 ]]; then
	SILENT=$9
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

#build=build
build=cmake-build-release

./$build/$prgm $N $START $END $COUNT $h $CORES $noX $SILENT
echo ""
echo "plotting for N = $N:"
$pth ./plotting/plot.py $N $show $noX

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
