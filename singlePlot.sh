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
	SILENT=$6
else
	SILENT=-
fi


if [[ "$OSTYPE" == "msys" ]]; then
    ./cmake-build-release/Bsc_ITP_MX.exe $N $J_START $J_END $J_COUNT $CORES $SILENT
    echo "plotting..."
    #python plot.py $N
else
    ./cmake-build-release/Bsc_ITP_MX $N $J_START $J_END $J_COUNT $CORES $SILENT
    echo "plotting..."
    python3 plot.py $N
fi

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
