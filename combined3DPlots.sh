#!/bin/bash

Jstart=0
Jend=5
Jcount=20
Tstart=0
Tend=5
Tcount=50000

noX=-X

start_time=$SECONDS


if [[ "$OSTYPE" == "msys" ]]; then
    for i in 6 8 10 12
    do
        python ./plotting/deleteData.py $i
        ./cmake-build-release/Bsc_ITP_MP.exe 3D $i $Jstart $Jend $Jcount $Tstart $Tend $Tcount -1 $noX silent
        echo "plotting..."
        python ./plotting/plot3D.py $i no-show
    done
else
    for i in 6 8 10 12 14 16
    do
        python3 ./plotting/deleteData.py $i
        ./cmake-build-release/Bsc_ITP_MP 3D $i $Jstart $Jend $Jcount $Tstart $Tend $Tcount 1 $noX silent
        echo "plotting..."
        python3 ./plotting/plot3D.py $i no-show
        echo ""
        echo ""
    done
fi

elapsed=$(( SECONDS - start_time ))
#elapsed_plots=$(( SECONDS - start_time_plots ))
echo ""
#echo "all plotting done, this took: $elapsed_plots seconds"
echo "all done, total elapsed time: $elapsed seconds"