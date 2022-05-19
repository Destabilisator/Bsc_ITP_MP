#!/bin/bash

Jstart=0
Jend=5
Jcount=20
Tstart=0
Tend=5
Tcount=50000

start_time=$SECONDS

if [[ "$OSTYPE" == "msys" ]]; then
    for i in 6 8 10
    do
        python ./plotting/deleteData.py $i
        ./cmake-build-release/Bsc_ITP_MX.exe 3D $i $Jstart $Jend $Jcount $Tstart $Tend $Tcount -1 silent
    done
    
    start_time_plots=$SECONDS
    for i in 6 8 10
    do
        python ./plotting/plot3D.py $i no-show
    done
else
    for i in 20
    do
        python3 ./plotting/deleteData.py $i
        ./cmake-build-release/Bsc_ITP_MX 3D $i $Jstart $Jend $Jcount $Tstart $Tend $Tcount 1 silent
    done
    
    start_time_plots=$SECONDS
    for i in 6 8 10
    do
        python3 ./plotting/plot3D.py $i no-show
    done
fi

elapsed=$(( SECONDS - start_time ))
elapsed_plots=$(( SECONDS - start_time_plots ))
echo ""
echo "all plotting done, this took: $elapsed_plots seconds"
echo "all done, total elapsed time: $elapsed seconds"