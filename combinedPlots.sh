#!/bin/bash

start=0
end=2.5

start_time=$SECONDS

if [[ "$OSTYPE" == "msys" ]]; then
    ./cmake-build-release/Bsc_ITP_MX.exe 6 $start $end 5000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX.exe 8 $start $end 5000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX.exe 10 $start $end 500 -1 silent
    ./cmake-build-release/Bsc_ITP_MX.exe 12 $start $end 50 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX.exe 14 $start $end 10 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX.exe 16 $start $end 5 -1 silent

    start_time_plots=$SECONDS
    python ./plotting/plotDeltaE.py
    python ./plotting/plotSpecificHeatJConst.py
    python ./plotting/plotSpecificHeatTConst.py
    python ./plotting/plotSusceptibilityJConst.py
    python ./plotting/plotSusceptibilityTConst.py
    python ./plotting/plotDispersion.py

else
    ./cmake-build-release/Bsc_ITP_MX 6 $start $end 1000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 8 $start $end 1000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 10 $start $end 500 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 12 $start $end 50 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX 14 $start $end 500 8 silent
    # ./cmake-build-release/Bsc_ITP_MX 16 $start $end 50 2 silent

    start_time_plots=$SECONDS
    python3 ./plotting/plotDeltaE.py
    python3 ./plotting/plotSpecificHeatTConst.py
    python3 ./plotting/plotSpecificHeatJConst.py
    python3 ./plotting/plotSusceptibilityTConst.py
    python3 ./plotting/plotSusceptibilityJConst.py
    python3 ./plotting/plotDispersion.py
fi

elapsed=$(( SECONDS - start_time ))
elapsed_plots=$(( SECONDS - start_time_plots ))
echo ""
echo "plotting done, this took: $elapsed_plots seconds"
echo "all done, total elapsed time: $elapsed seconds"
