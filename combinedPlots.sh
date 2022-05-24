#!/bin/bash

start=0
end=2.5

noX=-X

start_time=$SECONDS

if [[ "$OSTYPE" == "msys" ]]; then
    ./cmake-build-release/Bsc_ITP_MP.exe 6 $start $end 5000 -1 $noX silent
    ./cmake-build-release/Bsc_ITP_MP.exe 8 $start $end 5000 -1 $noX silent
    ./cmake-build-release/Bsc_ITP_MP.exe 10 $start $end 500 -1 $noX silent
    ./cmake-build-release/Bsc_ITP_MP.exe 12 $start $end 50 -1 $noX silent
    # ./cmake-build-release/Bsc_ITP_MP.exe 14 $start $end 10 -1 $noX silent
    # ./cmake-build-release/Bsc_ITP_MP.exe 16 $start $end 5 -1 $noX silent

    start_time_plots=$SECONDS
    python ./plotting/plotDeltaE.py
    python ./plotting/plotSpecificHeatJConst.py
    python ./plotting/plotSpecificHeatTConst.py
    python ./plotting/plotSusceptibilityJConst.py
    python ./plotting/plotSusceptibilityTConst.py
    python ./plotting/plotDispersion.py
    python ./plotting/plotSpinGap.py

else
    ./cmake-build-release/Bsc_ITP_MP 6 $start $end 50 -1 $noX silent # 10000
    ./cmake-build-release/Bsc_ITP_MP 8 $start $end 50 -1 $noX silent # 10000
    ./cmake-build-release/Bsc_ITP_MP 10 $start $end 50 -1 $noX silent # 5000
    ./cmake-build-release/Bsc_ITP_MP 12 $start $end 50 -1 $noX silent # 1000
    # ./cmake-build-release/Bsc_ITP_MP 14 $start $end 50 -1 $noX silent # 50
    # ./cmake-build-release/Bsc_ITP_MP 16 $start $end 50 -1 $noX silent #
    # ./cmake-build-release/Bsc_ITP_MP 18 $start $end 20 -1 $noX silent
    # ./cmake-build-release/Bsc_ITP_MP 20 $start $end 20 -1 $noX silent

    start_time_plots=$SECONDS
    # python3 ./plotting/plotDeltaE.py
    # python3 ./plotting/plotSpecificHeatTConst.py
    # python3 ./plotting/plotSpecificHeatJConst.py
    # python3 ./plotting/plotSusceptibilityTConst.py
    # python3 ./plotting/plotSusceptibilityJConst.py
    # python3 ./plotting/plotDispersion.py
    python3 ./plotting/plotSpinGap.py
fi

elapsed=$(( SECONDS - start_time ))
elapsed_plots=$(( SECONDS - start_time_plots ))
echo ""
echo "plotting done, this took: $elapsed_plots seconds"
echo "all done, total elapsed time: $elapsed seconds"
