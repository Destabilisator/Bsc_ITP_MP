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
    python plotDeltaE.py
    python plotSpecificHeatJConst.py
    python plotSpecificHeatTConst.py
    python plotSusceptibilityJConst.py
    python plotSusceptibilityTConst.py
    python plotDispersion.py

else
    ./cmake-build-release/Bsc_ITP_MX 6 $start $end 1000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 8 $start $end 1000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 10 $start $end 500 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 12 $start $end 50 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX 14 $start $end 500 8 silent
    # ./cmake-build-release/Bsc_ITP_MX 16 $start $end 50 2 silent

    start_time_plots=$SECONDS
    python3 plotDeltaE.py
    python3 plotSpecificHeatTConst.py
    python3 plotSpecificHeatJConst.py
    python3 plotSusceptibilityTConst.py
    python3 plotSusceptibilityJConst.py
    python3 plotDispersion.py
fi

elapsed=$(( SECONDS - start_time ))
elapsed_plots=$(( SECONDS - start_time_plots ))
echo ""
echo "plotting done, this took: $elapsed_plots seconds"
echo "all done, total elapsed time: $elapsed seconds"
