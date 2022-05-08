#!/bin/bash

start_time=$SECONDS

if [[ "$OSTYPE" == "msys" ]]; then
    # ./cmake-build-release/Bsc_ITP_MX.exe 6 0 2 5000 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX.exe 8 0 2 5000 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX.exe 10 0 2 5000 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX.exe 12 0 2 500 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX.exe 14 0 2 50 -1 silent

    # python plotDeltaE.py
    # python plotSpecificHeatJConst.py
    # python plotSpecificHeatTConst.py
    # python plotSusceptibilityJConst.py
    # python plotSusceptibilityTConst.py
    python plotDispersion.py
else
    ./cmake-build-release/Bsc_ITP_MX 6 0 2 5000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 8 0 2 5000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 10 0 2 5000 -1 silent
    ./cmake-build-release/Bsc_ITP_MX 12 0 2 500 -1 silent
    # ./cmake-build-release/Bsc_ITP_MX 14 0 2 50 -1 silent
    
    python3 plotDeltaE.py
    python3 plotSpecificHeatTConst.py
    python3 plotSpecificHeatJConst.py
    python3 plotSusceptibilityTConst.py
    python3 plotSusceptibilityJConst.py
    python3 plotDispersion.py
fi

elapsed=$(( SECONDS - start_time ))
echo "plotting done, this took: $elapsed seconds"
