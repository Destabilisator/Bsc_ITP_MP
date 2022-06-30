#!/bin/bash

start_time=$SECONDS
#build=build
build=cmake-build-release

if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

./$build/$prgm && echo "" && echo ""
$pyth plotting/plotBenchmarking.py

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
