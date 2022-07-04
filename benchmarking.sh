#!/bin/bash

start_time=$SECONDS
build=build
#build=cmake-build-release

if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

# delete old data
$pth ./plotting/deleteData.py -1 BRT
$pth ./plotting/deleteData.py -1 BMU

# generate new data
./$build/$prgm && echo "" && echo ""

# plot data
$pth plotting/plotBenchmarking.py

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
