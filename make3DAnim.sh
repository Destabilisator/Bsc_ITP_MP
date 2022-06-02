#!/bin/bash

Jstart=0.01
Jend=5
Jcount=20
Tstart=0.01
Tend=5
Tcount=50000

N=$1
h_START=$2
h_END=$3
h_STEP=$4

start_time=$SECONDS

if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

#build=build
build=cmake-build-release

for h in $(seq $h_START $h_STEP $h_END); do
    $pth ./plotting/deleteData.py $N
    ./$build/$prgm 3D $N $Jstart $Jend $Jcount $Tstart $Tend $Tcount $h -1 $noX silent
    echo "plotting..."
    $pth ./plotting/plot3D.py $N no-show C noX
    echo ""
    echo ""
done

$pth ./plotting/makeAnim.py $N

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
