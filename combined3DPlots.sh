#!/bin/bash

Jstart=0.01
Jend=5
Jcount=20
Tstart=0.01
Tend=5
Tcount=50000

noX=-X

start_time=$SECONDS

if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

for i in 6 8 10 12; do
    $pth ./plotting/deleteData.py $i
    ./cmake-build-release/$prgm 3D $i $Jstart $Jend $Jcount $Tstart $Tend $Tcount -1 $noX silent
    echo "plotting..."
    $pth ./plotting/plot3D.py $i no-show
    echo ""
    echo ""
done

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"