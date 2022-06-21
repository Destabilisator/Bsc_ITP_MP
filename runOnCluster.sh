#!/bin/bash

# args: low or high
if [[ "$#" -ge 1 ]]; then
	regime=$1
else
	regime=low
fi

# data
start=0.01
end=2.5
h=0.0

#plot
plotstart=25
plotend=$end

noX=-X
#noX=noX

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

for N in 6 8 10 12 14; do # 6 8 10 12 14
    # $pth plotting/deleteData.py $N SG
	$pth plotting/deleteData.py $N EE
	./$build/$prgm $N $start $end 50 $h 1 $noX silent && echo ""
done

$pth plotting/plotSpinGapQT.py low && echo "" && echo ""
$pth plotting/plotDeltaEQT.py low && echo "" && echo ""

# for N in 16 18 20 22 24 26 28 30 32; do
#     $pth plotting/deleteData.py $N SG
# 	$pth plotting/deleteData.py $N EE
# 	./$build/$prgm $N $start $end 50 $h 1 $noX silent && echo ""
# done

elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
