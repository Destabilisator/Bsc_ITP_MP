#!/bin/bash

# args: low or high
if [[ "$#" -ge 1 ]]; then
	regime=$1
else
	regime=low
fi

# data
start=0.01
end=2.0
h=0.0

#plot
plotstart=25
plotend=$end

noX=-X
#noX=noX

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

for N in 6 8 10 12 14; do # 
    # $pth plotting/deleteData.py $N SG
	# $pth plotting/deleteData.py $N EE
	./$build/$prgm $N $start $end 50 $h 10 $noX silent && echo ""
done

# # $pth plotting/plotSpinGapQT.py low && echo "" && echo ""
# # $pth plotting/plotDeltaEQT.py low && echo "" && echo ""

# for N in 16 18 20 22 24 26; do # 28 30 32
#     # $pth plotting/deleteData.py $N SG
# 	# $pth plotting/deleteData.py $N EE
# 	./$build/$prgm $N $start $end 50 $h 5 $noX silent && echo ""
# done


# N=14
# $pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
# ./$build/$prgm $N $start $end 50 $h 25 $noX silent && echo ""

N=16
$pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
./$build/$prgm $N $start $end 50 $h 25 $noX silent && echo ""

N=18
$pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
./$build/$prgm $N $start $end 50 $h 25 $noX silent && echo ""

N=20
$pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
./$build/$prgm $N $start $end 50 $h 5 $noX silent && echo ""

N=22
$pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
./$build/$prgm $N $start $end 50 $h 5 $noX silent && echo ""

N=24
$pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
./$build/$prgm $N $start $end 50 $h 2 $noX silent && echo ""

N=26
$pth plotting/deleteData.py $N SG
# # $pth plotting/deleteData.py $N EE
./$build/$prgm $N $start $end 50 $h 1 $noX silent && echo ""

# N=28
# $pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
# ./$build/$prgm $N $start $end 50 $h 1 $noX silent && echo ""

# N=30
# $pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
# ./$build/$prgm $N $start $end 50 $h 1 $noX silent && echo ""

# N=32
# $pth plotting/deleteData.py $N SG
# $pth plotting/deleteData.py $N EE
# ./$build/$prgm $N $start $end 50 $h 1 $noX silent && echo ""


elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"
