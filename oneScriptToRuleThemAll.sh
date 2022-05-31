#!/bin/bash

start_time=$SECONDS

# 2D
START=0.01
END=5

# 3D
J_START=0
J_END=5

T_START=0
T_END=5
T_COUNT=50000

# general
CORES=-1
noX=-X
show=no-show
SILENT=silent

# OS
if [[ "$OSTYPE" == "msys" ]]; then
	prgm=Bsc_ITP_MP.exe
	pth=python
else
	prgm=Bsc_ITP_MP
	pth=python3
fi

for N in 6 8 10 12 14 16 18 20 22 24 26 28 30 32; do
    $pth ./plotting/deleteData.py $N
done
echo ""
echo ""

# setting the pc on fire
./cmake-build-release/$prgm 6 $START $END 50000 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 6 $J_START $J_END 50 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 8 $START $END 50000 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 8 $J_START $J_END 50 $T_START 5$T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 10 $START $END 50000 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 10 $J_START $J_END 50 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 12 $START $END 5000 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 12 $J_START $J_END 50 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 14 $START $END 500 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 14 $J_START $J_END 50 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 16 $START $END 500 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 16 $J_START $J_END 50 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 18 $START $END 50 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 18 $J_START $J_END 50 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 20 $START $END 30 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 20 $J_START $J_END 30 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 24 $START $END 20 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 24 $J_START $J_END 20 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 28 $START $END 20 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 28 $J_START $J_END 20 $T_START $T_END $T_COUNT -1 -X silent && echo ""

./cmake-build-release/$prgm 32 $START $END 20 -1 -X silent && echo ""
./cmake-build-release/$prgm 3D 32 $J_START $J_END 20 $T_START $T_END $T_COUNT -1 -X silent && echo ""

# pc has been set on fire
elapsed=$(( SECONDS - start_time ))
echo "all done, total elapsed time: $elapsed seconds"