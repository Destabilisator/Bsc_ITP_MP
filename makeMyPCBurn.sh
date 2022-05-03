#!/bin/bash

if [[ $# == 4 ]]; then
	N=$1
	J_START=$2
	J_END=$3
	J_COUNT=$4
else
	N=4
	J_START=1
	J_END=2
	J_COUNT=30
fi

#case "$OSTYPE" in
#  solaris*) echo "SOLARIS" ;;
#  darwin*)  echo "OSX" ;; 
#  linux*)   echo "LINUX" ;;
#  bsd*)     echo "BSD" ;;
#  msys*)    echo "WINDOWS" ;;
#  cygwin*)  echo "ALSO WINDOWS" ;;
#  *)        echo "unknown: $OSTYPE" ;;
#esac


if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    ./cmake-build-release/Bsc_ITP_MX $N $J_START $J_END $J_COUNT
elif [[ "$OSTYPE" == "msys" ]]; then
    ./cmake-build-release/Bsc_ITP_MX.exe $N $J_START $J_END $J_COUNT
else
    echo "What kind of place is this?!"
	exit 1;
fi

echo "plotting..."

python plot.py