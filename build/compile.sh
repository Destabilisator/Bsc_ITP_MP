#!/bin/bash

# configure cmake
device_name=$(hostname)
start_time_configure=$(date +%s%3N)
echo "configuring cmake..."

cmake -DCMAKE_BUILD_TYPE=Release ..

current=$(date +%s%3N)
elapsed_configure="$(($current-$start_time_configure))"
echo "done, this took $elapsed_configure milliseconds"

# compile program
start_time_compile=$(date +%s%3N)
echo "compiling..."

make

current=$(date +%s%3N)
elapsed_compile="$(($current-$start_time_compile))"
elapsed_total="$(($current-$start_time_configure))"
echo "done, this took $elapsed_compile milliseconds"

# write times to history.txt
echo "" && echo "all done"
echo "configure time: $elapsed_configure milliseconds, compile time: $elapsed_compile milliseconds"
echo "total time: $elapsed_total milliseconds"

output="$device_name: $elapsed_configure - $elapsed_compile - $elapsed_total"
echo $output >> "history.txt"