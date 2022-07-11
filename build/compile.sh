#!/bin/bash

start_time_configure=$SECONDS
echo "configuring cmake..."

cmake -DCMAKE_BUILD_TYPE=Release ..

elapsed_configure=$(( SECONDS - start_time_configure ))
echo "done, this took $elapsed_configure seconds"
start_time_compile=$SECONDS
echo "compiling..."

make

elapsed_compile=$(( SECONDS - start_time_compile ))
elapsed_total=$(( SECONDS - start_time_configure ))
echo "done, this took $elapsed_compile seconds"
echo "" && echo "all done"
echo "configure time: $elapsed_configure seconds, compile time: $elapsed_compile seconds"
echo "total time: $elapsed_total seconds"
