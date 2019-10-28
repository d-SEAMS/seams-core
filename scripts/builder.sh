#!/usr/bin/env bash

topLevel=$(git rev-parse --show-toplevel)
cd $topLevel
if [[ "$PWD" =~ seams-core ]]; then
    rm -rf build
    mkdir build
    cd build
    # cmake .. -DCMAKE_CXX_FLAGS="-pg -fsanitize=address " -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
    cmake .. -DCMAKE_BUILD_TYPE=Debug
    # cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j4
    gdb --args yodaStruct --script ../lua_inputs/transition_diff.lua -f ../lua_inputs/parameter.txt -c ../lua_inputs/config.yml
fi
