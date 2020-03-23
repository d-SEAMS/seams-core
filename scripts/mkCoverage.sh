#!/usr/bin/env bash

topLevel=$(git rev-parse --show-toplevel)
cd $topLevel
if [[ "$PWD" =~ seams-core ]]; then
	rm -rf shellBuild
	mkdir shellBuild
	cd shellBuild
	# cmake .. -DCMAKE_CXX_FLAGS="-pg -fsanitize=address " -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg
	cmake .. -DCMAKE_BUILD_TYPE=Debug
	# cmake .. -DCMAKE_BUILD_TYPE=Release
	make testCoverage -j4
fi
