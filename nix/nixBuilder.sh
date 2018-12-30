source $stdenv/setup
mkdir pureBuild && cd pureBuild
cmake -DCMAKE_BUILD_TYPE=Release
make
