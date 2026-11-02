"""Write the CUDA source string of src/gpu_batch.cpp to a file a host
compiler can include (tests/gpu_kernels_host.cpp)."""
import sys

src = open(sys.argv[1], encoding="utf-8").read()
start = 'const char *kKernels = R"CUDA('
a = src.index(start) + len(start)
b = src.index(')CUDA";', a)
with open(sys.argv[2], "w", encoding="utf-8") as out:
    out.write(src[a:b])
