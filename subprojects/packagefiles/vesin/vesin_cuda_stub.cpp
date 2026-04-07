// Stub for CUDA functions when building without gpu-lite/CUDA support.
// vesin.cpp unconditionally references these symbols; this stub satisfies
// the linker while making CUDA calls return an error at runtime.

#include "vesin_cuda.hpp"
#include <stdexcept>

namespace vesin {
namespace cuda {

void free_neighbors(VesinNeighborList&) {
    // No-op: nothing to free without CUDA allocations
}

void neighbors(
    const double (*)[3],
    size_t,
    const double [3][3],
    const bool*,
    VesinOptions,
    VesinNeighborList&
) {
    throw std::runtime_error(
        "vesin: CUDA neighbor list requested but built without gpu-lite. "
        "Use VesinCPU device or rebuild with CUDA support.");
}

CudaNeighborListExtras* get_cuda_extras(VesinNeighborList*) {
    return nullptr;
}

} // namespace cuda
} // namespace vesin
