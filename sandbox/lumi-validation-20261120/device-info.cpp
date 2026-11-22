#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
int main() {
    int count = 0;
    if (hipGetDeviceCount(&count) != hipSuccess) return 1;
    for (int device = 0; device < count; ++device) {
        char bus[64]{};
        if (hipDeviceGetPCIBusId(bus, sizeof(bus), device) != hipSuccess) return 2;
        std::printf("rank=%s visible_count=%d device=%d pci_bus=%s\n",
                    std::getenv("SLURM_LOCALID"), count, device, bus);
    }
}
