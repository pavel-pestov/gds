#include "ocl/device.h"
#include "utils/logger.h"

using namespace gds;

int main(int argc, char * argv[])
{
    Logger logger;
    LOG("GDS online");
    std::vector<ocl::Device> devices = ocl::Device::get_devices();
    for (auto& device : devices)
    {
        device.print_info();
    }
    return 0;
}
