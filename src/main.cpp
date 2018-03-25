#include "core/device.h"
#include "utils/logger.h"

int main(int argc, char * argv[])
{
    Logger logger;
    LOG("GDS online");
    std::vector<core::Device> devices = core::Device::get_devices();
    for (auto& device : devices)
    {
        device.print_info();
    }
    return 0;
}
