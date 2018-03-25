#include "device.h"

namespace core {

Device::Device()
{
    _device.reset(new cl::Device());
}

}
