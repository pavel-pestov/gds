#ifndef _DEVICE_H_
#define _DEVICE_H_

#include <CL/cl2.hpp>
#include <memory>

namespace core {

class Device
{
public:
    Device();

private:
    std::unique_ptr<cl::Device> _device;
};

}

#endif
