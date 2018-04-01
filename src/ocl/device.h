#ifndef _DEVICE_H_
#define _DEVICE_H_

#include <CL/cl2.hpp>
#include <memory>

namespace gds {

namespace ocl {

class Device
{
public:
    Device(const cl_device_id &device);
    void print_info();
    static std::vector<Device> get_devices();
private:
    std::shared_ptr<cl::Device> _device;
};

}

}

#endif // _DEVICE_H_
