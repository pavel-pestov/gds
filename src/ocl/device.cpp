#include "device.h"

#include <CL/cl.h>

#include "../utils/logger.h"
#include "../utils/utils.h"

namespace gds {

namespace ocl {

Device::Device(const cl_device_id &device)
{
    _device = std::make_shared<cl::Device>(device);
}

template<class T> std::string get_device_info(cl_device_id device, cl_device_info param_name)
{
    std::vector<char> data(1024);
    size_t out_size;
    cl_int error = clGetDeviceInfo(device, param_name, data.size(), static_cast<void*>(data.data()), &out_size);
    if (error != CL_SUCCESS)
    {
        LOG("can't get device info");
        abort();
    }
    std::string out;
    for (int i = 0; i < out_size / sizeof(T); ++i)
    {
        out += to_string(reinterpret_cast<T*>(data.data())[i]);
        if (i < out_size / sizeof(T) - 1)
            out += ", ";
    }
    return out;
}

template<> std::string get_device_info<char>(cl_device_id device, cl_device_info param_name)
{
    std::vector<char> data(1024);
    size_t out_size;
    cl_int error = clGetDeviceInfo(device, param_name, data.size(), static_cast<void*>(data.data()), &out_size);
    if (error != CL_SUCCESS)
    {
        LOG("can't get device info");
        abort();
    }
    return std::string(&data[0]);
}

void Device::print_info()
{
    LOG("");
    LOG(std::string("CL_DEVICE_NAME = ") + get_device_info<char>(_device->get(),CL_DEVICE_NAME));
    LOG("");
    LOG(std::string("CL_DEVICE_VERSION = ") + get_device_info<char>(_device->get(),CL_DEVICE_VERSION));
    LOG(std::string("CL_DEVICE_OPENCL_C_VERSION = ") + get_device_info<char>(_device->get(),CL_DEVICE_OPENCL_C_VERSION));
    LOG(std::string("CL_DEVICE_MAX_COMPUTE_UNITS = ") + to_string(get_device_info<cl_uint>(_device->get(),CL_DEVICE_MAX_COMPUTE_UNITS)));
    LOG(std::string("CL_DEVICE_MAX_CLOCK_FREQUENCY = ") + to_string(get_device_info<cl_uint>(_device->get(),CL_DEVICE_MAX_CLOCK_FREQUENCY)));
    LOG(std::string("CL_DEVICE_MAX_MEM_ALLOC_SIZE = ") + to_string(get_device_info<cl_ulong>(_device->get(),CL_DEVICE_MAX_MEM_ALLOC_SIZE)));
    LOG(std::string("CL_DEVICE_GLOBAL_MEM_SIZE = ") + to_string(get_device_info<cl_ulong>(_device->get(),CL_DEVICE_GLOBAL_MEM_SIZE)));
    LOG(std::string("CL_DEVICE_AVAILABLE = ") + to_string(get_device_info<cl_bool>(_device->get(),CL_DEVICE_AVAILABLE)));
}

std::vector<Device> Device::get_devices()
{
    std::vector<cl_platform_id> platforms;
    cl_uint num_platforms = 0;
    while (num_platforms == platforms.size())
    {
        platforms.resize(platforms.size() + 1);
        cl_int error = clGetPlatformIDs(platforms.size(), platforms.data(), &num_platforms);
        if (error != CL_SUCCESS)
        {
            LOG("can't get platform");
            abort();
        }
    }
    platforms.resize(num_platforms);

    std::vector<Device> out;

    for (const auto& platform : platforms)
    {
        std::vector<cl_device_id> devices;
        cl_uint num_devices = 0;
        while (num_devices == devices.size()) {
            devices.resize(devices.size() + 1);
            cl_int error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, devices.size(), devices.data(), &num_devices);
            if (error != CL_SUCCESS)
            {
                LOG("can't get device");
                abort();
            }
        }
        devices.resize(num_devices);
        for (const auto& device : devices)
        {
            out.emplace_back(device);
        }
    }
    return out;
}

}

}
