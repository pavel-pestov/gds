#include "logger.h"

#include <iostream>

namespace gds {

Logger* Logger::_this = nullptr;

Logger::Logger()
{
    if (_this)
    {
        std::cerr << "can't create second logger" << std::endl;
        abort();
    }
    _this = this;
}

Logger::~Logger()
{
    std::lock_guard<std::mutex> lock(_mutex);
    _this = nullptr;
}

Logger& Logger::instance()
{
    if (!_this)
    {
        std::cerr << "logger not created" << std::endl;
        abort();
    }
    return *_this;
}

void Logger::log(const std::string& text)
{
    std::lock_guard<std::mutex> lock(_mutex);
    std::cout << text << std::endl;
}

}
