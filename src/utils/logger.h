#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <string>
#include <mutex>

#define LOG(text) gds::Logger::instance().log(text)

namespace gds {

class Logger
{
public:
    Logger();
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
    ~Logger();

    static Logger& instance();
    void log(const std::string& text);
private:
    static Logger* _this;
    std::mutex _mutex;
};

}

#endif // _LOGGER_H_
