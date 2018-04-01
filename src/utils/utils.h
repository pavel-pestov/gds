#ifndef _UTILS_H_
#define _UTILS_H_

#include <string>
#include <sstream>

namespace gds {

inline std::string to_string(const char* x)
{
    return std::string(x);
}

inline std::string to_string(const std::string& x)
{
    return x;
}

template<typename T> inline std::string to_string(const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

}

#endif // _UTILS_H_
