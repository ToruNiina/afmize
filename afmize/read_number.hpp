#ifndef AFMIZE_READ_NUMBER_HPP
#define AFMIZE_READ_NUMBER_HPP
#include <string>

namespace afmize
{

template<typename T>
T read_number(const std::string&);

template<>
inline double read_number<double>(const std::string& str)
{
    return std::stod(str);
}

template<>
float read_number<float>(const std::string& str)
{
    return std::stof(str);
}

} // afmize
#endif// AFMIZE_READ_NUMBER_HPP
