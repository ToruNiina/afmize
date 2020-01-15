#ifndef AFMIZE_INPUT_UTILITY_HPP
#define AFMIZE_INPUT_UTILITY_HPP
#include <type_traits>
#include <utility>
#include <stdexcept>
#include <extlib/toml11/toml.hpp>

namespace afmize
{

template<typename T>
std::enable_if_t<std::is_floating_point<T>::value, T>
read_as_angstrom(const toml::value& v, const std::string& name)
{
    using namespace std::literals::string_literals;

    if(v.is_floating()) {return v.as_floating();}
    if(v.is_integer())  {return v.as_integer();}
    if(!v.is_string())
    {
        throw std::runtime_error("`"s + name + "` has invalid type `"s +
            toml::stringize(v.type()) + "` (none of string, floating, "
            "or integer)."s);
    }

    const std::string str = v.as_string();
    toml::detail::location<std::string> loc(name, str);
    const auto result = toml::detail::parse_floating(loc);
    if(!result)
    {
        throw std::runtime_error(toml::format_error("[error] `"s + str +
            "` is not a valid floating-unit pair."s, v, "here"));
    }
    const auto val = result.unwrap().first;
    const std::string unit(loc.iter(), loc.end());

    if(unit == "pm")
    {
        return static_cast<T>(val * 0.01);
    }
    if(unit == "angst" || unit == "angstrom" || unit == u8"Å" || unit == u8"Å")
    {
        // left is U+212B "angstrom",
        // right is U+00C5 "latin capical letter A with ring above".
        return static_cast<T>(val);
    }
    else if(unit == "nm")
    {
        return static_cast<T>(val * 10.0);
    }
    else if(unit == "um" || unit == u8"μm")
    {
        return static_cast<T>(val * 10'000.0);
    }
    else if(unit == "mm")
    {
        return static_cast<T>(val * 10'000'000.0);
    }

    throw std::runtime_error(toml::format_error("[error] unknown length unit `"s
        + unit + "` appeared."s, v, "here"));
}

} // afmize
#endif// AFMIZE_UTILITY_HPP
