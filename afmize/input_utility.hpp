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
read_as_angstrom(const toml::value& v)
{
    using namespace std::literals::string_literals;

    const std::vector<std::string> hints{ // for error message.
            "\"1.0nm\"                : 1.0 nano meter.",
            "\"1.0Å\" or \"1.0angstrom\": 1.0 angstrom.",
            "\"1.0pm\"                : 1.0 pico meter.",
            "1.0 (floating)         : 1.0 angstrom.",
            "1   (integer)          : 1.0 angstrom."
        };

    if(v.is_floating()) {return v.as_floating();}
    if(v.is_integer())  {return v.as_integer();}
    if(!v.is_string())
    {
        throw std::runtime_error(toml::format_error("input has invalid type " +
            toml::stringize(v.type()) + ", none of string, floating, or integer",
            v, "here, only the following inputs are allowed", hints));
    }

    const std::string str = v.as_string();
    toml::detail::location<std::string> loc("internal region", str);
    const auto result = toml::detail::parse_floating(loc);
    if(!result)
    {
        throw std::runtime_error(toml::format_error("[error] `"s + str +
            "` is not a valid floating-unit pair."s,
            v, "here, only the following inputs are allowed", hints));
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
        + unit + "` appeared."s,
        v, "here, only the following inputs are allowed", hints));

}

} // afmize
#endif// AFMIZE_UTILITY_HPP
