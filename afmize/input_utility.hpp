#ifndef AFMIZE_INPUT_UTILITY_HPP
#define AFMIZE_INPUT_UTILITY_HPP
#include <type_traits>
#include <utility>
#include <extlib/toml/toml/toml.hpp>

namespace afmize
{

template<typename T>
decltype(toml::get<T>(std::declval<toml::value>()))
get(const toml::table& tab, const std::string& key, const std::string loc)
{
    // to show the error message
    try{return toml::get<T>(tab.at(key));}
    catch(std::out_of_range)
    {
        using namespace std::literals::string_literals;
        throw std::out_of_range(
                "key `"s + key + "` is not found in a table "s + loc);
    }
}

inline toml::value
find(const toml::table& tab, const std::string& key, const std::string loc)
{
    // to show the error message
    try{return tab.at(key);}
    catch(std::out_of_range)
    {
        using namespace std::literals::string_literals;
        throw std::out_of_range(
                "key `"s + key + "` is not found in a table "s + loc);
    }
}

template<typename T>
std::enable_if_t<std::is_floating_point<T>::value, T>
read_as_angstrom(const toml::value& v, const std::string& name)
{
    using namespace std::literals::string_literals;

    switch(v.which())
    {
        case toml::value::string_tag:
        {
            std::string num, unit;
            {
                const auto& str = toml::get<std::string>(v);
                auto iter  = str.begin();
                auto token = toml::detail::lex_float::invoke(iter, str.end());
                auto last  = std::remove(token->begin(), token->end(), '_');
                token->erase(last, token->end());
                num  = *token;
                unit = std::string(iter, str.end());
            }
            const double val = std::stod(num);
            if(unit.front() == '_'){unit.erase(unit.begin());}

            if(unit == "pm")
            {
                return static_cast<T>(val * 0.01);
            }
            if(unit == "angst" || unit == u8"Å" || unit == u8"Å")
            {
                // middle is U+212B "angstrom",
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
            else
            {
                throw std::runtime_error("unknown length unit `" + unit + "` appeared");
            }
        }
        case toml::value::integer_tag:
        {
            return static_cast<T>(toml::get<std::int64_t>(v));
        }
        case toml::value::float_tag:
        {
            return toml::get<T>(v);
        }
        default:
        {
            throw std::runtime_error("`"s + name + "` has "s +
                    "invalid type (none of string, float, or integer.)"s);
        }
    }
}

} // afmize
#endif// AFMIZE_UTILITY_HPP
