#ifndef AFMIZE_READER_BASE_HPP
#define AFMIZE_READER_BASE_HPP
#include <afmize/shapes.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

namespace afmize
{
template<typename realT>
class reader_base
{
  public:
    using snapshot_type   = std::vector<sphere<realT>>;
    using trajectory_type = std::vector<snapshot_type>;

    struct no_more_model{};

    reader_base(std::string nm): filename(std::move(nm)){}
    virtual ~reader_base() = default;

    virtual bool is_eof() = 0;

    virtual trajectory_type read_trajectory() = 0;
    virtual snapshot_type   read_snapshot()   = 0;

  protected:

    void highlight_columns(std::ostream& os, const std::string& line,
                            const std::size_t padding, const std::size_t width)
    {
        os << "> " << line << '\n';
        os << "  ";
        for(std::size_t i=0; i<padding; ++i) {os << ' ';}
        for(std::size_t i=0; i<width;   ++i) {os << '^';}
        os << '\n';
    }

    std::string get_substr(const std::string& line, std::size_t ln,
        std::string::size_type b, std::string::size_type w) noexcept
    {
        try
        {
            return line.substr(b, w);
        }
        catch(std::out_of_range)
        {
            std::cerr << "incomplete line found at line " << ln << " in file "
                      << filename << ".\n";
            std::cerr << '\t' << line << '\n';
            std::cerr << '\t' << std::string('^', line.size()) << '\n';
            std::cerr << "check your file.\n";
            std::exit(EXIT_FAILURE);
        }
        catch(std::bad_alloc)
        {
            std::cerr << "memory error detected.\nyou need to buy more RAM.\n";
            std::exit(EXIT_FAILURE);
        }
        catch(...)
        {
            std::cerr << "unexpected error detected while reading file `"
                      << filename << "`.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    static constexpr auto mes_float_format_err =
        "couldn't convert to floating-point value.\n"
        "If it looks ok, locale settings may have a problem. "
        "check your locale.\n";

    static constexpr auto mes_float_range_err =
        "the value exceeds the range of floating-point.\n";

    std::string filename;
};

} // afmize
#endif// AFMIZE_READER_BASE_HPP
