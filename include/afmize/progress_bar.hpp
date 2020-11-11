#ifndef AFMIZE_PROGRESS_BAR_HPP
#define AFMIZE_PROGRESS_BAR_HPP

#include <array>
#include <chrono>
#include <sstream>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

/* This code is originally implemented for other project named Mjolnir,       *
 * by the same author, Toru Niina (c) 2018. destributed under the MIT License.*/

namespace afmize
{

template<std::size_t Width>
class progress_bar
{
    //XXX: requires UTF-8. TODO: consider setting locale
    static constexpr auto full  = u8"█"; // U+2588 Full block
    static constexpr auto l7    = u8"▉"; // U+2589 Left seven eighths block
    static constexpr auto l6    = u8"▊"; // U+258A Left three quarters block
    static constexpr auto l5    = u8"▋"; // U+258B Left five eighths block
    static constexpr auto l4    = u8"▌"; // U+258C Left half block
    static constexpr auto l3    = u8"▍"; // U+258D Left three eighths block
    static constexpr auto l2    = u8"▎"; // U+258E Left one quarter block
    static constexpr auto l1    = u8"▏"; // U+258F Left one eighth block

  public:

    explicit progress_bar(std::size_t tot)
        : total_(tot), r_total_(1.0 / tot), start_(std::chrono::system_clock::now())
    {}
    ~progress_bar() = default;
    progress_bar(progress_bar const&) = default;
    progress_bar(progress_bar &&)     = default;
    progress_bar& operator=(progress_bar const&) = default;
    progress_bar& operator=(progress_bar &&)     = default;

    std::string format(std::size_t count)
    {
        if(count == 0){start_ = std::chrono::system_clock::now();}
        std::ostringstream os;

        //XXX: requires UTF-8.
        constexpr auto full  = u8"█"; // U+2588 Full block
        constexpr auto l7    = u8"▉"; // U+2589 Left seven eighths block
        constexpr auto l6    = u8"▊"; // U+258A Left three quarters block
        constexpr auto l5    = u8"▋"; // U+258B Left five eighths block
        constexpr auto l4    = u8"▌"; // U+258C Left half block
        constexpr auto l3    = u8"▍"; // U+258D Left three eighths block
        constexpr auto l2    = u8"▎"; // U+258E Left one quarter block
        constexpr auto l1    = u8"▏"; // U+258F Left one eighth block

        const double ratio = (count == total_) ? 1.0 :
            std::max(0.0, std::min(1.0, count * this->r_total_));

        std::array<char, 8> buf;
        buf.fill('\0');
        std::snprintf(buf.data(), 8, "%5.1f%%|", ratio * 100.0);
        os << '\r' << buf.data();

        const std::size_t filled = std::floor(ratio*Width);
        for(std::size_t i=0; i<filled; ++i)
        {
            os << full;
        }
        if(Width > filled)
        {
            switch(static_cast<std::size_t>(ratio * Width * 8) % 8)
            {
                case 0: {os << ' '; break;}
                case 1: {os << l1;  break;}
                case 2: {os << l2;  break;}
                case 3: {os << l3;  break;}
                case 4: {os << l4;  break;}
                case 5: {os << l5;  break;}
                case 6: {os << l6;  break;}
                case 7: {os << l7;  break;}
            }
            for(std::size_t i=1; i<Width - filled; ++i)
            {
                os << ' ';
            }
        }
        os << '|';

        if(ratio == 0.0)
        {
            os << std::flush;
            return os.str();
        }

        const auto current  = std::chrono::system_clock::now();
        const auto residual = std::chrono::duration_cast<std::chrono::milliseconds>(
            (current - start_) * (1.0 - ratio) / ratio).count() * 0.001;

        buf.fill('\0');
        if(residual < 60.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual);
            os << buf.data() << " seconds remaining ";
        }
        else if(residual < 60.0 * 60.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167);
            os << buf.data() << " minutes remaining ";
        }
        else if(residual < 60.0 * 60.0 * 24.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167 * 0.0167);
            os << buf.data() << " hours remaining   ";
        }
        else if(residual < 60.0 * 60.0 * 24.0 * 99.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167 * 0.0167 * 0.0417);
            os << buf.data() << " days remaining    ";
        }
        else
        {
            os << " over 100 days remaining";
        }
        return os.str();
    }

    void reset_total(std::size_t t)
    {
        this->total_   = t;
        this->r_total_ = 1.0 / t;
        this->start_   = std::chrono::system_clock::now();
    }

    std::size_t total() const noexcept {return total_;}

  private:

    std::size_t total_;
    double      r_total_;
    std::chrono::system_clock::time_point start_;
};


} // afmize
#endif// AFMIZE_PROGRESS_BAR_HPP
