#ifndef AFMIZE_COLOR_MAP_HPP
#define AFMIZE_COLOR_MAP_HPP
#include <extlib/pnm/pnm.hpp>

namespace afmize
{

// this implementation is made based on a image found in an issue #6032 of
// github:matplotlib/matplotlib repository which explains how the color changes
// when value changes.
// URL: https://github.com/matplotlib/matplotlib/issues/6032
template<typename Real, typename Pixel = pnm::rgb_pixel>
Pixel color_afmhot(const Real value, const Real min, const Real max)
{
    // assuming Pixel(R, G, B) and R, G, B are 8-bit integer
    assert(min <= value && value <= max);
    const Real sqrt2        = std::sqrt(Real(2.0));
    const Real sqrt2_plus_1 = sqrt2 + Real(1.0);

    const Real v = ((value - min) / (max - min)) * sqrt2_plus_1;

    if(v < Real(0.5))
    {
        const std::uint8_t R = std::floor(v * 256);
        return Pixel(R, 0, 0);
    }
    else if(v < Real(0.5) * sqrt2_plus_1)
    {
        const Real u = (v - Real(0.5)) / sqrt2;
        if     (u  < Real(0.0)) {return Pixel(127, 0,   0);}
        else if(Real(0.5) <= u) {return Pixel(255, 127, 0);}

        const std::uint8_t R = std::floor((Real(0.5) + u) * 256);
        const std::uint8_t G = std::floor(u               * 256);
        return Pixel(R, G, 0);
    }
    else if(v < Real(0.5) + sqrt2)
    {
        const Real u = (v - Real(0.5) * sqrt2_plus_1) / sqrt2;
        if     (u  < Real(0.0)) {return Pixel(255, 127,   0);}
        else if(Real(0.5) <= u) {return Pixel(255, 255, 127);}
        const std::uint8_t G = std::floor((Real(0.5) + u) * 256);
        const std::uint8_t B = std::floor(u               * 256);
        return Pixel(255, G, B);
    }
    else
    {
        const Real u = v - Real(0.5) - sqrt2;
        if     (u < Real(0.0)) {return Pixel(255, 255, 127);}
        else if(Real(0.5) < u) {return Pixel(255, 255, 255);}
        const std::uint8_t B = std::floor((Real(0.5) + u * 256));
        return Pixel(255, 255, B);
    }
}


} // afmize
#endif// AFMIZE_COLOR_MAP_HPP
