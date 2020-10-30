#ifndef AFMIZE_OUTPUT_UTILITY_HPP
#define AFMIZE_OUTPUT_UTILITY_HPP
#include "system.hpp"
#include "image.hpp"
#include <pnm/pnm.hpp>
#include <string>
#include <ostream>

namespace afmize
{

template<typename Real>
void write_xyz(const std::string& basename, const system<Real>& sys)
{
    std::ofstream ofs(basename + ".xyz", std::ios_base::app);
    ofs << sys.particles.size() << "\n\n";
    for(std::size_t i=0; i<sys.particles.size(); ++i)
    {
        const auto pos = sys.particles.at(i).center;
        ofs << std::setw(6) << sys.names.at(i) << ' '
            << std::setprecision(6) << std::setw(10) << pos[0]
            << std::setprecision(6) << std::setw(10) << pos[1]
            << std::setprecision(6) << std::setw(10) << pos[2]
            << '\n';
    }
    return ;
}

template<typename Real>
void write_tsv(const std::string& basename, const image<Real>& img)
{
    using namespace std::literals::string_literals;
    std::ofstream ofs(out + ".tsv"s);

    for(std::size_t y=0; y<img.y_pixel(); ++y)
    {
        for(std::size_t x=0; x<img.x_pixel(); ++x)
        {
            ofs << "\t" << std::setprecision(6) << std::setw(10) << img.at(x, y);
        }
        ofs << "\n";
    }
    return;
}

template<typename Real>
void write_ppm(const std::string& out, const image<Real>& img,
               const std::pair<Real, Real> height_range)
{
    using namespace std::literals::string_literals;
    const auto min_elem = height_range.first;
    const auto max_elem = height_range.second;

    pnm::image<pnm::rgb_pixel> ppm(img.x_pixel(), img.y_pixel());
    for(std::size_t i=0; i<img.x_pixel() * img.y_pixel(); ++i)
    {
        ppm.raw_access(i) =
            color_afmhot<Real, pnm::rgb_pixel>(img[i], min_elem, max_elem);
    }

    // origin of the image is upper left, but the physical origin is lower left.
    // to make it consistent, it simply reverses the image.
    pnm::image<pnm::rgb_pixel> reversed(ppm.x_size(), ppm.y_size());
    for(std::size_t i = 0; i<ppm.y_size(); ++i)
    {
        reversed[ppm.y_size() - i - 1] = ppm[i];
    }

    pnm::write(out + ".ppm"s, reversed, pnm::format::ascii);
    return;
}

template<typename Real>
void write_ppm(const std::string& out, const image<Real>& img)
{
    const auto minmax = std::minmax_element(img.begin(), img.end());
    const auto min_elem = *minmax.first;
    const auto max_elem = *minmax.second;

    write_ppm(img, out, std::make_pair(min_elem, max_elem));
    return;
}

} // afmize
#endif// AFMIZE_OUTPUT_UTILITY_HPP
