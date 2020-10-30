#ifndef AFMIZE_OUTPUT_UTILITY_HPP
#define AFMIZE_OUTPUT_UTILITY_HPP
#include "system.hpp"
#include "image.hpp"
#include "colormap.hpp"
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
    std::ofstream ofs(basename + ".tsv"s);

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

    write_ppm(out, img, std::make_pair(min_elem, max_elem));
    return;
}

template<typename Real>
void write_json(const std::string& out, const image<Real>& img, const system<Real>& sys)
{
    using namespace std::literals::string_literals;

    std::ofstream ofs(out + ".json");
    ofs << "{\n\t\"resolution\":{\"x\":" << sys.stage_info.x_resolution()
        << ", \"y\":" << sys.stage_info.y_resolution()
        << ", \"z\":" << sys.stage_info.z_resolution() << "},\n";
    ofs << "\t\"height\":[\n";
    for(std::size_t j=0; j<img.y_pixel(); ++j)
    {
        for(std::size_t i=0; i<img.x_pixel(); ++i)
        {
            const auto p = sys.stage_info.position_at(i, j);
            ofs << "\t\t{\"x\":" << p[0] << ", \"y\":" << p[1] << ", \"z\":" << p[2] << "}";
            if(j == img.y_pixel()-1 && i == img.x_pixel()-1)
            {
                ofs << "\n";
            }
            else {ofs << ",\n";}
        }
    }
    ofs << "\t]\n}";
    return;
}

template<typename Real>
void write_csv(const std::string& out, const image<Real>& img)
{
    using namespace std::literals::string_literals;

    std::ofstream ofs(out + ".csv");
    for(std::size_t j=0; j<img.y_pixel(); ++j)
    {
        for(std::size_t i=0; i<img.x_pixel(); ++i)
        {
            ofs << img(i, j);
            if(i+1 != img.x_pixel()) {ofs << ',';}
        }
        ofs << '\n';
    }
    return;
}

template<typename Real>
void write_svg(const std::string& out, const image<Real>& img, const system<Real>& sys,
               const Real scale_bar, const std::pair<Real, Real> height_range)
{
    using namespace std::literals::string_literals;
    std::ofstream svg(out + ".svg");
    if(!svg.good())
    {
        throw std::runtime_error("file open error: " + out);
    }

    const auto img_width  = img.x_pixel() * sys.stage_info.x_resolution();
    const auto img_height = img.y_pixel() * sys.stage_info.y_resolution();

    svg << "<svg width=\"" << img_width << "\" height=\"" << img_height << "\">\n";

    const auto minv = height_range.first;
    const auto maxv = height_range.second;

    for(std::size_t yi=0; yi<img.y_pixel(); ++yi)
    {
        // origin of the image is upper left, but the physical origin is lower left.
        // it makes the image physically "correct" (here, correct means top-view)
        const auto y_pos = (img.y_pixel() - yi - 1) * sys.stage_info.y_resolution();
        for(std::size_t xi=0; xi<img.x_pixel(); ++xi)
        {
            const auto x_pos = xi * sys.stage_info.x_resolution();
            const auto color = color_afmhot<Real>(img.at(xi, yi), minv, maxv);
            svg << "<rect x=\""   << x_pos << "\" y=\"" << y_pos
                << "\" width=\""  << sys.stage_info.x_resolution()
                << "\" height=\"" << sys.stage_info.y_resolution() << "\" style=\""
                << "fill:rgb(" << static_cast<int>(color.red)   << ','
                               << static_cast<int>(color.green) << ','
                               << static_cast<int>(color.blue)
                << ");stroke:none\"/>\n";
        }
    }

    if(scale_bar != Real(0.0))
    {
        const auto sb_width  = scale_bar;
        const auto sb_height = std::max(sys.stage_info.y_resolution() * 0.5, img_height * 0.01);

        const auto buf_x = img_width  * 0.05;
        const auto buf_y = img_height * 0.05;

        const auto sb_left = img_width  - buf_x - sb_width;
        const auto sb_up   = img_height - buf_y - sb_height;

        // scale bar
        svg << "<rect x=\""   << sb_left  << "\" y=\""      << sb_up
            << "\" width=\""  << sb_width << "\" height=\"" << sb_height
            << "\" style=\"fill:white;stroke:none\"/>\n";
    }
    svg << "</svg>\n";
    return;
}

template<typename Real>
void write_svg(const std::string& out, const image<Real>& img, const system<Real>& sys,
               const Real scale_bar)
{
    const auto minmax = std::minmax_element(img.begin(), img.end());
    const auto minv = *minmax.first;
    const auto maxv = *minmax.second;
    write_svg(out, img, sys, scale_bar, std::make_pair(minv, maxv));
    return;
}



} // afmize
#endif// AFMIZE_OUTPUT_UTILITY_HPP
