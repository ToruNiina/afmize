#ifndef AFMIZE_STAGE_HPP
#define AFMIZE_STAGE_HPP
#include <afmize/parameter.hpp>
#include <afmize/collision.hpp>
#include <afmize/system.hpp>
#include <vector>
#include <cassert>

namespace afmize
{


//
// Stage is a map from (x, y) to z.
//
template<typename Real>
struct stage
{
    using container      = std::vector<Real>;
    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    stage(Real x_res, Real y_res, Real z_res,
          std::pair<Real, Real> x_range, std::pair<Real, Real> y_range)
        : x_reso(x_res), y_reso(y_res), z_reso(z_res),
          x_lower(x_range.first), x_upper(x_range.second),
          y_lower(y_range.first), y_upper(y_range.second),
          x_pixels_(std::floor((x_range.second - x_range.first) / x_reso)),
          y_pixels_(std::floor((y_range.second - y_range.first) / y_reso)),
          heights(x_pixels_ * y_pixels_, 0)
    {}

    Real& operator[](std::size_t i)       noexcept {return heights[i];}
    Real  operator[](std::size_t i) const noexcept {return heights[i];}
    Real& at(std::size_t i)       {return heights.at(i);}
    Real  at(std::size_t i) const {return heights.at(i);}

    Real& operator()(std::size_t x, std::size_t y)       noexcept {return heights[y*x_pixels_+x];}
    Real  operator()(std::size_t x, std::size_t y) const noexcept {return heights[y*x_pixels_+x];}
    Real& at(std::size_t x, std::size_t y)       {return heights.at(y*x_pixels_+x);}
    Real  at(std::size_t x, std::size_t y) const {return heights.at(y*x_pixels_+x);}

    iterator        begin()       noexcept {return heights.begin();}
    iterator        end()         noexcept {return heights.end();}
    const_iterator  begin() const noexcept {return heights.begin();}
    const_iterator  end()   const noexcept {return heights.end();}
    const_iterator cbegin() const noexcept {return heights.cbegin();}
    const_iterator cend()   const noexcept {return heights.cend();}

    // for z, returns current height (at first it's zero)
    mave::vector<Real, 3> position_at(std::size_t x, std::size_t y) const noexcept
    {
        return mave::vector<Real, 3>{x_lower + x_reso * Real(x + 0.5),
                                     y_lower + y_reso * Real(y + 0.5),
                                     heights[y * x_pixels_ + x]
        };
    }

    Real x_resolution() const noexcept {return this->x_reso;}
    Real y_resolution() const noexcept {return this->y_reso;}
    Real z_resolution() const noexcept {return this->z_reso;}

    std::pair<Real, Real> x_range() const noexcept {return std::make_pair(x_lower, x_upper);}
    std::pair<Real, Real> y_range() const noexcept {return std::make_pair(y_lower, y_upper);}

    std::size_t x_pixel() const noexcept {return x_pixels_;}
    std::size_t y_pixel() const noexcept {return y_pixels_;}

    std::vector<Real> const& get_container() const noexcept {return heights;}

  private:

    Real x_reso, y_reso, z_reso;
    Real x_lower, x_upper, y_lower, y_upper;
    std::size_t x_pixels_, y_pixels_;
    container heights;
};

} // afmize
#endif// AFMIZE_STAGE_HPP
