#ifndef AFMIZE_STAGE_HPP
#define AFMIZE_STAGE_HPP
#include <mave/mave/mave.hpp>
#include "image.hpp"
#include <vector>
#include <cassert>

namespace afmize
{

//
// stage information such as resolution and range
//
template<typename Real>
struct stage
{
    using container      = std::vector<Real>;
    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    stage(Real x_res, Real y_res, Real z_res,
          std::pair<Real, Real> x_range, std::pair<Real, Real> y_range)
        : x_reso_(x_res), y_reso_(y_res), z_reso_(z_res),
          x_rreso_(1.0 / x_res), y_rreso_(1.0 / y_res),
          x_lower_(x_range.first), x_upper_(x_range.second),
          y_lower_(y_range.first), y_upper_(y_range.second),
          x_pixels_(std::floor((x_range.second - x_range.first) / x_reso_)),
          y_pixels_(std::floor((y_range.second - y_range.first) / y_reso_))
    {}

    // for z, returns current height (at first it's zero)
    mave::vector<Real, 3> position_at(std::size_t x, std::size_t y) const noexcept
    {
        return mave::vector<Real, 3>{x_lower_ + x_reso_ * Real(x + 0.5),
                                     y_lower_ + y_reso_ * Real(y + 0.5), 0.0};
    }

    std::pair<std::int64_t, std::int64_t>
    pixel_at(const mave::vector<Real, 3>& pos) const noexcept
    {
        return std::make_pair(
            static_cast<std::int64_t>(std::floor((pos[0] - x_lower_) * x_rreso_)),
            static_cast<std::int64_t>(std::floor((pos[1] - y_lower_) * y_rreso_)));
    }

    image<Real> create_image() const
    {
        return image<Real>(x_pixels_, y_pixels_);
    }

    Real x_resolution() const noexcept {return this->x_reso_;}
    Real y_resolution() const noexcept {return this->y_reso_;}
    Real z_resolution() const noexcept {return this->z_reso_;}

    std::pair<Real, Real> x_range() const noexcept {return std::make_pair(x_lower_, x_upper_);}
    std::pair<Real, Real> y_range() const noexcept {return std::make_pair(y_lower_, y_upper_);}

    std::size_t x_pixel() const noexcept {return x_pixels_;}
    std::size_t y_pixel() const noexcept {return y_pixels_;}

  private:

    Real x_reso_, y_reso_, z_reso_;
    Real x_rreso_, y_rreso_;
    Real x_lower_, x_upper_, y_lower_, y_upper_;
    std::size_t x_pixels_, y_pixels_;
};

} // afmize
#endif// AFMIZE_STAGE_HPP
