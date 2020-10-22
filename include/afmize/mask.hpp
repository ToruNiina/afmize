#ifndef AFMIZE_MASK_HPP
#define AFMIZE_MASK_HPP
#include <cstdint>
#include <vector>
#include <iterator>

// mask for image (stage). return true if masked.
namespace afmize
{

template<typename Real>
struct mask_nothing
{
    explicit mask_nothing(const stage<Real>& stage_info)
        : x_lower_(0), y_lower_(0),
          x_upper_(stage_info.x_pixel()), y_upper_(stage_info.y_pixel())
    {}

    constexpr bool operator()(const std::size_t, const std::size_t) const noexcept
    {
        return false; // nothing is masked.
    }

    std::size_t lower_bounding_x() const noexcept {return x_lower_;}
    std::size_t lower_bounding_y() const noexcept {return y_lower_;}
    std::size_t upper_bounding_x() const noexcept {return x_upper_;}
    std::size_t upper_bounding_y() const noexcept {return y_upper_;}

  private:

    std::size_t x_lower_;
    std::size_t x_upper_;
    std::size_t y_lower_;
    std::size_t y_upper_;
};

template<typename Real>
struct mask_by_rectangle
{
    mask_by_rectangle(const stage<Real>& stage_info,
            const Real x_lower, const Real x_upper,
            const Real y_lower, const Real y_upper)
    {
        const Real rreso_x = 1.0 / stage_info.x_resolution();
        const Real rreso_y = 1.0 / stage_info.y_resolution();

        const Real stage_lw_x = stage_info.x_range().first;
        const Real stage_lw_y = stage_info.y_range().first;
        const Real stage_up_x = stage_info.x_range().second;
        const Real stage_up_y = stage_info.y_range().second;

        this->x_lower_ = std::max<std::int64_t>(0,               std::floor((x_lower - stage_lw_x) * rreso_x));
        this->y_lower_ = std::max<std::int64_t>(0,               std::floor((y_lower - stage_lw_y) * rreso_y));
        this->x_upper_ = std::min<std::int64_t>(stage.x_pixel(), std::ceil ((stage_up_x - x_upper) * rreso_x));
        this->y_upper_ = std::min<std::int64_t>(stage.y_pixel(), std::ceil ((stage_up_y - y_upper) * rreso_y));
    }

    bool operator()(const std::size_t x, const std::size_t y) const noexcept
    {
        return (x_lower_ <= x && x <= x_upper_) &&
               (y_lower_ <= y && y <= y_upper_);
    }

    std::size_t lower_bounding_x() const noexcept {return x_lower_;}
    std::size_t lower_bounding_y() const noexcept {return y_lower_;}
    std::size_t upper_bounding_x() const noexcept {return x_upper_;}
    std::size_t upper_bounding_y() const noexcept {return y_upper_;}

  private:

    std::size_t x_lower_;
    std::size_t y_lower_;
    std::size_t x_upper_;
    std::size_t y_upper_;
};

template<typename Real>
struct mask_by_circle
{
    mask_by_circle(const stage<Real>& stage_info,
                   const Real x_center, const Real y_center, const Real radius)
    : x_reso_ (stage_info.x_resolution()),  y_reso_ (stage_info.y_resolution()),
      x_lower_(stage_info.x_range().first), x_upper_(stage_info.x_range().second),
      y_lower_(stage_info.y_range().first), y_upper_(stage_info.y_range().second),
      x_center_(x_center), y_center_(y_center), radius_(radius), radius_sq_(radius * radius)
    {
        const Real x_lower = x_center - radius;
        const Real y_lower = y_center - radius;
        const Real x_upper = x_center + radius;
        const Real y_upper = y_center + radius;

        const Real rreso_x = 1.0 / stage_info.x_resolution();
        const Real rreso_y = 1.0 / stage_info.y_resolution();

        const Real stage_lw_x = stage_info.x_range().first;
        const Real stage_lw_y = stage_info.y_range().first;
        const Real stage_up_x = stage_info.x_range().second;
        const Real stage_up_y = stage_info.y_range().second;

        this->bb_x_lower_ = std::max<std::int64_t>(0,               std::floor((x_lower - stage_lw_x) * rreso_x));
        this->bb_y_lower_ = std::max<std::int64_t>(0,               std::floor((y_lower - stage_lw_y) * rreso_y));
        this->bb_x_upper_ = std::min<std::int64_t>(stage.x_pixel(), std::ceil ((stage_up_x - x_upper) * rreso_x));
        this->bb_y_upper_ = std::min<std::int64_t>(stage.y_pixel(), std::ceil ((stage_up_y - y_upper) * rreso_y));
    }

    bool operator()(const std::size_t x, const std::size_t y) const noexcept
    {
        if((x < bb_x_lower_ || bb_x_upper_ < x) ||
           (y < bb_y_lower_ || bb_y_upper_ < y))
        {
            return false;
        }

        // check pixel(rectangle)-circle overlap

        const Real pixel_min_x = x_lower_ + x_reso_ *  x;
        const Real pixel_max_x = x_lower_ + x_reso_ * (x+1);
        const Real pixel_min_y = y_lower_ + y_reso_ *  y;
        const Real pixel_max_y = y_lower_ + y_reso_ * (y+1);

        const Real nearest_x = std::max(pixel_min_x, std::min(x_center_, pixel_max_x));
        const Real nearest_y = std::max(pixel_min_y, std::min(y_center_, pixel_max_y));

        const Real dist_sq = (nearest_x - x_center_) * (nearest_x - x_center_) +
                             (nearest_y - y_center_) * (nearest_y - y_center_);

        return dist_sq <= radius_sq_;
    }

    std::size_t lower_bounding_x() const noexcept {return bb_x_lower_;}
    std::size_t lower_bounding_y() const noexcept {return bb_y_lower_;}
    std::size_t upper_bounding_x() const noexcept {return bb_x_upper_;}
    std::size_t upper_bounding_y() const noexcept {return bb_y_upper_;}

  private:

    std::size_t bb_x_lower_; // bounding box
    std::size_t bb_y_lower_;
    std::size_t bb_x_upper_;
    std::size_t bb_y_upper_;
    Real x_reso_,   y_reso_;
    Real x_lower_,  x_upper_, y_lower_, y_upper_;
    Real x_center_, y_center_, radius_, radius_sq_;
};

} // afmize
#endif// AFMIZE_MASK_HPP
