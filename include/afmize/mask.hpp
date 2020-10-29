#ifndef AFMIZE_MASK_HPP
#define AFMIZE_MASK_HPP
#include <cstdint>
#include <vector>
#include <iterator>
#include <numeric>

// mask for image (stage). return true if masked.
namespace afmize
{

template<typename Real>
struct mask_nothing
{
    explicit mask_nothing(const stage<Real>& stage_info, const system<Real>&)
        : x_lower_(0), y_lower_(0),
          x_upper_(stage_info.x_pixel()), y_upper_(stage_info.y_pixel())
    {}

    constexpr bool operator()(const std::size_t, const std::size_t) const noexcept
    {
        return false; // nothing is masked.
    }

    std::size_t lower_bounding_x() const noexcept {return x_lower_;}
    std::size_t lower_bounding_y() const noexcept {return y_lower_;}
    std::size_t upper_bounding_x() const noexcept {return x_upper_ - 1;}
    std::size_t upper_bounding_y() const noexcept {return y_upper_ - 1;}

  private:

    std::size_t x_lower_;
    std::size_t x_upper_;
    std::size_t y_lower_;
    std::size_t y_upper_;
};

template<typename Real>
struct mask_by_rectangle
{
    // here, to make the area uniform, it uses sphere first
    mask_by_rectangle(const stage<Real>& stage_info, const system<Real>& mol)
        : mask_by_rectangle(stage_info, sphere<Real>{mol.bounding_radius,
            std::accumulate(
                mol.particles.begin(), mol.particles.end(), mave::vector<Real, 3>(0,0,0),
                [](const mave::vector<Real, 3>& s, const sphere<Real>& p) noexcept {
                    return s + p.center;
                }) / static_cast<Real>(mol.particles.size())
            })
    {}

    mask_by_rectangle(const stage<Real>& stage_info, const sphere<Real>& mol)
    {
        const std::int64_t hw_x = std::ceil(mol.radius / stage_info.x_resolution());
        const std::int64_t hw_y = std::ceil(mol.radius / stage_info.x_resolution());

        const Real stage_lw_x = stage_info.x_range().first;
        const Real stage_lw_y = stage_info.y_range().first;

        const std::int64_t xth = std::ceil((mol.center[0] - stage_lw_x) / stage_info.x_resolution());
        const std::int64_t yth = std::ceil((mol.center[1] - stage_lw_y) / stage_info.y_resolution());
        // both ends are included
        this->x_lower_ = std::max<std::int64_t>(0, xth - hw_x);
        this->x_upper_ = std::min<std::int64_t>(stage_info.x_pixel()-1, xth + hw_x);
        this->y_lower_ = std::max<std::int64_t>(0, yth - hw_y);
        this->y_upper_ = std::min<std::int64_t>(stage_info.y_pixel()-1, yth + hw_y);
    }

    bool operator()(const std::size_t x, const std::size_t y) const noexcept
    {
        // mask outside
        return (x < x_lower_ || x_upper_ < x) ||
               (y < y_lower_ || y_upper_ < y);
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

// template<typename Real>
// struct mask_by_circle
// {
//     mask_by_circle(const stage<Real>& stage_info, const system<Real>& mol)
//         : mask_by_rectangle(stage_info,
//             make_bounding_sphere_centered_at_geometric_center(mol.particles))
//     {}
// 
//     mask_by_circle(const stage<Real>& stage_info, const sphere<Real>& sph)
//         : mask_by_rectangle(stage_info, sph.center[0], sph.center[1], sph.radius)
//     {}
// 
//     bool operator()(const std::size_t x, const std::size_t y) const noexcept
//     {
//         if((x < bb_x_lower_ || bb_x_upper_ < x) ||
//            (y < bb_y_lower_ || bb_y_upper_ < y))
//         {
//             return false;
//         }
// 
//         // check pixel(rectangle)-circle overlap
// 
//         const Real pixel_min_x = x_lower_ + x_reso_ *  x;
//         const Real pixel_max_x = x_lower_ + x_reso_ * (x+1);
//         const Real pixel_min_y = y_lower_ + y_reso_ *  y;
//         const Real pixel_max_y = y_lower_ + y_reso_ * (y+1);
// 
//         const Real nearest_x = std::max(pixel_min_x, std::min(x_center_, pixel_max_x));
//         const Real nearest_y = std::max(pixel_min_y, std::min(y_center_, pixel_max_y));
// 
//         const Real dist_sq = (nearest_x - x_center_) * (nearest_x - x_center_) +
//                              (nearest_y - y_center_) * (nearest_y - y_center_);
// 
//         return dist_sq <= radius_sq_;
//     }
// 
//     std::size_t lower_bounding_x() const noexcept {return bb_x_lower_;}
//     std::size_t lower_bounding_y() const noexcept {return bb_y_lower_;}
//     std::size_t upper_bounding_x() const noexcept {return bb_x_upper_;}
//     std::size_t upper_bounding_y() const noexcept {return bb_y_upper_;}
// 
//   private:
// 
//     std::size_t bb_x_lower_; // bounding box
//     std::size_t bb_y_lower_;
//     std::size_t bb_x_upper_;
//     std::size_t bb_y_upper_;
//     Real x_reso_,   y_reso_;
//     Real x_lower_,  x_upper_, y_lower_, y_upper_;
//     Real x_center_, y_center_, radius_, radius_sq_;
// };

} // afmize
#endif// AFMIZE_MASK_HPP
