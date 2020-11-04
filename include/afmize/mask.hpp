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
    explicit mask_nothing(const system<Real>& sys)
        : x_lower_(0), y_lower_(0),
          x_upper_(sys.stage_info.x_pixel()), y_upper_(sys.stage_info.y_pixel())
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
    mask_by_rectangle(const system<Real>& mol)
        : mask_by_rectangle(mol, sphere<Real>{mol.bounding_radius,
            std::accumulate(
                mol.particles.begin(), mol.particles.end(), mave::vector<Real, 3>(0,0,0),
                [](const mave::vector<Real, 3>& s, const sphere<Real>& p) noexcept {
                    return s + p.center;
                }) / static_cast<Real>(mol.particles.size())
            })
    {}

    mask_by_rectangle(const system<Real>& mol, const sphere<Real>& bs)
    {
        const std::int64_t hw_x = std::ceil(bs.radius / mol.stage_info.x_resolution());
        const std::int64_t hw_y = std::ceil(bs.radius / mol.stage_info.x_resolution());

        const Real stage_lw_x = mol.stage_info.x_range().first;
        const Real stage_lw_y = mol.stage_info.y_range().first;

        const std::int64_t xth = std::ceil((bs.center[0] - stage_lw_x) / mol.stage_info.x_resolution());
        const std::int64_t yth = std::ceil((bs.center[1] - stage_lw_y) / mol.stage_info.y_resolution());
        // both ends are included
        this->x_lower_ = std::max<std::int64_t>(0, xth - hw_x);
        this->x_upper_ = std::min<std::int64_t>(mol.stage_info.x_pixel()-1, xth + hw_x);
        this->y_lower_ = std::max<std::int64_t>(0, yth - hw_y);
        this->y_upper_ = std::min<std::int64_t>(mol.stage_info.y_pixel()-1, yth + hw_y);
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

} // afmize
#endif// AFMIZE_MASK_HPP
