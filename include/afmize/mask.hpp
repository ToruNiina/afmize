#ifndef AFMIZE_MASK_HPP
#define AFMIZE_MASK_HPP
#include <cstdint>
#include <vector>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>

// mask for image (stage). return true if masked.
namespace afmize
{

template<typename Real>
struct mask_nothing
{
    explicit mask_nothing(const system<Real>& sys)
        : pixel_x_(sys.stage_info.x_pixel()), pixel_y_(sys.stage_info.y_pixel())
    {}
    mask_nothing(const system<Real>& sys,
                      const std::size_t, const std::size_t,
                      const std::size_t, const std::size_t)
        : pixel_x_(sys.stage_info.x_pixel()), pixel_y_(sys.stage_info.y_pixel())
    {}

    Real operator()(const image<Real>& img,
                    const std::size_t x, const std::size_t y) const noexcept
    {
        return img(x, y);
    }

    std::size_t size() const noexcept {return pixel_x_ * pixel_y_;}

    std::size_t pixel_x() const noexcept {return pixel_x_;}
    std::size_t pixel_y() const noexcept {return pixel_y_;}

    std::size_t lower_bounding_x() const noexcept {return 0;}
    std::size_t lower_bounding_y() const noexcept {return 0;}
    std::size_t upper_bounding_x() const noexcept {return pixel_x_;}
    std::size_t upper_bounding_y() const noexcept {return pixel_y_;}

  private:

    std::size_t pixel_x_;
    std::size_t pixel_y_;
};

template<typename Real>
struct mask_by_rectangle
{
    // here, to make the area uniform, it uses sphere first
    explicit mask_by_rectangle(const system<Real>& mol)
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
        this->x_upper_ = std::min<std::int64_t>(mol.stage_info.x_pixel(), xth + hw_x + 1);
        this->y_lower_ = std::max<std::int64_t>(0, yth - hw_y);
        this->y_upper_ = std::min<std::int64_t>(mol.stage_info.y_pixel(), yth + hw_y + 1);

        this->pixel_x_ = x_upper_ - x_lower_;
        this->pixel_y_ = y_upper_ - y_lower_;
    }

    // XXX pixel at upper is not included, but lower is included.
    mask_by_rectangle(const system<Real>&,
                      const std::size_t x_lower, const std::size_t x_size,
                      const std::size_t y_lower, const std::size_t y_size)
        : pixel_x_(x_size),           pixel_y_(y_size),
          x_lower_(x_lower),          y_lower_(y_lower),
          x_upper_(x_lower + x_size), y_upper_(y_lower + y_size)
    {}

    std::size_t size() const noexcept {return pixel_x_ * pixel_y_;}

    Real operator()(const image<Real>& img,
                    const std::size_t x, const std::size_t y) const
    {
        if(x_upper_ <= x + x_lower_)
        {
            throw std::out_of_range("afmize::mask: condition x + x_lower (" +
                    std::to_string(x) + " + " + std::to_string(x_lower_) + ") <= x_upper (" +
                    std::to_string(x_upper_) + ") is not satisfied.");
        }
        if(y_upper_ <= y + y_lower_)
        {
            throw std::out_of_range("afmize::mask: condition y + y_lower (" +
                    std::to_string(y) + " + " + std::to_string(y_lower_) + ") <= y_upper (" +
                    std::to_string(y_upper_) + ") is not satisfied.");

        }
        return img(x + x_lower_, y + y_lower_);
    }

    std::size_t pixel_x() const noexcept {return pixel_x_;}
    std::size_t pixel_y() const noexcept {return pixel_y_;}

    std::size_t lower_bounding_x() const noexcept {return x_lower_;}
    std::size_t lower_bounding_y() const noexcept {return y_lower_;}
    std::size_t upper_bounding_x() const noexcept {return x_upper_;}
    std::size_t upper_bounding_y() const noexcept {return y_upper_;}

  private:

    std::size_t pixel_x_;
    std::size_t pixel_y_;
    std::size_t x_lower_;
    std::size_t y_lower_;
    std::size_t x_upper_;
    std::size_t y_upper_;
};

} // afmize
#endif// AFMIZE_MASK_HPP
