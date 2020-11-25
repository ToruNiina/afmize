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
    mask_nothing(const system<Real>& sys,
                 const std::size_t, const std::size_t,
                 const std::size_t, const std::size_t)
        : pixel_x_(sys.stage_info.x_pixel()), pixel_y_(sys.stage_info.y_pixel())
    {}

    explicit mask_nothing(const image<Real>& img)
        : pixel_x_(img.x_pixel()), pixel_y_(img.y_pixel())
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
    explicit mask_by_rectangle(const image<Real>& img) // nonzero
    {
        x_lower_ = std::numeric_limits<std::size_t>::max();
        x_upper_ = 0;
        y_lower_ = std::numeric_limits<std::size_t>::max();
        y_upper_ = 0;
        for(std::size_t j=0; j<img.y_pixel(); ++j)
        {
            for(std::size_t i=0; i<img.x_pixel(); ++i)
            {
                if(img.at(i, j) != 0.0)
                {
                    x_upper_ = std::max(x_upper_, i+1);
                    y_upper_ = std::max(y_upper_, j+1);
                    x_lower_ = std::min(x_lower_, i);
                    y_lower_ = std::min(y_lower_, j);
                }
            }
        }
        pixel_x_ = x_upper_ - x_lower_;
        pixel_y_ = y_upper_ - y_lower_;
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
