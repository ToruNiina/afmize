#ifndef AFMIZE_IMAGE_HPP
#define AFMIZE_IMAGE_HPP
#include <vector>
#include <cassert>

namespace afmize
{

template<typename Real>
struct image
{
    using container      = std::vector<Real>;
    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    image(): x_pixels_(0), y_pixels_(0) {}
    image(const std::size_t x, const std::size_t y)
        : x_pixels_(x), y_pixels_(y), heights(x_pixels_ * y_pixels_, 0)
    {}

    Real& operator[](std::size_t i)       noexcept {return heights[i];}
    Real  operator[](std::size_t i) const noexcept {return heights[i];}
    Real& at(std::size_t i)       {return heights.at(i);}
    Real  at(std::size_t i) const {return heights.at(i);}

    Real& operator()(std::size_t x, std::size_t y)       noexcept
    {
        return heights[y * x_pixels_ + x];
    }
    Real  operator()(std::size_t x, std::size_t y) const noexcept
    {
        return heights[y * x_pixels_ + x];
    }
    Real& operator()(std::pair<std::size_t, std::size_t> xy)
    {
        return this->operator()(xy.first, xy.second);
    }
    Real  operator()(std::pair<std::size_t, std::size_t> xy) const
    {
        return this->operator()(xy.first, xy.second);
    }

    Real& at(std::size_t x, std::size_t y)
    {
        return heights.at(y * x_pixels_ + x);
    }
    Real  at(std::size_t x, std::size_t y) const
    {
        return heights.at(y * x_pixels_ + x);
    }

    Real& at(std::pair<std::size_t, std::size_t> xy)
    {
        return this->at(xy.first, xy.second);
    }
    Real  at(std::pair<std::size_t, std::size_t> xy) const
    {
        return this->at(xy.first, xy.second);
    }

    iterator        begin()       noexcept {return heights.begin();}
    iterator        end()         noexcept {return heights.end();}
    const_iterator  begin() const noexcept {return heights.begin();}
    const_iterator  end()   const noexcept {return heights.end();}
    const_iterator cbegin() const noexcept {return heights.cbegin();}
    const_iterator cend()   const noexcept {return heights.cend();}

    void resize(const std::size_t x, const std::size_t y)
    {
        heights.resize(x * y, 0.0);
        return;
    }

    std::size_t x_pixel() const noexcept {return x_pixels_;}
    std::size_t y_pixel() const noexcept {return y_pixels_;}

    std::vector<Real> const& get_container() const noexcept {return heights;}

  private:

    std::size_t x_pixels_, y_pixels_;
    container   heights;
};


} // afmize
#endif// AFMIZE_IMAGE_HPP
