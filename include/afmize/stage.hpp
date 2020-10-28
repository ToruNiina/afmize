#ifndef AFMIZE_STAGE_HPP
#define AFMIZE_STAGE_HPP
#include <afmize/parameter.hpp>
#include <afmize/collision.hpp>
#include <afmize/system.hpp>
#include <vector>
#include <cassert>

namespace afmize
{

template<typename Iterator>
struct range
{
    range(Iterator f, Iterator l): first(f), last(l) {}

    Iterator begin() const noexcept {return first;}
    Iterator end()   const noexcept {return last;}

  private:
    Iterator first, last;
};

template<typename Real>
struct cell_list
{
    struct list_elem
    {
        std::size_t particle_idx;
        std::size_t cell_idx;
        Real        z_coordinate;
    };

    cell_list(const Real grid_width_x, const Real grid_width_y,
              const std::pair<Real, Real>& lw, const std::pair<Real, Real>& up)
        : grid_width_x_(grid_width_x), grid_width_y_(grid_width_y),
          rwx_(1.0 / grid_width_x_), rwy_(1.0 / grid_width_y_),
          lower_(lw), upper_(up)
    {
        const auto width = up - lw;
        this->x_size_ = std::ceil(width[0] / grid_width_x);
        this->y_size_ = std::ceil(width[1] / grid_width_y);

        // It does not use any periodic boundary, so cells locating at the
        // boundary does not have 8 neighbors. So we make an always-empty cell
        // to represent pseudo-neighbors of such cells. This is the N+1-th cell.
        // So there are (x * y + 1) cells.
        //
        // 0        1    2,3             N,N+1
        // |        |     |              |
        // v        v     v              v
        // |i1|i2|i3|i4|i5|i6|i7|i8| ... |
        // '--------'-----'--------'...--'
        //   cell 0    1      3 ...
        //
        this->offsets_.resize(x_size_ * y_size_ + 2, 0);
        this->neighbors_.resize(x_size_ * y_size_);

        for(std::int64_t y=0; y<std::int64_t(y_size_); ++y)
        {
            for(std::int64_t x=0; x<std::int64_t(x_size_); ++x)
            {
                auto& neighbor = neighbors_.at(y * x_size_ + x);
                neighbor[0] = check_idx(x-1, y-1);
                neighbor[1] = check_idx(x  , y-1);
                neighbor[2] = check_idx(x+1, y-1);
                neighbor[3] = check_idx(x-1, y  );
                neighbor[4] = check_idx(x  , y  ); // includes itself
                neighbor[5] = check_idx(x+1, y  );
                neighbor[6] = check_idx(x-1, y+1);
                neighbor[7] = check_idx(x  , y+1);
                neighbor[8] = check_idx(x+1, y+1);
            }
        }
    }

    void construct(const std::vector<sphere<Real>>& particles)
    {
        this->list_.resize(particles.size());
        for(std::size_t i=0; i<particles.size(); ++i)
        {
            const auto pos = particles.at(i).center;
            list_.at(i) = list_elem{i, this->calc_index(pos), pos[2]};
        }

        // cell_idx: 0 -> N
        // z_coord:  large -> small
        std::sort(list_.begin(), list_.end(),
            [](const auto& lhs, const auto& rhs) -> bool {
                return (lhs.cell_idx == rhs.cell_idx) ? lhs.z_coordinate > rhs.z_coordinate :
                        lhs.cell_idx <  rhs.cell_idx;
            });

        // 0        1    2,3             N,N+1
        // |        |     |              |
        // v        v     v              v
        // |i1|i2|i3|i4|i5|i6|i7|i8| ... |
        // '--------'-----'--------'...--'
        //   cell 0    1      3 ...
        assert(offsets_.size() == x_size_ * y_size_ + 2);
        offsets_.at(0) = 0; // starting point
        offsets_.at(offsets_.size() - 2)  = list_.size(); // upper limit
        offsets_.at(offsets_.size() - 1)  = list_.size(); // empty N+1-th cell

        std::size_t current_offset = 0;
        std::size_t current_cell = 0;
        for(std::size_t i=0; i<list_.size(); ++i)
        {
            const auto& elem = list_.at(i);

            while(current_cell < elem.cell_idx)
            {
                offsets_.at(current_cell+1) = offsets_.at(current_cell) + current_offset;
                current_offset = 0;
                current_cell += 1;
            }
            current_offset += 1;
        }
        assert(offsets_.at(particles.size()) == list_.size());

        return ;
    }

    range<typename std::vector<list_elem>::const_iterator>
    cell_at(const mave::vector<Real, 3>& xy) const
    {
        const auto idx = calc_index(xy);
        return range<typename std::vector<list_elem>::const_iterator>{
            list_.begin() + offsets_.at(idx), list_.begin() + offsets_.at(idx + 1)};
    }

    std::array<std::size_t, 9> const&
    neighbors(const std::size_t i, const std::size_t j) const
    {
        return neighbors_.at(j * x_size_ + i);
    }

  private:

    std::size_t check_idx(std::int64_t x, std::int64_t y)
    {
        if(x < 0 || std::int64_t(x_size_) <= x ||
           y < 0 || std::int64_t(y_size_) <= y)
        {
            return x_size_ * y_size_; // the last one is always empty
        }
        return y * x_size_ + x;
    }

    std::size_t calc_index(const mave::vector<Real, 3>& v) const noexcept
    {
        const auto xi = static_cast<std::size_t>(std::floor((v[0] - lower_[0]) * rwx_));
        const auto yi = static_cast<std::size_t>(std::floor((v[1] - lower_[1]) * rwy_));
        return yi * x_size_ + xi;
    }

  private:

    Real grid_width_x_, grid_width_y_;
    Real rwx_, rwy_;
    std::size_t x_size_, y_size_;
    std::pair<Real, Real> lower_, upper_;

    std::vector<std::array<std::size_t, 8>> neighbors_;
    std::vector<std::size_t>                  offsets_;
    std::vector<list_elem>                       list_;
};



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
