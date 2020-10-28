#ifndef AFMIZE_CELL_LIST_HPP
#define AFMIZE_CELL_LIST_HPP
#include <afmize/shapes.hpp>
#include <afmize/collision.hpp>
#include <vector>
#include <cassert>
#include <iostream>

namespace afmize
{

template<typename Iterator>
struct range
{
    range(Iterator f, Iterator l): first(f), last(l) {}

    Iterator begin() const noexcept {return first;}
    Iterator end()   const noexcept {return last;}

    bool empty() {return first == last;}

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

    cell_list() : is_dirty_(true) {}

    void initialize(const Real stage_reso_x, const Real stage_reso_y,
                    const std::vector<sphere<Real>>& mol,
                    const std::size_t fine_factor = 2)
    {
        if( ! this->is_dirty_)
        {
            throw std::runtime_error("initialized twice");
        }

        this->is_dirty_ = false;
        this->max_radius_ = 0.0;
        for(const auto& sph : mol)
        {
            max_radius_ = std::max(max_radius_, sph.radius);
        }

        this->stage_reso_x_ = stage_reso_x;
        this->stage_reso_y_ = stage_reso_y;
        // Since the size of typical protein is in the order of 1~10 nm, most
        // of the pixels in an AFM image does not contain any particle. Thus
        // constructing a cell list that has the same grid size as the AFM
        // image might not be so efficient.
        //     Here we construct a minimal cell list that can contain molecule
        // regardless of the relative configuration (but we consider the
        // molecule is rigid-body).

        const auto bounding_sphere =
            make_bounding_sphere_centered_at_geometric_center(mol);

        // a square of which edge length is `bounding_sphere.radius * 2` is
        // enough large to containe the whole molecule.
        // The number of grids that can wrap the square is,

        const std::size_t grids_x = 2 * std::ceil(bounding_sphere.radius / stage_reso_x_);
        const std::size_t grids_y = 2 * std::ceil(bounding_sphere.radius / stage_reso_y_);

        // make grid a bit fine-grained to speed up more

        this->x_size_       = grids_x * fine_factor; // number of cells in x
        this->y_size_       = grids_y * fine_factor;
        this->cell_width_x_ = bounding_sphere.radius * 2 / x_size_;
        this->cell_width_y_ = bounding_sphere.radius * 2 / y_size_;
        this->rcw_x_        = 1.0 / cell_width_x_;
        this->rcw_y_        = 1.0 / cell_width_y_;
        this->region_x_     = bounding_sphere.radius * 2;
        this->region_y_     = bounding_sphere.radius * 2;

        // construct data container
        //
        // .--------- N+1 elems ---------.
        // 0        1     2, 3           N, N+1 (N-th element is always empty)
        // |        |     |--'           |---'
        // v        v     v              v
        // |i1|i2|i3|i4|i5|i6|i7|i8| ... |
        // '--------'-----'--------'...--'
        //   cell 0    1  |   3 ...
        //                +- 2nd cell is empty
        //
        // We construct an offset list that indicates which range corresponds to
        // the cell. Here, we add an always-empty cell as the last element to
        // represent out-of-bound region. So here we need to allocate N+1 +1
        // offsets.
        this->offsets_.resize(x_size_ * y_size_ + 2, 0);
        this->list_.resize(mol.size());
    }

    // construct cell list.
    // system must be the same (the orientation and position can change).
    void construct(const std::vector<sphere<Real>>& mol,
                   const aabb<Real>& bounding_box /* of mol */)
    {
        assert(!is_dirty_);
        assert(list_.size() == mol.size());

        // XXX here it assumes system::bounding_box is up to date.
        // It wraps a bit larger region,
        this->lower_ = bounding_box.lower;
        this->upper_[0] = this->lower_[0] + this->region_x_;
        this->upper_[1] = this->lower_[1] + this->region_y_;
        this->upper_[2] = bounding_box.upper[2];

        // calc cell indices
        for(std::size_t i=0; i<mol.size(); ++i)
        {
            const auto pos = mol[i].center;
            list_[i] = list_elem{i, this->calc_index(pos), pos[2]};
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
        offsets_.at(offsets_.size() - 1)  = list_.size(); // empty N+1-th cell

        std::size_t current_offset = 0;
        std::size_t current_cell   = 0;
        for(std::size_t i=0; i<list_.size(); ++i)
        {
            const auto& elem = list_.at(i);

            while(current_cell < elem.cell_idx)
            {
                offsets_.at(current_cell+1) = offsets_.at(current_cell) + current_offset;
                current_offset = 0;
                current_cell  += 1;
            }
            current_offset += 1;
        }
        while(current_cell + 1 < offsets_.size())
        {
            // no particles left. all the resting cells are empty.
            // don't forget the last one.
            offsets_.at(current_cell+1) = offsets_.at(current_cell) + current_offset;
            current_offset = 0;
            current_cell += 1;
        }
        assert(offsets_.at(offsets_.size() - 1) == list_.size());
        assert(offsets_.at(offsets_.size() - 2) == list_.size());

        return ;
    }

    range<typename std::vector<list_elem>::const_iterator>
    cell_at(const mave::vector<Real, 3>& xy) const
    {
        const auto x = static_cast<std::int32_t>(std::floor((xy[0] - lower_[0]) * rcw_x_));
        const auto y = static_cast<std::int32_t>(std::floor((xy[1] - lower_[1]) * rcw_y_));
        return cell_at(check_idx(x, y));
    }
    range<typename std::vector<list_elem>::const_iterator>
    cell_at(const std::size_t idx) const
    {
        return range<typename std::vector<list_elem>::const_iterator>{
                list_.begin() + offsets_.at(idx),
                list_.begin() + offsets_.at(idx + 1)
            };
    }

    // here, we don't care about the z coordinate of the probe.
    // sphere is not "open". It might collides with the particles above
    // while observation.
    void overwrapping_cells(const sphere<Real>& probe, std::vector<std::size_t>& out) const
    {
        out.clear();
        const std::int32_t ctr_x = std::floor((probe.center[0] - lower_[0]) * rcw_x_);
        const std::int32_t ctr_y = std::floor((probe.center[1] - lower_[1]) * rcw_y_);
        const std::int32_t ofs_x = (probe.radius + max_radius_) * rcw_x_;
        const std::int32_t ofs_y = (probe.radius + max_radius_) * rcw_y_;

        const auto nil = x_size_ * y_size_;
        for(std::int32_t y = ctr_y - ofs_y; y <= ctr_y + ofs_y; ++y)
        {
        for(std::int32_t x = ctr_x - ofs_x; x <= ctr_x + ofs_x; ++x)
        {
            const auto idx = check_idx(x, y);
            if(idx != nil)
            {
                out.push_back(idx);
            }
        }
        }
        return;
    }

    void overwrapping_cells(const circular_frustum<Real>& probe,
                            std::vector<std::size_t>& out) const
    {
        out.clear();

        if(this->upper_[2] < probe.apex[2])
        {
            return;
        }
        const auto r = probe.radius + max_radius_ +
            std::tan(probe.angle) * (this->upper_[2] - probe.apex[2]);

        const std::int32_t ctr_x = std::floor((probe.apex[0] - lower_[0]) * rcw_x_);
        const std::int32_t ctr_y = std::floor((probe.apex[1] - lower_[1]) * rcw_y_);
        const std::int32_t ofs_x = r * rcw_x_;
        const std::int32_t ofs_y = r * rcw_y_;

        const auto nil = x_size_ * y_size_;
        for(std::int32_t y = ctr_y - ofs_y; y <= ctr_y + ofs_y; ++y)
        {
        for(std::int32_t x = ctr_x - ofs_x; x <= ctr_x + ofs_x; ++x)
        {
            const auto idx = check_idx(x, y);
            if(idx != nil)
            {
                out.push_back(idx);
            }
        }
        }
        return;
    }

    void diagnosis(const std::vector<sphere<Real>>& mol) const
    {
        for(std::int32_t y=0; y<y_size_; ++y)
        {
            for(std::int32_t x=0; x<x_size_; ++x)
            {
                const auto elems = cell_at(check_idx(x, y));
                if( ! std::is_sorted(elems.begin(), elems.end(),
                            [](const auto& lhs, const auto& rhs) -> bool {
                                return lhs.z_coordinate > rhs.z_coordinate;
                            }))
                {
                    throw std::runtime_error("z coord is not sorted");
                }
                for(const auto& elem : elems)
                {
                    const auto p = mol.at(elem.particle_idx);

                    const auto x_ = std::floor((p.center[0] - lower_[0]) * rcw_x_);
                    const auto y_ = std::floor((p.center[1] - lower_[1]) * rcw_y_);

                    if(x != x_ || y != y_)
                    {
                        throw std::runtime_error("particle position inconsistent:"
                            " cell = (" + std::to_string(x) + ", " +
                                          std::to_string(y) + ")," +
                            " particle = (" + std::to_string(x_) + ", " +
                                              std::to_string(y_) + ").");
                    }
                }
            }
        }
        return ;
    }

  private:

    std::size_t check_idx(std::int32_t x, std::int32_t y) const
    {
        if(x < 0 || std::int32_t(x_size_) <= x ||
           y < 0 || std::int32_t(y_size_) <= y)
        {
            return x_size_ * y_size_; // the last one is always empty
        }
        return y * x_size_ + x;
    }

    std::size_t calc_index(const mave::vector<Real, 3>& v) const noexcept
    {
        const auto xi = static_cast<std::size_t>(std::floor((v[0] - lower_[0]) * rcw_x_));
        const auto yi = static_cast<std::size_t>(std::floor((v[1] - lower_[1]) * rcw_y_));
        assert(0 <= xi && xi < x_size_ && 0 <= yi && yi < y_size_);
        return yi * x_size_ + xi;
    }

  private:

    bool is_dirty_;
    Real max_radius_; // in system
    Real stage_reso_x_, stage_reso_y_; // width of a image pixel
    Real cell_width_x_, cell_width_y_; // width of a cell
    Real rcw_x_, rcw_y_;               // reciprocal cell width
    Real region_x_, region_y_;         // width of the whole region (fixed)
    std::size_t x_size_, y_size_;

    mave::vector<Real, 3> lower_, upper_;
    std::vector<std::size_t> offsets_;
    std::vector<list_elem>   list_;
};

} // afmize
#endif // AFMIZE_CELL_LIST_HPP
