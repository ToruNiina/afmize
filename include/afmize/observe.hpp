#ifndef AFMIZE_OBSERVE_HPP
#define AFMIZE_OBSERVE_HPP
#include "cell_list.hpp"
#include "collision.hpp"
#include "shapes.hpp"
#include "system.hpp"
#include "stage.hpp"
#include <map>

namespace afmize
{

template<typename Real>
Real discretize(const Real x, const Real resolution, const Real min_val)
{
    return std::round((x - min_val) / resolution) * resolution + min_val;
}

template<typename Real>
Real collide_at(const system<Real>& sys, const default_probe<Real>& probe,
                const Real bottom)
{
    // the point at which the probe collides with the stage
    Real height = bottom + probe.radius;
    for(const auto& sph : sys.particles)
    {
        const Real dz_sph = collision_z(sphere<Real>{
                probe.radius, probe.apex
            }, sph);
        const Real dz_frs = collision_z(circular_frustum<Real>{
                probe.angle, probe.radius, probe.apex
            }, sph);
        if(!std::isnan(dz_frs))
        {
            height = std::max(height, probe.apex[2] + dz_frs);
        }
        if(!std::isnan(dz_sph))
        {
            height = std::max(height, probe.apex[2] + dz_sph);
        }
    }
    // subtract probe radius to set the ground as 0
    return height - probe.radius;
}

template<typename Real>
Real smooth_at(const system<Real>& sys,
               const mave::vector<Real, 3>& pos, const Real bottom,
               const Real gamma, const Real sigma_x, const Real sigma_y)
{
    const Real rgamma    = 1.0 / gamma;
    const Real rsigma_x    = 1.0 / sigma_x;
    const Real rsigma_x_sq = rsigma_x * rsigma_x;
    const Real rsigma_y    = 1.0 / sigma_y;
    const Real rsigma_y_sq = rsigma_y * rsigma_y;

    Real expsum = std::exp(-bottom * rgamma);
    for(const auto& p : sys.particles)
    {
        const auto dr = p.center - pos;

        const auto dx_sq = dr[0] * dr[0];
        const auto dy_sq = dr[1] * dr[1];

        expsum += std::exp(-(dx_sq * rsigma_x_sq + dy_sq * rsigma_y_sq) +
                            (p.radius + p.center[2] - bottom) * rgamma);
    }
    return gamma * std::log(expsum);
}

template<typename Real>
struct ObserverBase
{
    virtual ~ObserverBase() = default;
    virtual image<Real> const& observe(const system<Real>&) = 0;
    virtual bool update_probe(const std::map<std::string, Real>&) = 0;
    virtual std::map<std::string, Real> get_probe() const = 0;
    virtual void print_probe(std::ostream&) const = 0;
    virtual image<Real> const&  get_image() const = 0;
};

template<typename Real, bool Descritize>
struct RigidObserver: public ObserverBase<Real>
{
    explicit RigidObserver(default_probe<Real> p)
        : probe(std::move(p))
    {}
    ~RigidObserver() override = default;

    // here we assume the stage locates z == 0.
    image<Real> const& observe(const system<Real>& sys) override
    {
        img_.resize(sys.stage_info.x_pixel(), sys.stage_info.y_pixel());

        // to skip pixels, we first calculate aabb of the probe at the
        const Real max_frustum_radius = probe.radius +
            std::tan(probe.angle) * (sys.bounding_box.upper[2] - probe.radius);

        const Real initial_z = sys.bounding_box.upper[2] + probe.radius;
        for(std::size_t j=0; j<sys.stage_info.y_pixel(); ++j)
        {
        for(std::size_t i=0; i<sys.stage_info.x_pixel(); ++i)
        {
            probe.apex    = sys.stage_info.position_at(i, j);
            probe.apex[2] = initial_z;

            aabb<Real> probe_aabb;
            probe_aabb.upper[0] = probe.apex[0] + max_frustum_radius;
            probe_aabb.upper[1] = probe.apex[1] + max_frustum_radius;
            probe_aabb.upper[2] = sys.bounding_box.upper[2];

            probe_aabb.lower[0] = probe.apex[0] - max_frustum_radius;
            probe_aabb.lower[1] = probe.apex[1] - max_frustum_radius;
            probe_aabb.lower[2] = 0.0;

            if( ! collides_with(probe_aabb, sys.bounding_box))
            {
                img_(i, j) = 0.0;
                continue;
            }

//             const auto min_height = collide_at(sys, probe, 0.0);


            // screen particles and get initial guess of the (minimum) height
            Real min_height = 0.0;
            for(const auto& elem : sys.cells.cell_at(probe.apex))
            {
                const auto& particle = sys.particles.at(elem.particle_idx);
                const auto  dz_sph   = collision_z(
                        sphere<Real>{probe.radius, probe.apex}, particle);
                if(!std::isnan(dz_sph))
                {
                    min_height = std::max(min_height, probe.apex[2] + dz_sph - probe.radius);
                    break;
                }
            }
            // Now we have initial guess. Next we check additional collision.
            probe.apex[2] = min_height + probe.radius;

            // additional collision with circular frustum
            sys.cells.overwrapping_cells(
                circular_frustum<Real>{probe.angle, probe.radius, probe.apex},
                this->index_buffer_);

            for(const auto cell_idx : index_buffer_)
            {
                for(const auto& elem : sys.cells.cell_at(cell_idx))
                {
                    if(elem.z_coordinate + sys.max_radius < min_height)
                    {
                        // elements are sorted by its z coordinate.
                        break;
                    }
                    const auto& particle = sys.particles.at(elem.particle_idx);

                    const auto dz = collision_z(
                            circular_frustum<Real>{probe.angle, probe.radius, probe.apex},
                            particle);

                    if(!std::isnan(dz)) // optional requires C++17 and it's 14...
                    {
                        // at the tip, we always have a (half) sphere.
                        // the top of the tip is always lower than the apex of
                        // circular frustum in probe.radius.
                        min_height = std::max(min_height, probe.apex[2] + dz - probe.radius);
                    }
                }
            }

            // additional collision with tip sphere
            probe.apex[2] = initial_z;
            sys.cells.overwrapping_cells(sphere<Real>{probe.radius, probe.apex},
                                         this->index_buffer_);
            for(const auto cell_idx : index_buffer_)
            {
                for(const auto& elem : sys.cells.cell_at(cell_idx))
                {
                    if(elem.z_coordinate + sys.max_radius < min_height)
                    {
                        break; // elements are sorted by its z coordinate.
                    }
                    const auto& particle = sys.particles.at(elem.particle_idx);

                    const auto dz_sph = collision_z(
                            sphere<Real>{probe.radius, probe.apex}, particle);

                    if(!std::isnan(dz_sph))
                    {
                        min_height = std::max(min_height, probe.apex[2] + dz_sph - probe.radius);
                    }
                }
            }

            if (Descritize)
            {
                img_(i, j) = afmize::discretize(min_height, sys.stage_info.z_resolution(), Real(0));
            }
            else
            {
                img_(i, j) = min_height;
            }
        }
        }
        return img_;
    }

    bool update_probe(const std::map<std::string, Real>& attr) override
    {
        bool is_updated = false;
        if(0.0 < attr.at("radius"))
        {
            probe.radius = attr.at("radius");
            is_updated   = true;
        }
        if(0.0 < attr.at("angle") && attr.at("angle") < 3.1416 * 0.5)
        {
            probe.angle = attr.at("angle");
            is_updated  = true;
        }
        return is_updated;
    }

    void print_probe(std::ostream& os) const override
    {
        os << probe.radius * 0.1 << " [nm] " << probe.angle * 180.0 / 3.1416 << " [rad]";
    }

    std::map<std::string, Real> get_probe() const override
    {
        return std::map<std::string, Real>{
            {"radius", probe.radius},
            {"angle",  probe.angle }
        };
    }

    image<Real> const&  get_image() const override {return img_;}

    default_probe<Real> probe;

  private:
    image<Real> img_;
    std::vector<std::size_t> index_buffer_; // avoid allocation
};

template<typename Real, bool Descritize>
struct SmoothObserver: public ObserverBase<Real>
{
    explicit SmoothObserver(const Real sigma,
            const Real gamma = 1.0 /*angst.*/, const Real cutoff_sigma = 5.0)
        : sigma_(sigma), gamma_(gamma), rgamma_(1.0 / gamma), cutoff_(cutoff_sigma)
    {}
    ~SmoothObserver() override = default;

    // here we assume the stage locates z == 0.
    image<Real> const& observe(const system<Real>& sys) override
    {
        img_.resize(sys.stage_info.x_pixel(), sys.stage_info.y_pixel());

        const auto rsigma    = 1.0 / sigma_;
        const auto rsigma_sq = rsigma * rsigma;
        const auto cutoff_sq = (sigma_ * cutoff_) * (sigma_ * cutoff_);

        std::vector<std::size_t> cells;
        cells.reserve(25);

        for(std::size_t j=0; j<sys.stage_info.y_pixel(); ++j)
        {
        for(std::size_t i=0; i<sys.stage_info.x_pixel(); ++i)
        {
            const auto pos = sys.stage_info.position_at(i, j);

            Real expsum = 1.0; // exp(z_0 = 0.0)

            sys.cells.overwrapping_cells(sigma_ * cutoff_, sigma_ * cutoff_, i, j, cells);
            for(const auto& cell_idx : cells)
            {
                for(const auto& elem : sys.cells.cell_at(cell_idx))
                {
                    const auto& p = sys.particles.at(elem.particle_idx);
                    const auto dr = p.center - pos;
                    const auto dx_sq = dr[0] * dr[0];
                    const auto dy_sq = dr[1] * dr[1];

                    if(cutoff_sq < dx_sq || cutoff_sq < dy_sq) {continue;}

                    expsum += std::exp(-(dx_sq * rsigma_sq + dy_sq * rsigma_sq) +
                                        (p.radius + p.center[2]) * rgamma_);
                }
            }
            const auto height = gamma_ * std::log(expsum);

            if (Descritize)
            {
                img_(i, j) = afmize::discretize(height, sys.stage_info.z_resolution(), Real(0));
            }
            else
            {
                img_(i, j) = height;
            }
        }
        }
        return img_;
    }

    bool update_probe(const std::map<std::string, Real>& attr) override
    {
        bool is_updated = false;
        if(0.0 < attr.at("sigma"))
        {
            sigma_ = attr.at("sigma");
            is_updated = true;
        }
        return is_updated;
    }
    void print_probe(std::ostream& os) const override
    {
        os << sigma_ * 0.1 << " [nm]";
    }
    std::map<std::string, Real> get_probe() const override
    {
        return std::map<std::string, Real>{
            {"sigma", sigma_},
        };
    }

    image<Real> const&  get_image() const override {return img_;}

  private:

    Real sigma_;
    Real gamma_;
    Real rgamma_;
    Real cutoff_;

    image<Real> img_;
    std::vector<std::size_t> index_buffer_; // avoid allocation
};


} // afmize
#endif// AFMIZE_OBSERVE_HPP
