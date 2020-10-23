#ifndef AFMIZE_OBSERVE_HPP
#define AFMIZE_OBSERVE_HPP
#include "shape.hpp"
#include "system.hpp"
#include "stage.hpp"

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

template<typename Real, bool Descritize = true>
struct rigid_observer
{
    void operator()(stage<Real>& stg, const systme<Real>& sys,
                    const Real bottom)
    {
        for(std::size_t j=0; j<stg.y_pixel(); ++j)
        {
            for(std::size_t i=0; i<stg.x_pixel(); ++i)
            {
                probe.apex    = stg.position_at(i, j);
                probe.apex[2] = initial_z;
                const auto height = afmize::collide_at(sys, probe, bottom);
                if (Descritize)
                {
                    stg(i, j) = afmize::discretize(height, stg.z_resolution(), bottom);
                }
                else
                {
                    stg(i, j) = height;
                }
            }
        }
        return;
    }
    default_probe<Real> probe;
};


} // afmize
#endif// AFMIZE_OBSERVE_HPP
