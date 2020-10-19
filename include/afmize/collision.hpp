#ifndef AFMIZE_COLLISION_HPP
#define AFMIZE_COLLISION_HPP
#include <afmize/shapes.hpp>
#include <limits>

namespace afmize
{

// if collides, returns distance in Z axis.
// otherwise,   returns NaN.
template<typename Real>
Real collision_z(const sphere<Real>& probe, const sphere<Real>& target) noexcept
{
    const auto dr = probe.center - target.center;
    const auto rr = probe.radius + target.radius;
    const auto zz = probe.center[2] - target.center[2];
    const auto D  = zz * zz - mave::length_sq(dr) + rr * rr;
    if(D < 0) {return std::numeric_limits<Real>::quiet_NaN();}

    return -zz + std::copysign(std::sqrt(D), zz);
}

template<typename Real>
Real collision_z(const circular_frustum<Real>& probe,
                 const sphere<Real>& target) noexcept
{
    const auto square  = [](const Real x) noexcept -> Real {return x * x;};
    const auto dist_xy = std::sqrt(square(probe.apex[0] - target.center[0]) +
                                   square(probe.apex[1] - target.center[1]));

    const auto cos_theta = std::cos(probe.angle);
    const auto threshold = probe.radius + cos_theta * target.radius;
    if(threshold < dist_xy)
    {
        if(probe.angle == 0.0) {return std::numeric_limits<Real>::quiet_NaN();}
        // collides at the lateral surface of circular frustum
        return target.center[2] -
            ((dist_xy - probe.radius - target.radius * cos_theta) /
              std::tan(probe.angle) - target.radius * std::sin(probe.angle))
            - probe.apex[2];
    }
    else if(probe.radius < dist_xy)
    {
        // collides at the edge of circular frustum
        return target.center[2] + std::sqrt(
                square(target.radius) - square(dist_xy - probe.radius)) -
                probe.apex[2];
    }
    else // trivial case.
    {
        return target.center[2] + target.radius - probe.apex[2];
    }
}

// check collision

template<typename Real>
bool collides_with(const sphere<Real>& probe, const sphere<Real>& target) noexcept
{
    const auto square = [](const Real x) noexcept -> Real {return x * x;};
    return mave::length_sq(probe.center - target.center) <
           square(probe.radius + target.radius);
}

template<typename Real>
bool collides_with(const circular_frustum<Real>& probe,
                   const sphere<Real>& target) noexcept
{
    const auto square  = [](const Real x) noexcept -> Real {return x * x;};
    const auto dist_xy = std::sqrt(square(probe.apex[0] - target.center[0]) +
                                   square(probe.apex[1] - target.center[1]));
    const auto dz = probe.apex[2] - target.center[2];

    const auto cos_theta = std::cos(probe.angle);
    const auto threshold = probe.radius + cos_theta * target.radius;
    if(threshold < dist_xy)
    {
        // collides at the lateral surface of circular frustum
        return cos_theta * (dist_xy - probe.radius + dz * std::tan(probe.angle))
                < target.radius;
    }
    else if(probe.radius < dist_xy)
    {
        // collides at the top of circular frustum
        if(dz < 0) {return true;}
        return (dz*dz + square(dist_xy - probe.radius)) < square(target.radius);
    }
    else // trivial case
    {
        if(dz < 0) {return true;}
        return dz < target.radius;
    }
}

} // afmize
#endif// AFMIZE_COLLISION_HPP
