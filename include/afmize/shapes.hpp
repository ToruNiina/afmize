#ifndef AFMIZE_SHAPES_HPP
#define AFMIZE_SHAPES_HPP
#include <mave/mave/mave.hpp>
#include <vector>

namespace afmize
{

template<typename Real>
struct sphere
{
    Real radius;
    mave::vector<Real, 3> center;
};

//          angle
// \       |  /
//  \      | /
//   \_____|/
//       '--'
//       radius
//
template<typename Real>
struct circular_frustum
{
    Real angle;
    Real radius;
    mave::vector<Real, 3> apex;
};

// default probe shape is sphere + frustum.
// at the apex of frustum, a sphere is attached and their radii are same.
template<typename Real>
struct default_probe
{
    Real angle;  // radian
    Real radius;
    mave::vector<Real, 3> apex;
};

// ----------------------------------------------------------------------------
// bounding box

template<typename Real>
struct aabb
{
    mave::vector<Real, 3> upper;
    mave::vector<Real, 3> lower;
};

template<typename Real>
aabb<Real> make_aabb(const sphere<Real>& sph) noexcept
{
    return aabb<Real>{
        sph.center + mave::vector<Real, 3>{sph.radius, sph.radius, sph.radius},
        sph.center - mave::vector<Real, 3>{sph.radius, sph.radius, sph.radius}
    };
}

template<typename Real>
aabb<Real> merge_aabb(const aabb<Real>& lhs, const aabb<Real>& rhs) noexcept
{
    return aabb<Real>{
        mave::max(lhs.upper, rhs.upper), mave::min(lhs.lower, rhs.lower)
    };
}

template<typename Real, typename Alloc>
aabb<Real> make_bounding_box(
        const std::vector<sphere<Real>, Alloc>& particles)
{
    if(particles.empty())
    {
        return aabb<Real>{mave::vector<Real, 3>{0.0, 0.0, 0.0},
                          mave::vector<Real, 3>{0.0, 0.0, 0.0}};
    }

    aabb<Real> bb = make_aabb(particles.front());
    for(const auto& p : particles)
    {
        bb = merge_aabb(bb, make_aabb(p));
    }
    return bb;
}

// ----------------------------------------------------------------------------
// "bounding" sphere.
// XXX Note that the center is fixed at the geometric center. It might decrease
// the runtime efficiency, but makes the calculation simpler. It is easy to
// calculate the geometric center from a configuration.

template<typename Real, typename Alloc>
sphere<Real> make_bounding_sphere_centered_at_geometric_center(
        const std::vector<sphere<Real>, Alloc>& particles)
{
    sphere<Real> sph;
    sph.center = mave::vector<Real, 3>(0,0,0);
    for(const auto& p : particles)
    {
        sph.center += p.center;
    }
    sph.center /= static_cast<Real>(particles.size());

    Real max_distance = 0.0;
    for(const auto& p : particles)
    {
        const auto dist = mave::length(p.center - sph.center) + p.radius;
        max_distance = std::max(max_distance, dist);
    }
    sph.radius = max_distance;
    return sph;
}

} // afmize
#endif // AFMIZE_SHAPES_HPP
