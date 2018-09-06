#ifndef AFMIZE_SHAPES_HPP
#define AFMIZE_SHAPES_HPP
#include <extlib/mave/mave/mave.hpp>

namespace afmize
{

template<typename Real>
struct sphere
{
    Real radius;
    mave::vector<Real, 3> center;
};

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
    Real angle;
    Real radius;
    mave::vector<Real, 3> apex;
};

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


} // afmize
#endif // AFMIZE_SHAPES_HPP
