#ifndef AFMIZE_SHAPES_HPP
#define AFMIZE_SHAPES_HPP
#include <extlib/mave/mave/mave.hpp>

namespace afmize
{

template<typename realT>
struct sphere
{
    realT radius;
    mave::vector<realT, 3> center;
};

template<typename realT>
struct circular_frustum
{
    realT angle;
    realT radius;
    mave::vector<realT, 3> apex;
};

// default probe shape is sphere + frustum.
// at the apex of frustum, a sphere is attached and their radii are same.
template<typename realT>
struct default_probe
{
    realT angle;
    realT radius;
    mave::vector<realT, 3> apex;
};

} // afmize
#endif // AFMIZE_SHAPES_HPP
