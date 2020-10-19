#ifndef AFMIZE_SYSTEM_HPP
#define AFMIZE_SYSTEM_HPP
#include <afmize/shapes.hpp>
#include <vector>

namespace afmize
{

template<typename Real>
struct system
{
    system(std::vector<sphere<Real>> mol)
        : particles(std::move(mol))
    {
        this->bounding_box = make_aabb(particles.front());
        for(std::size_t i=1; i<particles.size(); ++i)
        {
            this->bounding_box =
                merge_aabb(this->bounding_box, make_aabb(particles[i]));
        }
    }
    ~system() = default;

    aabb<Real>                bounding_box;
    std::vector<sphere<Real>> particles;
};




} // afmize
#endif // AFMIZE_SYSTEM_HPP
