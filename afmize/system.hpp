#ifndef AFMIZE_SYSTEM_HPP
#define AFMIZE_SYSTEM_HPP
#include <afmize/shapes.hpp>

namespace afmize
{

template<typename Real>
struct system
{
    system(std::vector<sphere<Real>> mol)
        : particles(std::move(mol))
    {
        this->shape = make_aabb(particles.front());
        for(std::size_t i=1; i<particles.size(); ++i)
        {
            this->shape = merge_aabb(this->shape, particles[i]);
        }
    }
    ~system() = default;

    aabb<Real>                shape;
    std::vector<sphere<Real>> particles;
};




} // afmize
#endif // AFMIZE_SYSTEM_HPP
