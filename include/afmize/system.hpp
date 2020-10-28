#ifndef AFMIZE_SYSTEM_HPP
#define AFMIZE_SYSTEM_HPP
#include <afmize/shapes.hpp>
#include <vector>

namespace afmize
{

template<typename Real>
struct system
{
    explicit system(std::pair<std::vector<std::string>, std::vector<sphere<Real>>> mol)
        : max_radius(0), particles(std::move(mol.second)), names(std::move(mol.first))
    {
        max_radius = particles.front().radius;
        this->bounding_box = make_aabb(particles.front());
        for(std::size_t i=1; i<particles.size(); ++i)
        {
            max_radius = std::max(max_radius, particles[i].radius);
            this->bounding_box =
                merge_aabb(this->bounding_box, make_aabb(particles[i]));
        }
    }
    ~system() = default;

    Real                      max_radius;
    aabb<Real>                bounding_box;
    std::vector<sphere<Real>> particles;
    std::vector<std::string>  names;
};

} // afmize
#endif // AFMIZE_SYSTEM_HPP
