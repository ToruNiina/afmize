#ifndef AFMIZE_NOISE_HPP
#define AFMIZE_NOISE_HPP
#include "image.hpp"
#include <random>

namespace afmize
{

template<typename Real, typename RNG>
void apply_noise(image<Real>& img, RNG& rng, const Real sigma = 3.0 /*angstrom*/,
                 const Real stage_position = 0.0)
{
    std::normal_distribution<Real> nrm(0.0, sigma);
    for(auto& pxl : img)
    {
        pxl = std::max(stage_position, pxl + nrm(rng));
    }
    return ;
}

} // afmize
#endif// AFMIZE_NOISE_HPP
