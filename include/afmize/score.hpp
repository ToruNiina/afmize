#ifndef AFMIZE_SCORE_HPP
#define AFMIZE_SCORE_HPP
#include "stage.hpp"
#include "mask.hpp"

namespace afmize
{

template<typename Real, typename Mask>
Real score_cosine_similarity(
        const stage<Real>& lhs, const stage<Real>& rhs, const Mask& mask)
{
    Real numer  = 0;
    Real denom1 = 0;
    Real denom2 = 0;

    for(std::size_t y=mask.lower_bounding_y(); y <= mask.upper_bounding_y(); ++y)
    {
        for(std::size_t x=mask.lower_bounding_x(); x <= mask.upper_bounding_x(); ++x)
        {
            // in a bounding rectangle, still there could be non-covered pixels.
            if(mask(x, y)) {continue;}

            const auto l = lhs(x, y);
            const auto r = rhs(x, y);

            numer  += l * r;
            denom1 += l * l;
            denom2 += r * r;
        }
    }
    return numer / (denom1 * denom2);
}

template<typename Real, typename Mask>
Real score_root_mean_squared_deviation(
        const stage<Real>& lhs, const stage<Real>& rhs, const Mask& mask)
{
    std::uint64_t N = 0;
    Real sd = 0.0;

    for(std::size_t y=mask.lower_bounding_y(); y <= mask.upper_bounding_y(); ++y)
    {
        for(std::size_t x=mask.lower_bounding_x(); x <= mask.upper_bounding_x(); ++x)
        {
            // in a bounding rectangle, still there could be non-covered pixels.
            if(mask(x, y)) {continue;}

            const auto l = lhs(x, y);
            const auto r = rhs(x, y);

            sd += (l - r) * (l - r);
            N  += 1;
        }
    }
    return std::sqrt(sd / static_cast<Real>(N));
}

} // afmize
#endif// AFMIZE_MASK_HPP
