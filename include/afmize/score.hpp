#ifndef AFMIZE_SCORE_HPP
#define AFMIZE_SCORE_HPP
#include "stage.hpp"
#include "mask.hpp"

namespace afmize
{

template<typename Real>
struct score_negative_cosine_similarity_t
{
    Real k; // modulating coefficient

    template<typename Mask>
    Real operator()(const stage<Real>& lhs, const stage<Real>& rhs, const Mask& mask) const
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

                const auto l = lhs.at(x, y);
                const auto r = rhs.at(x, y);

                numer  += l * r;
                denom1 += l * l;
                denom2 += r * r;
            }
        }
        return k * (1.0 - numer / (denom1 * denom2));
    }
};

template<typename Real>
struct score_root_mean_squared_deviation_t
{
    Real k; // modulating coefficient

    template<typename Mask>
    Real operator()(const stage<Real>& lhs, const stage<Real>& rhs, const Mask& mask) const
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
        return k * std::sqrt(sd / static_cast<Real>(N));
    }
};

} // afmize
#endif// AFMIZE_MASK_HPP
