#ifndef AFMIZE_SCORE_HPP
#define AFMIZE_SCORE_HPP
#include "image.hpp"
#include "mask.hpp"

namespace afmize
{

template<typename Real, typename Mask>
struct ScoreBase
{
    virtual ~ScoreBase() = default;
    virtual Real calc(const image<Real>&, const image<Real>&, const Mask&) const = 0;
};

template<typename Real, typename Mask>
struct NegativeCosineSimilarity: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit NegativeCosineSimilarity(const Real k_): k(k_) {}
    ~NegativeCosineSimilarity() override = default;

    Real calc(const image<Real>& lhs, const image<Real>& rhs, const Mask& mask) const override
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
        return k * (1.0 - numer / std::sqrt(denom1 * denom2));
    }
};

template<typename Real, typename Mask>
struct RootMeanSquareDeviation: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit RootMeanSquareDeviation(const Real k_): k(k_) {}
    ~RootMeanSquareDeviation() override = default;

    Real calc(const image<Real>& lhs, const image<Real>& rhs, const Mask& mask) const override
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
