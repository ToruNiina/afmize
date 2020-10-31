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
    virtual Real calc(const system<Real>&, const std::unique_ptr<ObserverBase<Real>>&,
                      const image<Real>&, const Mask&) const = 0;
};

template<typename Real, typename Mask>
struct NegativeCosineSimilarity: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit NegativeCosineSimilarity(const Real k_): k(k_) {}
    ~NegativeCosineSimilarity() override = default;

    Real calc(const system<Real>& sys, const std::unique_ptr<ObserverBase<Real>>& obs,
              const image<Real>& target, const Mask& mask) const override
    {
        Real numer  = 0;
        Real denom1 = 0;
        Real denom2 = 0;

        const auto& img = obs->observe(sys);

        for(std::size_t y=mask.lower_bounding_y(); y <= mask.upper_bounding_y(); ++y)
        {
            for(std::size_t x=mask.lower_bounding_x(); x <= mask.upper_bounding_x(); ++x)
            {
                // in a bounding rectangle, still there could be non-covered pixels.
                if(mask(x, y)) {continue;}

                const auto l = img.at(x, y);
                const auto r = target.at(x, y);

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

    Real calc(const system<Real>& sys, const std::unique_ptr<ObserverBase<Real>>& obs,
              const image<Real>& target, const Mask& mask) const override
    {
        std::uint64_t N = 0;
        Real sd = 0.0;

        const auto& img = obs->observe(sys);
        for(std::size_t y=mask.lower_bounding_y(); y <= mask.upper_bounding_y(); ++y)
        {
            for(std::size_t x=mask.lower_bounding_x(); x <= mask.upper_bounding_x(); ++x)
            {
                // in a bounding rectangle, still there could be non-covered pixels.
                if(mask(x, y)) {continue;}

                const auto l = img(x, y);
                const auto r = target(x, y);

                sd += (l - r) * (l - r);
                N  += 1;
            }
        }
        return k * std::sqrt(sd / static_cast<Real>(N));
    }
};

template<typename Real, typename Mask>
struct SumOfDifference: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit SumOfDifference(const Real k_): k(k_) {}
    ~SumOfDifference() override = default;

    Real calc(const system<Real>& sys, const std::unique_ptr<ObserverBase<Real>>& obs,
              const image<Real>& target, const Mask& mask) const override
    {
        Real s = 0.0;

        const auto& img = obs->observe(sys);
        for(std::size_t y=mask.lower_bounding_y(); y <= mask.upper_bounding_y(); ++y)
        {
            for(std::size_t x=mask.lower_bounding_x(); x <= mask.upper_bounding_x(); ++x)
            {
                // in a bounding rectangle, still there could be non-covered pixels.
                if(mask(x, y)) {continue;}

                const auto l = img(x, y);
                const auto r = target(x, y);

                s += std::abs(l - r);
            }
        }
        return k * s;
    }
};

} // afmize
#endif// AFMIZE_MASK_HPP
