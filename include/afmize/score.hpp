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
    virtual Real calc(const system<Real>& sys,
              const image<Real>& img,    const Mask& mask,
              const image<Real>& target, const Mask& target_mask) const = 0;
};

template<typename Real, typename Mask, bool UseZeroPixel>
struct NegativeCosineSimilarity: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit NegativeCosineSimilarity(const Real k_): k(k_) {}
    ~NegativeCosineSimilarity() override = default;

    Real calc(const system<Real>& sys,
              const image<Real>& img,    const Mask& mask,
              const image<Real>& target, const Mask& target_mask) const override
    {
        assert(mask.size() == target_mask.size());

        Real numer  = 0;
        Real denom1 = 0;
        Real denom2 = 0;

        for(std::size_t y=0; y < mask.pixel_y(); ++y)
        {
            for(std::size_t x=0; x < mask.pixel_x(); ++x)
            {
                if( ! UseZeroPixel && mask(img, x, y) == 0)
                {
                    continue;
                }
                const auto l =        mask(img,    x, y);
                const auto r = target_mask(target, x, y);

                numer  += l * r;
                denom1 += l * l;
                denom2 += r * r;
            }
        }
        if(denom1 == 0 || denom2 == 0)
        {
            return std::numeric_limits<Real>::infinity();
        }
        return k * (1.0 - numer / std::sqrt(denom1 * denom2));
    }
};

template<typename Real, typename Mask, bool UseZeroPixel>
struct NegativeCorrelation: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit NegativeCorrelation(const Real k_): k(k_) {}
    ~NegativeCorrelation() override = default;

    Real calc(const system<Real>& sys,
              const image<Real>& img,    const Mask& mask,
              const image<Real>& target, const Mask& target_mask) const override
    {
        assert(mask.size() == target_mask.size());

        std::size_t N = 0;
        Real mean1 = 0;
        Real mean2 = 0;
        for(std::size_t y=0; y < mask.pixel_y(); ++y)
        {
            for(std::size_t x=0; x < mask.pixel_x(); ++x)
            {
                if( ! UseZeroPixel && mask(img, x, y) == 0)
                {
                    continue;
                }
                const auto l =        mask(img,    x, y);
                const auto r = target_mask(target, x, y);

                mean1 += l;
                mean2 += r;
                N     += 1;
            }
        }
        if(N == 0 || mean1 == 0.0 || mean2 == 0.0)
        {
            return std::numeric_limits<Real>::infinity();
        }
        mean1 /= static_cast<Real>(N);
        mean2 /= static_cast<Real>(N);

        Real numer  = 0;
        Real denom1 = 0;
        Real denom2 = 0;

        for(std::size_t y=0; y < mask.pixel_y(); ++y)
        {
            for(std::size_t x=0; x < mask.pixel_x(); ++x)
            {
                if( ! UseZeroPixel && mask(img, x, y) == 0)
                {
                    continue;
                }
                const auto l =        mask(img,    x, y);
                const auto r = target_mask(target, x, y);

                numer  += (l - mean1) * (r - mean2);
                denom1 += (l - mean1) * (l - mean1);
                denom2 += (r - mean2) * (r - mean2);
            }
        }
        return k * (1.0 - numer / std::sqrt(denom1 * denom2));
    }
};


template<typename Real, typename Mask, bool UseZeroPixel>
struct RootMeanSquareDeviation: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit RootMeanSquareDeviation(const Real k_): k(k_) {}
    ~RootMeanSquareDeviation() override = default;

    Real calc(const system<Real>& sys,
              const image<Real>& img,    const Mask& mask,
              const image<Real>& target, const Mask& target_mask) const override
    {
        assert(mask.size() == target_mask.size());

        std::uint64_t N = 0;
        Real sd = 0.0;

        for(std::size_t y=0; y < mask.pixel_y(); ++y)
        {
            for(std::size_t x=0; x < mask.pixel_x(); ++x)
            {
                if( ! UseZeroPixel && mask(img, x, y) == 0)
                {
                    continue;
                }
                const auto l =        mask(img,    x, y);
                const auto r = target_mask(target, x, y);

                sd += (l - r) * (l - r);
                N  += 1;
            }
        }
        if(N == 0)
        {
            return std::numeric_limits<Real>::infinity();
        }
        return k * std::sqrt(sd / static_cast<Real>(N));
    }
};

template<typename Real, typename Mask, bool UseZeroPixel>
struct SumOfDifference: public ScoreBase<Real, Mask>
{
    Real k; // modulating coefficient

    explicit SumOfDifference(const Real k_): k(k_) {}
    ~SumOfDifference() override = default;

    Real calc(const system<Real>& sys,
              const image<Real>& img,    const Mask& mask,
              const image<Real>& target, const Mask& target_mask) const override
    {
        assert(mask.size() == target_mask.size());

        Real s = 0.0;

        for(std::size_t y=0; y < mask.pixel_y(); ++y)
        {
            for(std::size_t x=0; x < mask.pixel_x(); ++x)
            {
                if( ! UseZeroPixel && mask(img, x, y) == 0)
                {
                    continue;
                }

                const auto l =        mask(img,    x, y);
                const auto r = target_mask(target, x, y);

                s += std::abs(l - r);
            }
        }
        return k * s;
    }
};

template<typename Real, typename Mask>
struct TopographicalPenalty: public ScoreBase<Real, Mask>
{
    Real penalty;   // penalty per particle in the forbidden region
    Real reward;    // reward per particle in the favorable region
    Real thickness; // thickness of the favorable region

    TopographicalPenalty(const Real p, const Real r, const Real th)
        : penalty(p), reward(r), thickness(th)
    {}
    ~TopographicalPenalty() override = default;

    Real calc(const system<Real>& sys,
              const image<Real>& img,    const Mask& mask,
              const image<Real>& target, const Mask& target_mask) const override
    {
        assert(mask.size() == target_mask.size());

        Real score = 0.0;
        for(std::size_t y=0; y < mask.pixel_y(); ++y)
        {
            for(std::size_t x=0; x < mask.pixel_x(); ++x)
            {
                const auto threshold_penalty = target_mask(target, x, y);
                const auto threshold_reward  = target_mask(target, x, y) - thickness;

                for(const auto& p : sys.particles) // speedup
                {
                    if(threshold_penalty < p.center[2] + p.radius)
                    {
                        // higher energy
                        score += penalty;
                    }
                    else if(threshold_reward < p.center[2] + p.radius)
                    {
                        // lower energy
                        score -= reward;
                    }
                }
            }
        }
        return score;
    }
};

template<typename Real, typename Mask>
struct PixelPenalty: public ScoreBase<Real, Mask>
{
    Real penalty;   // penalty per pixel in the forbidden region
    Real reward;    // reward per pixel in the favorable region
    Real thickness; // thickness of the favorable region

    PixelPenalty(const Real p, const Real r, const Real th)
        : penalty(p), reward(r), thickness(th)
    {}
    ~PixelPenalty() override = default;

    Real calc(const system<Real>& sys,
              const image<Real>& img,    const Mask& mask,
              const image<Real>& target, const Mask& target_mask) const override
    {
        assert(mask.size() == target_mask.size());

        Real score = 0.0;
        for(std::size_t y=0; y < mask.pixel_y(); ++y)
        {
            for(std::size_t x=0; x < mask.pixel_x(); ++x)
            {
                const auto threshold_penalty = target_mask(target, x, y);
                const auto threshold_reward  = target_mask(target, x, y) - thickness;

                if(threshold_penalty < mask(img, x, y))
                {
                    // higher energy
                    score += penalty;
                }
                else if(threshold_reward < mask(img, x, y))
                {
                    // lower energy
                    score -= reward;
                }
            }
        }
        return score;
    }
};

} // afmize
#endif// AFMIZE_MASK_HPP
