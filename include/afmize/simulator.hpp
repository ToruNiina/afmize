#ifndef AFMIZE_SIMULATOR_HPP
#define AFMIZE_SIMULATOR_HPP
#include "system.hpp"
#include "stage.hpp"

namespace afmize
{
struct SimulatorBase
{
    virtual ~SimulatorBase() = default;
    virtual void run() = 0;
};

template<typename Real, typename Score>
struct SimulatedAnnealingSimulator : public SimulatorBase<Real>
{
    SimulatedAnnealingSimulator(stage<Real> stg, system<Real> sys,
                                rigid_observer<Real> obs, Score score)
        : stg_(std::move(stg)), sys_(std::move(sys)), obs_(std::move(obs)),
        score_(std::move(score))
    {}

    ~SimulatedAnnealingSimulator() override = default;

    void run() override
    {
        // TODO
        return ;
    }

  private:
    stage<Real>          stg_;
    system<Real>         sys_;
    rigid_observer<Real> obs_;
    Score              score_;
};

} // afmize
#endif// AFMIZE_SIMULATOR_HPP
