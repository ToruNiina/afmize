#ifndef AFMIZE_SIMULATOR_BASE_HPP
#define AFMIZE_SIMULATOR_BASE_HPP
#include "system.hpp"
#include "image.hpp"

namespace afmize
{

template<typename Real>
struct SimulatorBase
{
    virtual ~SimulatorBase() = default;
    virtual void run() = 0;
    virtual bool run(const std::size_t) = 0;
    virtual bool step() = 0;
    virtual std::size_t total_step() const noexcept = 0;
    virtual std::size_t current_step() const noexcept = 0;

    virtual system<Real> const& current_state() const noexcept = 0;
    virtual system<Real>&       current_state()       noexcept = 0;
    virtual image<Real> const&  current_image() const noexcept = 0;
};

} // afmize
#endif// AFMIZE_SIMULATOR_HPP
