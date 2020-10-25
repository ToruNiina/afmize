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
    virtual void run(const std::size_t) = 0;
};

template<typename Real>
struct LinearSchedule
{
    LinearSchedule(
            const Real init, const Real last, const std::size_t total_step)
        : init_(init), last_(last), total_step_(total_step)
    {}

    Real operator()(const std::size_t tstep) const noexcept
    {
        const auto t = static_cast<Real>(tstep) / static_cast<Real>(total_step_);
        return t * (last_ - init_) + init_;
    }

  private:
    Real init_;
    Real last_;
    std::size_t total_step_;
};

template<typename Real>
struct ExponentialSchedule
{
    ExponentialSchedule(
            const Real init, const Real last, const std::size_t total_step)
        : init_(init), last_(last), coef_(std::log(last / init)),
          total_step_(total_step)
    {}

    Real operator()(const std::size_t tstep) const noexcept
    {
        const auto t = static_cast<Real>(tstep) / static_cast<Real>(total_step_);
        return init_ * std::exp(coef_ * t);
    }

  private:
    Real init_;
    Real last_;
    Real coef_;
    std::size_t total_step_;
};

// Brownian Movement under Simulated Annealing with given score function and mask.
template<typename Real, typename Score, typename Mask, typename Schedule>
struct SimulatedAnnealingSimulator : public SimulatorBase<Real>
{
    SimulatedAnnealingSimulator(const std::size_t total_step,
            const std::uint32_t seed, stage<Real> ref,
            stage<Real> stg, system<Real> sys, rigid_observer<Real> obs,
            Score score, Mask mask, Schedule schedule)
        : step_(0), total_step_(total_step),
          reference_(std::move(ref)),
          rng_(seed),
          nrm_(0.0, 1.0),
          stg_(std::move(stg)),
          sys_(std::move(sys)),
          obs_(std::move(obs)),
          score_(std::move(score)),
          schedule_(std::move(schedule))
    {
        current_enregy_ = score_(reference_, stage_, Mask(stg_, sys_));
    }
    ~SimulatedAnnealingSimulator() override = default;

    void run() override
    {
        while(step_ < total_step_)
        {
            this->step();
        }
    }
    void run(const std::size_t steps) override
    {
        const auto until = std::min(this->step_ + steps, this->total_step_);
        while(step_ < until)
        {
            this->step();
        }
        return ;
    }

    void step()
    {
        const auto temperature = schedule_(step_);
        const auto beta = 1.0 / temperature;

        // translation

        const auto dx = nrm_(rng_) * sigma_dx_;
        const auto dy = nrm_(rng_) * sigma_dy_;
        this->try_translation(mave::vector<Real, 3>(dx, dy, 0.0), beta);

        // rotation
        mave::matrix<Real, 3, 3> mat;

        // around x
        const auto rot_x = (this->generate_01() * 2.0 - 1.0) * max_drotx_;
        mat.zero();
        mat(0, 0) = 1.0;
        mat(1, 1) =  std::cos(rot_x);
        mat(1, 2) = -std::sin(rot_x);
        mat(2, 1) =  std::sin(rot_x);
        mat(2, 2) =  std::cos(rot_x);
        this->try_rotation(mat, beta);

        // rot around y
        const auto rot_y = (this->generate_01() * 2.0 - 1.0) * max_droty_;
        mat.zero();
        mat(0, 0) =  std::cos(rot_x);
        mat(0, 2) =  std::sin(rot_x);
        mat(1, 1) = 1.0;
        mat(2, 0) = -std::sin(rot_x);
        mat(2, 2) =  std::cos(rot_x);
        this->try_rotation(mat, beta);

        // rot around z
        const auto rot_z = (this->generate_01() * 2.0 - 1.0) * max_drotz_;
        mat.zero();
        mat(0, 0) =  std::cos(rot_x);
        mat(0, 1) = -std::sin(rot_x);
        mat(1, 0) =  std::sin(rot_x);
        mat(1, 1) =  std::cos(rot_x);
        mat(2, 2) = 1.0;
        this->try_rotation(mat, beta);

        this->step_ += 1;
        return ;
    }

    rigid_observer<Real>& observer() noexcept {return obs_;}

    std::size_t current_t() const noexcept {return step_;}

    Real& sigma_dx()        noexcept {return sigma_dx_;}
    Real  sigma_dx()  const noexcept {return sigma_dx_;}

    Real& sigma_dy()        noexcept {return sigma_dy_;}
    Real  sigma_dy()  const noexcept {return sigma_dy_;}

    Real& max_drotx()       noexcept {return max_drotx_;}
    Real  max_drotx() const noexcept {return max_drotx_;}

    Real& max_droty()       noexcept {return max_droty_;}
    Real  max_droty() const noexcept {return max_droty_;}

    Real& max_drotz()       noexcept {return max_drotz_;}
    Real  max_drotz() const noexcept {return max_drotz_;}

  private:

    void try_translation(const mave::vector<Real, 3>& dr, const Real beta)
    {
        next_ = sys_;

        // translation does not change the shape of bounding box.
        next_.bounding_box.upper += dr;
        next_.bounding_box.lower += dr;

        // avoid particle from escaping observation stage
        if(next_.bounding_box.upper[0] < stage.lower()[0] ||
           next_.bounding_box.upper[1] < stage.lower()[1] ||
           stage.upper()[0] < next_.bounding_box.lower[0] ||
           stage.upper()[1] < next_.bounding_box.lower[1])
        {
            return ;
        }

        // apply the movement
        for(auto& p : this->next_.particles)
        {
            p.center += dr;
        }

        // generate pseudo AFM image
        // bottom == 0 because the bottom object is aligned to 0.0.
        obs_(stg_, next_, 0.0);

        // calculate score depending on score function
        const auto energy = score_(reference_, stg_, Mask(stg_, sys_));
        const auto deltaE = energy - current_energy_;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            sys_ = next_;
            current_energy_ = energy;;
        }
        return;
    }

    void try_rotation(const mave::matrix<Real, 3, 3>& dr, const Real beta)
    {
        next_ = sys_;

        // apply the movement
        for(auto& p : this->next_.particles)
        {
            p.center += dr;
        }

        next_.bounding_box = make_bounding_box(next_.particles);

        // avoid particle from escaping observation stage
        if(next_.bounding_box.upper[0] < stage.lower()[0] ||
           next_.bounding_box.upper[1] < stage.lower()[1] ||
           stage.upper()[0] < next_.bounding_box.lower[0] ||
           stage.upper()[1] < next_.bounding_box.lower[1])
        {
            return ;
        }

        // align the bottom to the xy plane (z=0.0)
        for(auto& p : this->next_.particles)
        {
            p.center[2] -= next_.bounding_box.lower[2];
        }
        next_.bounding_box.upper[2] -= next_.bounding_box.lower[2];
        next_.bounding_box.lower[2]  = 0.0;

        // generate pseudo AFM image
        // bottom == 0 because the bottom object is aligned to 0.0.
        obs_(stg_, next_, 0.0);

        // calculate score depending on score function
        const auto energy = score_(reference_, stg_, Mask(stg_, sys_));
        const auto deltaE = energy - current_energy_;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            sys_ = next_;
            current_energy_ = energy;;
        }
        return;
    }

    Real generate_01() noexcept
    {
        return std::generate_canonical<Real, std::numeric_limits<Real>::digits>(rng_);
    }

  private:
    std::size_t         step_;
    std::size_t   total_step_;
    Real      current_energy_;

    Real sigma_dx_;
    Real sigma_dy_;
    Real max_drotx_;
    Real max_droty_;
    Real max_drotz_;

    std::mt19937                   rng_;
    std::normal_distribution<Real> nrm_;
    stage<Real>                    stg_;
    stage<Real>              reference_;
    system<Real>                   sys_;
    system<Real>                  next_;
    rigid_observer<Real>           obs_;
    Score                        score_;
    Schedule                  schedule_;
};

} // afmize
#endif// AFMIZE_SIMULATOR_HPP
