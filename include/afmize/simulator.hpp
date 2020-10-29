#ifndef AFMIZE_SIMULATOR_HPP
#define AFMIZE_SIMULATOR_HPP
#include "progress_bar.hpp"
#include "observe.hpp"
#include "shapes.hpp"
#include "system.hpp"
#include "stage.hpp"
#include "score.hpp"
#include <random>
#include <iostream>

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
    virtual stage<Real> const&  current_image() const noexcept = 0;
};

template<typename Real>
struct ScheduleBase
{
    virtual ~ScheduleBase() = default;
    virtual Real temperature(const std::size_t tstep) const noexcept = 0;
};

template<typename Real>
struct LinearSchedule : public ScheduleBase<Real>
{
    LinearSchedule(
            const Real init, const Real last, const std::size_t total_step)
        : init_(init), last_(last), total_step_(total_step)
    {}
    ~LinearSchedule() override = default;

    Real temperature(const std::size_t tstep) const noexcept override
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
struct ExponentialSchedule: public ScheduleBase<Real>
{
    ExponentialSchedule(
            const Real init, const Real last, const std::size_t total_step)
        : init_(init), last_(last), coef_(std::log(last / init)),
          total_step_(total_step)
    {}
    ~ExponentialSchedule() override = default;

    Real temperature(const std::size_t tstep) const noexcept override
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
template<typename Real, typename Mask>
struct SimulatedAnnealingSimulator : public SimulatorBase<Real>
{
    SimulatedAnnealingSimulator(const std::size_t total_step,
            const std::size_t   save,
            const std::uint32_t seed,
            const Real sigma_x,  const Real sigma_y,  const Real sigma_z,
            const Real max_rotx, const Real max_roty, const Real max_rotz,
            const Real max_dprobe_radius, const Real max_dprobe_angle,
            stage<Real>         ref,
            stage<Real>         stg,
            system<Real>        sys,
            std::unique_ptr<ObserverBase<Real>> obs,
            std::unique_ptr<ScoreBase<Real, Mask>> score,
            std::unique_ptr<ScheduleBase<Real>> schedule)
        : step_(0),
          total_step_(total_step),
          save_step_(save),
          sigma_dx_(sigma_x),
          sigma_dy_(sigma_y),
          sigma_dz_(sigma_z),
          max_drotx_(max_rotx),
          max_droty_(max_roty),
          max_drotz_(max_rotz),
          max_dprobe_radius_(max_dprobe_radius),
          max_dprobe_angle_(max_dprobe_angle),
          rng_(seed),
          nrm_(0.0, 1.0),
          stg_(stg),
          tmp_stg_(std::move(stg)),
          reference_(std::move(ref)),
          sys_(sys),
          next_(sys),
          obs_(std::move(obs)),
          score_(std::move(score)),
          schedule_(std::move(schedule)),
          bar_(total_step)
    {
        // FIXME
        // - split stage into image and resolution info
        // - move resolution information to system
        sys_.cells.initialize(stg_.x_resolution(), stg_.y_resolution(),
                              sys.particles);
        sys_.cells.construct(sys_.particles, sys_.bounding_box);
        obs_->observe(stg_, sys_);
        current_energy_ = score_->calc(reference_, stg_, Mask(stg_, sys_));
    }
    ~SimulatedAnnealingSimulator() override = default;

    void run() override
    {
        while(step_ < total_step_)
        {
            this->step();
        }
    }
    bool run(const std::size_t steps) override
    {
        const auto until = std::min(this->step_ + steps, this->total_step_);
        while(step_ < until)
        {
            this->step();
        }
        return step_ < total_step_;
    }

    bool step() override
    {
        if(this->step_ % save_step_ == 0)
        {
            const auto p = obs_->get_probe();
            std::cerr << bar_.format(this->step_);
            std::cout << "{ \"step\":"   << this->step_
                      << ", \"energy\":" << this->current_energy_
                      << ", \"probe radius [nm]\":" << p.radius * 0.1
                      << ", \"probe angle [deg]\":" << p.angle * 180.0 / 3.1416
                      << "}," << std::endl;
        }

        const auto temperature = schedule_->temperature(step_);
        const auto beta = 1.0 / temperature;

        // translation

        const auto dx = nrm_(rng_) * sigma_dx_;
        const auto dy = nrm_(rng_) * sigma_dy_;
//         const auto dz = nrm_(rng_) * sigma_dz_;
//         std::cerr << "dxy = " << dx << ", " << dy << std::endl;
        this->try_translation(mave::vector<Real, 3>(dx, dy, 0.0), beta);

        // rotation
        mave::matrix<Real, 3, 3> mat;

        // around x
        const auto rot_x = (this->generate_01() * 2.0 - 1.0) * max_drotx_;
        const auto cos_x = std::cos(rot_x);
        const auto sin_x = std::sin(rot_x);
//         std::cerr << "x rotation = " << rot_x << std::endl;
        mat.zero();
        mat(0, 0) = 1.0;
        mat(1, 1) =  cos_x;
        mat(1, 2) = -sin_x;
        mat(2, 1) =  sin_x;
        mat(2, 2) =  cos_x;
        this->try_rotation(mat, beta);

        // rot around y
        const auto rot_y = (this->generate_01() * 2.0 - 1.0) * max_droty_;
        const auto cos_y = std::cos(rot_y);
        const auto sin_y = std::sin(rot_y);
//         std::cerr << "y rotation = " << rot_y << std::endl;
        mat.zero();
        mat(0, 0) =  cos_y;
        mat(0, 2) =  sin_y;
        mat(1, 1) = 1.0;
        mat(2, 0) = -sin_y;
        mat(2, 2) =  cos_y;
        this->try_rotation(mat, beta);

        // rot around z
        const auto rot_z = (this->generate_01() * 2.0 - 1.0) * max_drotz_;
        const auto cos_z = std::cos(rot_z);
        const auto sin_z = std::sin(rot_z);
//         std::cerr << "z rotation = " << rot_z << std::endl;
        mat.zero();
        mat(0, 0) =  cos_z;
        mat(0, 1) = -sin_z;
        mat(1, 0) =  sin_z;
        mat(1, 1) =  cos_z;
        mat(2, 2) = 1.0;
        this->try_rotation(mat, beta);

        const auto drad = (this->generate_01() * 2.0 - 1.0) * max_dprobe_radius_;
        const auto dang = (this->generate_01() * 2.0 - 1.0) * max_dprobe_angle_;

        this->try_probe_change(drad, dang, beta);

        this->step_ += 1;
        return step_ < total_step_;
    }

    system<Real> const& current_state() const noexcept override {return sys_;}
    system<Real>&       current_state()       noexcept override {return sys_;}

    stage<Real> const& current_image() const noexcept override
    {
        return stg_;
    }

    std::unique_ptr<ObserverBase<Real>>& observer() noexcept {return obs_;}

    std::size_t total_step()   const noexcept override {return total_step_;}
    std::size_t current_step() const noexcept override {return step_;}

    Real& sigma_dx()        noexcept {return sigma_dx_;}
    Real  sigma_dx()  const noexcept {return sigma_dx_;}

    Real& sigma_dy()        noexcept {return sigma_dy_;}
    Real  sigma_dy()  const noexcept {return sigma_dy_;}

    Real& sigma_dz()        noexcept {return sigma_dz_;}
    Real  sigma_dz()  const noexcept {return sigma_dz_;}

    Real& max_drotx()       noexcept {return max_drotx_;}
    Real  max_drotx() const noexcept {return max_drotx_;}

    Real& max_droty()       noexcept {return max_droty_;}
    Real  max_droty() const noexcept {return max_droty_;}

    Real& max_drotz()       noexcept {return max_drotz_;}
    Real  max_drotz() const noexcept {return max_drotz_;}

  private:

    void try_translation(mave::vector<Real, 3> dr, const Real beta)
    {
        next_ = sys_;

        // translation does not change the shape of bounding box.
        next_.bounding_box.upper += dr;
        next_.bounding_box.lower += dr;

        // avoid particle from escaping observation stage.
        // To fix the number of pixels when calculating score, it does not
        // allow particle sticks out of the stage region.
        if(next_.bounding_box.lower[0] < tmp_stg_.x_range().first  ||
           next_.bounding_box.lower[1] < tmp_stg_.y_range().first  ||
           tmp_stg_.x_range().second < next_.bounding_box.upper[0] ||
           tmp_stg_.y_range().second < next_.bounding_box.upper[1])
        {
            return ;
        }

//         if(next_.bounding_box.lower[2] < 0.0)
        {
            const auto offset = next_.bounding_box.lower[2];
            dr[2]                       -= offset;
            next_.bounding_box.upper[2] -= offset;
            next_.bounding_box.lower[2] -= offset;
        }

        // apply the movement
        for(auto& p : this->next_.particles)
        {
            p.center += dr;
        }

        // update cell list
        next_.cells.construct(next_.particles, next_.bounding_box);

        // generate pseudo AFM image
        // bottom == 0 because the bottom object is aligned to 0.0.
        obs_->observe(tmp_stg_, next_);

        // calculate score depending on score function
        const auto energy = score_->calc(reference_, tmp_stg_, Mask(tmp_stg_, sys_));
        const auto deltaE = energy - current_energy_;

//         std::cerr << "beta = " << beta << ", dE = " << deltaE
//                   << ", Enext = " << energy << ", Ecurr = " << current_energy_
//                   << ", prob = " << std::exp(-deltaE * beta) << std::endl;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            sys_ = next_;
            current_energy_ = energy;
            stg_ = tmp_stg_;
        }
        return;
    }

    void try_rotation(const mave::matrix<Real, 3, 3>& rot, const Real beta)
    {
        next_ = sys_;

        // apply the movement.
        // 1. move center to origin
        // 2. apply rotation matrix
        // 3. move back to the original center

        mave::vector<Real, 3> com(0.0, 0.0, 0.0);
        for(const auto& p : this->next_.particles)
        {
            com += p.center;
        }
        com /= static_cast<Real>(next_.particles.size());

        mave::matrix<Real, 4, 4> translation1;
        translation1.zero();
        translation1(0, 0) = 1.0;
        translation1(1, 1) = 1.0;
        translation1(2, 2) = 1.0;
        translation1(3, 3) = 1.0;
        translation1(0, 3) = -com[0];
        translation1(1, 3) = -com[1];
        translation1(2, 3) = -com[2];

        mave::matrix<Real, 4, 4> rotation;
        rotation.zero();
        rotation(0, 0) = rot(0, 0);
        rotation(0, 1) = rot(0, 1);
        rotation(0, 2) = rot(0, 2);
        rotation(1, 0) = rot(1, 0);
        rotation(1, 1) = rot(1, 1);
        rotation(1, 2) = rot(1, 2);
        rotation(2, 0) = rot(2, 0);
        rotation(2, 1) = rot(2, 1);
        rotation(2, 2) = rot(2, 2);
        rotation(3, 3) = 1.0;

        mave::matrix<Real, 4, 4> translation2;
        translation2.zero();
        translation2(0, 0) = 1.0;
        translation2(1, 1) = 1.0;
        translation2(2, 2) = 1.0;
        translation2(3, 3) = 1.0;
        translation2(0, 3) = com[0];
        translation2(1, 3) = com[1];
        translation2(2, 3) = com[2];

        const mave::matrix<Real, 4, 4> matrix = translation2 * rotation * translation1;
        for(auto& p : this->next_.particles)
        {
            mave::vector<Real, 4> r(p.center[0], p.center[1], p.center[2], 1.0);
            r = matrix * r;
            p.center[0] = r[0];
            p.center[1] = r[1];
            p.center[2] = r[2];
        }

        next_.bounding_box = make_bounding_box(next_.particles);

        // avoid particle from escaping observation stage
        if(next_.bounding_box.lower[0] < tmp_stg_.x_range().first  ||
           next_.bounding_box.lower[1] < tmp_stg_.y_range().first  ||
           tmp_stg_.x_range().second < next_.bounding_box.upper[0] ||
           tmp_stg_.y_range().second < next_.bounding_box.upper[1])
        {
            return ;
        }

//         if(next_.bounding_box.lower[2] < 0.0)
        {
            // align the bottom to the xy plane (z=0.0)
            for(auto& p : this->next_.particles)
            {
                p.center[2] -= next_.bounding_box.lower[2];
            }
            next_.bounding_box.upper[2] -= next_.bounding_box.lower[2];
            next_.bounding_box.lower[2]  = 0.0;
        }

        // update cell list
        next_.cells.construct(next_.particles, next_.bounding_box);

        // generate pseudo AFM image
        // bottom == 0 because the bottom object is aligned to 0.0.
        obs_->observe(tmp_stg_, next_);

        // calculate score depending on score function
        const auto energy = score_->calc(reference_, tmp_stg_, Mask(tmp_stg_, sys_));
        const auto deltaE = energy - current_energy_;

//         std::cerr << "beta = " << beta << ", dE = " << deltaE
//                   << ", Enext = " << energy << ", Ecurr = " << current_energy_
//                   << ", prob = " << std::exp(-deltaE * beta) << std::endl;

//         std::cerr << "trial energy = " << energy << std::endl;
//         std::cerr << "trial deltaE = " << deltaE << std::endl;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            sys_ = next_;
            current_energy_ = energy;
            stg_ = tmp_stg_;
        }
        return;
    }

    void try_probe_change(const Real drad, const Real dang, const Real beta)
    {
        const auto current_probe = obs_->get_probe();
        auto next_probe = current_probe;

        next_probe.radius += drad;
        if(next_probe.radius < 0)
        {
            next_probe.radius = current_probe.radius;
        }
        next_probe.angle  += dang;
        if(next_probe.angle < 0)
        {
            next_probe.angle = 0;
        }
        else if(3.14159265 * 0.5 <= next_probe.angle)
        {
            next_probe.angle = 3.14159265 * 0.5;
        }
        obs_->update_probe(next_probe);

        // generate pseudo AFM image
        // bottom == 0 because the bottom object is aligned to 0.0.
        obs_->observe(tmp_stg_, sys_);

        // calculate score depending on score function
        const auto energy = score_->calc(reference_, tmp_stg_, Mask(tmp_stg_, sys_));
        const auto deltaE = energy - current_energy_;

//         std::cerr << "beta = " << beta << ", dE = " << deltaE
//                   << ", Enext = " << energy << ", Ecurr = " << current_energy_
//                   << ", prob = " << std::exp(-deltaE * beta) << std::endl;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            // system is not updated
            current_energy_ = energy;
            stg_ = tmp_stg_;
        }
        else
        {
            // if energy increases, go back to the original probe
            obs_->update_probe(current_probe);
        }
        return;
    }


    Real generate_01() noexcept
    {
        return std::generate_canonical<Real, std::numeric_limits<Real>::digits>(rng_);
    }

  private:
    std::size_t         step_;
    std::size_t    save_step_;
    std::size_t   total_step_;
    Real      current_energy_;

    Real sigma_dx_;
    Real sigma_dy_;
    Real sigma_dz_;
    Real max_drotx_;
    Real max_droty_;
    Real max_drotz_;

    Real max_dprobe_radius_;
    Real max_dprobe_angle_;

    std::mt19937                             rng_;
    std::normal_distribution<Real>           nrm_;
    stage<Real>                              stg_;
    stage<Real>                          tmp_stg_;
    stage<Real>                        reference_;
    system<Real>                             sys_;
    system<Real>                            next_;
    std::unique_ptr<ObserverBase<Real>>      obs_;
    std::unique_ptr<ScoreBase<Real, Mask>> score_;
    std::unique_ptr<ScheduleBase<Real>> schedule_;
    afmize::progress_bar<70>                 bar_;
};

} // afmize
#endif// AFMIZE_SIMULATOR_HPP
