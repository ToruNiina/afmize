#ifndef AFMIZE_ANNEALING_SIMULATOR_HPP
#define AFMIZE_ANNEALING_SIMULATOR_HPP
#include "simulator_base.hpp"
#include "progress_bar.hpp"
#include "output_utility.hpp"
#include "observe.hpp"
#include "shapes.hpp"
#include "system.hpp"
#include "stage.hpp"
#include "score.hpp"
#include <map>
#include <random>
#include <iostream>

namespace afmize
{

template<typename Real>
struct ScheduleBase
{
    virtual ~ScheduleBase() = default;
    virtual Real temperature(const std::size_t tstep) const noexcept = 0;
    virtual Real init_temperature() const noexcept = 0;
    virtual Real last_temperature() const noexcept = 0;
    virtual void set_init_temperature(const Real) = 0;
    virtual void set_last_temperature(const Real) = 0;
};

template<typename Real>
struct LinearSchedule final: public ScheduleBase<Real>
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

    Real init_temperature() const noexcept override {return init_;}
    Real last_temperature() const noexcept override {return last_;}

    void set_init_temperature(const Real init) override
    {
        this->init_ = init;
    }
    void set_last_temperature(const Real last) override
    {
        this->last_ = last;
    }

  private:
    Real init_;
    Real last_;
    std::size_t total_step_;
};

template<typename Real>
struct ExponentialSchedule final: public ScheduleBase<Real>
{
    ExponentialSchedule(
            const Real init, const Real last, const std::size_t total_step)
        : init_(init), last_(last), coef_(std::log(last / init)),
          total_step_(total_step)
    {
        if(last_ <= 0.0)
        {
            last_ = 1e-12;
        }
    }
    ~ExponentialSchedule() override = default;

    Real temperature(const std::size_t tstep) const noexcept override
    {
        const auto t = static_cast<Real>(tstep) / static_cast<Real>(total_step_);
        return init_ * std::exp(coef_ * t);
    }

    Real init_temperature() const noexcept override {return init_;}
    Real last_temperature() const noexcept override {return last_;}

    void set_init_temperature(const Real init) override
    {
        this->init_ = init;
        if(this->init_ <= 0)
        {
            this->init_ = 1e-12;
        }
        coef_ = std::log(last_ / init_);
    }
    void set_last_temperature(const Real last) override
    {
        this->last_ = last;
        if(this->last_ <= 0)
        {
            this->last_ = 1e-12;
        }
        coef_ = std::log(last_ / init_);
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
            const std::map<std::string, Real>& max_dprobe,
            image<Real>         ref,
            system<Real>        sys,
            std::unique_ptr<ObserverBase<Real>> obs,
            std::unique_ptr<ScoreBase<Real, Mask>> score,
            std::unique_ptr<ScheduleBase<Real>> schedule,
            const std::string& out)
        : step_(0),
          total_step_(total_step),
          save_step_(save),
          sigma_dx_(sigma_x),
          sigma_dy_(sigma_y),
          sigma_dz_(sigma_z),
          max_drotx_(max_rotx),
          max_droty_(max_roty),
          max_drotz_(max_rotz),
          dprobe_(max_dprobe),
          rng_(seed),
          nrm_(0.0, 1.0),
          reference_(std::move(ref)),
          sys_(sys),
          next_(sys),
          obs_(std::move(obs)),
          score_(std::move(score)),
          schedule_(std::move(schedule)),
          bar_(total_step),
          output_basename_(out)
    {
        sys_.cells.initialize(sys.stage_info.x_resolution(),
                              sys.stage_info.y_resolution(),
                              sys.particles);
        sys_.cells.construct(sys_.particles, sys_.bounding_box);
        this->img_ = obs_->observe(sys_);
        current_energy_ = score_->calc(sys_, img_, Mask(img_), reference_, Mask(img_));

        // clear output content
        {
            std::ofstream trj(output_basename_ + ".xyz");
        }
        {
            std::ofstream ene(output_basename_ + ".log");
            ene << "# step energy probe_shape\n";
        }

        if(schedule_->init_temperature() < 0.0)
        {
            this->warm_up();
        }
    }
    ~SimulatedAnnealingSimulator() override = default;

    void run() override
    {
        while(this->step()) {}
        return ;
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
        if(this->step_ >= total_step_)
        {
            // --------------------------------------------------------------------
            // output final result

            afmize::write_xyz(output_basename_, this->current_state());
            afmize::write_ppm(output_basename_ + "_result", this->current_image());
            afmize::write_tsv(output_basename_ + "_result", this->current_image());

            std::ofstream ene(output_basename_ + ".log", std::ios::app);
            ene << this->step_
                << " " << this->current_energy_;
            obs_->print_probe(ene);
            ene << "\n";
            std::cerr << bar_.format(this->step_);

            return false;
        }

        if(this->step_ % save_step_ == 0)
        {
            const auto fname = output_basename_ + "_" + std::to_string(this->current_step());
            afmize::write_xyz(output_basename_, this->current_state());
            afmize::write_ppm(fname,           this->current_image());
            afmize::write_tsv(fname,           this->current_image());

            std::ofstream ene(output_basename_ + ".log", std::ios::app);
            ene << this->step_
                << " " << this->current_energy_;
            obs_->print_probe(ene);
            ene << "\n";

            std::cerr << bar_.format(this->step_);
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

        // try to change probe
        this->try_probe_change(beta);

        this->step_ += 1;
        return true;
    }

    system<Real> const& current_state() const noexcept override {return sys_;}
    system<Real>&       current_state()       noexcept override {return sys_;}

    image<Real> const& current_image() const noexcept override
    {
        return img_;
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

        this->translate(dr, next_);

        // avoid particle from escaping observation stage.
        // To fix the number of pixels when calculating score, it does not
        // allow particle sticks out of the stage region.
        if(next_.bounding_box.lower[0] < sys_.stage_info.x_range().first  ||
           next_.bounding_box.lower[1] < sys_.stage_info.y_range().first  ||
           sys_.stage_info.x_range().second < next_.bounding_box.upper[0] ||
           sys_.stage_info.y_range().second < next_.bounding_box.upper[1])
        {
            return ;
        }

        // calculate score depending on score function
        const auto& next_img = obs_->observe(next_);
        const auto energy = score_->calc(next_, obs_->get_image(), Mask(next_img), reference_, Mask(next_img));
        const auto deltaE = energy - current_energy_;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            sys_ = next_;
            current_energy_ = energy;
            img_ = obs_->get_image();
        }
        return;
    }

    void try_rotation(const mave::matrix<Real, 3, 3>& rot, const Real beta)
    {
        next_ = sys_;

        this->rotate(rot, next_);

        // avoid particle from escaping observation stage
        if(next_.bounding_box.lower[0] < sys_.stage_info.x_range().first  ||
           next_.bounding_box.lower[1] < sys_.stage_info.y_range().first  ||
           sys_.stage_info.x_range().second < next_.bounding_box.upper[0] ||
           sys_.stage_info.y_range().second < next_.bounding_box.upper[1])
        {
            return ;
        }

        // calculate score depending on score function
        const auto& next_img = obs_->observe(next_);
        const auto energy = score_->calc(next_, obs_->get_image(), Mask(next_img),
                                                reference_, Mask(next_img));
        const auto deltaE = energy - current_energy_;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            sys_ = next_;
            current_energy_ = energy;
            img_ = obs_->get_image();
        }
        return;
    }

    void try_probe_change(const Real beta)
    {
        const auto prev_state = obs_->get_probe();
        auto probe_state = prev_state;
        for(auto& kv : probe_state)
        {
            kv.second += (this->generate_01() * 2.0 - 1.0) * dprobe_.at(kv.first);
        }
        if( ! obs_->update_probe(probe_state))
        {
            return ;
        }

        // calculate score depending on score function
        const auto& next_img = obs_->observe(next_);
        const auto energy = score_->calc(next_, obs_->get_image(), Mask(next_img), reference_, Mask(next_img));
        const auto deltaE = energy - current_energy_;

//         std::cerr << "beta = " << beta << ", dE = " << deltaE
//                   << ", Enext = " << energy << ", Ecurr = " << current_energy_
//                   << ", prob = " << std::exp(-deltaE * beta) << std::endl;

        if(deltaE <= 0.0 || this->generate_01() < std::exp(-deltaE * beta))
        {
            // system is not updated
            current_energy_ = energy;
            img_ = obs_->get_image();
        }
        else
        {
            // if energy increases, go back to the original probe
            obs_->update_probe(probe_state);
        }
        return;
    }

    Real generate_01() noexcept
    {
        return std::generate_canonical<Real, std::numeric_limits<Real>::digits>(rng_);
    }

    void translate(mave::vector<Real, 3> dr, system<Real>& target)
    {
        // translation does not change the shape of bounding box.
        target.bounding_box.upper += dr;
        target.bounding_box.lower += dr;

        if(target.bounding_box.lower[2] < 0.0)
        {
            const auto offset = target.bounding_box.lower[2];
            dr[2]                        -= offset;
            target.bounding_box.upper[2] -= offset;
            target.bounding_box.lower[2] -= offset;
        }

        // apply the movement
        for(auto& p : target.particles)
        {
            p.center += dr;
        }

        // update cell list
        target.cells.construct(target.particles, target.bounding_box);
        return ;
    }
    void rotate(const mave::matrix<Real, 3, 3>& rot, system<Real>& target)
    {
        // 1. move center to origin
        // 2. apply rotation matrix
        // 3. move back to the original center

        mave::vector<Real, 3> com(0.0, 0.0, 0.0);
        for(const auto& p : target.particles)
        {
            com += p.center;
        }
        com /= static_cast<Real>(target.particles.size());

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
        for(auto& p : target.particles)
        {
            mave::vector<Real, 4> r(p.center[0], p.center[1], p.center[2], 1.0);
            r = matrix * r;
            p.center[0] = r[0];
            p.center[1] = r[1];
            p.center[2] = r[2];
        }
        target.bounding_box = make_bounding_box(target.particles);

        if(target.bounding_box.lower[2] < 0.0)
        {
            // align the bottom to the xy plane (z=0.0)
            for(auto& p : target.particles)
            {
                p.center[2] -= target.bounding_box.lower[2];
            }
            target.bounding_box.upper[2] -= target.bounding_box.lower[2];
            target.bounding_box.lower[2]  = 0.0;
        }
        // update cell list
        target.cells.construct(target.particles, target.bounding_box);
        return;
    }

    void warm_up()
    {
        constexpr std::size_t N = 100;
        std::vector<Real> energy_differences;
        energy_differences.reserve(N);

        const auto initial_configuration = sys_;

        Real current_energy = this->current_energy_;
        std::bernoulli_distribution half_half(0.5);

        std::cerr << "collecting warm up information ..." << std::endl;
        afmize::progress_bar<70> bar(N);
        for(std::size_t i=0; i<N; ++i)
        {
            std::cerr << bar.format(i);
            next_ = sys_;

            const auto dx = nrm_(rng_) * sigma_dx_;
            const auto dy = nrm_(rng_) * sigma_dy_;

            translate(mave::vector<Real, 3>(dx, dy, 0.0), next_);

            // rot around x
            const auto rot_x = (this->generate_01() * 2.0 - 1.0) * max_drotx_;
            const auto cos_x = std::cos(rot_x);
            const auto sin_x = std::sin(rot_x);
            mave::matrix<Real, 3, 3> mat_x;
            mat_x.zero();
            mat_x(0, 0) = 1.0;
            mat_x(1, 1) =  cos_x;
            mat_x(1, 2) = -sin_x;
            mat_x(2, 1) =  sin_x;
            mat_x(2, 2) =  cos_x;

            // rot around y
            const auto rot_y = (this->generate_01() * 2.0 - 1.0) * max_droty_;
            const auto cos_y = std::cos(rot_y);
            const auto sin_y = std::sin(rot_y);
            mave::matrix<Real, 3, 3> mat_y;
            mat_y.zero();
            mat_y(0, 0) =  cos_y;
            mat_y(0, 2) =  sin_y;
            mat_y(1, 1) = 1.0;
            mat_y(2, 0) = -sin_y;
            mat_y(2, 2) =  cos_y;

            // rot around z
            const auto rot_z = (this->generate_01() * 2.0 - 1.0) * max_drotz_;
            const auto cos_z = std::cos(rot_z);
            const auto sin_z = std::sin(rot_z);
            mave::matrix<Real, 3, 3> mat_z;
            mat_z.zero();
            mat_z(0, 0) =  cos_z;
            mat_z(0, 1) = -sin_z;
            mat_z(1, 0) =  sin_z;
            mat_z(1, 1) =  cos_z;
            mat_z(2, 2) = 1.0;

            rotate(mat_z * mat_y * mat_x, next_);


            const auto& next_img = obs_->observe(next_);
            const auto energy = score_->calc(next_, next_img  , Mask(next_img),
                                                    reference_, Mask(next_img));
            const auto deltaE = energy - current_energy;

            if(0.0 < deltaE)
            {
                energy_differences.push_back(deltaE);
            }

            if(50 <= N && energy_differences.size() <= 10)
            {
                if(half_half(rng_))
                {
                    sys_ = next_;
                    current_energy = energy;
                }
            }
        }
        std::cerr << bar.format(N) << std::endl;

        this->sys_ = initial_configuration;

        std::sort(energy_differences.begin(), energy_differences.end());
        const auto initial_temperature = energy_differences.at(energy_differences.size() / 2) / std::log(2.0);
        this->schedule_->set_init_temperature(initial_temperature);
        return ;
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

    std::map<std::string, Real> dprobe_;

    std::mt19937                             rng_;
    std::normal_distribution<Real>           nrm_;
    image<Real>                              img_;
    image<Real>                        reference_;
    system<Real>                             sys_;
    system<Real>                            next_;
    std::unique_ptr<ObserverBase<Real>>      obs_;
    std::unique_ptr<ScoreBase<Real, Mask>> score_;
    std::unique_ptr<ScheduleBase<Real>> schedule_;
    afmize::progress_bar<70>                 bar_;
    std::string                  output_basename_;
};

} // afmize
#endif// AFMIZE_SIMULATOR_HPP
