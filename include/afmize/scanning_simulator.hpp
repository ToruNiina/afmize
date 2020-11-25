#ifndef AFMIZE_SCANNING_SIMULATOR_HPP
#define AFMIZE_SCANNING_SIMULATOR_HPP
#include "simulator_base.hpp"
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

template<typename Real, typename Mask>
struct ScanningSimulator : public SimulatorBase<Real>
{
    constexpr static Real pi = 3.14159265;

    struct location
    {
        mave::matrix<Real, 3, 3> rot;
        std::size_t x_offset;
        std::size_t y_offset;
        std::size_t z_offset;
    };

    ScanningSimulator(const std::size_t num_div, const std::size_t num_save,
        Real dz_,
        image<Real> ref, system<Real> sys,
        std::unique_ptr<ObserverBase<Real>>    obs,
        std::unique_ptr<ScoreBase<Real, Mask>> score,
        std::string out)
        : step_(0), num_save_(num_save), next_output_(0.0), doutput_percent_(1.0),
          num_div_(num_div),
          dtheta_(2 * 3.14159265 / num_div),
          dz_(dz_),
          max_height_(*std::max_element(ref.begin(), ref.end())),
          img_(ref.x_pixel(), ref.y_pixel()),
          reference_(std::move(ref)),
          sys_(sys),
          init_(sys),
          obs_(std::move(obs)),
          score_(std::move(score)),
          bar_(0),
          output_basename_(std::move(out))
    {
        sys_.cells.initialize(sys.stage_info.x_resolution(),
                              sys.stage_info.y_resolution(),
                              sys.particles);
        sys_.cells.construct(sys_.particles, sys_.bounding_box);
        init_.cells.initialize(sys.stage_info.x_resolution(),
                              sys.stage_info.y_resolution(),
                              sys.particles);
        init_.cells.construct(init_.particles, init_.bounding_box);
        this->center_[0] = 0.0;
        this->center_[1] = 0.0;
        this->center_[2] = 0.0;
        for(const auto& p : this->sys_.particles)
        {
            this->center_ += p.center;
        }
        this->center_ /= static_cast<Real>(sys_.particles.size());

        // --------------------------------------------------------------------
        // cache scanning angles
        if(num_div_ % 4 != 0)
        {
            throw std::runtime_error("Scanning: number of divisions should be "
                                     "a multiple of 4.");
        }
        const auto drot_y = 2 * pi / num_div_;

        // +z
        this->axes_rot_.push_back(mave::matrix<Real, 3, 3>( 1.0,  0.0,  0.0,
                                                            0.0,  1.0,  0.0,
                                                            0.0,  0.0,  1.0));
        // upper halves
        for(std::size_t i=1; i<num_div_/4; ++i)
        {
            const auto cos_y = std::cos(drot_y * i);
            const auto sin_y = std::sin(drot_y * i);
            const mave::matrix<Real, 3, 3> rot_y1( cos_y, 0.0, sin_y,
                                                     0.0, 1.0,   0.0,
                                                  -sin_y, 0.0, cos_y);
            const auto drot_z = (2 * pi) / (i * 4);
            for(std::size_t j=0; j<i*4; ++j)
            {
                const auto cos_z = std::cos(drot_z * j);
                const auto sin_z = std::sin(drot_z * j);

                const mave::matrix<Real, 3, 3> rot_z(cos_z, -sin_z, 0.0,
                                                     sin_z,  cos_z, 0.0,
                                                       0.0,    0.0, 1.0);
                axes_rot_.push_back(rot_z * rot_y1);
            }
        }
        // on Equator
        {
            const auto cos_y = 0.0;
            const auto sin_y = 1.0;
            const mave::matrix<Real, 3, 3> rot_y( cos_y, 0.0, sin_y,
                                                    0.0, 1.0,   0.0,
                                                 -sin_y, 0.0, cos_y);
            const auto drot_z = 2 * pi / num_div_;
            for(std::size_t j=0; j<num_div_; ++j)
            {
                const auto cos_z = std::cos(drot_z * j);
                const auto sin_z = std::sin(drot_z * j);

                const mave::matrix<Real, 3, 3> rot_z(cos_z, -sin_z, 0.0,
                                                     sin_z,  cos_z, 0.0,
                                                       0.0,    0.0, 1.0);
                axes_rot_.push_back(rot_z * rot_y);
            }
        }
        // -z
        this->axes_rot_.push_back(mave::matrix<Real, 3, 3>(-1.0,  0.0,  0.0,
                                                            0.0, -1.0,  0.0,
                                                            0.0,  0.0, -1.0));
        // lower halves
        for(std::size_t i=1; i<num_div_/4; ++i)
        {
            const auto cos_y = std::cos(drot_y * i);
            const auto sin_y = std::sin(drot_y * i);
            const mave::matrix<Real, 3, 3> rot_y2(-cos_y, 0.0,  sin_y,
                                                     0.0, 1.0,    0.0,
                                                  -sin_y, 0.0, -cos_y);

            const auto drot_z = (2 * pi) / (i * 4);
            for(std::size_t j=0; j<i*4; ++j)
            {
                const auto cos_z = std::cos(drot_z * j);
                const auto sin_z = std::sin(drot_z * j);

                const mave::matrix<Real, 3, 3> rot_z(cos_z, -sin_z, 0.0,
                                                     sin_z,  cos_z, 0.0,
                                                       0.0,    0.0, 1.0);

                axes_rot_.push_back(rot_z * rot_y2);
            }
        }

        if(axes_rot_.size() > 10000)
        {
            doutput_percent_ = 0.1;
        }
        else
        {
            doutput_percent_ = 1.0;
        }

        std::cout << axes_rot_.size() << " rotational configuration will be searched" << std::endl;

        // setup progress bar
        bar_.reset_total(axes_rot_.size());
    }
    ~ScanningSimulator() override = default;

    void run() override
    {
        while(this->step()) {}
        std::cerr << bar_.format(this->step_);
    }
    bool run(const std::size_t steps) override
    {
        const auto until = std::min(this->step_ + steps, this->axes_rot_.size());
        while(step_ < until)
        {
            this->step();
        }
        std::cerr << bar_.format(this->step_);
        return step_ < this->axes_rot_.size();
    }

    bool step() override
    {
        std::cerr << bar_.format(this->step_);
        if(this->step_ == this->axes_rot_.size())
        {
            this->output_status();
            return false;
        }

        mave::vector<Real, 3> vtx1(0.0, 0.0, 1.0);
        mave::vector<Real, 3> vtx2(0.0, 0.0, 0.0);
        mave::vector<Real, 3> vtx3(1.0, 0.0, 0.0);

        const auto n = axes_rot_.at(this->step_) *
                       mave::vector<Real, 3>(0.0, 0.0, 1.0);
        const auto dtheta = 2 * pi / num_div_;
        for(std::size_t i=0; i<num_div_; ++i)
        {
            sys_ = init_;
            const auto theta = dtheta * i;
            const auto sin_t = std::sin(theta);
            const auto cos_t = std::cos(theta);
            mave::matrix<Real, 3, 3> rot;
            rot.zero();
            rot(0, 0) = n[0] * n[0] * (1.0 - cos_t) +        cos_t;
            rot(0, 1) = n[0] * n[1] * (1.0 - cos_t) - n[2] * sin_t;
            rot(0, 2) = n[0] * n[2] * (1.0 - cos_t) + n[1] * sin_t;

            rot(1, 0) = n[1] * n[0] * (1.0 - cos_t) + n[2] * sin_t;
            rot(1, 1) = n[1] * n[1] * (1.0 - cos_t) +        cos_t;
            rot(1, 2) = n[1] * n[2] * (1.0 - cos_t) - n[0] * sin_t;

            rot(2, 0) = n[2] * n[0] * (1.0 - cos_t) - n[1] * sin_t;
            rot(2, 1) = n[2] * n[1] * (1.0 - cos_t) + n[0] * sin_t;
            rot(2, 2) = n[2] * n[2] * (1.0 - cos_t) +        cos_t;

            const auto mat = rot * axes_rot_.at(this->step_);

            for(auto& p : this->sys_.particles)
            {
                p.center -= this->center_;
                p.center  = mat * p.center;
                p.center += this->center_;
            }

            sys_.bounding_box = make_bounding_box(sys_.particles);
            // align the bottom to the xy plane (z=0.0)
            for(auto& p : this->sys_.particles)
            {
                p.center[2] -= sys_.bounding_box.lower[2];
            }
            sys_.bounding_box.upper[2] -= sys_.bounding_box.lower[2];
            sys_.bounding_box.lower[2]  = 0.0;

            sys_.cells.construct(sys_.particles, sys_.bounding_box);

            const auto v1 = mat * vtx1;
            const auto v2 = mat * vtx2;
            const auto v3 = mat * vtx3;

            this->scan_translation(location{mat, 0, 0});
        }
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

    std::size_t total_step()   const noexcept override {return this->axes_rot_.size();}
    std::size_t current_step() const noexcept override {return step_;}

  private:

    void scan_translation(location loc)
    {
        const std::size_t z_len = 1 + std::ceil(
                std::max(0.0, this->max_height_ - sys_.bounding_box.upper[2]) /
                this->dz_);

        for(std::size_t z_ofs=0; z_ofs < z_len; ++z_ofs)
        {
            loc.z_offset = z_ofs;
            auto img = obs_->observe(sys_);

            const Mask mask(img);
            const std::size_t x_rem = reference_.x_pixel() - mask.pixel_x();
            const std::size_t y_rem = reference_.y_pixel() - mask.pixel_y();

            for(std::size_t y_ofs=0; y_ofs < y_rem; ++y_ofs)
            {
                loc.y_offset = y_ofs;
                for(std::size_t x_ofs=0; x_ofs < x_rem; ++x_ofs)
                {
                    loc.x_offset = x_ofs;

                    const Mask target_mask(sys_, x_ofs, mask.pixel_x(),
                                                 y_ofs, mask.pixel_y());
                    const auto penalty = this->score_->calc(sys_, img, mask, reference_, target_mask);

                    const auto found = std::lower_bound(high_score_.begin(), high_score_.end(),
                            penalty, [](const auto& lhs, const Real& p) {return lhs.second < p;});

                    if(high_score_.size() < num_save_ || found != high_score_.end())
                    {
                        high_score_.insert(found, std::make_pair(loc, penalty));
                        if(num_save_ < high_score_.size())
                        {
                            high_score_.pop_back();
                        }
                    }
                }
            }
            for(auto& p : sys_.particles)
            {
                p.center[2] += this->dz_;
            }
        }
        return;
    }
    void output_status()
    {
        std::cerr << "\nwriting best " << high_score_.size() << " conformations ... ";
        std::ofstream trj(output_basename_ + ".xyz");
        std::ofstream ene(output_basename_ + ".log");
        ene << "# idx energy\n";

        const auto width = std::to_string(high_score_.size()).size();

        std::size_t idx = 0;
        for(const auto& best : high_score_)
        {
            sys_ = init_;
            const Mask mask(sys_);

            const auto& loc = best.first;
            const auto& rot = loc.rot;

            for(auto& p : this->sys_.particles)
            {
                p.center -= this->center_;
                p.center  = rot * p.center;
                p.center += this->center_;
            }
            sys_.bounding_box = make_bounding_box(sys_.particles);

            const auto trans = mave::vector<Real, 3>(
                sys_.stage_info.x_resolution() * (static_cast<std::int64_t>(loc.x_offset) - static_cast<std::int64_t>(mask.lower_bounding_x())),
                sys_.stage_info.y_resolution() * (static_cast<std::int64_t>(loc.y_offset) - static_cast<std::int64_t>(mask.lower_bounding_y())),
               -sys_.bounding_box.lower[2] + this->dz_ * loc.z_offset);

            // align the bottom to the xy plane (z=0.0)
            for(auto& p : this->sys_.particles)
            {
                p.center += trans;
            }
            sys_.bounding_box.upper += trans;
            sys_.bounding_box.lower += trans;

            sys_.cells.construct(sys_.particles, sys_.bounding_box);

            const auto& img = obs_->observe(sys_);
            std::ostringstream oss;
            oss << output_basename_ << "_" << std::setw(width) << std::setfill('0') << idx;

            afmize::write_svg(oss.str(), img, sys_, 0.0);
            afmize::write_tsv(oss.str(), img);
            afmize::write_xyz(output_basename_, this->sys_);

            ene << idx++ << ' ' << best.second << '\n';
        }
        std::cerr << "done.\n";
    }

  private:
    std::size_t step_;
    std::size_t num_save_;
    Real        next_output_;
    Real        doutput_percent_;

    std::size_t num_div_;
    Real        dtheta_;
    std::vector<mave::matrix<Real, 3, 3>> axes_rot_;
    mave::vector<Real, 3> center_;

    Real dz_;
    Real max_height_;

    // pairof{pairof{axis rotation, rotation around axis}, score}
    std::vector<std::pair<location, Real>> high_score_;

    image<Real>                              img_;
    image<Real>                        reference_;
    system<Real>                             sys_;
    system<Real>                            init_;
    std::unique_ptr<ObserverBase<Real>>      obs_;
    std::unique_ptr<ScoreBase<Real, Mask>> score_;
    afmize::progress_bar<70>                 bar_;
    std::string                  output_basename_;
};

} // afmize
#endif// AFMIZE_SIMULATOR_HPP
