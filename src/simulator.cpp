#include <afmize/annealing_simulator.hpp>
#include <afmize/scanning_simulator.hpp>
#include <afmize/parameter.hpp>
#include <afmize/progress_bar.hpp>
#include <afmize/pdb_reader.hpp>
#include <afmize/xyz_reader.hpp>
#include <afmize/input_utility.hpp>
#include <afmize/output_utility.hpp>
#include <afmize/noise.hpp>
#include <afmize/colormap.hpp>
#include <random>

namespace afmize
{

template<typename Real>
std::unique_ptr<reader_base<Real>>
open_file(const std::string& fname)
{
    if(fname.substr(fname.size() - 4, 4) == ".pdb")
    {
        return std::make_unique<pdb_reader<Real>>(fname);
    }
    else if(fname.substr(fname.size() - 6, 6) == ".movie")
    {
        return std::make_unique<pdb_reader<Real>>(fname);
    }
    else if(fname.substr(fname.size() - 4, 4) == ".xyz")
    {
        return std::make_unique<xyz_reader<Real>>(fname);
    }
    else
    {
        throw std::runtime_error("afmize supports only pdb or xyz");
    }
}

template<typename Real>
std::unique_ptr<ObserverBase<Real>> read_observation_method(const toml::value& config)
{
    constexpr Real pi         = Real(3.1415926535897932384626);
    constexpr Real deg_to_rad = pi / Real(180.0);

    const auto& p = toml::find(config, "probe");
    default_probe<Real> probe;
    probe.angle  = toml::find<Real>(p, "angle") * deg_to_rad;
    probe.radius = read_as_angstrom<Real>(toml::find(p, "radius"));

    // If z resolution is set, that means height is descritized in an image.
    const auto z_descritized = config.at("stage").at("resolution").contains("z");
    if(z_descritized)
    {
        return std::make_unique<RigidObserver<Real, true>>(std::move(probe));
    }
    else
    {
        return std::make_unique<RigidObserver<Real, false>>(std::move(probe));
    }
}

template<typename Real, typename Mask>
std::unique_ptr<ScoreBase<Real, Mask>>
read_score_function(const toml::value& config)
{
    const auto& score = config.at("score");

    auto method = toml::find<std::string>(score, "method");
    std::transform(method.begin(), method.end(), method.begin(),
            [](const char c) -> char {return std::tolower(c);});

    if(method == "cosine similarity" || method == "cosinesimilarity")
    {
        const auto use_zero = toml::find<bool>(score, "use_zero_pixel_in_model");
        const auto k = toml::find<Real>(score, "k");
        if(use_zero)
        {
            return std::make_unique<NegativeCosineSimilarity<Real, Mask, true>>(k);
        }
        else
        {
            return std::make_unique<NegativeCosineSimilarity<Real, Mask, false>>(k);
        }
    }
    else if(method == "correlation")
    {
        const auto use_zero = toml::find<bool>(score, "use_zero_pixel_in_model");
        const auto k = toml::find<Real>(score, "k");
        if(use_zero)
        {
            return std::make_unique<NegativeCorrelation<Real, Mask, true>>(k);
        }
        else
        {
            return std::make_unique<NegativeCorrelation<Real, Mask, false>>(k);
        }
    }
    else if(method == "rmsd")
    {
        const auto use_zero = toml::find<bool>(score, "use_zero_pixel_in_model");
        const auto k = toml::find<Real>(score, "k");
        if(use_zero)
        {
            return std::make_unique<RootMeanSquareDeviation<Real, Mask, true>>(k);
        }
        else
        {
            return std::make_unique<RootMeanSquareDeviation<Real, Mask, false>>(k);
        }
    }
    else if(method == "topographical penalty")
    {
        const auto penalty   = toml::find<Real>(score, "penalty");
        const auto reward    = toml::find<Real>(score, "reward");
        const auto thickness = read_as_angstrom<Real>(toml::find(score, "thickness"));
        return std::make_unique<TopographicalPenalty<Real, Mask>>(
                penalty, reward, thickness);
    }
    else if(method == "pixel penalty")
    {
        const auto penalty   = toml::find<Real>(score, "penalty");
        const auto reward    = toml::find<Real>(score, "reward");
        const auto thickness = read_as_angstrom<Real>(toml::find(score, "thickness"));
        return std::make_unique<PixelPenalty<Real, Mask>>(penalty, reward, thickness);
    }
    else
    {
        throw std::runtime_error("unknown score function: " + method);
    }
}

template<typename Real>
std::unique_ptr<ScheduleBase<Real>>
read_temperature_schedule(const toml::value& config)
{
    const auto& sim = toml::find(config, "simulator");
    const auto method = toml::find<std::string>(sim, "schedule", "method");
    if(method == "linear")
    {
        const auto total_t = toml::find<std::size_t>(sim, "steps");
        const auto initial = toml::find<Real>(sim, "schedule", "initial");
        const auto final   = toml::find<Real>(sim, "schedule", "final");

        return std::make_unique<LinearSchedule<Real>>(initial, final, total_t);
    }
    else if(method == "exponential")
    {
        const auto total_t = toml::find<std::size_t>(sim, "steps");
        const auto initial = toml::find<Real>(sim, "schedule", "initial");
        const auto final   = toml::find<Real>(sim, "schedule", "final");

        return std::make_unique<ExponentialSchedule<Real>>(initial, final, total_t);
    }
    else
    {
        throw std::runtime_error("unknown schedule method: " + method);
    }
}

template<typename Real>
image<Real> read_reference_image(const toml::value& config, const stage<Real>& stg)
{
    const auto refname = toml::find<std::string>(config, "image", "reference");
    if(refname.substr(refname.size() - 4, 4) == ".tsv")
    {
        const auto x = toml::find<std::size_t>(config, "image", "pixels", "x");
        const auto y = toml::find<std::size_t>(config, "image", "pixels", "y");
        image<Real> img(x, y);

        std::ifstream ifs(refname);
        if(!ifs.good())
        {
            throw std::runtime_error("file open error: " + refname);
        }
        for(std::size_t y=0; y<img.y_pixel(); ++y)
        {
            for(std::size_t x=0; x<img.x_pixel(); ++x)
            {
                ifs >> img.at(x, y);
            }
        }
        return img;
    }
    else
    {
        throw std::runtime_error("[error] image format not supported: " + refname);
    }
}

template<typename Real>
system<Real> read_system(const toml::value& config)
{
    const auto reader = open_file<Real>(toml::find<std::string>(config, "initial", "input"));

    // image lower coordinate
    const auto lower_x = read_as_angstrom<Real>(toml::find(config, "image", "lower", "x"));
    const auto lower_y = read_as_angstrom<Real>(toml::find(config, "image", "lower", "y"));

    // read stage setting
    // [stage]
    // pixels.x     = 80
    // pixels.y     = 80
    // resolution.x = "1.0nm"
    // resolution.y = "1.0nm"
    const auto pixel_x = toml::find<std::size_t>(config, "stage", "pixels", "x");
    const auto pixel_y = toml::find<std::size_t>(config, "stage", "pixels", "y");

    const auto& resolution = toml::find(config, "stage", "resolution");
    const auto reso_x = read_as_angstrom<Real>(toml::find(resolution, "x"));
    const auto reso_y = read_as_angstrom<Real>(toml::find(resolution, "y"));
    const auto reso_z = read_as_angstrom<Real>(
            toml::find_or(resolution, "z", toml::value(0.0)));

    stage<Real> stg(reso_x, reso_y, reso_z,
        std::make_pair(lower_x, lower_x + reso_x * pixel_x),
        std::make_pair(lower_y, lower_y + reso_y * pixel_y));

    return system<Real>(reader->read_snapshot(), std::move(stg));
}

template<typename Real, typename Mask>
std::unique_ptr<SimulatedAnnealingSimulator<Real, Mask>>
read_simulated_annealing_simulator(const toml::value& config, system<Real> init)
{
    constexpr Real pi         = Real(3.1415926535897932384626);
    constexpr Real deg_to_rad = pi / Real(180.0);
    const auto& sim = toml::find(config, "simulator");
    const auto  out = toml::find<std::string>(sim, "output");

    const auto seed  = toml::find<std::uint32_t>(sim, "seed");
    const auto steps = toml::find<std::size_t  >(sim, "steps");
    const auto save  = toml::find<std::size_t  >(sim, "save");

    const auto sx = afmize::read_as_angstrom<Real>(toml::find(sim, "sigma_x"));
    const auto sy = afmize::read_as_angstrom<Real>(toml::find(sim, "sigma_y"));
    const auto sz = afmize::read_as_angstrom<Real>(toml::find(sim, "sigma_z"));
    const auto rx = toml::find<Real>(sim, "maxrot_x") * deg_to_rad;
    const auto ry = toml::find<Real>(sim, "maxrot_y") * deg_to_rad;
    const auto rz = toml::find<Real>(sim, "maxrot_z") * deg_to_rad;

    const auto pr = afmize::read_as_angstrom<Real>(toml::find(sim, "max_dradius"));
    const auto pa = toml::find<Real>(sim, "max_dangle") * deg_to_rad;

    auto ref = read_reference_image(config, init.stage_info);

    return std::make_unique<SimulatedAnnealingSimulator<Real, Mask>>(
        steps, save, seed, sx, sy, sz, rx, ry, rz, pr, pa,
        std::move(ref),
        std::move(init),
        read_observation_method<Real>(config),
        read_score_function<Real, Mask>(config),
        read_temperature_schedule<Real>(config),
        out);
}

template<typename Real, typename Mask>
std::unique_ptr<ScanningSimulator<Real, Mask>>
read_scanning_simulator(const toml::value& config, system<Real> init)
{
    const auto out = toml::find<std::string>(config, "simulator", "output");

    const auto& sim = toml::find(config, "simulator");
    const auto num_div  = toml::find<std::size_t>(sim, "num_division");
    const auto num_save = toml::find<std::size_t>(sim, "save");
    const auto dz       = read_as_angstrom<Real>(toml::find(sim, "dz"));

    auto ref = read_reference_image(config, init.stage_info);

    return std::make_unique<ScanningSimulator<Real, Mask>>(
        num_div, num_save, dz,
        std::move(ref),
        std::move(init),
        read_observation_method<Real>(config),
        read_score_function<Real, Mask>(config),
        out);
}

template<typename Real>
std::unique_ptr<SimulatorBase<Real>>
read_simulator(const toml::value& config, system<Real> init)
{
    const auto& sim = toml::find(config, "simulator");
    const auto algo = toml::find<std::string>(sim, "method");
    if(algo == "SA" || algo == "SimulatedAnnealing")
    {
        const auto mask  = toml::find<std::string>(config, "score", "mask");
        if(mask == "rectangular")
        {
            return read_simulated_annealing_simulator<Real, mask_by_rectangle<Real>>(
                    config, std::move(init));
        }
        else if(mask == "none")
        {
            return read_simulated_annealing_simulator<Real, mask_nothing<Real>>(
                    config, std::move(init));
        }
        else
        {
            throw std::runtime_error("unknown mask: " + mask);
        }
    }
    else if(algo == "Scanning")
    {
        return read_scanning_simulator<Real, mask_by_rectangle<Real>>(
                config, std::move(init));
    }
    else
    {
        throw std::runtime_error("unknown simulation algorhtm: " + algo);
    }
}

} // afmize

int main(int argc, char** argv)
{
    using Real = double;

    if(argc != 2)
    {
        std::cerr << "usage: afmize <input_file.toml>\n";
        return 1;
    }

    const std::string config_file(argv[1]);
    if(config_file.substr(config_file.size() - 5, 5) != ".toml")
    {
        std::cerr << "afmize requires toml file as an input\n";
        return 1;
    }

    const auto config = toml::parse(config_file);

    // -----------------------------------------------------------------------
    // update global parameter, radii

    if(config.contains("radii"))
    {
        const auto& radii = toml::find(config, "radii");
        for(const auto& kv : toml::find_or<toml::table>(radii, "atom", toml::table{}))
        {
            afmize::parameter<Real>::radius_atom[kv.first] =
                afmize::read_as_angstrom<Real>(kv.second);
        }
        for(const auto& res : toml::find_or<toml::table>(radii, "residue", toml::table{}))
        {
            for(const auto& atm : res.second.as_table())
            {
                afmize::parameter<Real>::radius_residue[res.first][atm.first] =
                    afmize::read_as_angstrom<Real>(atm.second);
            }
        }
    }

    for(const auto& kv : afmize::parameter<Real>::radius_atom)
    {
        std::cerr << "-- radius of ATOM:" << std::setw(4) << kv.first
                  << " = " << kv.second << '\n';
    }
    for(const auto& kv1 : afmize::parameter<Real>::radius_residue)
    {
        for(const auto& kv2 : kv1.second)
        {
            std::cerr << "-- radius of ATOM:" << std::setw(4) << kv2.first
                      << " in RES:" << std::setw(5) << kv1.first
                      << " = " << kv2.second << '\n';
        }
    }

    afmize::system<Real> sys = afmize::read_system<Real>(config);
    // align the bottom to z == 0
    if(sys.bounding_box.lower[2] != 0)
    {
        const mave::vector<Real, 3> offset{0, 0, sys.bounding_box.lower[2]};
        for(auto& p : sys.particles)
        {
            p.center -= offset;
        }
        sys.bounding_box.lower -= offset;
        sys.bounding_box.upper -= offset;
    }

    auto sim = afmize::read_simulator(config, std::move(sys));
    sim->run();

    return 0;
}
