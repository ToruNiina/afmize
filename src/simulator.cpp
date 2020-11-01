#include <afmize/simulator.hpp>
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
std::unique_ptr<ObserverBase<Real>> read_observation_method(const toml::value& sim)
{
    constexpr Real pi         = Real(3.1415926535897932384626);
    constexpr Real deg_to_rad = pi / Real(180.0);

    default_probe<Real> probe;
    probe.angle  = toml::find<Real>(sim, "probe", "angle") * deg_to_rad;
    probe.radius = read_as_angstrom<Real>(toml::find(sim, "probe", "radius"));

    const auto z_descritized = toml::find<bool>(sim, "image", "descritized");
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
read_score_function(const toml::value& sim)
{
    const auto method = toml::find<std::string>(sim, "score", "method");
    if(method == "correlation")
    {
        const auto k = toml::find<Real>(sim, "score", "k");
        return std::make_unique<NegativeCosineSimilarity<Real, Mask>>(k);
    }
    else if(method == "RMSD")
    {
        const auto k = toml::find<Real>(sim, "score", "k");
        return std::make_unique<RootMeanSquareDeviation<Real, Mask>>(k);
    }
    else if(method == "topographical penalty")
    {
        const auto penalty   = toml::find<Real>(sim, "score", "penalty");
        const auto reward    = toml::find<Real>(sim, "score", "reward");
        const auto thickness = read_as_angstrom<Real>(toml::find(sim, "score", "thickness"));
        return std::make_unique<TopographicalPenalty<Real, Mask>>(
                penalty, reward, thickness);
    }
    else if(method == "pixel penalty")
    {
        const auto penalty   = toml::find<Real>(sim, "score", "penalty");
        const auto reward    = toml::find<Real>(sim, "score", "reward");
        const auto thickness = read_as_angstrom<Real>(toml::find(sim, "score", "thickness"));
        return std::make_unique<PixelPenalty<Real, Mask>>(
                penalty, reward, thickness);
    }
    else
    {
        throw std::runtime_error("unknown score function: " + method);
    }
}

template<typename Real>
std::unique_ptr<ScheduleBase<Real>>
read_temperature_schedule(const toml::value& sim)
{
    const auto method = toml::find<std::string>(sim, "schedule", "method");
    if(method == "linear")
    {
        const auto total_t = toml::find<std::size_t>(sim, "algorithm", "steps");
        const auto initial = toml::find<Real>(sim, "schedule", "initial");
        const auto final   = toml::find<Real>(sim, "schedule", "final");

        return std::make_unique<LinearSchedule<Real>>(initial, final, total_t);
    }
    else if(method == "exponential")
    {
        const auto total_t = toml::find<std::size_t>(sim, "algorithm", "steps");
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
image<Real> read_reference_image(const toml::value& sim, const stage<Real>& stg)
{
    const auto refname = toml::find<std::string>(sim, "image", "reference");
    image<Real> img;

    if(refname.substr(refname.size() - 4, 4) == ".pdb" ||
       refname.substr(refname.size() - 4, 4) == ".xyz")
    {
        // generate reference image
        system<Real> answer(open_file<Real>(refname)->read_snapshot(), stg);

        if(answer.bounding_box.lower[2] != 0)
        {
            const mave::vector<Real, 3> offset{0, 0, answer.bounding_box.lower[2]};
            for(auto& p : answer.particles)
            {
                p.center -= offset;
            }
            answer.bounding_box.lower -= offset;
            answer.bounding_box.upper -= offset;
        }

        answer.cells.initialize(stg.x_resolution(), stg.y_resolution(),
                                answer.particles);
        answer.cells.construct(answer.particles, answer.bounding_box);

        auto obs = read_observation_method<Real>(sim);

        // XXX
        default_probe<Real> p;
        p.radius = 50.0; // 5 nm
        p.angle  = 20.0 * 3.1416 / 180.0;
        obs->update_probe(p);

        img = obs->observe(answer);

        const Real sigma = read_as_angstrom<Real>(toml::find(sim, "image", "noise"));
        std::random_device dev{};
        std::mt19937 mt(dev());
        apply_noise(img, mt, sigma);

        write_tsv("reference", img);
        write_ppm("reference", img);
        write_xyz("reference", answer);
    }

    return img;
}

template<typename Real, typename Mask>
std::unique_ptr<SimulatedAnnealingSimulator<Real, Mask>>
read_simulated_annealing_simulator(const toml::value& sim, system<Real> init)
{
    constexpr Real pi         = Real(3.1415926535897932384626);
    constexpr Real deg_to_rad = pi / Real(180.0);

    const auto steps = toml::find<std::size_t>(sim, "algorithm", "steps");
    const auto save_step = toml::find<std::size_t>(sim, "algorithm", "output");
    const auto seed = toml::find<std::uint32_t>(sim, "seed");

    const auto sx = afmize::read_as_angstrom<Real>(toml::find(sim, "algorithm", "sigma_x"));
    const auto sy = afmize::read_as_angstrom<Real>(toml::find(sim, "algorithm", "sigma_y"));
    const auto sz = afmize::read_as_angstrom<Real>(toml::find(sim, "algorithm", "sigma_z"));
    const auto rx = toml::find<Real>(sim, "algorithm", "maxrot_x") * deg_to_rad;
    const auto ry = toml::find<Real>(sim, "algorithm", "maxrot_y") * deg_to_rad;
    const auto rz = toml::find<Real>(sim, "algorithm", "maxrot_z") * deg_to_rad;

    const auto pr = afmize::read_as_angstrom<Real>(toml::find(sim, "algorithm", "max_dradius"));
    const auto pa = toml::find<Real>(sim, "algorithm", "max_dangle") * deg_to_rad;

    auto ref = read_reference_image(sim, init.stage_info);

    return std::make_unique<SimulatedAnnealingSimulator<Real, Mask>>(
        steps, save_step, seed, sx, sy, sz, rx, ry, rz, pr, pa,
        std::move(ref),
        std::move(init),
        read_observation_method<Real>(sim),
        read_score_function<Real, Mask>(sim),
        read_temperature_schedule<Real>(sim));
}

template<typename Real>
std::unique_ptr<SimulatorBase<Real>>
read_simulator(const toml::value& config, system<Real> init)
{
    const auto& sim = toml::find(config, "simulator");
    const auto algo = toml::find<std::string>(sim, "algorithm", "method");
    if(algo == "SA" || algo == "SimulatedAnnealing")
    {
        const auto mask  = toml::find<std::string>(sim, "image", "mask");
        const auto score = toml::find<std::string>(sim, "score", "method");
        if(mask == "rectangular")
        {
            return read_simulated_annealing_simulator<Real, mask_by_rectangle<Real>>(
                    sim, std::move(init));
        }
        else if(mask == "none")
        {
            return read_simulated_annealing_simulator<Real, mask_nothing<Real>>(
                    sim, std::move(init));
        }
        else
        {
            throw std::runtime_error("unknown mask: " + mask);
        }
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

    if(config.as_table().count("radii") == 1)
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

    // read input output files
    // [file]
    // input = "input.pdb"
    // output.basename = "output"
    const auto& file   = toml::find(config, "file");
    const auto  reader = afmize::open_file<Real>(toml::find<std::string>(file, "input"));
    const auto  output_basename = toml::find<std::string>(file, "output", "basename");

    // read stage setting
    //
    // [stage.resolution]
    // x = "1.0nm"
    // y = "1.0nm"
    // z = "0.64angstrom"
    // [stage.range]
    // x = ["0.0nm", "100.0nm"]
    // y = ["0.0nm", "100.0nm"]
    const auto& resolution = toml::find(config, "stage", "resolution");
    const auto& range      = toml::find(config, "stage", "range");
    const auto range_x = toml::find<std::array<toml::value, 2>>(range, "x");
    const auto range_y = toml::find<std::array<toml::value, 2>>(range, "y");
    afmize::stage<Real> stg(
        afmize::read_as_angstrom<Real>(toml::find(resolution, "x")),
        afmize::read_as_angstrom<Real>(toml::find(resolution, "y")),
        afmize::read_as_angstrom<Real>(toml::find(resolution, "z")),
        std::make_pair(afmize::read_as_angstrom<Real>(range_x[0]),
                       afmize::read_as_angstrom<Real>(range_x[1])),
        std::make_pair(afmize::read_as_angstrom<Real>(range_y[0]),
                       afmize::read_as_angstrom<Real>(range_y[1]))
        );

    // simulator
    //
    // [simulator]
    // seed              = 123456789
    // algorithm.method  = "SA"
    // algorithm.steps   = 1_000_000
    // algorithm.output  = 100
    // schedule.method   = "linear" | "exponential" # (required if "SA")
    // schedule.initial  = 10.0
    // schedule.final    =  0.0
    // probe.radius      = "1.0nm"
    // probe.angle       = 10.0 # degree
    // image.reference   = "pdb file" | "asd file"
    // image.descritized = false
    // image.mask        = "rectangular" | "none"
    // score.method      = "correlation" | "RMSD"
    // score.k           = 10.0


    afmize::system<Real> sys(reader->read_snapshot(), std::move(stg));
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

    {
        // clear file content
        std::ofstream ofs(output_basename + ".xyz");
    }
    {
        // write the first state
        afmize::write_xyz(output_basename,        sim->current_state());
        afmize::write_ppm(output_basename + "_0", sim->current_image());
        afmize::write_tsv(output_basename + "_0", sim->current_image());
    }

    const auto save_step = toml::find<std::size_t>(config, "simulator", "algorithm", "output");
    while(sim->step())
    {
        if(sim->current_step() % save_step == 0)
        {
            const auto fname = output_basename + "_"+ std::to_string(sim->current_step());
            afmize::write_ppm(fname, sim->current_image());
            afmize::write_tsv(fname, sim->current_image());
            afmize::write_xyz(output_basename, sim->current_state());
        }
    }
    {
        const auto fname = output_basename + "_"+ std::to_string(sim->current_step());
        afmize::write_ppm(fname, sim->current_image());
        afmize::write_tsv(fname, sim->current_image());
        afmize::write_xyz(output_basename, sim->current_state());
    }
    return 0;
}
