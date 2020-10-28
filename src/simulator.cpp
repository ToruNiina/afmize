#include <afmize/simulator.hpp>
#include <afmize/parameter.hpp>
#include <afmize/progress_bar.hpp>
#include <afmize/pdb_reader.hpp>
#include <afmize/xyz_reader.hpp>
#include <afmize/input_utility.hpp>
#include <afmize/colormap.hpp>

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
void write_xyz(const std::string& basename, const system<Real>& sys)
{
    std::ofstream ofs(basename + ".xyz", std::ios_base::app);
    ofs << sys.particles.size() << "\n\n";
    for(std::size_t i=0; i<sys.particles.size(); ++i)
    {
        const auto pos = sys.particles.at(i).center;
        ofs << std::setw(6) << sys.names.at(i) << ' '
            << std::setprecision(6) << std::setw(10) << pos[0]
            << std::setprecision(6) << std::setw(10) << pos[1]
            << std::setprecision(6) << std::setw(10) << pos[2]
            << '\n';
    }
    return ;
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
stage<Real> read_reference_image(const toml::value& sim, const stage<Real>& stg)
{
    stage<Real>  stage(stg);
    const auto refname = toml::find<std::string>(sim, "image", "reference");

    if(refname.substr(refname.size() - 4, 4) == ".pdb")
    {
        // generate reference image
        system<Real> answer(open_file<Real>(refname)->read_snapshot());

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

        answer.cells.initialize(stage.x_resolution(), stage.y_resolution(),
                                answer.particles);
        answer.cells.construct(answer.particles, answer.bounding_box);

        const auto obs = read_observation_method<Real>(sim);

        obs->observe(stage, answer);
    }
    return stage;
}

template<typename Real, typename Mask>
std::unique_ptr<SimulatedAnnealingSimulator<Real, Mask>>
read_simulated_annealing_simulator(const toml::value& sim, stage<Real> stg, system<Real> init)
{
    const auto steps = toml::find<std::size_t>(sim, "algorithm", "steps");
    const auto save_step = toml::find<std::size_t>(sim, "algorithm", "output");
    const auto seed = toml::find<std::uint32_t>(sim, "seed");

    auto ref = read_reference_image(sim, stg);

    return std::make_unique<SimulatedAnnealingSimulator<Real, Mask>>(
        steps, save_step, seed,
        std::move(ref),
        std::move(stg),
        std::move(init),
        read_observation_method<Real>(sim),
        read_score_function<Real, Mask>(sim),
        read_temperature_schedule<Real>(sim));
}

template<typename Real>
std::unique_ptr<SimulatorBase<Real>> read_simulator(const toml::value& config,
        stage<Real> stg, system<Real> init)
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
                    sim, std::move(stg), std::move(init));
        }
        else if(mask == "none")
        {
            return read_simulated_annealing_simulator<Real, mask_nothing<Real>>(
                    sim, std::move(stg), std::move(init));
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

template<typename Real>
void write_tsv(const std::string& out, const stage<Real>& stg)
{
    using namespace std::literals::string_literals;
    std::ofstream ofs(out);

    for(std::size_t y=0; y<stg.y_pixel(); ++y)
    {
        for(std::size_t x=0; x<stg.x_pixel(); ++x)
        {
            ofs << "\t" << std::setprecision(6) << std::setw(10) << stg.at(x, y);
        }
        ofs << "\n";
    }
    return;
}


template<typename Real>
void write_ppm(const stage<Real>& stg, const std::string& out,
               const std::pair<Real, Real> height_range)
{
    using namespace std::literals::string_literals;
    const auto min_elem = height_range.first;
    const auto max_elem = height_range.second;

    pnm::image<pnm::rgb_pixel> ppm(stg.x_pixel(), stg.y_pixel());
    for(std::size_t i=0; i<stg.x_pixel() * stg.y_pixel(); ++i)
    {
        ppm.raw_access(i) =
            color_afmhot<Real, pnm::rgb_pixel>(stg[i], min_elem, max_elem);
    }

    // origin of the image is upper left, but the physical origin is lower left.
    // to make it consistent, it simply reverses the image.
    pnm::image<pnm::rgb_pixel> reversed(ppm.x_size(), ppm.y_size());
    for(std::size_t i = 0; i<ppm.y_size(); ++i)
    {
        reversed[ppm.y_size() - i - 1] = ppm[i];
    }

    pnm::write(out + ".ppm"s, reversed, pnm::format::ascii);
    return;
}

template<typename Real>
void write_ppm(const stage<Real>& stg, const std::string& out)
{
    const auto minmax = std::minmax_element(stg.begin(), stg.end());
    const auto min_elem = *minmax.first;
    const auto max_elem = *minmax.second;

    write_ppm(stg, out, std::make_pair(min_elem, max_elem));
    return;
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

    auto sim = afmize::read_simulator(config, std::move(stg),
            afmize::system<Real>(reader->read_snapshot()));

    const auto save_step = toml::find<std::size_t>(config, "simulator", "algorithm", "output");
    while(sim->step())
    {
        if(sim->current_step() % save_step == 0)
        {
            const auto fname = output_basename + "_"+ std::to_string(sim->current_step());
            afmize::write_ppm(sim->current_image(), fname);
            afmize::write_xyz(fname, sim->current_state());
            afmize::write_tsv(fname + ".tsv", sim->current_image());
        }
    }
    afmize::write_xyz(output_basename, sim->current_state());
    return 0;
}
