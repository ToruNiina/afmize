#define TOML11_COLORIZE_ERROR_MESSAGE
#include <afmize/stage.hpp>
#include <afmize/system.hpp>
#include <afmize/colormap.hpp>
#include <afmize/xyz_reader.hpp>
#include <afmize/pdb_reader.hpp>
#include <afmize/input_utility.hpp>
#include <afmize/output_utility.hpp>
#include <afmize/progress_bar.hpp>
#include <afmize/observe.hpp>
#include <afmize/noise.hpp>
#include <limits>
#include <string>

namespace afmize
{

template<typename Real>
std::unique_ptr<afmize::reader_base<Real>>
open_file(const std::string& fname)
{
    if(fname.substr(fname.size() - 4, 4) == ".pdb")
    {
        return std::make_unique<afmize::pdb_reader<Real>>(fname);
    }
    else if(fname.substr(fname.size() - 6, 6) == ".movie")
    {
        return std::make_unique<afmize::pdb_reader<Real>>(fname);
    }
    else if(fname.substr(fname.size() - 4, 4) == ".xyz")
    {
        return std::make_unique<afmize::xyz_reader<Real>>(fname);
    }
    else
    {
        throw std::runtime_error("afmize supports only pdb or xyz");
    }
}

struct output_format_flags
{
    output_format_flags() : csv(false), tsv(false), json(false), ppm(false), svg(false) {}

    static output_format_flags
    from_inputs(std::vector<std::string> inputs)
    {
        output_format_flags self;
        for(const auto& format : inputs)
        {
            if(format == "csv")
            {
                set_flag(self.csv, format);
            }
            if(format == "tsv")
            {
                set_flag(self.tsv, format);
            }
            else if(format == "json")
            {
                set_flag(self.json, format);
            }
            else if(format == "ppm")
            {
                set_flag(self.ppm, format);
            }
            else if(format == "svg")
            {
                set_flag(self.svg, format);
            }
            else
            {
                std::cerr << "warning: file.output.formats contains an invalid flag \""
                          << format << "\"\n";
            }
        }
        return self;
    }

    static void
    set_flag(bool& flag, const std::string& input)
    {
        if(flag)
        {
            std::cerr << "warning: file.output.formats contains duplicated flags \""
                      << input << "\"\n";
        }
        else
        {
            flag = true;
        }
    }

    bool csv;
    bool tsv;
    bool json;
    bool ppm;
    bool svg;
};


} // afmize

int main(int argc, char** argv)
{
    using namespace std::literals::string_literals;

    using Real = double;
    constexpr Real pi             = Real(3.1415926535897932384626);
    constexpr Real deg_to_rad     = pi / Real(180.0);

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

    std::cerr << "reading config file..." << std::endl;

    const auto config = toml::parse(config_file);

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
    const auto& file   = toml::find(config, "file");
    const auto  reader = afmize::open_file<Real>(
            toml::find<std::string>(file, "input")
        );
    const auto& output = toml::find(file, "output");
    const auto  output_basename = toml::find<std::string>(output, "basename");
    const auto  output_formats = afmize::output_format_flags::from_inputs(
            toml::find<std::vector<std::string>>(output, "formats"));

    std::cerr << "-- " << reader->size() << " snapshots are found\n";

    // image size information
    const auto& resolution = toml::find(config, "resolution");
    const auto& range      = toml::find(config, "range");
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

    const auto scale_bar_length = afmize::read_as_angstrom<Real>(
        toml::find_or<toml::value>(config, "scale_bar", toml::value{{"length", toml::value(0.0)}}).at("length"));

    // probe size information
    const auto& probe_tab  = toml::find(config, "probe");
    const auto& probe_size = toml::find(probe_tab, "size");
    afmize::default_probe<Real> probe{
            toml::find<Real>(probe_size, "angle") * deg_to_rad,
            afmize::read_as_angstrom<Real>(toml::find(probe_size, "radius")),
            mave::vector<Real, 3>{0, 0, 0}
        };

    // stage information ...
    const toml::value nan_v(std::numeric_limits<Real>::quiet_NaN());
    const toml::value stage_tab = toml::find_or<toml::value>(config, "stage", toml::value(toml::table{}));
    const bool stage_align      = toml::find_or<bool>(stage_tab, "align", false);
    const Real stage_position   = [&]() -> Real {
        if(stage_tab.contains("position"))
        {
            return afmize::read_as_angstrom<Real>(toml::find(stage_tab, "position"));
        }
        else
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
    }();

    if(stage_align && std::isnan(stage_position))
    {
        std::cerr << toml::format_error("[error] "
                "`stage.align` requires `stage.position`.", stage_tab,
                "in this table", {"In order to align system on the stage, "
                "you need to input the z-coordinate position of your stage."
                }) << std::endl;
        return EXIT_FAILURE;
    }

    const Real noise_sigma = [&]() -> Real {
        if(config.contains("noise"))
        {
            return afmize::read_as_angstrom<Real>(toml::find(config, "noise"));
        }
        else
        {
            return 0.0;
        }
    }();
    std::random_device dev;
    std::mt19937 rng(dev());

    // color range information ...
    const auto& cmap = toml::find_or(config, "colormap", toml::value{});
    const auto cmap_min = afmize::read_as_angstrom<Real>(
        toml::find_or(cmap, "min", nan_v));
    const auto cmap_max = afmize::read_as_angstrom<Real>(
        toml::find_or(cmap, "max", nan_v));

    // image generation method
    const auto method = toml::find_or(config, "method", "rigid");
    if(method != "rigid" && method != "smooth")
    {
        std::cerr << toml::format_error("[error] unknown method: " + method,
            config, "here", {"expected one of the followings.",
                "- \"rigid\" : well-known rigid body collidion based method (default)",
                "- \"smooth\": smooth function approximation [T. Niina et al., JCTC (2020)]"
            }) << std::endl;

        return EXIT_FAILURE;
    }

    const auto sigma_x_v =
        config.contains("sigma_x") ? toml::find(config, "sigma_x") :
        config.contains("sigma")   ? toml::find(config, "sigma")   : nan_v;
    const auto sigma_y_v =
        config.contains("sigma_y") ? toml::find(config, "sigma_y") :
        config.contains("sigma")   ? toml::find(config, "sigma")   : nan_v;

    if(config.contains("sigma") && config.contains("sigma_x"))
    {
        std::cout << toml::format_error(
                "when \"sigma_x\" is defined, \"sigma\" will be ignored.",
                config.at("sigma_x"), "sigma for x is defined here",
                config.at("sigma"), "the default sigma will be ignored");
    }
    if(config.contains("sigma") && config.contains("sigma_y"))
    {
        std::cout << toml::format_error(
                "when \"sigma_y\" is defined, \"sigma\" will be ignored.",
                config.at("sigma_y"), "sigma for y is defined here",
                config.at("sigma"), "the default sigma will be ignored");
    }

    const auto sigma_x = afmize::read_as_angstrom<Real>(sigma_x_v);
    const auto sigma_y = afmize::read_as_angstrom<Real>(sigma_y_v);
    const auto gamma = afmize::read_as_angstrom<Real>(
        toml::find_or(config, "gamma", nan_v));

    afmize::progress_bar<70> bar(
        (reader->size() == 1) ? stg.x_pixel() * stg.y_pixel() : reader->size()
        );

    std::cerr << "done. creating image..." << std::endl;
    auto img = stg.create_image();
    try
    {
        std::size_t index = 0;
        while(!reader->is_eof())
        {
            if(reader->size() != 1)
            {
                std::cerr << bar.format(index);
            }

            const auto sys = [=](auto sys) -> afmize::system<Real> {
                if(stage_align) // align the lower edge of bounding box to stage
                {
                    assert(!std::isnan(stage_position));
                    const auto dz = sys.bounding_box.lower[2] - stage_position;
                    for(auto& p : sys.particles)
                    {
                        p.center[2] -= dz;
                    }
                }
                return sys;
            }(afmize::system<Real>(reader->read_snapshot(), stg));

            const Real bottom = std::isnan(stage_position) ?
                                sys.bounding_box.lower[2] : stage_position;

            const Real initial_z = sys.bounding_box.upper[2] + probe.radius;

            for(std::size_t j=0; j<img.y_pixel(); ++j)
            {
                for(std::size_t i=0; i<img.x_pixel(); ++i)
                {
                    probe.apex    = sys.stage_info.position_at(i, j);
                    probe.apex[2] = initial_z;

                    if(method == "rigid")
                    {
                        if(sys.stage_info.z_resolution() != 0)
                        {
                            img(i, j) = afmize::discretize(
                                afmize::collide_at(sys, probe, bottom),
                                sys.stage_info.z_resolution(),
                                bottom);
                        }
                        else
                        {
                            img(i, j) = afmize::collide_at(sys, probe, bottom);
                        }
                    }
                    else if(method == "smooth")
                    {
                        img(i, j) = afmize::smooth_at(sys, probe.apex, bottom,
                                                      gamma, sigma_x, sigma_y);
                    }
                    else
                    {
                        assert(false); // never reach here
                    }

                    if(reader->size() == 1)
                    {
                        std::cerr << bar.format(j * stg.x_pixel() + i);
                    }
                }
            }

            if(noise_sigma != 0.0)
            {
                afmize::apply_noise(img, rng, noise_sigma, stage_position);
            }

            if(reader->size() == 1)
            {
                std::cerr << bar.format(stg.y_pixel() * stg.x_pixel());
            }

            std::string outname(output_basename);
            if(reader->size() != 1)
            {
                std::ostringstream oss;
                oss << '_' << std::setfill('0')
                    << std::setw(std::to_string(reader->size()).size()) << index;
                outname += oss.str();
            }

            if(output_formats.csv)
            {
                afmize::write_csv (outname, img);
            }
            if(output_formats.tsv)
            {
                afmize::write_tsv (outname, img);
            }

            if(output_formats.json)
            {
                afmize::write_json(outname, img, sys);
            }

            if(std::isnan(cmap_min) || std::isnan(cmap_max))
            {
                if(output_formats.ppm)
                {
                    afmize::write_ppm (outname, img);
                }
                if(output_formats.svg)
                {
                    afmize::write_svg (outname, img, sys, scale_bar_length);
                }
            }
            else
            {
                if(output_formats.ppm)
                {
                    afmize::write_ppm(outname, img, std::make_pair(cmap_min, cmap_max));
                }
                if(output_formats.svg)
                {
                    afmize::write_svg(outname, img, sys, scale_bar_length,
                                      std::make_pair(cmap_min, cmap_max));
                }
            }
            ++index;
        }
    }
    catch(afmize::reader_base<Real>::no_more_model)
    {
        ; // do nothing
    }
    std::cerr << "done." << std::endl;
    return 0;
}
