#include <afmize/stage.hpp>
#include <afmize/system.hpp>
#include <afmize/colormap.hpp>
#include <afmize/xyz_reader.hpp>
#include <afmize/pdb_reader.hpp>
#include <afmize/input_utility.hpp>
#include <afmize/progress_bar.hpp>
#include <extlib/pnm/pnm.hpp>
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

template<typename Real>
Real descretize(const Real x, const Real resolution, const Real min_val)
{
    return std::round((x - min_val) / resolution) * resolution + min_val;
}

template<typename Real>
Real collide_at(const system<Real>& sys, const default_probe<Real>& probe,
                const Real bottom)
{
    // the point at which the probe collides with the stage
    Real height = bottom + probe.radius;
    for(const auto& sph : sys.particles)
    {
        const Real dz_sph = collision_z(sphere<Real>{
                probe.radius, probe.apex
            }, sph);
        const Real dz_frs = collision_z(circular_frustum<Real>{
                probe.angle, probe.radius, probe.apex
            }, sph);

        if(!std::isnan(dz_frs))
        {
            height = std::max(height, probe.apex[2] + dz_frs);
        }
        if(!std::isnan(dz_sph))
        {
            height = std::max(height, probe.apex[2] + dz_sph);
        }
    }
    // subtract probe radius to set the ground as 0
    return height - probe.radius;
}

template<typename Real>
void write_json(const stage<Real>& stg, const std::string& out)
{
    using namespace std::literals::string_literals;

    std::ofstream ofs(out + ".json");
    ofs << "{\n\t\"resolution\":{\"x\":" << stg.x_resolution()
        << ", \"y\":" << stg.y_resolution()
        << ", \"z\":" << stg.z_resolution() << "},\n";
    ofs << "\t\"height\":[\n";
    for(std::size_t j=0; j<stg.y_pixel(); ++j)
    {
        for(std::size_t i=0; i<stg.x_pixel(); ++i)
        {
            const auto p = stg.position_at(i, j);
            ofs << "\t\t{\"x\":" << p[0] << ", \"y\":" << p[1] << ", \"z\":" << p[2] << "}";
            if(j == stg.y_pixel()-1 && i == stg.x_pixel()-1) {ofs << "\n";}
            else {ofs << ",\n";}
        }
    }
    ofs << "\t]\n}";
    return;
}

template<typename Real>
void write_csv(const stage<Real>& stg, const std::string& out)
{
    using namespace std::literals::string_literals;

    std::ofstream ofs(out + ".csv");
    for(std::size_t j=0; j<stg.y_pixel(); ++j)
    {
        for(std::size_t i=0; i<stg.x_pixel(); ++i)
        {
            ofs << stg(i, j) << ' ';
        }
        ofs << '\n';
    }
    return;
}

template<typename Real>
void write_ppm(const stage<Real>& stg, const std::string& out)
{
    using namespace std::literals::string_literals;
    const auto minmax = std::minmax_element(stg.begin(), stg.end());
    const auto min_elem = *minmax.first;
    const auto max_elem = *minmax.second;

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
    pnm::write(out + ".ppm"s, reversed, pnm::format::binary);
    return;
}
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

    if(config.count("radii") == 1)
    {
        const auto& radii = toml::get<toml::table>(config.at("radii"));
        if(radii.count("atom") == 1)
        {
            for(const auto& kv : toml::get<toml::table>(radii.at("atom")))
            {
                afmize::parameter<Real>::radius_atom[kv.first] =
                    afmize::read_as_angstrom<Real>(kv.second, "radii.atom");
            }
        }
        if(radii.count("residue") == 1)
        {
            for(const auto& res : toml::get<toml::table>(radii.at("residue")))
            {
                for(const auto& atm : toml::get<toml::table>(res.second))
                {
                    afmize::parameter<Real>::radius_residue[res.first][atm.first] =
                        afmize::read_as_angstrom<Real>(atm.second,
                                "radii."s + res.first + "."s + atm.first);
                }
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
    const auto& file   = afmize::get<toml::table>(config, "file", "root");
    const auto  reader = afmize::open_file<Real>(
            afmize::get<std::string>(file, "input", "[file]")
        );
    const auto output = afmize::get<std::string>(file, "output", "[file]");

    std::cerr << "-- " << reader->size() << " snapshots are found\n";

    // image size information
    const auto& resolution = afmize::get<toml::table>(config, "resolution", "root");
    const auto& range      = afmize::get<toml::table>(config, "range", "root");
    auto range_x = afmize::get<std::pair<toml::value, toml::value>>(range, "x", "[range]");
    auto range_y = afmize::get<std::pair<toml::value, toml::value>>(range, "y", "[range]");
    afmize::stage<Real> stg(
            afmize::read_as_angstrom<Real>(
                afmize::find(resolution, "x", "[resolution]"), "resolution.x"),
            afmize::read_as_angstrom<Real>(
                afmize::find(resolution, "y", "[resolution]"), "resolution.y"),
            afmize::read_as_angstrom<Real>(
                afmize::find(resolution, "z", "[resolution]"), "resolution.z"),
            std::make_pair(
                afmize::read_as_angstrom<Real>(range_x.first,  "range_x"),
                afmize::read_as_angstrom<Real>(range_x.second, "range_x")),
            std::make_pair(
                afmize::read_as_angstrom<Real>(range_y.first,  "range_y"),
                afmize::read_as_angstrom<Real>(range_y.second, "range_y"))
        );

    // probe size information
    const auto& probe_tab  = afmize::get<toml::table>(config, "probe", "root");
    const auto& probe_size = afmize::get<toml::table>(probe_tab, "size", "[probe]");
    afmize::default_probe<Real> probe{
            afmize::get<Real>(probe_size, "angle",  "[probe.size]") * deg_to_rad,
            afmize::read_as_angstrom<Real>(
                afmize::find(probe_size, "radius", "[probe.size]"),
                "probe.size.radius"),
            mave::vector<Real, 3>{0, 0, 0}
        };

    // stage information ...
    const bool stage_align    = [](const toml::table& root) -> bool {
        if(root.count("stage")  == 0) {return false;}
        const auto& stage = afmize::get<toml::table>(root, "stage", "root");
        if(stage.count("align") == 0) {return false;}
        return toml::get<bool>(stage.at("align"));
    }(config);

    const Real stage_position = [](const toml::table& root) -> Real {
        if(root.count("stage")  == 0)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        const auto& stage = afmize::get<toml::table>(root, "stage", "root");
        if(stage.count("position") == 0)
        {
            return std::numeric_limits<Real>::quiet_NaN();
        }
        return afmize::read_as_angstrom<Real>(stage.at("position"), "stage.position");
    }(config);

    if(stage_align && std::isnan(stage_position))
    {
        std::cerr << "stage.align == true but stage.position is not defined!" << std::endl;
        std::cerr << "In order to align system on the stage, you need to input\n"
                     " the z-coordinate position of your stage." << std::endl;
        return EXIT_FAILURE;
    }

    afmize::progress_bar<70> bar(
        (reader->size() == 1) ? stg.x_pixel() * stg.y_pixel() : reader->size()
        );

    std::cerr << "done. creating image..." << std::endl;
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
            }(afmize::system<Real>(reader->read_snapshot()));

            const Real bottom = std::isnan(stage_position) ?
                                sys.bounding_box.lower[2] : stage_position;

            const Real initial_z = sys.bounding_box.upper[2] + probe.radius;

            for(std::size_t j=0; j<stg.y_pixel(); ++j)
            {
                for(std::size_t i=0; i<stg.x_pixel(); ++i)
                {
                    probe.apex    = stg.position_at(i, j);
                    probe.apex[2] = initial_z;

                    stg(i, j) = afmize::descretize(
                        afmize::collide_at(sys, probe, bottom),
                        stg.z_resolution(),
                        bottom);

                    if(reader->size() == 1)
                    {
                        std::cerr << bar.format(j * stg.x_pixel() + i);
                    }
                }
            }
            if(reader->size() == 1)
            {
                std::cerr << bar.format(stg.y_pixel() * stg.x_pixel());
            }

            std::string outname(output);
            if(reader->size() != 1)
            {
                std::ostringstream oss;
                oss << '_' << std::setfill('0')
                    << std::setw(std::to_string(reader->size()).size()) << index;
                outname += oss.str();
            }

            afmize::write_ppm (stg, output);
            afmize::write_csv (stg, output);
            afmize::write_json(stg, output);
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
