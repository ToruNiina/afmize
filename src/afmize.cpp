#include <afmize/stage.hpp>
#include <afmize/system.hpp>
#include <afmize/colormap.hpp>
#include <afmize/xyz_reader.hpp>
#include <afmize/pdb_reader.hpp>
#include <extlib/toml/toml/toml.hpp>
#include <extlib/pnm/pnm.hpp>
#include <string>

namespace afmize
{

template<typename T>
decltype(toml::get<T>(std::declval<toml::value>()))
get(const toml::table& tab, const std::string& key, const std::string loc)
{
    // to show the error message
    try{return toml::get<T>(tab.at(key));}
    catch(std::out_of_range)
    {
        using namespace std::literals::string_literals;
        throw std::out_of_range(
                "key `"s + key + "` is not found in a table "s + loc);
    }
}

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
Real collide_at(const system<Real>& sys, const default_probe<Real>& probe)
{
    Real height = sys.bounding_box.lower[2];
    for(const auto& sph : sys.particles)
    {
        const Real dz_sph = collision_z(sphere<Real>{
                probe.radius, probe.apex
            }, sph);
        const Real dz_frs = collision_z(circular_frustum<Real>{
                probe.angle, probe.radius, probe.apex
            }, sph);

        height = std::max(height, probe.apex[2] + dz_frs);
        if(!std::isnan(dz_sph))
        {
            height = std::max(height, probe.apex[2] + dz_sph);
        }
    }
    return height;
}

template<typename Real>
void write_json(const stage<Real>& stg, const std::string& out)
{
    using namespace std::literals::string_literals;

    std::ofstream ofs(out + ".json");
    ofs << "{\"height\":[\n";
    for(std::size_t j=0; j<stg.y_pixel(); ++j)
    {
        for(std::size_t i=0; i<stg.x_pixel(); ++i)
        {
            const auto p = stg.position_at(i, j);
            ofs << "\t{\"x\":" << p[0] << ", \"y\":" << p[1] << ", \"z\":" << p[2] << "}";
            if(j == stg.y_pixel()-1 && i == stg.x_pixel()-1) {ofs << "\n";}
            else {ofs << ",\n";}
        }
    }
    ofs << "]}";
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

    // origin of image is upper left, but the physical origin is lower left.
    // to workaround this, it simply reverse the image.
    pnm::image<pnm::rgb_pixel> reversed(stg.x_pixel(), stg.y_pixel());
    std::copy(ppm.rbegin(), ppm.rend(), reversed.begin());

    pnm::write(out + ".ppm"s, reversed, pnm::format::binary);
    return;
}
} // afmize

int main(int argc, char** argv)
{
    using Real = double;
    constexpr Real pi             = Real(3.1415926535897932384626);
    constexpr Real deg_to_rad     = pi / Real(180.0);
    constexpr Real nm_to_angstrom = 10;

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

    std::cout << "reading config file..." << std::endl;

    const auto config = toml::parse(config_file);

    if(config.count("radii") == 1)
    {
        for(const auto& kv : afmize::get<toml::table>(config, "radii", "root"))
        {
            const auto atom_name   = kv.first;
            const auto atom_radius = toml::get<Real>(kv.second);
            afmize::parameter<Real>::radius[atom_name] = atom_radius;
        }
    }

    // read input output files
    const auto& file   = afmize::get<toml::table>(config, "file", "root");
    const auto  reader = afmize::open_file<Real>(
            afmize::get<std::string>(file, "input", "[file]")
        );
    const auto output = afmize::get<std::string>(file, "output", "[file]");

    // image size information
    const auto& resolution = afmize::get<toml::table>(config, "resolution", "root");
    const auto& range      = afmize::get<toml::table>(config, "range", "root");
    auto range_x = afmize::get<std::pair<Real, Real>>(range, "x", "[range]");
    auto range_y = afmize::get<std::pair<Real, Real>>(range, "y", "[range]");
    afmize::stage<Real> stg(
            afmize::get<Real>(resolution, "x", "[resolution]") * nm_to_angstrom,
            afmize::get<Real>(resolution, "y", "[resolution]") * nm_to_angstrom,
            afmize::get<Real>(resolution, "x", "[resolution]") * nm_to_angstrom,
            std::make_pair(range_x.first * nm_to_angstrom, range_x.second * nm_to_angstrom),
            std::make_pair(range_y.first * nm_to_angstrom, range_y.second * nm_to_angstrom)
        );

    // probe size information
    const auto& probe_tab  = afmize::get<toml::table>(config, "probe", "root");
    const auto& probe_size = afmize::get<toml::table>(probe_tab, "size", "[probe]");
    afmize::default_probe<Real> probe{
            afmize::get<Real>(probe_size, "angle",  "[probe.size]") * deg_to_rad,
            afmize::get<Real>(probe_size, "radius", "[probe.size]") * nm_to_angstrom,
            mave::vector<Real, 3>{0, 0, 0}
        };

    std::cout << "done." << std::endl;

    try
    {
        while(!reader->is_eof())
        {
            const afmize::system<Real> sys(reader->read_snapshot());
            const Real initial_z = sys.bounding_box.upper[2] + probe.radius;

            for(std::size_t j=0; j<stg.y_pixel(); ++j)
            {
                for(std::size_t i=0; i<stg.x_pixel(); ++i)
                {
                    probe.apex    = stg.position_at(i, j);
                    probe.apex[2] = initial_z;

                    stg(i, j) = afmize::descretize(
                        afmize::collide_at(sys, probe),
                        stg.z_resolution(),
                        sys.bounding_box.lower[2]);
                }
            }

            afmize::write_ppm (stg, output);
            afmize::write_csv (stg, output);
            afmize::write_json(stg, output);
        }
    }
    catch(afmize::reader_base<Real>::no_more_model)
    {
        ; // do nothing
    }
    return 0;
}