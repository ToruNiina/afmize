#ifndef AFMIZE_READ_XYZ_HPP
#define AFMIZE_READ_XYZ_HPP
#include <afmize/parameter.hpp>
#include <afmize/read_number.hpp>
#include <afmize/reader_base.hpp>
#include <sstream>

namespace afmize
{

template<typename realT>
class xyz_reader final : public reader_base<realT>
{
  public:

    using base_type       = reader_base<realT>;
    using snapshot_type   = typename base_type::snapshot_type;
    using trajectory_type = typename base_type::trajectory_type;

    xyz_reader(const std::string& fname)
        : base_type(fname), ln(0), xyz(fname)
    {
        if(!xyz.good()) {throw std::runtime_error("file open error: " + fname);}
    }
    ~xyz_reader() override = default;

    bool is_eof() override {xyz.peek(); return xyz.eof();}

    trajectory_type read_trajectory() override
    {
        using namespace std::literals::string_literals;

        this->xyz.seekg(0, std::ios_base::beg);
        this->ln = 0;

        trajectory_type traj;
        while(this->is_eof())
        {
            try
            {
                traj.push_back(this->read_snapshot());
                xyz.peek();
            }
            catch(typename base_type::no_more_model)
            {
                break;
            }
        }
        return traj;
    }

    snapshot_type read_snapshot() override
    {
        if(this->is_eof())
        {
            throw typename base_type::no_more_model{};
        }

        const auto n_elem = [this] {
            std::string l;
            std::getline(this->xyz, l);
            this->ln += 1;
            return l;
        }();
        std::size_t N = 0;
        try
        {
            N = std::stoull(n_elem);
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format in number of element at line "
                      << ln << ".\n";
            std::cerr << "> " << n_elem << '\n';
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "too many atoms in one snapshot at line"
                      << ln << ".\n";
            std::cerr << "> " << n_elem << '\n';
            std::exit(EXIT_FAILURE);
        }

        {
            std::string comment;
            std::getline(this->xyz, comment);
        }

        snapshot_type snapshot;
        snapshot.reserve(N);
        while(this->is_eof() && snapshot.size() < N)
        {
            std::string line;
            std::getline(this->xyz, line);
            snapshot.push_back(this->read_atom(line));
        }
        if(snapshot.size() != N)
        {
            throw std::runtime_error("a snapshot in the xyz file (" +
                    this->filename + std::string(") lacks particle"));
        }
        return snapshot;
    }

  private:

    sphere<realT> read_atom(const std::string& line)
    {
        std::istringstream iss(line);
        std::string name, x, y, z;
        iss >> name >> x >> y >> z;

        sphere<realT> particle;

        try
        {
            particle.radius = parameter<realT>::radius.at(name);
        }
        catch(std::out_of_range)
        {
            std::cerr << "unknown atom name found at line " << ln << ".\n";
            this->highlight_columns(std::cerr, line, 0, name.size());
            std::exit(EXIT_FAILURE);
        }

        try
        {
            particle.center[0] = read_number<realT>(x);
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            std::cerr << "> " << line << '\n';
            std::cerr << base_type::mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            std::cerr << "> " << line << '\n';
            std::cerr << base_type::mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        try
        {
            particle.center[1] = read_number<realT>(x);
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            std::cerr << "> " << line << '\n';
            std::cerr << base_type::mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            std::cerr << "> " << line << '\n';
            std::cerr << base_type::mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        try
        {
            particle.center[2] = read_number<realT>(x);
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            std::cerr << "> " << line << '\n';
            std::cerr << base_type::mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            std::cerr << "> " << line << '\n';
            std::cerr << base_type::mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        return particle;
    }

  private:

    std::size_t   ln;
    std::ifstream xyz;
};

} // afmize
#endif // AFMIZE_READ_xyz_HPP
