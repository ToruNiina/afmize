#ifndef AFMIZE_READ_PDB_HPP
#define AFMIZE_READ_PDB_HPP
#include <afmize/parameter.hpp>
#include <afmize/read_number.hpp>
#include <afmize/shapes.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

namespace afmize
{

template<typename realT>
class pdb_reader
{
  public:

    using snapshot = std::vector<sphere<realT>>;

    pdb_reader(const std::string& fname) : model_found(false), ln(0), pdb(fname)
    {
        if(!pdb.good()) {throw std::runtime_error("file open error: " + fname);}

        while(!pdb.eof())
        {
            std::string line;
            std::getline(this->pdb, line);
            try
            {
                if(line.substr(0, 5) == "MODEL")
                {
                    this->model_found = true;
                    break;
                }
            }
            catch(...)
            {
                continue;
            }
            this->pdb.peek();
        }
        this->pdb.seekg(0, std::ios_base::beg);
    }
    ~pdb_reader() = default;

    bool is_eof() {pdb.peek(); return pdb.eof();}

    std::vector<snapshot> read_all()
    {
        this->pdb.seekg(0, std::ios_base::beg);
        this->ln = 0;

        std::vector<snapshot> models;
        while(!pdb.eof())
        {
            try
            {
                models.push_back(this->read_model());
                pdb.peek();
            }
            catch(std::runtime_error re)
            {
                if(re.what() != "no more MODEL in PDB file")
                {
                    throw re;
                }
                break;
            }
        }
        return models;
    }

    snapshot read_model()
    {
        if(model_found)
        {
            while(!pdb.eof())
            {
                const auto line = [this] {
                    std::string l;
                    std::getline(this->pdb, l);
                    this->ln += 1;
                    return l;
                }();

                const auto header = get_substr(line, 0, 5);
                if(header == "MODEL")
                {
                    pdb.peek();
                    break;
                }
            }
            if(pdb.eof())
            {
                throw std::runtime_error("no more MODEL in PDB file");
            }
        }

        snapshot snap;
        while(!pdb.eof())
        {
            const auto line = [this] {
                std::string l;
                std::getline(this->pdb, l);
                this->ln += 1;
                return l;
            }();

            const auto header = get_substr(line, 0, 6);
            if(model_found && (header == "ENDMDL" || header == "MODEL "))
            {
                break;
            }
            if(header != "ATOM  " && header != "HETATM")
            {
                continue;
            }
            snap.push_back(read_atom(line));
            this->pdb.peek();
        }

        pdb.peek();
        return snap;
    }

  private:

    sphere<realT> read_atom(const std::string& line)
    {
        const auto header = line.substr(0, 6);
        if(header != "ATOM  " && header != "HETATM")
        {
            throw std::invalid_argument("internal error: invalid line");
        }

        sphere<realT> particle;

        // according to wwPDB 3.3
        //> Alignment of one-letter atom name such as C starts at column 14,
        //> while two-letter atom name such as FE starts at column 13.
        try
        {
            if(!std::isupper(line.at(12)))
            {
                particle.radius =
                    parameter<realT>::radius.at(this->get_substr(line, 13, 1));
            }
            else
            {
                particle.radius =
                    parameter<realT>::radius.at(this->get_substr(line, 12, 2));
            }
        }
        catch(std::out_of_range)
        {
            std::cerr << "unknown atom name found at line " << ln << ".\n";
            highlight_columns(std::cerr, line, 12, 4);
            std::cerr << "according to wwPDB 3.3,\n> one-letter atom name "
                         "such as C starts at column 14,\n> while two-letter "
                         "atom name such as FE starts at column 13.\n";
            std::cerr << "see also:\nhttp://www.wwpdb.org/documentation/"
                         "file-format-content/format33/sect9.html#ATOM\n";
            if(line.at(12) == 'H')
            {
                std::cerr << "WARNING: It looks like a hydrogen. "
                          << "continue reading...";
                particle.radius = parameter<realT>::radius.at("H");
                std::cerr << "\n\n";
            }
            else
            {
                std::exit(EXIT_FAILURE);
            }
        }

        try
        {
            particle.center[0] =
                read_number<realT>(this->get_substr(line, 30, 8));
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            highlight_columns(std::cerr, line, 30, 8);
            std::cerr << mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            highlight_columns(std::cerr, line, 30, 8);
            std::cerr << mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        try
        {
            particle.center[1] =
                read_number<realT>(this->get_substr(line, 38, 8));
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            highlight_columns(std::cerr, line, 38, 8);
            std::cerr << mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            highlight_columns(std::cerr, line, 38, 8);
            std::cerr << mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        try
        {
            particle.center[2] =
                read_number<realT>(this->get_substr(line, 46, 8));
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            highlight_columns(std::cerr, line, 46, 8);
            std::cerr << mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            highlight_columns(std::cerr, line, 46, 8);
            std::cerr << mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        return particle;
    }

    void highlight_columns(std::ostream& os, const std::string& line,
                            const std::size_t padding, const std::size_t width)
    {
        os << "> " << line << '\n';
        os << "  ";
        for(std::size_t i=0; i<padding; ++i) {os << ' ';}
        for(std::size_t i=0; i<width;   ++i) {os << '^';}
        os << '\n';
    }

    std::string get_substr(const std::string& line,
        std::string::size_type b, std::string::size_type w) noexcept
    {
        try
        {
            return line.substr(b, w);
        }
        catch(std::out_of_range)
        {
            std::cerr << "incomplete ATOM line found at line " << ln << ".\n";
            std::cerr << '\t' << line << '\n';
            std::cerr << '\t' << std::string('^', line.size()) << '\n';
            std::cerr << "check your pdb file.\n";
            std::exit(EXIT_FAILURE);
        }
        catch(std::bad_alloc)
        {
            std::cerr << "memory error detected.\nyou need to buy more RAM.\n";
            std::exit(EXIT_FAILURE);
        }
        catch(...)
        {
            std::cerr << "unexpected error detected while reading pdb file.\n";
            std::exit(EXIT_FAILURE);
        }
    }

  private:

    static constexpr auto mes_float_format_err =
        "couldn't convert to floating-point value.\n"
        "If it looks ok, locale settings may have a problem. "
        "check your locale.\n";

    static constexpr auto mes_float_range_err =
        "the value exceeds the range of floating-point.\n";

    bool          model_found;
    std::size_t   ln; // line number
    std::ifstream pdb;
};

} // afmize
#endif // AFMIZE_READ_PDB_HPP
