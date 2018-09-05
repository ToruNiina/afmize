#ifndef AFMIZE_READ_PDB_HPP
#define AFMIZE_READ_PDB_HPP
#include <afmize/parameter.hpp>
#include <afmize/read_number.hpp>
#include <afmize/reader_base.hpp>

namespace afmize
{

template<typename realT>
class pdb_reader final : public reader_base<realT>
{
  public:

    using base_type       = reader_base<realT>;
    using snapshot_type   = typename base_type::snapshot_type;
    using trajectory_type = typename base_type::trajectory_type;

    pdb_reader(const std::string& fname, const bool read_hetatms = false)
        : base_type(fname), model_found(false), read_HETATMs(read_hetatms),
          ln(0), pdb(fname)
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
    ~pdb_reader() override = default;

    bool is_eof() override {pdb.peek(); return pdb.eof();}

    trajectory_type read_trajectory() override
    {
        this->pdb.seekg(0, std::ios_base::beg);
        this->ln = 0;

        std::vector<snapshot_type> models;
        while(!pdb.eof())
        {
            try
            {
                models.push_back(this->read_snapshot());
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

    snapshot_type read_snapshot() override
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

                const auto header = this->get_substr(line, ln, 0, 5);
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

        snapshot_type snap;
        while(!pdb.eof())
        {
            const auto line = [this] {
                std::string l;
                std::getline(this->pdb, l);
                this->ln += 1;
                return l;
            }();

            const auto header = this->get_substr(line, ln, 0, 6);
            if(model_found && (header == "ENDMDL" || header == "MODEL "))
            {
                break;
            }
            if(header == "ATOM  " || (this->read_HETATMs && header == "HETATM"))
            {
                snap.push_back(read_atom(line));
            }
            this->pdb.peek();
        }

        pdb.peek();
        return snap;
    }

  private:

    sphere<realT> read_atom(const std::string& line)
    {
        const auto header = line.substr(0, 6);
        if(!(header == "ATOM  " || (this->read_HETATMs && header == "HETATM")))
        {
            throw std::invalid_argument("internal error: invalid line: " + line);
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
                    parameter<realT>::radius.at(this->get_substr(line, ln, 13, 1));
            }
            else
            {
                particle.radius =
                    parameter<realT>::radius.at(this->get_substr(line, ln, 12, 2));
            }
        }
        catch(std::out_of_range)
        {
            std::cerr << "unknown atom name found at line " << ln << ".\n";
            this->highlight_columns(std::cerr, line, 12, 4);
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
                read_number<realT>(this->get_substr(line, ln, 30, 8));
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            this->highlight_columns(std::cerr, line, 30, 8);
            std::cerr << base_type::mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            this->highlight_columns(std::cerr, line, 30, 8);
            std::cerr << base_type::mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        try
        {
            particle.center[1] =
                read_number<realT>(this->get_substr(line, ln, 38, 8));
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            this->highlight_columns(std::cerr, line, 38, 8);
            std::cerr << base_type::mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            this->highlight_columns(std::cerr, line, 38, 8);
            std::cerr << base_type::mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        try
        {
            particle.center[2] =
                read_number<realT>(this->get_substr(line, ln, 46, 8));
        }
        catch(std::invalid_argument)
        {
            std::cerr << "invalid format at line " << ln << ".\n";
            this->highlight_columns(std::cerr, line, 46, 8);
            std::cerr << base_type::mes_float_format_err;
            std::exit(EXIT_FAILURE);
        }
        catch(std::out_of_range)
        {
            std::cerr << "invalid value at line" << ln << ".\n";
            this->highlight_columns(std::cerr, line, 46, 8);
            std::cerr << base_type::mes_float_range_err;
            std::exit(EXIT_FAILURE);
        }

        return particle;
    }

  private:


    bool          model_found;
    bool          read_HETATMs;
    std::size_t   ln; // line number
    std::ifstream pdb;
};

} // afmize
#endif // AFMIZE_READ_PDB_HPP
