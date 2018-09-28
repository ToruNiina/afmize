#ifndef AFMIZE_PARAMETER_HPP
#define AFMIZE_PARAMETER_HPP
#include <string>
#include <map>

namespace afmize
{

template<typename realT>
struct parameter
{
    static std::map<std::string, realT> radius_atom;
    static std::map<std::string, std::map<std::string, realT>> radius_residue;
};

// Van der Waals radii. unit is Angstrom. overwritable from config file.
// references:
// A. Bondi (1964). "van der Waals Volumes and Radii".
//     J. Phys. Chem. 68: 441. doi:10.1021/j100785a001.
// M. Mantina; A.C. Chamberlin; R. Valero; C.J. Cramer; D.G. Truhlar (2009).
//     "Consistent van der Waals Radii for the Whole Main Group".
//     J. Phys. Chem. A. 113 (19): 5806–12. doi:10.1021/jp8111556.
template<typename realT>
std::map<std::string, realT> parameter<realT>::radius_atom = {
    { "H", 1.20},
    {"HE", 1.40},
    { "B", 1.92},
    { "C", 1.70},
    { "N", 1.55},
    { "O", 1.52},
    { "F", 1.47},
    {"NE", 1.54},
    {"NA", 2.27},
    {"MG", 1.73},
    {"AL", 1.84},
    {"SI", 2.10},
    { "P", 1.80},
    { "S", 1.80},
    {"CL", 1.75},
    {"AR", 1.88},
    { "K", 2.75},
    {"CA", 2.31},
};

template<typename realT>
std::map<std::string /* residue */, std::map<std::string /* atom */, realT>>
parameter<realT>::radius_residue;

} // afmize
#endif// AFMIZE_PARAMETER_HPP
