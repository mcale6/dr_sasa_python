#pragma once
#include "common.hpp"

// Core constants
inline constexpr float DEFAULT_PROBE_RADIUS = 1.4f;

// Standard SASA reference values
inline const std::map<std::string, float> STANDARD_SASA_VALUES = {
    {"ALA", 92.40211452},
    {"GLN", 184.71688},
    {"TRP", 235.3483229},
    {"SER", 122.353095},
    {"ARG", 219.941475},
    {"ILE", 158.07277},
    {"ASN", 146.06073},
    {"ASP", 154.71124},
    {"HIS", 194.617376},
    {"MET", 191.547244},
    {"LYS", 201.689792},
    {"LEU", 169.04496452},
    {"THR", 145.526463},
    {"PHE", 201.7065},
    {"TYR", 209.07715},
    {"GLU", 160.482561},
    {"CYS", 130.28853},
    {"PRO", 126.877028},
    {"GLY", 80.23533},
    {"VAL", 138.90233},
    {"A", 163.5912717},
    {"C", 162.86302775},
    {"T", 173.03036075},
    {"G", 164.253836},
    {"DA", 163.5912717},
    {"DC", 162.86302775},
    {"DT", 173.03036075},
    {"DG", 164.253836}
};

// Classification constants
inline const std::vector<std::string> BACKBONE_ATOMS = {
    "C", "CA", "N", "O", "P", "OP1", "OP2", "O3'", "O5'", "C3'", "C4'", "C5'"
};

inline const std::vector<std::string> POLAR_BACKBONE = {
    "N", "O", "OP1", "OP2", "O3'", "O5'"
};

inline const std::vector<std::string> HYD_BACKBONE = {
    "C", "CA", "C3'", "C4'", "C5'"
};

inline const std::vector<std::string> POLAR_SIDECHAIN = {
    "NZ", "OG", "OG1", "ND1", "NE2", "OE1", "OE2", "NE", "NH1", "NH2",
    "OD1", "OD2", "OH", "SE", "SG", "N1", "N2", "N3", "N4", "N6", "N7", "N9",
    "O2", "O4", "O6"
};

inline const std::vector<std::string> HYD_SIDECHAIN = {
    "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3",
    "CZ", "CZ2", "CZ3", "CH2", "C2", "C4", "C5", "C6", "C8"
};

inline const std::map<std::string, std::vector<std::string>> MAJOR_GROOVE = {
    {"A", {"N6", "N7", "C8", "C5"}},
    {"G", {"O6", "N7", "C8", "C5"}},
    {"C", {"N4", "C5", "C6"}},
    {"T", {"O4", "C5", "C6", "C7"}},
    {"DA", {"N6", "N7", "C8", "C5"}},
    {"DG", {"O6", "N7", "C8", "C5"}},
    {"DC", {"N4", "C5", "C6"}},
    {"DT", {"O4", "C5", "C6", "C7"}}
};

inline const std::map<std::string, std::vector<std::string>> MINOR_GROOVE = {
    {"A", {"N1", "C2", "N3"}},
    {"G", {"N2", "N3", "C2"}},
    {"C", {"O2", "N3", "C2"}},
    {"T", {"O2", "N3", "C2"}},
    {"DA", {"N1", "C2", "N3"}},
    {"DG", {"N2", "N3", "C2"}},
    {"DC", {"O2", "N3", "C2"}},
    {"DT", {"O2", "N3", "C2"}}
};
