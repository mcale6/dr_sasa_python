#include "common.hpp"
#include "atom_bindings.hpp"
#include "exceptions.hpp"
#include "constants.hpp"
#include "utils.hpp"

namespace py = pybind11;

// Base class for all SASA calculators
class SASACalculator {
protected:
    VDWcontainer vdw_radii_;
    float probe_radius_;
    int compute_mode_;

public:
    SASACalculator(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0)
        : vdw_radii_(""), probe_radius_(probe_radius), compute_mode_(compute_mode) {
        vdw_radii_.GenPoints();
    }
};

// Simple SASA Calculator
class SimpleSASA : public SASACalculator {
public:
    using SASACalculator::SASACalculator;

    py::dict calculate(const std::string& pdb_file,
                      bool print_output = false,
                      const std::string& output_name = "output");

    py::dict calculate_from_atoms(std::vector<atom_struct> atoms,
                                bool print_output = false,
                                const std::string& output_name = "output");
};

// Generic SASA Calculator
class GenericSASA : public SASACalculator {
public:
    using SASACalculator::SASACalculator;

    py::dict calculate(const std::string& pdb_file,
                      std::vector<std::vector<std::string>>& chains,
                      bool include_matrix = true,
                      bool print_output = false,
                      const std::string& output_name = "output");

    py::dict calculate_from_atoms(std::vector<atom_struct> atoms,
                                std::vector<std::vector<std::string>>& chains,
                                bool include_matrix = true,
                                bool print_output = false,
                                const std::string& output_name = "output");
};

// Decoupled SASA Calculator
class DecoupledSASA : public SASACalculator {
public:
    using SASACalculator::SASACalculator;

    py::dict calculate(const std::string& pdb_file,
                      std::vector<std::vector<std::string>>& chains,
                      bool include_matrix = true,
                      bool print_output = false,
                      const std::string& output_name = "output");

    py::dict calculate_from_atoms(std::vector<atom_struct> atoms,
                                std::vector<std::vector<std::string>>& chains,
                                bool include_matrix = true,
                                bool print_output = false,
                                const std::string& output_name = "output");
};
