#include "src/common.hpp"
#include "src/atom_bindings.hpp"
#include "src/sasa_calculators.hpp"
#include "src/utils.hpp"

PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = R"pbdoc(
        DR_SASA Python Bindings
        ----------------------
        
        This module provides Python bindings for the DR_SASA library,
        a fast and accurate solvent accessible surface area calculator.
    )pbdoc";
    
    py::register_exception<SASAError>(m, "SASAError");  // Register custom exceptions
    bind_atom_struct(m); // Bind core data structures
    
    py::class_<SimpleSASA>(m, "SimpleSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = DEFAULT_PROBE_RADIUS,
             py::arg("compute_mode") = 0)
        .def("calculate", &SimpleSASA::calculate,
             py::arg("pdb_file"),
             py::arg("print_output") = false,
             py::arg("output_name") = "simple",
             R"pbdoc(
                Calculate SASA from PDB file.

                Args:
                    pdb_file (str): Path to PDB file
                    print_output (bool): Whether to generate printed output
                    output_name (str): Base name for output files

                Returns:
                    dict: Results including SASA values and printed output if requested
             )pbdoc")
        .def("calculate_from_atoms", &SimpleSASA::calculate_from_atoms,
             py::arg("atoms"),
             py::arg("print_output") = false,
             py::arg("output_name") = "simple",
             "Calculate SASA from list of atom_struct objects");

    py::class_<GenericSASA>(m, "GenericSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = DEFAULT_PROBE_RADIUS,
             py::arg("compute_mode") = 0)
        .def("calculate", &GenericSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains") = std::vector<std::vector<std::string>>(),
             py::arg("include_matrix") = true,
             py::arg("print_output") = false,
             py::arg("output_name") = "generic",
             R"pbdoc(
                Calculate SASA with chain-based analysis.

                Args:
                    pdb_file (str): Path to PDB file
                    chains (List[List[str]]): Chain selections
                    include_matrix (bool): Whether to include matrix analysis
                    print_output (bool): Whether to generate printed output
                    output_name (str): Base name for output files

                Returns:
                    dict: Results including SASA values and printed output if requested
             )pbdoc")
        .def("calculate_from_atoms", &GenericSASA::calculate_from_atoms,
             py::arg("atoms"),
             py::arg("chains") = std::vector<std::vector<std::string>>(),
             py::arg("include_matrix") = true,
             py::arg("print_output") = false,
             py::arg("output_name") = "generic",
             "Calculate SASA from list of atom_struct objects");

    py::class_<DecoupledSASA>(m, "DecoupledSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = DEFAULT_PROBE_RADIUS,
             py::arg("compute_mode") = 0)
        .def("calculate", &DecoupledSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains") = std::vector<std::vector<std::string>>(),
             py::arg("include_matrix") = true,
             py::arg("print_output") = false,
             py::arg("output_name") = "decoupled",
             R"pbdoc(
                Calculate decoupled SASA analysis.

                Args:
                    pdb_file (str): Path to PDB file
                    chains (List[List[str]]): Chain selections
                    include_matrix (bool): Whether to include matrix analysis
                    print_output (bool): Whether to generate printed output
                    output_name (str): Base name for output files

                Returns:
                    dict: Results including SASA values and printed output if requested
             )pbdoc")
        .def("calculate_from_atoms", &DecoupledSASA::calculate_from_atoms,
             py::arg("atoms"),
             py::arg("chains") = std::vector<std::vector<std::string>>(),
             py::arg("include_matrix") = true,
             py::arg("print_output") = false,
             py::arg("output_name") = "decoupled",
             "Calculate SASA from list of atom_struct objects");

    // Module attributes
    m.attr("__version__") = "0.5.0";
    m.attr("__author__") = "Original: Ribeiro J., Ríos-Vera C., Melo F., Schüller A.";
    m.attr("DEFAULT_PROBE_RADIUS") = DEFAULT_PROBE_RADIUS;
}