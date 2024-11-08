#include "src/common.hpp"
#include "src/atom_bindings.hpp"
#include "src/simple_sasa.hpp"
#include "src/generic_sasa.hpp"
#include "src/decoupled_sasa.hpp"
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
             py::arg("pdb_file"));

    py::class_<GenericSASA>(m, "GenericSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = DEFAULT_PROBE_RADIUS,
             py::arg("compute_mode") = 0)
        .def("calculate", &GenericSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains") = std::vector<std::vector<std::string>>(),
             py::arg("include_matrix") = true);

    py::class_<DecoupledSASA>(m, "DecoupledSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = DEFAULT_PROBE_RADIUS,
             py::arg("compute_mode") = 0)
        .def("calculate", &DecoupledSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains") = std::vector<std::vector<std::string>>(),
             py::arg("include_matrix") = true);

    // Convenience functions
    m.def("calculate_simple_sasa", 
        [](const std::string& pdb_file, float probe_radius = DEFAULT_PROBE_RADIUS) {
            SimpleSASA calculator(probe_radius);
            return calculator.calculate(pdb_file);
        }, 
        py::arg("pdb_file"), 
        py::arg("probe_radius") = DEFAULT_PROBE_RADIUS);

    m.def("calculate_delta_sasa", 
        [](const std::string& pdb_file, 
           std::vector<std::vector<std::string>>& chains,  
           float probe_radius = DEFAULT_PROBE_RADIUS,
           bool include_matrix = true) {
            GenericSASA calculator(probe_radius);
            return calculator.calculate(pdb_file, chains, include_matrix);
        }, 
        py::arg("pdb_file"), 
        py::arg("chains"), 
        py::arg("probe_radius") = DEFAULT_PROBE_RADIUS,
        py::arg("include_matrix") = true);

    m.attr("__version__") = "0.5.0";
    m.attr("__author__") = "Original: Ribeiro J., Ríos-Vera C., Melo F., Schüller A.";
    m.attr("DEFAULT_PROBE_RADIUS") = DEFAULT_PROBE_RADIUS;
}