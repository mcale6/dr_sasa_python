#include "atom_bindings.hpp"
#include "utils.hpp"

namespace py = pybind11;

void bind_atom_struct(py::module& m) {
    py::class_<atom_struct>(m, "AtomStruct")
        // Constructor with default values
        .def(py::init<uint32_t, int, std::string, std::string, std::string, 
                     std::string, std::string, std::string, std::string,
                     float, float, float, std::string, float, float, std::string>(),
            py::arg("id") = 0,                    // Default ID
            py::arg("resi") = 0,                  // Default residue number
            py::arg("icode") = "",                // Default insertion code
            py::arg("name") = "",                 // Default atom name
            py::arg("resn") = "",                 // Default residue name
            py::arg("chain") = "",                // Default chain ID
            py::arg("element") = "",              // Default element
            py::arg("structure") = "",            // Default structure
            py::arg("mol_type") = "ATOM",         // Default molecule type
            py::arg("x") = 0.0f,                  // Default x coordinate
            py::arg("y") = 0.0f,                  // Default y coordinate
            py::arg("z") = 0.0f,                  // Default z coordinate
            py::arg("altloc") = "",               // Default alternate location
            py::arg("occupancy") = 1.0f,          // Default occupancy
            py::arg("tfactor") = 0.0f,            // Default temperature factor
            py::arg("charge") = "")               // Default charge
        
        // Essential properties for basic functionality
        .def_readwrite("ID", &atom_struct::ID)                    // Atom ID
        .def_readwrite("NAME", &atom_struct::NAME)                // Atom name (e.g., "CA")
        .def_readwrite("RESN", &atom_struct::RESN)                // Residue name (e.g., "ALA")
        .def_readwrite("CHAIN", &atom_struct::CHAIN)              // Chain identifier
        .def_readwrite("RESI", &atom_struct::RESI)                // Residue number
        .def_readwrite("iCODE", &atom_struct::iCODE)             // Insertion code
        .def_readwrite("SASA", &atom_struct::SASA)               // Solvent accessible surface area
        .def_readwrite("MOL_TYPE", &atom_struct::MOL_TYPE)       // Molecule type
        .def_readwrite("STRUCT_TYPE", &atom_struct::STRUCT_TYPE)  // Structure type
        .def_readwrite("ELEMENT", &atom_struct::ELEMENT)          // Chemical element
        .def_readwrite("HETATM", &atom_struct::HETATM)           // HETATM flag
        .def_readwrite("ACTIVE", &atom_struct::ACTIVE)           // Active flag
        
        // Properties for coordinates and measurements
        .def_property("COORDS",
            [](const atom_struct& a) {
                return py::make_tuple(a.COORDS[0], a.COORDS[1], a.COORDS[2]);
            },
            [](atom_struct& a, py::sequence coords) {
                if (py::len(coords) != 3) {
                    throw std::runtime_error("Coordinates must have exactly 3 elements");
                }
                std::vector<float> new_coords = {
                    coords[0].cast<float>(),
                    coords[1].cast<float>(),
                    coords[2].cast<float>()
                };
                a.COORDS = new_coords;
            })
        .def_property("COORDS",
            [](const atom_struct& a) {
                return py::array_t<float>({3}, a.COORDS.data());
            },
            [](atom_struct& a, py::array_t<float> array) {
                if (array.size() != 3) {
                    throw std::runtime_error("Coordinates must be length 3");
                }
                a.COORDS = numpy_to_vector(array);
            })
                .def_property_readonly("contacts",
            [](const atom_struct& a) {
                return conversion::atom_contacts_to_dict(a);
            })
        .def_property_readonly("interaction_partners",
            [](const atom_struct& a) {
                return py::array_t<uint32_t>(
                    a.INTERACTION_SASA_P.size(),
                    a.INTERACTION_SASA_P.data()
                );
            })
        .def_readwrite("VDW", &atom_struct::VDW)                 // Van der Waals radius
        .def_readwrite("AREA", &atom_struct::AREA)               // Surface area
        
        // Additional properties
        .def_readwrite("ALTLOC", &atom_struct::ALTLOC)           // Alternate location
        .def_readwrite("OCCUPANCY", &atom_struct::OCCUPANCY)     // Occupancy
        .def_readwrite("TFACTOR", &atom_struct::TFACTOR)         // Temperature factor
        .def_readwrite("CHARGE", &atom_struct::CHARGE)           // Charge
        
        // Methods
        .def("sID", &atom_struct::sID)       // Get string ID
        .def("rsID", &atom_struct::rsID)     // Get residue string ID
        .def("print", [](const atom_struct& self) {
            atom_struct copy = self;
            return copy.print();
        })
        
        // String representation
        .def("__repr__", [](const atom_struct& self) {
            atom_struct copy = self;
            return "<AtomStruct: " + copy.print() + ">";
        });
}