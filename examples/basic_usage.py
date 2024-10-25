import dr_sasa_py

def example_basic_sasa():
    """Basic SASA calculation example"""
    # Initialize DrSASA
    dr_sasa = dr_sasa_py.DrSASA()
    
    # Calculate SASA
    results = dr_sasa.calculate_sasa(
        pdb_file="tests/data/1bl0.pdb",
        keep_hetatm=True
    )
    
    print(f"Total SASA: {results['total_sasa']:.2f} Å²")
    print(f"Number of atoms: {len(results['atom_sasa'])}")
    
def example_chain_analysis():
    """Chain-based analysis example"""
    dr_sasa = dr_sasa_py.DrSASA(probe_radius=1.4)
    
    # Calculate delta SASA between chains
    results = dr_sasa.calculate_delta_sasa(
        pdb_file="tests/data/4ins.pdb",
        chains=[["A", "B"], ["C", "D"]]
    )
    
    print(f"Interface area: {results['total_interface_area']:.2f} Å²")
    print(f"Number of interface atoms: {len(results['interface_atoms'])}")

if __name__ == "__main__":
    print("Basic SASA calculation example:")
    example_basic_sasa()
    print("\nChain analysis example:")
    example_chain_analysis()