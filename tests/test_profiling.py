from pathlib import Path
import time
import sys
sys.path.append(str("/home/alessio/dr_sasa_python/build/lib"))
sys.path.append(str("/home/alessio/dr_sasa_python/bindings/python/utils"))
import dr_sasa_py

# Constants
TEST_DATA_DIR = Path(__file__).parent / "data"
TEST_FILES = {
    "small": "pred.pdb",    # Small protein complex
    #"medium": "1bl0.pdb",   # Medium-sized protein-DNA complex
    #"large": "6gwp.pdb"     # Large complex
}

def run_sasa_profiling():
    """
    Profile SASA calculations using both SimpleSASA and GenericSASA.
    Each calculator is tested with 3 different sized structures.
    """
    probe_radius = 1.4
    
    # Initialize calculators
    calculators = {
        "simple": dr_sasa_py.SimpleSASA(probe_radius=probe_radius),
        #"generic": dr_sasa_py.GenericSASA(probe_radius=probe_radius)
    }

    # Process each test file with both calculators
    for file_size, filename in TEST_FILES.items():
        pdb_path = TEST_DATA_DIR / filename
        if not pdb_path.exists():
            print(f"Test file not found: {pdb_path}")
            continue

        print(f"\nProcessing {file_size} structure: {filename}")
        
        for calc_name, calculator in calculators.items():
            print(f"\nUsing {calc_name} calculator:")
            
            # Time the calculation
            start_time = time.time()
            
            try:
                if calc_name == "simple":
                    results = calculator.calculate(str(pdb_path))
                else:
                    # For GenericSASA, also test chain interactions
                    results = calculator.calculate(
                        str(pdb_path),
                        chains=[["A"], ["B"]],
                        include_matrix=True
                    )
                
                end_time = time.time()
                duration = end_time - start_time
                
                # Print basic statistics
                n_atoms = len(results["atom_data"])
                n_residues = len(results["residue_data"])
                print(f"Time taken: {duration:.3f} seconds")
                print(f"Atoms processed: {n_atoms}")
                print(f"Residues processed: {n_residues}")
                print(f"Processing speed: {n_atoms/duration:.1f} atoms/second")
                
            except Exception as e:
                print(f"Error processing {filename} with {calc_name}: {str(e)}")

if __name__ == "__main__":
    run_sasa_profiling() #py-spy record -o profile.svg --native --format flamegraph  --rate 100 python tests/test_profiling.py