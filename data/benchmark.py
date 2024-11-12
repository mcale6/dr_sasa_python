import json
import glob
import time
import numpy as np
from pathlib import Path
import dr_sasa_py as sasa
from tqdm import tqdm

def parse_chain_spec(chain_spec):
    """Parse chain specification like 'A:BC' into two lists ['A'], ['B', 'C']"""
    if ':' not in chain_spec:
        return [], []
    parts = chain_spec.split(':')
    return list(parts[0]), list(parts[1])

def benchmark_dataset(dataset_json, pdb_dir):
    # Load dataset
    with open(dataset_json) as f:
        dataset = json.load(f)
    
    results = {}
    calculator = sasa.SimpleSASA()
    
    # Process each PDB
    for pdb_id, info in tqdm(dataset.items()):
        pdb_file = Path(pdb_dir) / f"{pdb_id}.pdb"
        if not pdb_file.exists():
            print(f"Missing PDB file: {pdb_file}")
            continue
            
        try:
            # Parse chain specification
            chain1, chain2 = parse_chain_spec(info['Interacting_chains'])
            
            # Time the calculation
            start_time = time.time()
            result = calculator.calculate(str(pdb_file))
            calc_time = time.time() - start_time
            
            # Get arrays
            sasa = result['sasa']
            bsa = result['bsa']
            rasa = result['rasa']
            
            # Calculate statistics
            total_bsa = np.sum(bsa)
            avg_rasa = np.mean(rasa)
            
            results[pdb_id] = {
                'calculation_time': calc_time,
                'total_bsa': float(total_bsa),
                'avg_rasa': float(avg_rasa),
                'reference_bsa': info.get('BSA', None),
                'error': None
            }
            
        except Exception as e:
            results[pdb_id] = {
                'error': str(e)
            }
    
    return results

if __name__ == "__main__":
    dataset_json = "dataset.json"
    pdb_dir = "pdbs"
    
    results = benchmark_dataset(dataset_json, pdb_dir)
    
    # Print summary
    print("\nBenchmark Summary:")
    print("-----------------")
    success = len([r for r in results.values() if not r.get('error')])
    total = len(results)
    print(f"Successfully processed: {success}/{total} structures")
    
    # Calculate statistics for successful runs
    times = [r['calculation_time'] for r in results.values() if not r.get('error')]
    if times:
        print(f"Average calculation time: {np.mean(times):.3f}s")
        print(f"Min calculation time: {np.min(times):.3f}s")
        print(f"Max calculation time: {np.max(times):.3f}s")
    
    # Compare with reference BSA values
    bsa_errors = []
    for pdb_id, result in results.items():
        if result.get('error'):
            continue
        if result['reference_bsa'] is not None:
            error = abs(result['total_bsa'] - result['reference_bsa'])
            bsa_errors.append(error)
            
    if bsa_errors:
        print(f"\nBSA Comparison:")
        print(f"Average absolute error: {np.mean(bsa_errors):.2f} Å²")
        print(f"Max absolute error: {np.max(bsa_errors):.2f} Å²")
    
    # Save detailed results
    with open("benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)