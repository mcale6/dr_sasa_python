import json
import time
import sys
from pathlib import Path
from typing import Dict, Any
import pandas as pd
from tqdm import tqdm
import argparse
sys.path.append(str("/home/alessio/dr_sasa_python/build/lib"))
import dr_sasa_py

def calculate_interface_metrics(result: Dict) -> float:
    return result["intra_bsa_matrix"]["atom_matrix"].sum()

def benchmark_dataset(dataset_json: str, pdb_dir: str) -> Dict[str, Any]:
    # Load reference dataset
    with open(dataset_json) as f:
        reference_data = json.load(f)
    
    results = {}
    
    # Initialize calculator
    calculator = dr_sasa_py.GenericSASA(probe_radius=1.4)
    
    # Process each structure
    for pdb_id, ref_info in tqdm(reference_data.items()):
        pdb_file = Path(pdb_dir) / f"{pdb_id}.pdb"
        if not pdb_file.exists():
            print(f"Missing PDB file: {pdb_file}")
            continue
            
        try:
            # Calculate SASA
            start_time = time.time()
            result = calculator.calculate(
                str(pdb_file),
                include_matrix=True
            )
            calc_time = time.time() - start_time
            
            # Calculate BSA
            calculated_bsa = calculate_interface_metrics(result)
            reference_bsa = ref_info['BSA']
            
            # Calculate differences
            absolute_diff = calculated_bsa - reference_bsa
            relative_diff = ((calculated_bsa + 1) - (reference_bsa + 1)) / (reference_bsa + 1)
            
            # Store results
            results[pdb_id] = {
                'calculation_time': calc_time,
                'calculated_bsa': calculated_bsa,
                'reference_bsa': reference_bsa,
                'absolute_diff': absolute_diff,
                'relative_diff': relative_diff,
                'error': None
            }
            
        except Exception as e:
            results[pdb_id] = {
                'error': str(e)
            }
            print(f"Error processing {pdb_id}: {str(e)}")
    
    return results

def analyze_results(results: Dict[str, Any]) -> pd.DataFrame:
    """Analyze benchmark results and create summary statistics."""
    records = []
    
    for pdb_id, result in results.items():
        if result.get('error'):
            continue
            
        records.append({
            'pdb_id': pdb_id,
            'calculation_time': result['calculation_time'],
            'calculated_bsa': result['calculated_bsa'],
            'reference_bsa': result['reference_bsa'],
            'absolute_diff': result['absolute_diff'],
            'relative_diff': result['relative_diff']
        })
    
    return pd.DataFrame(records)

def main():
    parser = argparse.ArgumentParser(description='Compare calculated BSA values with reference dataset')
    parser.add_argument('-dataset_json', default = "data/dataset.json", help='Path to reference dataset JSON file')
    parser.add_argument('-pdb_dir', default  ="data/PRODIGYdataset_fixed" , help='Directory containing PDB files')
    parser.add_argument('-output_dir', default = "data/benchmark_results", help='Directory for output files')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Run benchmark
    print("\nRunning BSA calculations...")
    results = benchmark_dataset(args.dataset_json, args.pdb_dir)
    
    # Analyze results
    print("\nAnalyzing results...")
    summary_df = analyze_results(results)
    
    # Save results
    summary_df.to_csv(output_dir / "bsa_comparison.csv", index=False)

    # Print summary statistics
    print("\nBenchmark Summary:")
    print("-" * 40)
    print(f"Structures processed: {len(summary_df)}")
    print(f"Average calculation time: {summary_df['calculation_time'].mean():.3f}s")
    
    print("\nBSA Comparison Statistics:")
    print(f"Mean absolute difference: {summary_df['absolute_diff'].mean():.2f}")
    print(f"Median absolute difference: {summary_df['absolute_diff'].median():.2f}")
    print(f"Mean relative difference: {summary_df['relative_diff'].mean():.2%}")
    print(f"Median relative difference: {summary_df['relative_diff'].median():.2%}")

if __name__ == "__main__":
    main()