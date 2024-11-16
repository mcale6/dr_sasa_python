import json
import time
import numpy as np
from typing import *
from pathlib import Path
import sys
from tqdm import tqdm
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
sys.path.append(str("/home/alessio/dr_sasa_python/bindings/python/utils"))
sys.path.append(str("/home/alessio/dr_sasa_python/build/lib"))
import dr_sasa_py
from utils import *

# Add the path to dr_sasa_py if needed
project_root = Path(__file__).parent.parent
sys.path.append(str(project_root / "build/lib"))

class SASAMetrics:
    """Calculate and store SASA-related metrics."""
    
    @staticmethod
    def calculate_bsa(sasa_free: float, sasa_complexed: float) -> float:
        """Calculate Buried Surface Area (BSA).
        
        BSA = ∆SASA = SASA_free - SASA_complexed
        """
        return sasa_free - sasa_complexed
    
    @staticmethod
    def calculate_absolute_diff(value1: float, value2: float) -> float:
        """Calculate absolute difference.
        
        Diff(SASA_i) = SASA_i^ToolX - SASA_i^dr_sasa
        """
        return value1 - value2
    
    @staticmethod
    def calculate_relative_diff(value1: float, value2: float) -> float:
        """Calculate relative difference.
        
        Diff_rel(SASA_i) = ((SASA_i^ToolX + 1) - (SASA_i^dr_sasa + 1))/(SASA_i^dr_sasa + 1)
        """
        return ((value1 + 1) - (value2 + 1))/(value2 + 1)
    
    @staticmethod
    def calculate_error_per_atom(value: float, reference: float, n_atoms: int) -> float:
        """Calculate error per atom.
        
        ε = |A-A_ref|/N
        """
        return abs(value - reference)/n_atoms

def calculate_interface_metrics(dfs: Dict[str, pd.DataFrame]) -> Dict[str, float]:
    """Calculate interface metrics from DataFrame representations.
    
    Returns dictionary with:
    - BSA for each chain
    - Interface area
    - Complex surface
    - Uncomplexed surface
    """
    metrics = {}
    
    # Chain A metrics
    metrics['chain_A_BSA'] = dfs['overlaps'][
        dfs['overlaps']["source_residue"].str.startswith("A_")
    ].overlap_area.sum()
    
    # Chain B metrics
    metrics['chain_B_BSA'] = dfs['overlaps'][
        dfs['overlaps']["source_residue"].str.startswith("B_")
    ].overlap_area.sum()
    
    # Interface area
    metrics['interface_area'] = dfs['overlaps'].overlap_area.sum()/2
    
    # Complex surface areas
    metrics['chain_A_complex_surface'] = dfs['atoms'][
        dfs['atoms'].chain == "A"
    ].sasa.sum()
    
    metrics['chain_B_complex_surface'] = dfs['atoms'][
        dfs['atoms'].chain == "B"
    ].sasa.sum()
    
    metrics['total_complex_surface'] = dfs['atoms'].sasa.sum()
    
    # Uncomplexed surface areas
    metrics['chain_A_uncomplex_surface'] = (
        metrics['chain_A_BSA'] + metrics['chain_A_complex_surface']
    )
    
    return metrics

def benchmark_dataset(dataset_json: str, 
                     pdb_dir: str, 
                     tool: str = "dr_sasa",
                     reference_tool: Optional[str] = None) -> Dict[str, Any]:
    """Run benchmark comparing SASA calculations.
    
    Args:
        dataset_json: Path to dataset JSON
        pdb_dir: Directory containing PDB files
        tool: Tool to benchmark
        reference_tool: Optional reference tool for comparison
        
    Returns:
        Dictionary containing benchmark results
    """
    # Load dataset
    with open(dataset_json) as f:
        dataset = json.load(f)
    
    results = {}
    metrics = SASAMetrics()
    
    # Initialize calculator
    calculator = dr_sasa_py.GenericSASA(probe_radius=1.4)
    
    # Process each structure
    for pdb_id, info in tqdm(dataset.items()):
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
            
            # Convert to DataFrames
            dfs = convert_to_dataframes(result)
            
            # Calculate metrics
            interface_metrics = calculate_interface_metrics(dfs)
            
            # Store results
            results[pdb_id] = {
                'calculation_time': calc_time,
                'metrics': interface_metrics,
                'atom_data': result['atom_data'],
                'residue_data': result['residue_data'],
                'error': None
            }
            
            # Compare with reference if provided
            if reference_tool:
                ref_result = calculate_with_reference(pdb_file, reference_tool)
                
                # Calculate differences
                for metric_name, metric_value in interface_metrics.items():
                    abs_diff = metrics.calculate_absolute_diff(
                        metric_value,
                        ref_result['metrics'][metric_name]
                    )
                    rel_diff = metrics.calculate_relative_diff(
                        metric_value,
                        ref_result['metrics'][metric_name]
                    )
                    error_per_atom = metrics.calculate_error_per_atom(
                        metric_value,
                        ref_result['metrics'][metric_name],
                        len(dfs['atoms'])
                    )
                    
                    results[pdb_id][f'{metric_name}_abs_diff'] = abs_diff
                    results[pdb_id][f'{metric_name}_rel_diff'] = rel_diff
                    results[pdb_id][f'{metric_name}_error_per_atom'] = error_per_atom
            
        except Exception as e:
            results[pdb_id] = {
                'error': str(e)
            }
            print(f"Error processing {pdb_id}: {str(e)}")
    
    return results

def analyze_results(results: Dict[str, Any]) -> pd.DataFrame:
    """Analyze benchmark results and create summary statistics."""
    summary = []
    
    # Process each structure
    for pdb_id, result in results.items():
        if result.get('error'):
            continue
            
        record = {
            'pdb_id': pdb_id,
            'calculation_time': result['calculation_time']
        }
        
        # Add all metrics
        for metric_name, metric_value in result['metrics'].items():
            record[metric_name] = metric_value
            
        # Add differences if present
        for key in result:
            if key.endswith(('_abs_diff', '_rel_diff', '_error_per_atom')):
                record[key] = result[key]
                
        summary.append(record)
    
    return pd.DataFrame(summary)


if __name__ == "__main__":
    # Configuration
    dataset_json = "data/dataset.json"
    pdb_dir = Path("/path/to/pdbs")
    output_dir = Path("benchmark_results")
    output_dir.mkdir(exist_ok=True)
    
    # Example usage of direct calculation
    test_pdb = "/path/to/test.pdb"
    calc = dr_sasa_py.GenericSASA(probe_radius=1.4)
    result = calc.calculate(str(test_pdb), include_matrix=True, print_output=True)
    
    # Convert to DataFrames for easy analysis
    dfs = convert_to_dataframes(result)
    
    # Access patterns
    chain_a_bsa = dfs['overlaps'][
        dfs['overlaps']["source_residue"].str.startswith("A_")
    ].overlap_area.sum()
    
    chain_a_sasa = dfs['atoms'][dfs['atoms'].chain == "A"].sasa.sum()
    chain_a_dsasa = dfs['residues'][dfs['residues'].chain == "A"].dsasa.sum()
    
    # Run benchmark
    results = benchmark_dataset(dataset_json, pdb_dir, "dr_sasa", reference_tool="ToolX")
    
    # Analyze results
    summary_df = analyze_results(results)
    
    # Create visualizations
    plot_results(summary_df, output_dir)
    
    # Save results
    summary_df.to_csv(output_dir / "benchmark_summary.csv")
    with open(output_dir / "full_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    # Print summary statistics
    print("\nBenchmark Summary:")
    print("-" * 40)
    print(f"Structures processed: {len(summary_df)}")
    print("\nTiming Statistics:")
    print(f"Average calculation time: {summary_df['calculation_time'].mean():.3f}s")
    print(f"Median calculation time: {summary_df['calculation_time'].median():.3f}s")
    
    if 'chain_A_BSA_error_per_atom' in summary_df.columns:
        print("\nError Statistics:")
        for col in summary_df.columns:
            if col.endswith('_error_per_atom'):
                print(f"{col}:")
                print(f"  Mean: {summary_df[col].mean():.4f}")
                print(f"  Median: {summary_df[col].median():.4f}")
                print(f"  Std: {summary_df[col].std():.4f}")
