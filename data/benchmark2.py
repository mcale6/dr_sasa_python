import json
from datetime import datetime
import time
from pathlib import Path
from typing import Dict, Any, List, Tuple
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import psutil
import os
import threading
import dr_sasa_python as sasa
import freesasa
from Bio.PDB import *
from Bio.PDB.SASA import ShrakeRupley
import subprocess
import tempfile

class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for NumPy types."""
    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                          np.int16, np.int32, np.int64, np.uint8,
                          np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

class CPUMonitor:
    def __init__(self):
        self.process = psutil.Process()
        self.measurements = []
        self._monitoring = False
    
    def start(self):
        """Start CPU monitoring in background."""
        self._monitoring = True
        self.measurements = []
        self._start_time = time.time()
        threading.Thread(target=self._monitor, daemon=True).start()
    
    def stop(self):
        """Stop CPU monitoring and return statistics."""
        self._monitoring = False
        if not self.measurements:
            return {}
            
        cpu_percents = [m['cpu_percent'] for m in self.measurements]
        thread_counts = [m['num_threads'] for m in self.measurements]
        memory_usage = [m['memory_percent'] for m in self.measurements]
        python_threads = [m['python_threads'] for m in self.measurements]
        openmp_threads = [m['openmp_threads'] for m in self.measurements]
        system_threads = [m['system_threads'] for m in self.measurements]
        
        return {
            'duration': time.time() - self._start_time,
            'stats': {
                'cpu_percent': {
                    'mean': np.mean(cpu_percents),
                    'max': np.max(cpu_percents),
                    'min': np.min(cpu_percents)
                },
                'threads': {
                    'mean': np.mean(thread_counts),
                    'max': np.max(thread_counts),
                    'min': np.min(thread_counts),
                    'python': {
                        'mean': np.mean(python_threads),
                        'max': np.max(python_threads)
                    },
                    'openmp': {
                        'mean': np.mean(openmp_threads),
                        'max': np.max(openmp_threads)
                    },
                    'system': {
                        'mean': np.mean(system_threads),
                        'max': np.max(system_threads)
                    }
                },
                'memory_percent': {
                    'mean': np.mean(memory_usage),
                    'max': np.max(memory_usage),
                    'min': np.min(memory_usage)
                },
                'hardware': {
                    'physical_cores': psutil.cpu_count(logical=False),
                    'logical_cores': psutil.cpu_count(),
                    'numa_nodes': len(psutil.numa_partitions()) if hasattr(psutil, 'numa_partitions') else 1
                }
            },
            'num_measurements': len(self.measurements)
        }

    def _monitor(self):
        """Background monitoring function."""
        while self._monitoring:
            try:
                # Get OpenMP thread count from environment
                openmp_threads = int(os.environ.get('OMP_NUM_THREADS', '0'))
                # Get Python thread count
                python_threads = threading.active_count()
                # Get total process threads
                total_threads = self.process.num_threads()
                # Calculate system threads (total minus Python and OpenMP)
                system_threads = max(0, total_threads - python_threads - openmp_threads)
                
                self.measurements.append({
                    'timestamp': time.time() - self._start_time,
                    'cpu_percent': self.process.cpu_percent(),
                    'num_threads': total_threads,
                    'python_threads': python_threads,
                    'openmp_threads': openmp_threads,
                    'system_threads': system_threads,
                    'memory_percent': self.process.memory_percent()
                })
                time.sleep(0.1)  # Sample every 100ms
            except Exception as e:
                print(f"Monitoring error: {str(e)}")
                pass

def calculate_original_dr_sasa(pdb_file: str, dr_sasa_exec: str) -> Tuple[np.ndarray, float]:
    """Calculate SASA using original C++ implementation."""
    # Create temporary output PDB file
    with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w+') as temp_output:
        start_time = time.time()
        
        # Run dr_sasa command
        try:
            cmd = [dr_sasa_exec, "-m", "0", "-i", pdb_file, "-o", temp_output.name]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"dr_sasa failed: {result.stderr}")
            
            # Read output PDB and extract B-factors (SASA values)
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', temp_output.name)
            
            # Extract SASA values from B-factor column
            atom_sasa = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom_sasa.append(atom.get_bfactor())
            
        except Exception as e:
            raise RuntimeError(f"Error running original dr_sasa: {str(e)}")
        
        calc_time = time.time() - start_time
        return np.array(atom_sasa), calc_time

def calculate_dr_sasa(pdb_file: str) -> Tuple[np.ndarray, float]:
    """Calculate SASA using dr_sasa."""
    calculator = sasa.SimpleSASA(probe_radius=1.4)
    
    start_time = time.time()
    result = calculator.calculate(str(pdb_file))
    calc_time = time.time() - start_time
    
    # Extract per-atom SASA values
    atom_sasa = np.asarray([result["atom_data"][str(i)]["sasa"] for i in list(result["atom_data"].keys())])
    
    return atom_sasa, calc_time

def calculate_freesasa(pdb_file: str) -> Tuple[np.ndarray, float]:
    """Calculate SASA using FreeSASA."""
    start_time = time.time()
    
    # Configure FreeSASA
    parameters = freesasa.Parameters(
        {'algorithm': freesasa.ShrakeRupley,
         'probe-radius': 1.4,
         'n-points': 100}
    )
    
    # Calculate SASA
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure, parameters)
    
    # Extract per-atom values
    atom_sasa = np.array([result.atomArea(i) for i in range(structure.nAtoms())])
    
    calc_time = time.time() - start_time
    return atom_sasa, calc_time

def calculate_biopython(pdb_file: str) -> Tuple[np.ndarray, float]:
    """Calculate SASA using Biopython."""
    start_time = time.time()
    
    # Parse structure
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    
    # Calculate SASA
    sr = ShrakeRupley()
    sr.compute(structure, level="A")
    
    # Extract per-atom values
    atom_sasa = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_sasa.append(atom.sasa)
    
    calc_time = time.time() - start_time
    return np.array(atom_sasa), calc_time

def calculate_differences(values1: np.ndarray, values2: np.ndarray) -> Dict:
    """Calculate statistical differences between two sets of atomic SASA values."""
    if len(values1) != len(values2):
        raise ValueError(f"Arrays have different lengths: {len(values1)} vs {len(values2)}")
        
    differences = values1 - values2
    squared_diff = np.square(differences)
    
    return {
        'mean_diff': float(np.mean(differences)),
        'std_diff': float(np.std(differences)),
        'rmsd': float(np.sqrt(np.mean(squared_diff))),
        'max_abs_diff': float(np.max(np.abs(differences))),
        'mean_squared_error': float(np.mean(squared_diff)),
        'correlation': float(np.corrcoef(values1, values2)[0,1])
    }

def analyze_structure(pdb_file: str, dr_sasa_exec: str) -> Dict:
    """Analyze a structure using different SASA calculation methods."""
    results = {}
    
    try:
        # Calculate SASA with each method
        dr_sasa_values, dr_sasa_time = calculate_dr_sasa(pdb_file)
        original_values, original_time = calculate_original_dr_sasa(pdb_file, dr_sasa_exec)
        freesasa_values, freesasa_time = calculate_freesasa(pdb_file)
        biopython_values, biopython_time = calculate_biopython(pdb_file)
        
        results['times'] = {
            'dr_sasa_py': dr_sasa_time,
            'dr_sasa_original': original_time,
            'freesasa': freesasa_time,
            'biopython': biopython_time
        }
        
        # Compare results
        results['comparisons'] = {
            'dr_sasa_py_vs_original': calculate_differences(dr_sasa_values, original_values),
            'dr_sasa_py_vs_freesasa': calculate_differences(dr_sasa_values, freesasa_values),
            'dr_sasa_py_vs_biopython': calculate_differences(dr_sasa_values, biopython_values),
            'original_vs_freesasa': calculate_differences(original_values, freesasa_values),
            'original_vs_biopython': calculate_differences(original_values, biopython_values),
            'freesasa_vs_biopython': calculate_differences(freesasa_values, biopython_values)
        }
        
        # Store number of atoms and raw values for detailed analysis
        results['n_atoms'] = len(dr_sasa_values)
        results['values'] = {
            'dr_sasa_py': dr_sasa_values.tolist(),
            'dr_sasa_original': original_values.tolist(),
            'freesasa': freesasa_values.tolist(),
            'biopython': biopython_values.tolist()
        }
        
    except Exception as e:
        results['error'] = str(e)
    
    return results

def create_comparison_plots(results: Dict, output_dir: Path):
    """Create detailed comparison plots between different implementations."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Scatter plot of SASA values
    plt.figure(figsize=(12, 8))
    for pdb_id, result in results.items():
        if 'error' not in result:
            plt.scatter(result['values']['dr_sasa_original'],
                       result['values']['dr_sasa_py'],
                       alpha=0.5, label=pdb_id)
    
    plt.plot([0, plt.xlim()[1]], [0, plt.xlim()[1]], 'k--')
    plt.xlabel('Original dr_sasa SASA values')
    plt.ylabel('Python dr_sasa SASA values')
    plt.title('Comparison of Original vs Python Implementation')
    plt.savefig(output_dir / 'implementation_comparison.png')
    plt.close()
    
    # Distribution of differences
    differences = []
    for result in results.values():
        if 'error' not in result:
            diff = np.array(result['values']['dr_sasa_py']) - \
                   np.array(result['values']['dr_sasa_original'])
            differences.extend(diff)
    
    plt.figure(figsize=(10, 6))
    sns.histplot(differences, bins=50)
    plt.xlabel('SASA Difference (Python - Original)')
    plt.ylabel('Count')
    plt.title('Distribution of SASA Differences')
    plt.savefig(output_dir / 'differences_distribution.png')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Compare SASA calculations across different tools')
    parser.add_argument('-dataset_json', default="data/dataset.json", help='Path to dataset JSON file')
    parser.add_argument('-pdb_dir', default="data/PRODIGYdataset", help='Directory containing PDB files')
    parser.add_argument('-output_dir', default="data/benchmark_results", help='Directory for output files')
    parser.add_argument('-dr_sasa_exec', default="/home/alessio/dr_sasa_python/dr_sasa_n/build/dr_sasa",  help='Path to original dr_sasa executable')
    parser.add_argument('--max_structures', type=int, default=None, help='Maximum number of structures to process (default: all)')
    
    args = parser.parse_args()
    

    # Load dataset
    with open(args.dataset_json) as f:
        dataset = json.load(f)
    
    # Initialize results storage
    all_results = {}
    timing_stats = {
        'dr_sasa_py': [],
        'dr_sasa_original': [],
        'freesasa': [],
        'biopython': []
    }
    comparison_stats = {
        'dr_sasa_py_vs_original': {'rmsd': [], 'correlation': [], 'max_diff': []},
        'dr_sasa_py_vs_freesasa': {'rmsd': [], 'correlation': [], 'max_diff': []},
        'dr_sasa_py_vs_biopython': {'rmsd': [], 'correlation': [], 'max_diff': []},
        'original_vs_freesasa': {'rmsd': [], 'correlation': [], 'max_diff': []},
        'original_vs_biopython': {'rmsd': [], 'correlation': [], 'max_diff': []},
        'freesasa_vs_biopython': {'rmsd': [], 'correlation': [], 'max_diff': []}
    }
    
    # Initialize CPU monitor
    cpu_monitor = CPUMonitor()
    cpu_monitor.start()
    
    # Create output directory
    current_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_dir = Path(f"{args.output_dir}_{current_time}")
    output_dir.mkdir(exist_ok=True)
        
    # Process structures
    pdb_list = list(dataset.keys())
    if args.max_structures:
        pdb_list = pdb_list[:args.max_structures]
    
    print(f"\nProcessing {len(pdb_list)} structures...")
    for pdb_id in tqdm(pdb_list):
        pdb_file = str(Path(args.pdb_dir) / f"{pdb_id}.pdb")
        if not Path(pdb_file).exists():
            print(f"Missing PDB file: {pdb_file}")
            continue
        
        try:
            # Analyze structure with all methods
            results = analyze_structure(pdb_file, args.dr_sasa_exec)
            all_results[pdb_id] = results
            
            if 'error' not in results:
                # Collect timing statistics
                for method, time in results['times'].items():
                    timing_stats[method].append(time)
                
                # Collect comparison statistics
                for comp_name, stats in results['comparisons'].items():
                    comparison_stats[comp_name]['rmsd'].append(stats['rmsd'])
                    comparison_stats[comp_name]['correlation'].append(stats['correlation'])
                    comparison_stats[comp_name]['max_diff'].append(stats['max_abs_diff'])
            else:
                print(f"Error processing {pdb_id}: {results['error']}")
        
        except Exception as e:
            print(f"Unexpected error processing {pdb_id}: {str(e)}")
            continue
    
    # Stop CPU monitoring
    cpu_stats = cpu_monitor.stop()
    
    # Print comprehensive results
    print("\nAnalysis Results")
    print("=" * 50)
    
    print("\nHardware Configuration:")
    print("-" * 40)
    print(f"Physical cores: {psutil.cpu_count(logical=False)}")
    print(f"Logical cores: {psutil.cpu_count()}")
    
    print("\nProcessing Summary:")
    print("-" * 40)
    successful = len([r for r in all_results.values() if 'error' not in r])
    failed = len(all_results) - successful
    print(f"Total structures processed: {len(all_results)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    
    print("\nTiming Statistics (seconds):")
    print("-" * 40)
    for method, times in timing_stats.items():
        if times:
            print(f"{method}:")
            print(f"  Mean: {np.mean(times):.3f} ±{np.std(times):.3f}")
            print(f"  Median: {np.median(times):.3f}")
            print(f"  Min: {np.min(times):.3f}")
            print(f"  Max: {np.max(times):.3f}")
    
    print("\nSASA Comparison Statistics:")
    print("-" * 40)
    for comp_name, stats in comparison_stats.items():
        if stats['rmsd']:
            print(f"\n{comp_name}:")
            print(f"  RMSD: {np.mean(stats['rmsd']):.6f} ±{np.std(stats['rmsd']):.6f}")
            print(f"  Correlation: {np.mean(stats['correlation']):.6f} ±{np.std(stats['correlation']):.6f}")
            print(f"  Max difference: {np.mean(stats['max_diff']):.6f} ±{np.std(stats['max_diff']):.6f}")
    
    print("\nCPU Usage Summary:")
    print("-" * 40)
    print(f"Average CPU usage: {cpu_stats['stats']['cpu_percent']['mean']:.1f}%")
    print(f"Peak CPU usage: {cpu_stats['stats']['cpu_percent']['max']:.1f}%")
    print(f"Average thread count: {cpu_stats['stats']['threads']['mean']:.1f}")
    print(f"Memory usage: {cpu_stats['stats']['memory_percent']['mean']:.1f}%")
    
    # Save detailed results
    print("\nSaving detailed results...")
    
    # Save raw results
    with open(output_dir / "full_results.json", "w") as f:
        json.dump(all_results, f, cls=NumpyEncoder, indent=2)
    
    # Create and save summary DataFrame
    summary_data = []
    for pdb_id, result in all_results.items():
        if 'error' not in result:
            record = {
                'pdb_id': pdb_id,
                'n_atoms': result['n_atoms']
            }
            
            # Add timing information
            for method, time in result['times'].items():
                record[f'{method}_time'] = time
            
            # Add comparison metrics
            for comp_name, stats in result['comparisons'].items():
                for metric in ['rmsd', 'correlation', 'mean_squared_error', 'max_abs_diff']:
                    record[f'{comp_name}_{metric}'] = stats[metric]
            
            summary_data.append(record)
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_dir / "comparison_summary.csv", index=False)
    
    # Create visualization plots
    print("\nGenerating visualization plots...")
    try:
        create_comparison_plots(all_results, output_dir)
    except Exception as e:
        print(f"Error creating plots: {str(e)}")
    
    print(f"\nAnalysis complete. Results saved in: {output_dir}")

if __name__ == "__main__":
    main()