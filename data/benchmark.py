import json
import time
import sys
from pathlib import Path
from typing import Dict, Any
import pandas as pd
from tqdm import tqdm
import argparse
import psutil
import os
import threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.append(str("/home/alessio/dr_sasa_python/build/lib"))
import dr_sasa_py

def calculate_interface_metrics(result: Dict) -> float:
    return result["intra_bsa_matrix"]["atom_matrix"].sum()

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
        
        return {
            'duration': time.time() - self._start_time,
            'cpu_percent': {
                'mean': np.mean(cpu_percents),
                'max': np.max(cpu_percents),
                'min': np.min(cpu_percents)
            },
            'threads': {
                'mean': np.mean(thread_counts),
                'max': np.max(thread_counts),
                'min': np.min(thread_counts)
            },
            'memory_percent': {
                'mean': np.mean(memory_usage),
                'max': np.max(memory_usage),
                'min': np.min(memory_usage)
            },
            'num_measurements': len(self.measurements)
        }
    
    def _monitor(self):
        """Background monitoring function."""
        while self._monitoring:
            try:
                self.measurements.append({
                    'timestamp': time.time() - self._start_time,
                    'cpu_percent': self.process.cpu_percent(),
                    'num_threads': self.process.num_threads(),
                    'memory_percent': self.process.memory_percent()
                })
                time.sleep(0.1)  # Sample every 100ms
            except:
                pass

def benchmark_dataset(dataset_json: str, pdb_dir: str) -> Dict[str, Any]:
    """Run benchmark comparing BSA calculations with improved CPU monitoring."""
    with open(dataset_json) as f:
        reference_data = json.load(f)
    
    results = {}
    calculator = dr_sasa_py.GenericSASA(probe_radius=1.4)
    
    # Initialize CPU monitor
    cpu_monitor = CPUMonitor()
    
    print("\nStarting CPU monitoring...")
    cpu_monitor.start()
    
    for pdb_id, ref_info in tqdm(reference_data.items()):
        pdb_file = Path(pdb_dir) / f"{pdb_id}.pdb"
        if "1ATN" in pdb_id:
            break
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
            reference_bsa = float(ref_info['BSA'])  # Ensure float type
            
            # Calculate differences
            absolute_diff = float(calculated_bsa - reference_bsa)
            relative_diff = float(((calculated_bsa + 1) - (reference_bsa + 1)) / (reference_bsa + 1))
            
            # Store results
            results[pdb_id] = {
                'calculation_time': float(calc_time),
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
    
    # Stop monitoring and get statistics
    cpu_stats = cpu_monitor.stop()
    results['cpu_statistics'] = cpu_stats
    
    return results

def print_detailed_cpu_summary(cpu_stats: Dict):
    """Print detailed CPU usage summary."""
    print("\nDetailed CPU Usage Summary:")
    print("-" * 40)
    print(f"Available CPU cores: {psutil.cpu_count()}")
    print(f"Available physical cores: {psutil.cpu_count(logical=False)}")
    
    if not cpu_stats:
        print("No CPU statistics available")
        return
        
    print(f"\nMonitoring duration: {cpu_stats['duration']:.2f} seconds")
    print(f"Number of measurements: {cpu_stats['num_measurements']}")
    
    print("\nCPU Usage:")
    print(f"  Average: {cpu_stats['cpu_percent']['mean']:.1f}%")
    print(f"  Maximum: {cpu_stats['cpu_percent']['max']:.1f}%")
    print(f"  Minimum: {cpu_stats['cpu_percent']['min']:.1f}%")
    
    print("\nThread Count:")
    print(f"  Average: {cpu_stats['threads']['mean']:.1f}")
    print(f"  Maximum: {cpu_stats['threads']['max']}")
    print(f"  Minimum: {cpu_stats['threads']['min']}")
    
    print("\nMemory Usage:")
    print(f"  Average: {cpu_stats['memory_percent']['mean']:.1f}%")
    print(f"  Maximum: {cpu_stats['memory_percent']['max']:.1f}%")
    print(f"  Minimum: {cpu_stats['memory_percent']['min']:.1f}%")

def analyze_results(results: Dict[str, Any]) -> pd.DataFrame:
    """Analyze benchmark results with CPU statistics."""
    records = []
    
    for pdb_id, result in results.items():
        if pdb_id == 'cpu_statistics' or result.get('error'):
            continue
            
        records.append({
            'pdb_id': pdb_id,
            'calculation_time': result['calculation_time'],
            'calculated_bsa': result['calculated_bsa'],
            'reference_bsa': result['reference_bsa'],
            'absolute_diff': result['absolute_diff'],
            'relative_diff': result['relative_diff'],
            'pre_calc_threads': result['cpu_usage']['pre_calculation']['num_threads'],
            'post_calc_threads': result['cpu_usage']['post_calculation']['num_threads'],
            'pre_calc_cpu_percent': result['cpu_usage']['pre_calculation']['cpu_percent'],
            'post_calc_cpu_percent': result['cpu_usage']['post_calculation']['cpu_percent'],
            'cpu_affinity': result['cpu_usage']['pre_calculation']['cpu_affinity']
        })
    
    return pd.DataFrame(records)

def main():
    parser = argparse.ArgumentParser(description='Compare calculated BSA values with reference dataset')
    parser.add_argument('-dataset_json', default = "data/dataset.json", help='Path to reference dataset JSON file')
    parser.add_argument('-pdb_dir', default  ="data/PRODIGYdataset_fixed" , help='Directory containing PDB files')
    parser.add_argument('-output_dir', default = "data/benchmark_results", help='Directory for output files')
    
   
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("\nRunning BSA calculations with CPU monitoring...")
    results = benchmark_dataset(args.dataset_json, args.pdb_dir)
    
    # Print summaries
    successful_results = {k: v for k, v in results.items() 
                        if k != 'cpu_statistics' and not v.get('error')}
    
    print("\nBenchmark Summary:")
    print("-" * 40)
    print(f"Structures processed: {len(successful_results)}")
    
    calc_times = [r['calculation_time'] for r in successful_results.values()]
    abs_diffs = [r['absolute_diff'] for r in successful_results.values()]
    rel_diffs = [r['relative_diff'] for r in successful_results.values()]
    
    print(f"Average calculation time: {np.mean(calc_times):.3f}s")
    
    print("\nBSA Comparison Statistics:")
    print(f"Mean absolute difference: {np.mean(abs_diffs):.2f}")
    print(f"Median absolute difference: {np.median(abs_diffs):.2f}")
    print(f"Mean relative difference: {np.mean(rel_diffs):.2%}")
    print(f"Median relative difference: {np.median(rel_diffs):.2%}")
    
    print_detailed_cpu_summary(results.get('cpu_statistics', {}))
    
    # Save results
    print("\nSaving results...")
    with open(output_dir / "full_results.json", "w") as f:
        json.dump(results, f, cls=NumpyEncoder, indent=2)
    
    # Save summary as CSV
    summary_data = []
    for pdb_id, result in successful_results.items():
        summary_data.append({
            'pdb_id': pdb_id,
            'calculation_time': result['calculation_time'],
            'calculated_bsa': result['calculated_bsa'],
            'reference_bsa': result['reference_bsa'],
            'absolute_diff': result['absolute_diff'],
            'relative_diff': result['relative_diff']
        })
    
    pd.DataFrame(summary_data).to_csv(output_dir / "bsa_comparison.csv", index=False)

if __name__ == "__main__":
    main()
