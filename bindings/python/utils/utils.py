import pandas as pd
import numpy as np
from typing import Dict, Any
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

def convert_to_dataframes(results: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    """Convert SASA calculation results to pandas DataFrames.
    
    Args:
        results: Dictionary containing SASA calculation results
        
    Returns:
        Dictionary of DataFrames:
            - 'atoms': Atom-level data with detailed ASA components
            - 'residues': Residue-level data
            - 'contacts': Contact data (if present)
            - 'overlaps': Overlap data (if present)
    """
    dfs = {}
    
    # Convert atom data to DataFrame with all ASA components
    atom_data = results['atom_data']
    atom_records = []
    
    for atom_id, atom_info in atom_data.items():
        record = {
            'id': int(atom_id),
            'name': atom_info['name'],
            'resname': atom_info['resname'],
            'chain': atom_info['chain'],
            'resid': atom_info['resid'],
            'struct_type': atom_info['struct_type'],
            'x': atom_info['coords'][0],
            'y': atom_info['coords'][1],
            'z': atom_info['coords'][2],
            'sphere_area': atom_info['sphere_area'],
            'sasa': atom_info['sasa'],
            'dsasa': atom_info.get('dsasa', None),  # New field
            'vdw': atom_info.get('vdw', None),  # New field
            'polar': atom_info['polar'],
            'charge': atom_info['charge'],
            
            'total_asa': atom_info.get('total_asa', None),
            'bb_asa': atom_info.get('bb_asa', None),
            'sc_asa': atom_info.get('sc_asa', None),
            'majorgroove_asa': atom_info.get('majorgroove_asa', None),
            'minorgroove_asa': atom_info.get('minorgroove_asa', None),
            'nogroove_asa': atom_info.get('nogroove_asa', None),
            
            'polar_asa': atom_info.get('polar_asa', None),
            'polar_bb_asa': atom_info.get('polar_bb_asa', None),
            'polar_sc_asa': atom_info.get('polar_sc_asa', None),
            'polar_majorgroove_asa': atom_info.get('polar_majorgroove_asa', None),
            'polar_minorgroove_asa': atom_info.get('polar_minorgroove_asa', None),
            'polar_nogroove_asa': atom_info.get('polar_nogroove_asa', None),
            
            'hyd_asa': atom_info.get('hyd_asa', None),
            'hyd_bb_asa': atom_info.get('hyd_bb_asa', None),
            'hyd_sc_asa': atom_info.get('hyd_sc_asa', None),
            'hyd_majorgroove_asa': atom_info.get('hyd_majorgroove_asa', None),
            'hyd_minorgroove_asa': atom_info.get('hyd_minorgroove_asa', None),
            'hyd_nogroove_asa': atom_info.get('hyd_nogroove_asa', None),
            
            'lig_asa': atom_info.get('lig_asa', None),
            'lig_polar_asa': atom_info.get('lig_polar_asa', None),
            'lig_hyd_asa': atom_info.get('lig_hyd_asa', None)
        }
        atom_records.append(record)
    
    dfs['atoms'] = pd.DataFrame(atom_records).set_index('id')

    # Convert residue data to DataFrame with potential NaN handling
    residue_records = []
    contact_records = []  # Separate DataFrame for contacts
    overlap_records = []  # Separate DataFrame for overlaps
    
    for res_idx, res in enumerate(results['residue_data']):
        # Basic residue record
        record = {
            'residue_id': f"{res['chain']}_{res['resname']}_{res['resid']}",
            'chain': res['chain'],
            'resname': res['resname'],
            'resid': res['resid'],
            'total_sasa': res['total_sasa'],
            'total_area': res['total_area'],
            'n_atoms': res['n_atoms'],
            'center_x': res['center'][0],
            'center_y': res['center'][1],
            'center_z': res['center'][2]
        }
        
        # Add optional fields if present
        for field in ['standard_sasa', 'dsasa']:
            if field in res:
                record[field] = res[field]
                
        residue_records.append(record)
        
        # Process contacts into separate DataFrame
        if res['contacts']:
            for contact_id, contact_info in res['contacts'].items():
                contact_record = {
                    'source_residue': record['residue_id'],
                    'contact_atom_id': contact_id,
                    'contact_area': contact_info['contact_area'],
                    'distance': contact_info['distance'],
                    'struct_type': contact_info['struct_type']
                }
                contact_records.append(contact_record)
        
        # Process overlaps into separate DataFrame
        if res['overlaps']:
            for overlap in res['overlaps']:
                overlap_record = {
                    'source_residue': record['residue_id'],
                    'atoms': ','.join(map(str, overlap['atoms'])),
                    'overlap_area': overlap['overlap_area'],
                    'normalized_area': overlap['normalized_area'],
                    'buried_area': overlap['buried_area']
                }
                overlap_records.append(overlap_record)
    
    # Create main residue DataFrame
    dfs['residues'] = pd.DataFrame(residue_records)
    dfs['residues'].set_index('residue_id', inplace=True)
    
    # Create contacts DataFrame if contacts exist
    if contact_records:
        dfs['contacts'] = pd.DataFrame(contact_records)
        
    # Create overlaps DataFrame if overlaps exist
    if overlap_records:
        dfs['overlaps'] = pd.DataFrame(overlap_records)
    
    # Add residue index if present
    if 'residue_index' in results:
        dfs['residue_index'] = pd.DataFrame.from_dict(
            results['residue_index'], 
            orient='index',
            columns=['index']
        )
    
    return dfs

def plot_results(summary_df: pd.DataFrame, output_dir: Path):
    """Create visualization plots for benchmark results."""
    # Set style
    plt.style.use('seaborn')
    
    # Plot distributions
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    
    # BSA distribution
    sns.histplot(data=summary_df, x='chain_A_BSA', ax=axes[0,0])
    axes[0,0].set_title('Distribution of Chain A BSA')
    
    # Interface area distribution
    sns.histplot(data=summary_df, x='interface_area', ax=axes[0,1])
    axes[0,1].set_title('Distribution of Interface Areas')
    
    # Calculation time vs structure size
    sns.scatterplot(
        data=summary_df,
        x='total_complex_surface',
        y='calculation_time',
        ax=axes[1,0]
    )
    axes[1,0].set_title('Calculation Time vs Structure Size')
    
    # If differences present, plot error distribution
    if 'chain_A_BSA_error_per_atom' in summary_df.columns:
        sns.boxplot(
            data=summary_df.melt(
                value_vars=[col for col in summary_df.columns if col.endswith('_error_per_atom')]
            ),
            x='variable',
            y='value',
            ax=axes[1,1]
        )
        axes[1,1].set_title('Error per Atom Distribution')
        axes[1,1].set_xticklabels(axes[1,1].get_xticklabels(), rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'benchmark_summary.png', dpi=300, bbox_inches='tight')
    plt.close()


def analyze_dataframes(dfs: Dict[str, pd.DataFrame]) -> Dict[str, Any]:
    analysis = {}
    
    # Atom-level analysis
    if 'atoms' in dfs:
        atom_df = dfs['atoms']
        analysis['atoms'] = {
            'total_atoms': len(atom_df),
            'average_sasa': atom_df['sasa'].mean(),
            'total_sasa': atom_df['sasa'].sum(),
            'polar_atoms': atom_df['polar'].sum(),
            'chains': atom_df['chain'].unique().tolist(),
            'residues_count': len(atom_df.groupby(['chain', 'resid']).count())
        }
    
    # Residue-level analysis
    if 'residues' in dfs:
        res_df = dfs['residues']
        analysis['residues'] = {
            'total_residues': len(res_df),
            'average_sasa': res_df['total_sasa'].mean(),
            'total_sasa': res_df['total_sasa'].sum(),
            'chains': res_df['chain'].unique().tolist(),
            'residue_types': res_df['resname'].value_counts().to_dict()
        }
        
        if 'dsasa' in res_df.columns:
            analysis['residues']['total_dsasa'] = res_df['dsasa'].sum()
            analysis['residues']['interface_residues'] = (res_df['dsasa'] > 0).sum()
    
    # Contact analysis
    if 'contacts' in dfs:
        contact_df = dfs['contacts']
        analysis['contacts'] = {
            'total_contacts': len(contact_df),
            'average_contact_area': contact_df['contact_area'].mean(),
            'total_contact_area': contact_df['contact_area'].sum(),
            'average_distance': contact_df['distance'].mean()
        }
    
    # Overlap analysis
    if 'overlaps' in dfs:
        overlap_df = dfs['overlaps']
        analysis['overlaps'] = {
            'total_overlaps': len(overlap_df),
            'average_overlap_area': overlap_df['overlap_area'].mean(),
            'total_buried_area': overlap_df['buried_area'].sum(),
            'average_normalized_area': overlap_df['normalized_area'].mean()
        }
    
    return analysis