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
            - 'atoms': Atom-level data with surface and contact information
            - 'residues': Residue-level data with detailed surface analysis
            - 'chains': Chain-level data with type and aggregate information
            - 'contacts': Contact information between atoms
            - 'overlaps': Overlap groups and their properties
    """
    dfs = {}
    
    # Convert atom data to DataFrame
    atom_records = []
    for atom_id, atom_info in results['atoms'].items():
        # Get contact and overlap data safely
        contacts = atom_info.get('contacts', {})
        nonbonded = contacts.get('nonbonded', {})
        overlap_groups = contacts.get('overlap_groups', [])
        
        # Calculate total contact area
        total_contact_area = sum(contact.get('area', 0.0) for contact in nonbonded.values())
        
        # Calculate total overlap area
        total_overlap_area = sum(overlap.get('area', 0.0) for overlap in overlap_groups)
        
        record = {
            'id': int(atom_id),
            # Basic Properties
            'name': atom_info['name'],
            'resname': atom_info['resname'],
            'chain': atom_info['chain'],
            'resid': atom_info['resid'],
            'index': atom_info['index'],
            'x': atom_info['coords'][0],
            'y': atom_info['coords'][1],
            'z': atom_info['coords'][2],
            
            # Surface Analysis
            'sphere_area': atom_info['surface']['sphere_area'],
            'sasa': atom_info['surface']['sasa'],
            'buried_area': atom_info['surface']['buried_area'],
            'contact_area': atom_info['surface']['contact_area'],
            'dsasa': atom_info['surface']['dsasa'],
            
            # Properties
            'vdw': atom_info['properties']['vdw'],
            'polar': atom_info['properties']['polar'],
            'charge': atom_info['properties']['charge'],
            'struct_type': atom_info['properties']['struct_type'],
            
            # Contact metrics
            'n_contacts': len(nonbonded),
            'n_overlaps': len(overlap_groups),
            'total_contact_area': total_contact_area,
            'total_overlap_area': total_overlap_area
        }
        atom_records.append(record)
    
    dfs['atoms'] = pd.DataFrame(atom_records).set_index('id')

    # Convert residue data to DataFrame
    residue_records = []
    for residue_id, res in results['residues'].items():
        contacts = res.get('contacts', {})
        overlaps = res.get('overlaps', [])
        
        # Calculate total contact and overlap areas for residue
        total_contact_area = sum(contact.get('contact_area', 0.0) for contact in contacts.values())
        total_overlap_area = sum(overlap.get('overlap_area', 0.0) for overlap in overlaps)
        
        record = {
            'residue_id': residue_id,
            # Identifiers
            'chain': res['identifiers']['chain'],
            'resname': res['identifiers']['name'],
            'resid': res['identifiers']['number'],
            
            # Structure
            'n_atoms': res['structure']['n_atoms'],
            'center_x': res['structure']['center'][0],
            'center_y': res['structure']['center'][1],
            'center_z': res['structure']['center'][2],
            
            # Surface Areas
            'total_sasa': res['surface']['total_sasa'],
            'total_area': res['surface']['total_area'],
            'standard_sasa': res['surface']['standard_sasa'],
            'dsasa': res['surface']['dsasa'],
            
            # Contact metrics
            'n_contacts': len(contacts),
            'n_overlaps': len(overlaps),
            'total_contact_area': total_contact_area,
            'total_overlap_area': total_overlap_area
        }
        residue_records.append(record)
    
    dfs['residues'] = pd.DataFrame(residue_records).set_index('residue_id')

    # Convert chain data to DataFrame with more information
    chain_records = []
    for chain_id, chain in results['chains'].items():
        # Get all residues for this chain from residues DataFrame
        chain_residues = [r for r in residue_records if r['chain'] == chain_id]
        
        record = {
            'chain_id': chain_id,
            'type': chain['type'],
            'n_residues': len(chain['residues']),
            'n_atoms': sum(r['n_atoms'] for r in chain_residues),
            'total_sasa': sum(r['total_sasa'] for r in chain_residues),
            'total_area': sum(r['total_area'] for r in chain_residues),
            'total_dsasa': sum(r['dsasa'] for r in chain_residues),
            'total_contacts': sum(r['n_contacts'] for r in chain_residues),
            'total_overlaps': sum(r['n_overlaps'] for r in chain_residues),
            'total_contact_area': sum(r['total_contact_area'] for r in chain_residues),
            'total_overlap_area': sum(r['total_overlap_area'] for r in chain_residues)
        }
        chain_records.append(record)
    
    dfs['chains'] = pd.DataFrame(chain_records).set_index('chain_id')
    
    # Create detailed contact network DataFrame
    contact_records = []
    for atom_id, atom_info in results['atoms'].items():
        contacts = atom_info.get('contacts', {}).get('nonbonded', {})
        for contact_id, contact_data in contacts.items():
            record = {
                'atom_id': int(atom_id),
                'contact_id': int(contact_id),
                'area': contact_data['area'],
                'distance': contact_data.get('distance', float('nan'))
            }
            contact_records.append(record)
    
    if contact_records:
        dfs['contacts'] = pd.DataFrame(contact_records)

    # Create detailed overlap groups DataFrame
    overlap_records = []
    for atom_id, atom_info in results['atoms'].items():
        overlap_groups = atom_info.get('contacts', {}).get('overlap_groups', [])
        for i, overlap in enumerate(overlap_groups):
            record = {
                'atom_id': int(atom_id),
                'overlap_index': i,
                'atoms': ','.join(map(str, overlap['atoms'])),
                'area': overlap['area'],
                'normalized_area': overlap['normalized_area'],
                'buried_area': overlap['buried_area']
            }
            overlap_records.append(record)
    
    if overlap_records:
        dfs['overlaps'] = pd.DataFrame(overlap_records)

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
