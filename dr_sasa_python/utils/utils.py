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
        
        # Get surface components safely
        surface = atom_info.get('surface', {})
        backbone = surface.get('backbone', {})
        sidechain = surface.get('sidechain', {})
        groove = surface.get('groove', {})
        ligand = surface.get('ligand', {})
        
        # Calculate total contact area
        total_contact_area = sum(contact.get('area', 0.0) for contact in nonbonded.values())
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
            
            # Basic Surface Analysis
            'sphere_area': surface['sphere_area'],
            'sasa': surface['sasa'],
            'buried_area': surface['buried_area'],
            'contact_area': surface['contact_area'],
            'dsasa': surface['dsasa'],
            'total_asa': surface['total_asa'],
            
            # Backbone Components
            'backbone_total': backbone.get('total', 0.0),
            'backbone_polar': backbone.get('polar', 0.0),
            'backbone_hydrophobic': backbone.get('hydrophobic', 0.0),
            
            # Sidechain Components
            'sidechain_total': sidechain.get('total', 0.0),
            'sidechain_polar': sidechain.get('polar', 0.0),
            'sidechain_hydrophobic': sidechain.get('hydrophobic', 0.0),
            
            # Groove Components
            'major_groove_total': groove.get('major', {}).get('total', 0.0),
            'major_groove_polar': groove.get('major', {}).get('polar', 0.0),
            'major_groove_hydrophobic': groove.get('major', {}).get('hydrophobic', 0.0),
            
            'minor_groove_total': groove.get('minor', {}).get('total', 0.0),
            'minor_groove_polar': groove.get('minor', {}).get('polar', 0.0),
            'minor_groove_hydrophobic': groove.get('minor', {}).get('hydrophobic', 0.0),
            
            'no_groove_total': groove.get('none', {}).get('total', 0.0),
            'no_groove_polar': groove.get('none', {}).get('polar', 0.0),
            'no_groove_hydrophobic': groove.get('none', {}).get('hydrophobic', 0.0),
            
            # Ligand Components
            'ligand_total': ligand.get('total', 0.0),
            'ligand_polar': ligand.get('polar', 0.0),
            'ligand_hydrophobic': ligand.get('hydrophobic', 0.0),
            
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

    # Convert residue data to DataFrame with similar component organization
    residue_records = []
    for residue_id, res in results['residues'].items():
        contacts = res.get('contacts', {})
        overlaps = res.get('overlaps', [])
        surface = res.get('surface', {})
        
        total_contact_area = sum(contact.get('contact_area', 0.0) for contact in contacts.values())
        total_overlap_area = sum(overlap.get('overlap_area', 0.0) for overlap in overlaps)
        
        record = {
            'residue_id': residue_id,
            'chain': res['identifiers']['chain'],
            'resname': res['identifiers']['name'],
            'resid': res['identifiers']['number'],
            
            'n_atoms': res['structure']['n_atoms'],
            'center_x': res['structure']['center'][0],
            'center_y': res['structure']['center'][1],
            'center_z': res['structure']['center'][2],
            
            'total_sasa': surface['total_sasa'],
            'total_area': surface['total_area'],
            'standard_sasa': surface['standard_sasa'],
            'dsasa': surface['dsasa'],
            
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

def plot_results(df_dict: Dict[str, pd.DataFrame], output_dir: Path):
    """Create comprehensive visualization plots for surface analysis results.
    
    Args:
        df_dict: Dictionary of DataFrames containing atoms, residues, chains data
        output_dir: Directory to save plot files
    """
    plt.style.use('seaborn')
    
    # 1. Surface Area Distribution Analysis
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    
    # Surface area distributions by type
    surface_data = df_dict['atoms'].melt(
        value_vars=['backbone_total', 'sidechain_total', 
                   'major_groove_total', 'minor_groove_total'],
        var_name='component', value_name='area'
    )
    sns.boxplot(data=surface_data, x='component', y='area', ax=axes[0,0])
    axes[0,0].set_title('Distribution of Surface Areas by Component')
    axes[0,0].set_xticklabels(axes[0,0].get_xticklabels(), rotation=45)
    
    # Polar vs Hydrophobic distribution
    polar_data = df_dict['atoms'].melt(
        value_vars=['backbone_polar', 'backbone_hydrophobic',
                   'sidechain_polar', 'sidechain_hydrophobic'],
        var_name='component', value_name='area'
    )
    sns.boxplot(data=polar_data, x='component', y='area', ax=axes[0,1])
    axes[0,1].set_title('Polar vs Hydrophobic Surface Distribution')
    axes[0,1].set_xticklabels(axes[0,1].get_xticklabels(), rotation=45)
    
    # SASA vs Buried Area scatter
    sns.scatterplot(
        data=df_dict['atoms'],
        x='sasa',
        y='buried_area',
        hue='struct_type',
        alpha=0.6,
        ax=axes[1,0]
    )
    axes[1,0].set_title('SASA vs Buried Area by Structure Type')
    
    # Contact analysis
    sns.scatterplot(
        data=df_dict['atoms'],
        x='n_contacts',
        y='total_contact_area',
        hue='struct_type',
        size='n_overlaps',
        alpha=0.6,
        ax=axes[1,1]
    )
    axes[1,1].set_title('Contact Analysis')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'surface_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Residue Level Analysis
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    
    # Residue SASA distribution by type
    sns.boxplot(
        data=df_dict['residues'],
        x='resname',
        y='total_sasa',
        ax=axes[0,0]
    )
    axes[0,0].set_title('SASA Distribution by Residue Type')
    axes[0,0].set_xticklabels(axes[0,0].get_xticklabels(), rotation=45)
    
    # Residue Contact Analysis
    sns.scatterplot(
        data=df_dict['residues'],
        x='total_sasa',
        y='total_contact_area',
        hue='chain',
        size='n_contacts',
        alpha=0.6,
        ax=axes[0,1]
    )
    axes[0,1].set_title('Residue Contact Analysis')
    
    # dSASA vs Standard SASA
    sns.scatterplot(
        data=df_dict['residues'],
        x='standard_sasa',
        y='dsasa',
        hue='resname',
        alpha=0.6,
        ax=axes[1,0]
    )
    axes[1,0].set_title('dSASA vs Standard SASA')
    
    # Residue burial analysis
    df_dict['residues']['burial_ratio'] = 1 - (df_dict['residues']['total_sasa'] / 
                                              df_dict['residues']['total_area'])
    sns.boxplot(
        data=df_dict['residues'],
        x='resname',
        y='burial_ratio',
        ax=axes[1,1]
    )
    axes[1,1].set_title('Residue Burial Ratio Distribution')
    axes[1,1].set_xticklabels(axes[1,1].get_xticklabels(), rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'residue_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Chain Level Analysis
    if len(df_dict['chains']) > 1:  # Only if multiple chains
        fig, axes = plt.subplots(2, 2, figsize=(16, 16))
        
        # Chain surface composition
        chain_surface = df_dict['chains'][['type', 'total_sasa', 'total_area', 'total_dsasa']]
        chain_surface_melted = chain_surface.melt(
            id_vars=['type'],
            var_name='metric',
            value_name='area'
        )
        sns.barplot(
            data=chain_surface_melted,
            x='type',
            y='area',
            hue='metric',
            ax=axes[0,0]
        )
        axes[0,0].set_title('Chain Surface Composition')
        
        # Contact distribution by chain
        sns.barplot(
            data=df_dict['chains'],
            x='type',
            y='total_contacts',
            ax=axes[0,1]
        )
        axes[0,1].set_title('Contacts by Chain Type')
        
        # Chain interaction network
        if 'contacts' in df_dict:
            contact_matrix = pd.pivot_table(
                df_dict['contacts'],
                values='area',
                index='atom_id',
                columns='contact_id',
                aggfunc='sum'
            ).fillna(0)
            sns.heatmap(contact_matrix, ax=axes[1,0], cmap='viridis')
            axes[1,0].set_title('Chain Contact Matrix')
        
        # Chain type distribution
        sns.countplot(
            data=df_dict['chains'],
            x='type',
            ax=axes[1,1]
        )
        axes[1,1].set_title('Distribution of Chain Types')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'chain_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()

    # 4. Optional: Groove Analysis for Nucleic Acids
    nucleic_data = df_dict['atoms'][
        df_dict['atoms']['struct_type'].isin(['DNA', 'RNA'])
    ]
    if not nucleic_data.empty:
        fig, axes = plt.subplots(2, 2, figsize=(16, 16))
        
        # Major vs Minor groove distribution
        groove_data = nucleic_data.melt(
            value_vars=['major_groove_total', 'minor_groove_total', 'no_groove_total'],
            var_name='groove',
            value_name='area'
        )
        sns.boxplot(data=groove_data, x='groove', y='area', ax=axes[0,0])
        axes[0,0].set_title('Major vs Minor Groove Surface Distribution')
        
        # Groove polarity analysis
        groove_polarity = nucleic_data.melt(
            value_vars=['major_groove_polar', 'major_groove_hydrophobic',
                       'minor_groove_polar', 'minor_groove_hydrophobic'],
            var_name='component',
            value_name='area'
        )
        sns.boxplot(data=groove_polarity, x='component', y='area', ax=axes[0,1])
        axes[0,1].set_title('Groove Polarity Analysis')
        
        # Groove contact analysis
        if 'contacts' in df_dict:
            groove_contacts = df_dict['contacts'][
                df_dict['contacts']['atom_id'].isin(nucleic_data.index)
            ]
            sns.scatterplot(
                data=groove_contacts,
                x='distance',
                y='area',
                ax=axes[1,0]
            )
            axes[1,0].set_title('Groove Contact Distance vs Area')
        
        # Base-specific groove analysis
        base_groove = nucleic_data.groupby('resname')[
            ['major_groove_total', 'minor_groove_total']
        ].mean()
        base_groove.plot(kind='bar', ax=axes[1,1])
        axes[1,1].set_title('Average Groove Area by Base')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'groove_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()