import pandas as pd
import numpy as np
from typing import Dict, Any


def convert_to_dataframes(results: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    dfs = {}
    atom_info = results['atoms']
    coords = np.array(atom_info['coordinates'])
    # Base atom dataframe
    base_df = pd.DataFrame({
        'structure': atom_info['structure'],
        'id': atom_info['ids'],
        'name': atom_info['names'],
        'residue': atom_info['residues'],
        'chain': atom_info['chains'],
        'element': atom_info['elements'],
        'struct_type': atom_info['struct_types'],
        'mol_type': atom_info['mol_types'],
        'altloc': atom_info['altlocs'],
        'icode': atom_info['icodes'],
        'residue_number': atom_info['residue_numbers'],
        'is_hetatm': atom_info['is_hetatm'],
        'is_active': atom_info['is_active'],
        'is_terminal': atom_info['is_terminal'],
        'atom_type': atom_info['atom_types'],
        'atom_type_40': atom_info['atom_types_40'],
        'int_num': atom_info['int_nums'],
        'ef_int_num': atom_info['ef_int_nums'],
        'dtype': atom_info['dtypes'],
        'polarity': atom_info['polarity'],
        'x': coords[:,0],
        'y': coords[:,1],
        'z': coords[:,2],
        'radius': atom_info['radii'],
        'radius2': atom_info['radii2'],
        'vdw': atom_info['vdw'],
        'occupancy': atom_info['occupancies'],
        'b_factor': atom_info['b_factors'],
        'charge': atom_info['charges']
    })
    
    # Surface area dataframe
    dfs['surface'] = pd.DataFrame({
        'area': results['surface']['area'],
        'sasa': results['surface']['sasa'],
        'buried_area': results['surface']['buried_area'],
        'fast_dsasa': results['surface']['fast_dsasa'],
        'accessible': results['surface']['accessible']
    })

    # Interaction data
    interaction_data = []
    for atom_idx, (interacting_atoms, sasa_atoms, contact_map) in enumerate(zip(
            results['interactions']['atoms'],
            results['interactions']['sasa_atoms'], 
            results['interactions']['contact_areas'])):
        
        # Add each interaction for this atom
        for interacting_atom in interacting_atoms:
            interaction_data.append({
                'atom_index': atom_idx,
                'interacting_atom': interacting_atom,
                'is_sasa_interaction': interacting_atom in sasa_atoms,
                'contact_area': contact_map.get(interacting_atom, 0.0)
            })
    
    if interaction_data:
        dfs['interactions'] = pd.DataFrame(interaction_data)

    # Matrix handling if present (using ChainPair keys)
    if 'matrices' in results:
        for chain_pair, matrix in results['matrices']['atom_matrix'].items():
            matrix_id = f"matrix_{chain_pair[0]}_{chain_pair[1]}"
            dfs[matrix_id] = pd.DataFrame(
                matrix,
                index=results['matrices']['row_atoms'][chain_pair],
                columns=results['matrices']['col_atoms'][chain_pair]
            )

    return dfs

def get_analysis_summary(dfs: Dict[str, pd.DataFrame]) -> str:
    atoms_df = dfs['atoms']
    
    summary = f"""Analysis Summary:

Structure Overview:
- Total atoms: {len(atoms_df):,}
- Molecular types: {', '.join(atoms_df['mol_type'].unique())}
- Chains: {', '.join(atoms_df['chain'].unique())}

Surface Analysis:
- Total area: {atoms_df['area'].sum():.1f} Å²
- Total SASA: {atoms_df['sasa'].sum():.1f} Å²
- Total buried area: {atoms_df['buried_area'].sum():.1f} Å²
- Mean accessibility: {atoms_df['accessible'].mean():.1%}

Energy Analysis:
- Total energy: {atoms_df['energy'].sum():.1f}
- Mean DD energy: {atoms_df['dd_energy'].mean():.2f}
- Mean BSA energy: {atoms_df['bsa_energy'].mean():.2f}

Interaction Statistics:
- Interacting atoms: {len(dfs['interactions']):,}
- Mean contact area: {dfs['interactions']['contact_area'].mean():.1f} Å²
- Total overlap groups: {len(dfs['overlaps']):,}
- Mean overlap area: {dfs['overlaps']['overlap_area'].mean():.1f} Å²

Surface Analysis by Molecular Type:
{atoms_df.groupby('mol_type').agg({
    'sasa': 'sum',
    'buried_area': 'sum',
    'accessible': 'mean'
}).round(1)}
"""
    return summary