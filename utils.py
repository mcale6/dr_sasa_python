import pandas as pd
import numpy as np
from typing import Dict, Any


def convert_to_dataframes(results: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    dfs = {}
    
    # Base atom dataframe
    atom_info = results['atoms']
    coords = np.array(atom_info['coordinates'])
    base_df = pd.DataFrame({
        'id': atom_info['ids'],
        'residue': atom_info['residues'],
        'chain': atom_info['chains'],
        'element': atom_info['elements'],
        'mol_type': atom_info['mol_types'],
        'residue_number': atom_info['residue_numbers'],
        'is_hetatm': atom_info['is_hetatm'],
        'atom_type': atom_info['atom_types'],
        'polarity': atom_info['polarity'],
        'x': coords[:,0],
        'y': coords[:,1],
        'z': coords[:,2],
        'radius': atom_info['radii'],
        'radius2': atom_info['radii2'],
        'vdw': atom_info['vdw'],
        'b_factor': atom_info['b_factors']
    })
    dfs['atoms'] = base_df
    
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
        
        for interacting_atom in interacting_atoms:
            interaction_data.append({
                'atom_index': atom_idx,
                'interacting_atom': interacting_atom,
                'is_sasa_interaction': interacting_atom in sasa_atoms,
                'contact_area': contact_map.get(interacting_atom, 0.0)
            })
    
    if interaction_data:
        dfs['interactions'] = pd.DataFrame(interaction_data)

    # Matrix handling with new format
    if 'matrices' in results:
        # Atom level matrix
        atom_matrix = results['matrices']['atom']
        dfs['atom_matrix'] = pd.DataFrame(
            atom_matrix['values'],
            index=atom_matrix['row_labels'],
            columns=atom_matrix['col_labels']
        )
        
        # Residue level matrix
        res_matrix = results['matrices']['residue']
        dfs['residue_matrix'] = pd.DataFrame(
            res_matrix['values'],
            index=res_matrix['row_labels'],
            columns=res_matrix['col_labels']
        )

    # Add residue interaction DataFrame
    if 'residue_interactions' in results:
        ri = results['residue_interactions']
        dfs['residue_interactions'] = pd.DataFrame({
            'residue1': ri['residue1'],
            'residue2': ri['residue2'],
            'buried_area': ri['buried_area'],
            'dsasa_res1': ri['dsasa_res1'],
            'dsasa_res2': ri['dsasa_res2']
        })
        
    # Surface summary if present
    if 'surface_summary' in results:
        summary_data = []
        for obj_id, summary in results['surface_summary'].items():
            summary_data.append({
                'object_id': obj_id,
                'complex_surface': summary['complex_surface'],
                'buried_by_other': summary['buried_by_other'],     # A<---B
                'buries_in_other': summary['buries_in_other'],     # A--->B
                'interface_area': summary['interface_area']
            })
        dfs['surface_summary'] = pd.DataFrame(summary_data)

    return dfs
