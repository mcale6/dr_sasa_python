import pandas as pd
import numpy as np
from typing import Dict, Any


def convert_to_dataframes(results: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    dfs = {}
    
    # Convert atom data to DataFrame
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
            'polar': atom_info['polar'],
            'charge': atom_info['charge']
        }
        atom_records.append(record)
    
    dfs['atoms'] = pd.DataFrame(atom_records).set_index('id')

    # Convert residue data to DataFrame
    residue_records = []
    for res in results['residue_data']:
        record = {
            'chain': res['chain'],
            'resname': res['resname'],
            'resid': res['resid'],
            'total_sasa': res['total_sasa'],
            'total_area': res['total_area'],
            'standard_sasa': res['standard_sasa'],
            'dsasa': res['dsasa'],
            'n_atoms': res['n_atoms'],
            'center_x': res['center'][0],
            'center_y': res['center'][1],
            'center_z': res['center'][2],
           # 'buried_area': res['buried_area']
        }
        
        # Add contacts if present
        if res['contacts']:
            for contact_id, contact_info in res['contacts'].items():
                record[f'contact_{contact_id}_area'] = contact_info['contact_area']
                record[f'contact_{contact_id}_distance'] = contact_info['distance']
                record[f'contact_{contact_id}_type'] = contact_info['struct_type']
        
        # Add overlaps if present
        if res['overlaps']:
            for i, overlap in enumerate(res['overlaps']):
                record[f'overlap_{i}_residues'] = ','.join(map(str, overlap['residues']))
                record[f'overlap_{i}_area'] = overlap['overlap_area']
                record[f'overlap_{i}_norm_area'] = overlap['normalized_area']
                record[f'overlap_{i}_buried'] = overlap['buried_area']
        
        residue_records.append(record)
    
    dfs['residues'] = pd.DataFrame(residue_records)

    # Handle matrix data if present
    if 'interaction_matrices' in results:
        for key, matrix_data in results['interaction_matrices'].items():
            matrix_df = pd.DataFrame(
                matrix_data['matrix'],
                index=matrix_data['row_labels'],
                columns=matrix_data['col_labels']
            )
            dfs[f'interaction_matrix_{key}'] = matrix_df

    if 'intra_matrices' in results:
        atom_matrix = results['intra_matrices']['atom_matrix']
        residue_matrix = results['intra_matrices']['residue_matrix']
        
        dfs['intra_atom_matrix'] = pd.DataFrame(
            atom_matrix,
            index=results['intra_matrices']['atom_labels'],
            columns=results['intra_matrices']['atom_labels']
        )
        
        dfs['intra_residue_matrix'] = pd.DataFrame(
            residue_matrix,
            index=results['intra_matrices']['residue_labels'],
            columns=results['intra_matrices']['residue_labels']
        )

    return dfs