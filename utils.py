import pandas as pd
import numpy as np
from typing import Dict, Any

def convert_to_dataframes(results: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    dfs = {}
    
    # 1. Atom-level information DataFrame
    atom_info = results['atom_info']
    coords = np.array(atom_info['coordinates']).reshape(-1, 3)
    
    dfs['atoms'] = pd.DataFrame({
        'name': atom_info['names'],
        'element': atom_info['elements'],
        'mol_type': atom_info['mol_types'],
        'x': coords[:, 0],
        'y': coords[:, 1],
        'z': coords[:, 2],
        'radius': atom_info['radii'],
        'occupancy': atom_info['occupancies'],
        'b_factor': atom_info['b_factors'],
        'charge': atom_info['charges'],
        'is_hetatm': atom_info['is_hetatm'],
        'atom_type': atom_info['atom_types'],
        'sasa': results['sasa']['values'],
        'delta_sasa': results['sasa']['delta'],
        'relative_sasa': results['sasa']['relative']
    })
    
    # 2. SASA by molecular type
    dfs['sasa_by_type'] = pd.DataFrame([
        {'mol_type': k, 'sasa': v}
        for k, v in results['sasa']['by_mol_type'].items()
    ])
    
    # 3. Interface analysis
    # 3.1 Interface by chain
    dfs['interface_by_chain'] = pd.DataFrame([
        {'chain': k, 'interface_area': v}
        for k, v in results['interface']['by_chain'].items()
    ])
    
    # 3.2 Interface atoms
    interface_atoms = []
    for chain, atom_indices in results['interface']['atoms'].items():
        for idx in atom_indices:
            interface_atoms.append({
                'chain': chain,
                'atom_index': idx,
                'atom_name': dfs['atoms'].iloc[idx]['name'],
                'mol_type': dfs['atoms'].iloc[idx]['mol_type'],
                'sasa': dfs['atoms'].iloc[idx]['sasa'],
                'delta_sasa': dfs['atoms'].iloc[idx]['delta_sasa']
            })
    dfs['interface_atoms'] = pd.DataFrame(interface_atoms)
    
    # 4. Molecular analysis
    # 4.1 Atoms and residues by type
    mol_composition = pd.DataFrame([
        {
            'mol_type': k,
            'atom_count': results['molecular']['atoms_by_type'][k],
            'residue_count': results['molecular']['residues_by_type'][k]
        }
        for k in results['molecular']['atoms_by_type'].keys()
    ])
    dfs['molecular_composition'] = mol_composition
    
    # 4.2 Type contacts
    type_contacts = []
    for (type1, type2), contact_value in results['molecular']['type_contacts'].items():
        type_contacts.append({
            'type1': type1,
            'type2': type2,
            'contacts': contact_value
        })
    dfs['type_contacts'] = pd.DataFrame(type_contacts)
    
    # 5. Interactions (if available)
    if 'interactions' in results:
        # 5.1 Interactions by molecular type
        interactions = []
        for (type1, type2), pairs in results['interactions']['by_type'].items():
            for atom1_idx, atom2_idx in pairs:
                interactions.append({
                    'type1': type1,
                    'type2': type2,
                    'atom1_idx': atom1_idx,
                    'atom2_idx': atom2_idx,
                    'atom1_name': dfs['atoms'].iloc[atom1_idx]['name'],
                    'atom2_name': dfs['atoms'].iloc[atom2_idx]['name'],
                    'energy': results['interactions']['energies'][len(interactions)]
                })
        dfs['interactions'] = pd.DataFrame(interactions)
    
    # Add summary statistics
    dfs['summary'] = pd.DataFrame([{
        'total_atoms': len(dfs['atoms']),
        'total_interface_area': results['interface']['total_area'],
        'buried_surface_area': results['interface']['buried_surface_area'],
        'total_sasa': dfs['atoms']['sasa'].sum(),
        'mean_sasa': dfs['atoms']['sasa'].mean(),
        'interface_atoms': len(dfs['interface_atoms']),
        'molecular_types': len(dfs['molecular_composition']),
        'total_interactions': len(dfs['interactions']) if 'interactions' in dfs else 0
    }])
    
    return dfs

def get_analysis_summary(dfs: Dict[str, pd.DataFrame]) -> str:
    summary = dfs['summary'].iloc[0]
    mol_comp = dfs['molecular_composition']
    
    analysis = f"""Analysis Summary:
    
Structure Overview:
- Total atoms: {summary['total_atoms']:,}
- Molecular types: {summary['molecular_types']}
- HETATM atoms: {dfs['atoms']['is_hetatm'].sum():,}

Molecular Composition:
{mol_comp.to_string(index=False)}

Surface Analysis:
- Total SASA: {summary['total_sasa']:.1f} Å²
- Mean atom SASA: {summary['mean_sasa']:.1f} Å²
- Total interface area: {summary['total_interface_area']:.1f} Å²
- Buried surface area: {summary['buried_surface_area']:.1f} Å²
- Interface atoms: {summary['interface_atoms']:,}

Interface by Chain:
{dfs['interface_by_chain'].to_string(index=False)}

Type Contacts:
{dfs['type_contacts'].to_string(index=False)}
"""
    
    if 'interactions' in dfs:
        analysis += f"\nInteractions:\n"
        analysis += f"- Total interactions: {summary['total_interactions']:,}\n"
        analysis += f"- Mean interaction energy: {dfs['interactions']['energy'].mean():.2f}\n"
    
    return analysis
