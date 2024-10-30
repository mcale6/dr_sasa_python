import pandas as pd
import numpy as np
from typing import Dict, Any

def convert_to_dataframes(results: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
   dfs = {}
   
   # Base atom dataframe
   atom_info = results['atom_info']
   coords = np.array(atom_info['coordinates'])
   base_df = pd.DataFrame({
       'name': atom_info['names'],
       'element': atom_info['elements'],
       'mol_type': atom_info['mol_types'],
       'x': coords[:,0], 
       'y': coords[:,1],
       'z': coords[:,2],
       'radius': atom_info['radii'],
       'occupancy': atom_info['occupancies'],
       'b_factor': atom_info['b_factors'],
       'charge': atom_info['charges'],
       'is_hetatm': atom_info['is_hetatm'],
       'atom_type': atom_info['atom_types']
   })
   
   # Add SASA data if present
   if 'sasa' in results:
       base_df['sasa'] = results['sasa']['values']
       base_df['delta_sasa'] = results['sasa']['delta']
       base_df['relative_sasa'] = results['sasa']['relative']
       
       dfs['sasa_by_type'] = pd.DataFrame([
           {'mol_type': k, 'sasa': v} 
           for k,v in results['sasa']['by_mol_type'].items()
       ])
       
   dfs['atoms'] = base_df

   # Matrices if present
   if 'matrices' in results:
       matrices = results['matrices']
       if 'inter_molecular' in matrices:
           inter = matrices['inter_molecular']
           for mol_type, matrix in inter['atomic'].items():
               dfs[f'inter_matrix_{mol_type[0]}_{mol_type[1]}'] = pd.DataFrame(
                   matrix,
                   index=inter['row_atoms'][mol_type],
                   columns=inter['col_atoms'][mol_type]
               )
               
       if 'intra_molecular' in matrices:
           intra = matrices['intra_molecular']
           n = len(intra['col_atoms'])
           dfs['intra_matrix'] = pd.DataFrame(
               np.array(intra['atomic']).reshape(n,n),
               index=intra['row_atoms'],
               columns=intra['col_atoms']
           )
   
   # Interface data if present
   if 'interface' in results:
       interface = results['interface']
       dfs['interface'] = pd.DataFrame([{
           'chain': chain,
           'atom_index': idx,
           'atom_name': base_df.iloc[idx]['name'],
           'sasa': base_df.iloc[idx].get('sasa', None),
           'delta_sasa': base_df.iloc[idx].get('delta_sasa', None)
       } for chain, indices in interface['atoms'].items()
         for idx in indices])
           
       dfs['interface_summary'] = pd.DataFrame([{
           'total_area': interface['total_area'],
           'buried_area': interface['buried_surface_area'],
           **{f'chain_{k}_area': v for k,v in interface['by_chain'].items()}
       }])

   return dfs

def get_analysis_summary(dfs: Dict[str, pd.DataFrame]) -> str:
    summary = dfs['summary'].iloc[0]
    mol_comp = dfs['molecular_composition']
    
    analysis = f"""Analysis Summary:
    
Structure Overview:
- Total atoms: {summary['total_atoms']:,}
- Molecular types: {summary['molecular_types']}

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
