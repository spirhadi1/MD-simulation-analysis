import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
import pandas as pd
import time
import argparse

# Argument parsing
parser = argparse.ArgumentParser(description='Compute dihedral angles from a trajectory.')
parser.add_argument('topology', type=str, help='Path to the topology file')
parser.add_argument('trajectory', type=str, help='Path to the trajectory file')
parser.add_argument('--output', type=str, default='dihedral_angles.pkl', help='Output file name (default: dihedral_angles.pkl)')
args = parser.parse_args()

start_time = time.time()

# Load the universe
u = mda.Universe(args.topology, args.trajectory)
residues = u.select_atoms('protein').residues

# Function to compute dihedral angles and convert to sin and cos
def compute_dihedral_angles(dihedral_selections, angle_name):
    angles_data = {}
    dihedral_angles = [selection for selection in dihedral_selections if selection is not None]
    
    if dihedral_angles:
        dihedral_results = Dihedral(dihedral_angles).run()
        for idx, selection in enumerate(dihedral_selections):
            if selection is not None and idx < dihedral_results.angles.shape[1]:
                resid = selection.atoms[0].residue.resid
                segid = selection.atoms[0].residue.segid
                resname = selection.atoms[0].residue.resname
                angles = dihedral_results.angles[:, idx]
                
                column_sin = f"{segid}:{resname}:{resid}:{angle_name}_sin"
                column_cos = f"{segid}:{resname}:{resid}:{angle_name}_cos"
                angles_data[column_sin] = np.sin(np.radians(angles))
                angles_data[column_cos] = np.cos(np.radians(angles))
    
    return angles_data

# Compute sin and cos for phi, psi, omega, and chi1 angles
phi_angles = compute_dihedral_angles([res.phi_selection() for res in residues], "Phi")
psi_angles = compute_dihedral_angles([res.psi_selection() for res in residues], "Psi")
omega_angles = compute_dihedral_angles([res.omega_selection() for res in residues], "Omega")
chi1_angles = compute_dihedral_angles([res.chi1_selection() for res in residues], "Chi1")

# Function to compute custom dihedrals
def compute_custom_dihedrals(selection_dict, angle_name):
    angles_data = {}
    dihedral_angles = []
    residue_indices = []
    
    for res in residues:
        resname = res.resname
        segid = res.segid
        atoms = selection_dict.get(resname)
        if atoms:
            selected_atoms = res.atoms.select_atoms(f'name {" ".join(atoms)}')
            if len(selected_atoms) == 4:
                dihedral_angles.append(selected_atoms)
                residue_indices.append((segid, resname, res.resid))
    
    if dihedral_angles:
        dihedral_results = Dihedral(dihedral_angles).run()
        for idx, (segid, resname, resid) in enumerate(residue_indices):
            if idx < dihedral_results.angles.shape[1]:
                angles = dihedral_results.angles[:, idx]
                column_sin = f"{segid}:{resname}:{resid}:{angle_name}_sin"
                column_cos = f"{segid}:{resname}:{resid}:{angle_name}_cos"
                angles_data[column_sin] = np.sin(np.radians(angles))
                angles_data[column_cos] = np.cos(np.radians(angles))
    
    return angles_data

# Dictionaries for chi2, chi3, chi4, and chi5 selections
chi2_dict = {
    'ARG': ['CA', 'CB', 'CG', 'CD'], 'ASN': ['CA', 'CB', 'CG', 'OD1'],
    'ASP': ['CA', 'CB', 'CG', 'OD1'], 'GLN': ['CA', 'CB', 'CG', 'CD'],
    'GLU': ['CA', 'CB', 'CG', 'CD'], 'HIS': ['CA', 'CB', 'CG', 'ND1'],
    'ILE': ['CA', 'CB', 'CG1', 'CD1'], 'LEU': ['CA', 'CB', 'CG', 'CD1'],
    'LYS': ['CA', 'CB', 'CG', 'CD'], 'MET': ['CA', 'CB', 'CG', 'SD'],
    'PHE': ['CA', 'CB', 'CG', 'CD1'], 'PRO': ['CA', 'CB', 'CG', 'CD'],
    'TRP': ['CA', 'CB', 'CG', 'CD1'], 'TYR': ['CA', 'CB', 'CG', 'CD1']
}

chi3_dict = {
    'ARG': ['CB', 'CG', 'CD', 'NE'], 'GLN': ['CB', 'CG', 'CD', 'OE1'],
    'GLU': ['CB', 'CG', 'CD', 'OE1'], 'LYS': ['CB', 'CG', 'CD', 'CE'],
    'MET': ['CB', 'CG', 'SD', 'CE']
}

chi4_dict = {
    'ARG': ['CG', 'CD', 'NE', 'CZ'], 'LYS': ['CG', 'CD', 'CE', 'NZ']
}

chi5_dict = {
    'ARG': ['CD', 'NE', 'CZ', 'NH1']
}

# Compute sin and cos for custom dihedrals
chi2_angles = compute_custom_dihedrals(chi2_dict, "Chi2")
chi3_angles = compute_custom_dihedrals(chi3_dict, "Chi3")
chi4_angles = compute_custom_dihedrals(chi4_dict, "Chi4")
chi5_angles = compute_custom_dihedrals(chi5_dict, "Chi5")

# Combine all angle data
all_angles = {**phi_angles, **psi_angles, **omega_angles, **chi1_angles,
              **chi2_angles, **chi3_angles, **chi4_angles, **chi5_angles}

# Convert to DataFrame
df = pd.DataFrame(all_angles)

# Save to file
df.to_pickle(args.output)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time:.2f} seconds")
print(f"Dihedral angles saved to {args.output}")

