#!/usr/bin/env python

import argparse
import pytraj as pt
import numpy as np
import pandas as pd

# Function to determine protein and ligand residue range
def get_protein_and_ligand_residues(top_file):
    top = pt.load_topology(top_file)
    
    # Define a list of standard amino acid residue names
    protein_residue_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                             'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                             'THR', 'TRP', 'TYR', 'VAL']
    
    # Define a list of common ion residue names
    ion_residue_names = ['Na+', 'Cl-', 'K+', 'Mg+', 'Ca+', 'Zn+', 'Fe+', 'Fe2+', 'Fe3+']
    
    # Identify protein residues
    protein_residues = [res for res in top.residues if res.name in protein_residue_names]
    
    # Identify ligand residues (non-protein, non-ion, and non-water)
    ligand_residues = [
        res for res in top.residues 
        if res.name not in protein_residue_names 
        and res.name not in ion_residue_names 
        and res.name != 'WAT'
    ]
    
    if not protein_residues:
        raise ValueError("No protein residues found in topology file.")
    
    # Use the 'index' attribute to get residue numbers
    first_res = protein_residues[0].index + 1  # Residue indices start from 0, so add 1
    last_res = protein_residues[-1].index + 1
    protein_range = f":{first_res}-{last_res}"
    
    ligand_range = ""
    if ligand_residues:
        ligand_indices = [str(res.index + 1) for res in ligand_residues]  # Add 1 to match residue numbers
        ligand_range = ",".join(ligand_indices)
    
    selection = protein_range if not ligand_range else f"{protein_range},{ligand_range}"
    return selection

# Function to calculate hydrogen bond frequency
def calculate_hbond_frequency(traj_file, top_file):
    selection = get_protein_and_ligand_residues(top_file)
    
    # Load every 10th frame of the trajectory
    traj = pt.load(traj_file, top=top_file, frame_indices=range(0, pt.iterload(traj_file, top=top_file).n_frames, 10))
    
    # Calculate hydrogen bonds
    data = pt.search_hbonds(traj, solvent_donor=f':WAT@O & around 5 {selection}', solvent_acceptor=f':WAT & around 5 {selection}')
    
    # Process hydrogen bond data
    hbond_matrix = data.values[4:].T.astype(float)  # Exclude the first four rows ('total_solute_hbonds', 'HB00000[UV]', 'HB00000[Bridge]', 'HB00000[ID]')
    hbond_matrix_binary = np.clip(hbond_matrix, 0, 1)  # Turn any number over 1 to 1
    hbond_frequency = np.sum(hbond_matrix_binary, axis=0) / len(traj)
    
    return data, hbond_frequency

# Function to save hydrogen bond frequency data to a CSV file
def save_hbond_frequency_to_csv(data, hbond_frequency, output_file):
    pairs = data.donor_acceptor[0:]
    df = pd.DataFrame([hbond_frequency], columns=pairs)
    df.to_csv(output_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Calculate hydrogen bond frequency for a solute in a trajectory considering interactions with water.")
    parser.add_argument("topology", help="Topology file")
    parser.add_argument("trajectory", help="Trajectory file")
    parser.add_argument("--output", default="hbond_frequency.csv", help="Output CSV filename (default: hbond_frequency.csv)")
    
    args = parser.parse_args()
    
    # Calculate hydrogen bond frequency
    data, hbond_frequency = calculate_hbond_frequency(args.trajectory, args.topology)
    
    # Save the hydrogen bond frequency data to a CSV file
    save_hbond_frequency_to_csv(data, hbond_frequency, args.output)
    print(f"Hydrogen bond frequency data saved to '{args.output}'.")

if __name__ == "__main__":
    main()
