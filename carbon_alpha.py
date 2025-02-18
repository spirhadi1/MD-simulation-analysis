import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pickle
import time
import argparse

def main():
    # Define the argument parser
    parser = argparse.ArgumentParser(description="Calculate pairwise distances between alpha carbons in a trajectory.")
    parser.add_argument("topology", help="Path to the topology file (e.g., PSF, PDB, etc.)")
    parser.add_argument("trajectory", help="Path to the trajectory file (e.g., DCD, NC, etc.)")
    parser.add_argument(
        "--output", 
        default="pairwise_ca_distances.pkl", 
        help="Output file name (default: pairwise_ca_distances.pkl)"
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Start timing
    start_time = time.time()

    # Load the universe
    print("Loading the universe...")
    u = mda.Universe(args.topology, args.trajectory)

    # Select all alpha carbons
    print("Selecting alpha carbons...")
    alpha_carbons = u.select_atoms('name CA')

    # Get residue information for column names
    print("Extracting residue information...")
    residue_info = [
        f"{res.segid}:{res.resname}:{res.resid}" for res in alpha_carbons.residues
    ]

    # Generate all unique residue pairs for column names
    print("Generating residue pairs...")
    n_residues = len(alpha_carbons)
    pairs = [(i, j) for i in range(n_residues) for j in range(i + 1, n_residues)]
    columns = [f"{residue_info[i]}:{residue_info[j]}" for i, j in pairs]

    # Precompute pair indices for NumPy indexing
    pair_indices = np.array(pairs)

    # Preload trajectory positions
    print("Preloading trajectory positions...")
    positions = np.array([ts.positions[alpha_carbons.indices] for ts in u.trajectory])

    # Preallocate distances array
    print("Allocating distance array...")
    total_frames = len(positions)
    distances = np.empty((total_frames, len(pairs)))

    # Calculate distances for all frames
    print("Calculating distances...")
    for frame_idx in range(total_frames):
        dist_matrix = distance_array(positions[frame_idx], positions[frame_idx])
        distances[frame_idx] = dist_matrix[pair_indices[:, 0], pair_indices[:, 1]]

    # Save distances and column titles to a pickle file
    print(f"Saving results to {args.output}...")
    with open(args.output, "wb") as f:
        pickle.dump({"distances": distances, "columns": columns}, f)

    # Print elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total elapsed time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()

