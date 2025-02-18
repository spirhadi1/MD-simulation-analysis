#Written by Somayeh Pirhadi, Schiffer Lab, 2023, UMass Chan Medical School
'''This code calculates the pairwise vdw and electrostatic interactions within a range of residues. It outputs 2 csv files called vdw-interactions.csv and elec-interactions.csv. The trajectory must have box information not removed. 
#provide topology, trajectory, first and last residue number as input (Usage: ./nonbonded.sh <topology_file> <trajectory_file> <start_res> <end_res>). for example: ./nonbonded.sh comp.prmtop 11md-2.dcd 1 5
#To analyze a single pdb structure, you need to prepare comp.prmtop and comp.inpcrd out of it : ./nonbonded.sh comp.prmtop comp.inpcrd 1 5'''

#!/bin/bash

# Check if all required arguments are provided
if [ $# -ne 4 ]; then
    echo "Usage: $0 <topology_file> <trajectory_file> <start_res> <end_res>"
    exit 1
fi

# Read command-line arguments
topology_file=$1
trajectory_file=$2
start_res=$3
end_res=$4

# Initialize arrays to hold average values
declare -A avg_elec_values
declare -A avg_vdw_values

# Generate cpptraj commands
cpptraj_commands=""
for ((i=$start_res; i<=$end_res-1; i++)); do
  for ((j=$i+1; j<=$end_res; j++)); do
    lie_outfile="lie-${i}-${j}.txt"
    cpptraj_commands+="lie :$i :$j out $lie_outfile\n"
  done
done

# Run cpptraj with generated commands
echo -e "parm $topology_file\ntrajin $trajectory_file\n$cpptraj_commands" | cpptraj

# Loop through residue pairs
# Loop through residue pairs
for ((i=$start_res; i<=$end_res-1; i++)); do
  for ((j=$i+1; j<=$end_res; j++)); do
    lie_outfile="lie-${i}-${j}.txt"
    
    # Calculate the averages for elec and vdw columns
    avg_elec=$(awk 'NR > 1 { elec_sum += $2 } END { if (NR > 1) print elec_sum/(NR-1); else print 0 }' "$lie_outfile")
    avg_vdw=$(awk 'NR > 1 { vdw_sum += $3 } END { if (NR > 1) print vdw_sum/(NR-1); else print 0 }' "$lie_outfile")
    
    # Store average values in arrays
    avg_elec_values["$i-$j"]=$avg_elec
    avg_vdw_values["$i-$j"]=$avg_vdw
    
    # Remove the intermediate file
    rm "$lie_outfile"
  done
done


# Print CSV header for electrostatic interactions
elec_header=""
for ((i=$start_res; i<=$end_res-1; i++)); do
  for ((j=$i+1; j<=$end_res; j++)); do
    elec_header+="e${i}-e${j},"
  done
done
echo "$elec_header" > elec-interactions.csv

# Print CSV header for VDW interactions
vdw_header=""
for ((i=$start_res; i<=$end_res-1; i++)); do
  for ((j=$i+1; j<=$end_res; j++)); do
    vdw_header+="v${i}-v${j},"
  done
done
echo "$vdw_header" > vdw-interactions.csv

# Print values for electrostatic interactions
elec_row=""
for ((i=$start_res; i<=$end_res-1; i++)); do
  for ((j=$i+1; j<=$end_res; j++)); do
    elec_row+="${avg_elec_values["$i-$j"]},"
  done
done
echo "$elec_row" >> elec-interactions.csv

# Print values for VDW interactions
vdw_row=""
for ((i=$start_res; i<=$end_res-1; i++)); do
  for ((j=$i+1; j<=$end_res; j++)); do
    vdw_row+="${avg_vdw_values["$i-$j"]},"
  done
done
echo "$vdw_row" >> vdw-interactions.csv
echo " Files vdw-interactions.csv and elec-interactions.csv have been written successfully"
echo ""

