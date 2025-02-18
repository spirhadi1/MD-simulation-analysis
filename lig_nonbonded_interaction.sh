#Written by Somayeh Pirhadi, Schiffer Lab, 2023, UMass Chan Medical School
'''The code calcultes the electrostatic and vdw interactions of a ligand with residues whithin 6 angstrom. The trajectory must have box information not removed.
It outputs two files called "ligand_elec.csv" and "ligand_vdw.csv". Usage: lig-interaction.sh <topology_file> <trajectory_file> <reference frame> <ligand_id>. Example: ./lig-interaction.sh comp.prmtop 11md-2.dcd 11md.rst7 199'''


#!/bin/bash

if [ $# -ne 4 ]; then
    echo "Usage: $0 <topology_file> <trajectory_file> <reference> <ligand_id>"
    exit 1
fi

# Read command-line arguments
topology_file=$1
trajectory_file=$2
reference=$3
ligand_id=$4

# Run the cpptraj command and store the output in a variable
cpptraj_output=$(cpptraj -p $topology_file -c $reference --resmask ":$ligand_id<:6")

# Initialize an empty array to store the numbers
numbers=()

# Iterate through each line of the cpptraj output
while IFS= read -r line; do
    # Skip the first line (header)
    if [ "$line" != "#Res  Name First  Last Natom #Orig  #Mol C I" ]; then
        # Use awk to extract the number from the first column and add it to the array
        number=$(echo "$line" | awk '{print $1}')
        numbers+=("$number")
    fi
done <<< "$cpptraj_output"

# Create a temporary input file for the lie commands
lie_commands_file="lie_commands.txt"
rm -f "$lie_commands_file"  # Remove if it already exists

# Generate the lie commands and write them to the input file
for number in "${numbers[@]}"; do
    lie_outfile="lie-LIG_to_res${number}.txt"  
    lie_command="lie energy$number :$ligand_id :$number out $lie_outfile"
    echo "$lie_command" >> "$lie_commands_file"
done

# Run the generated lie commands using cpptraj, importing the parm and traj
cpptraj -p $topology_file -y $trajectory_file -i "$lie_commands_file"

# Create CSV files for electrostatic and van der Waals averages
elec_averages_file="ligand_elec.csv"
vdw_averages_file="ligand_vdw.csv"

# Initialize arrays to store electrostatic and van der Waals average values
declare -A avg_elec_values
declare -A avg_vdw_values

# Extract electrostatic and van der Waals values and calculate averages
for number in "${numbers[@]}"; do
   lie_outfile="lie-LIG_to_res${number}.txt"
   
   # Calculate the averages for elec and vdw columns
   avg_elec=$(awk 'NR > 1 { elec_sum += $2 } END { if (NR > 1) print elec_sum/(NR-1); else print 0 }' "$lie_outfile")
   avg_vdw=$(awk 'NR > 1 { vdw_sum += $3 } END { if (NR > 1) print vdw_sum/(NR-1); else print 0 }' "$lie_outfile")
    
   # Store average values in arrays
   avg_elec_values["$number"]=$avg_elec
   avg_vdw_values["$number"]=$avg_vdw
    
   # Remove the intermediate file
   rm "$lie_outfile"
done

# Write the CSV files for electrostatic and van der Waals averages with formatted headers
elec_averages_file="ligand_elec.csv"
vdw_averages_file="ligand_vdw.csv"

# Format the headers

vdw_header=""
for number in "${numbers[@]}"; do
   vdw_header+="v${ligand_id}-v${number},"
done
echo "$vdw_header" > ligand_vdw.csv

elec_header=""
for number in "${numbers[@]}"; do
   elec_header+="e${ligand_id}-e${number},"
done
echo "$elec_header" > ligand_elec.csv

# Format the data values
elec_row=""
for number in "${numbers[@]}"; do
   elec_row+="${avg_elec_values["$number"]},"
done
echo "$elec_row" >> ligand_elec.csv

# Print values for VDW interactions
vdw_row=""
for number in "${numbers[@]}"; do
   vdw_row+="${avg_vdw_values["$number"]},"
done
echo "$vdw_row" >> ligand_vdw.csv

# Clean up: remove the temporary lie commands file
rm "$lie_commands_file"

echo "Electrostatic averages have been written to $elec_averages_file"
echo "Van der Waals averages have been written to $vdw_averages_file"

