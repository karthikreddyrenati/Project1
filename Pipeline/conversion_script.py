
# Combined Script
import numpy as np
import MDAnalysis as mda
import pandas as pd
import json


# Part 1: CRD to CSV Script
def read_crd_file(file_path):
    coordinates = []
    with open(file_path, 'r') as file:
        for line in file:
            values = line.split()
            # Assuming each line has a multiple of 3 values (x, y, z coordinates)
            for i in range(0, len(values), 3):
                if i+2 < len(values):  # Ensure there's a complete set of coordinates
                    x, y, z = float(values[i]), float(values[i+1]), float(values[i+2])
                    coordinates.append([x, y, z])
    return np.array(coordinates)

def save_as_csv(array, output_file_path):
    np.savetxt(output_file_path, array, delimiter=',', fmt='%f')  # Using '%f' for floating-point format

# Usage
file_path = '1h9t_md_restrain.crd'
coordinates = read_crd_file(file_path)

# Save to CSV
csv_file_path = '../Data/1h9t_md_restrain_crd.csv'
save_as_csv(coordinates, csv_file_path)
print(f"Coordinates saved to {csv_file_path}")


# Part 2: PDB to JSON Script

# Load the PDB file
u = mda.Universe('1h9t_md.pdb') # right here, please enter the actual file path and file

# Create a list of dictionaries, each containing data for one atom
data = []
for atom in u.atoms:
    atom_data = {
        'name': atom.name,
        'x': atom.position[0],
        'y': atom.position[1],
        'z': atom.position[2]
    }
    data.append(atom_data)

# Convert the list of dictionaries into a DataFrame
df = pd.DataFrame(data)

# Using defaultdict to store values
records = defaultdict(list)

for index, row in df.iterrows():
    records[row['name']].append({"X": row['x'], "Y": row['y'], "Z": row['z']})

# Formatting the output
for name, values in records.items():
    print(f'"{name}": {values}')

with open('1h9t_md_pdb.json', 'w') as file:
    json.dump(records, file, indent=4)




# Part 3: XYZ to JSON Script
def read_gz_xyz_file(file_path):
    # Read the file into a DataFrame, assuming whitespace delimiters
    # The compression='gzip' option tells pandas to decompress the file
    df = pd.read_csv(file_path, compression='gzip', delim_whitespace=True, names=['atoms', 'x', 'y', 'z'], skiprows=2)
    return df

def aggregate_coordinates(group):
    return [{'x': row['x'], 'y': row['y'], 'z': row['z']} for _, row in group.iterrows()]

def main():
    # Usage example
    file_path = '1h9t_traj.xyz.gz'
    df = read_gz_xyz_file(file_path)

    # Open a file to write
    json_file_path = '1h9t_traj_xyz.json'
    with open(json_file_path, 'w') as file:
        file.write('{')  # Start of the JSON object

        # Iterate through each group, process, and write to file
        for i, (atom, group) in enumerate(df.groupby('atoms')):
            if i > 0:  # Add a comma before the next group if it's not the first one
                file.write(',')
            json_string = json.dumps({atom: aggregate_coordinates(group)}, indent=4)
            file.write(json_string[1:-1])  # Write the group JSON without outer braces

        file.write('}')  # End of the JSON object

    print(f"Data saved to {json_file_path}")

if __name__ == '__main__':
    main()

