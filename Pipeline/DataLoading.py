import luigi
import pandas as pd
import numpy as np
import MDAnalysis as mda
import pandas as pd
from collections import defaultdict
import json

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

def load_pdb_file(file_path):
    u = mda.Universe(file_path)
    return u

def process_pdb_data(universe):
    data = []
    for atom in universe.atoms:
        atom_data = {
            'name': atom.name,
            'x': atom.position[0],
            'y': atom.position[1],
            'z': atom.position[2]
        }
        data.append(atom_data)
    df = pd.DataFrame(data)

    # Using defaultdict to store values
    records = defaultdict(list)
    for index, row in df.iterrows():
        records[row['name']].append({"X": row['x'], "Y": row['y'], "Z": row['z']})

    # Formatting the output as JSON
    return json.dumps(records, indent=4)

def save_as_csv(array, output_file_path):
    np.savetxt(output_file_path, array, delimiter=',', fmt='%f')

def read_gz_xyz_file(file_path):
    # Read the file into a DataFrame, assuming whitespace delimiters
    # The compression='gzip' option tells pandas to decompress the file
    df = pd.read_csv(file_path, compression='gzip', delim_whitespace=True, names=['atoms', 'x', 'y', 'z'], skiprows=2)
    return df

def aggregate_coordinates(group):
    return [{'x': row['x'], 'y': row['y'], 'z': row['z']} for _, row in group.iterrows()]


class LoadDataTask(luigi.Task):
    file_type = luigi.Parameter()
    def output(self):
        # Mention path for pipeline task output if needed
        if self.file_type == 'crd':
            return luigi.LocalTarget("Data/raw_data_crd.csv")
        elif self.file_type == "pdb":
            return luigi.LocalTarget("temp_processed_data.json")
        elif self.file_type == "xyz":
            return luigi.LocalTarget("Data/raw_data_xyz.json")


    def run(self):
        # Implement data loading logic here
        if self.file_type == 'crd': # crd file
            file_path = 'Data/1h9t_md_restrain.crd'
            coordinates = read_crd_file(file_path)
            save_as_csv(coordinates, self.output().path)
        elif self.file_type == 'pdb': # pdb file
            pdb_file_path = '1h9t_md.pdb'
            universe = load_pdb_file(pdb_file_path)
            processed_data = process_pdb_data(universe)
            temp_json_path = 'temp_processed_data.json'
            with open(temp_json_path, 'w') as file:
                file.write(processed_data)
            with self.output().open('w') as out_file:
                out_file.write(temp_json_path)
        elif self.file_type == 'xyz': # xyz file
            xyz_file_path = 'Data/1h9t_traj.xyz.gz'
            df = read_gz_xyz_file(xyz_file_path)
            with open(self.output().path, 'w') as file:
                file.write('{')  # Start of the JSON object

                for i, (atom, group) in enumerate(df.groupby('atoms')):
                    if i > 0:
                        file.write(',')
                    json_string = json.dumps({atom: aggregate_coordinates(group)}, indent=4)
                    file.write(json_string[1:-1])

                file.write('}')



