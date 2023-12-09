import luigi
import pandas as pd
import numpy as np
# import MDAnalysis as mda
import pandas as pd
from collections import defaultdict
import json

class LoadDataTask(luigi.Task):
    file_type = luigi.Parameter()
    def output(self):
        # # Mention path for pipeline task output if needed
        # if self.file_type == 'crd':
        #     return luigi.LocalTarget("Data/raw_data_crd.csv")
        # elif self.file_type == "pdb":
        #     return luigi.LocalTarget("temp_processed_data.json")
        # elif self.file_type == "xyz":
        #     return luigi.LocalTarget("Data/raw_data_xyz.json")
        pass


    def run(self):
        pass



