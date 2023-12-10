import luigi
import pandas as pd
import numpy as np
# import MDAnalysis as mda
import pandas as pd
from collections import defaultdict
import json

class LoadDataTask(luigi.Task):
    # file_type = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget("Data/dataload_stage_output.csv")


    def run(self):
        

        df = pd.read_csv('Data/1h9t_traj.xyz.gz', compression='gzip', delim_whitespace=True, names=['atoms', 'x', 'y', 'z'], skiprows=2)

        sample_df = df.head()

        sample_df.to_csv(self.output().path, index=False)