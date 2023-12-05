import luigi
from Pipeline.DataLoading import LoadDataTask
import pandas as pd
import json

class PreprocessDataTask(luigi.Task):
    def requires(self):
        # Requires previous task to complete
        return LoadDataTask()

    def output(self):
        return luigi.LocalTarget("Data/preprocessed_data.csv")

    def run(self):
        #Logic for pre-processing

        dataFrame = pd.read_csv("Data/raw_data.csv")
        dataFrame.dropna()
        
        dataFrame.to_csv(self.output().path, index=False)


