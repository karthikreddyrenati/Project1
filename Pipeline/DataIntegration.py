import luigi
from Pipeline.DataPreprocessing import PreprocessDataTask
import pandas as pd

class IntegrateDataTask(luigi.Task):
    def requires(self):
        return PreprocessDataTask()

    def output(self):
        return luigi.LocalTarget("Data/integrated_data.csv")

    def run(self):
        # Implement data integration logic here
        data = pd.read_csv("Data/preprocessed_data.csv")

        data.to_csv(self.output().path, index=False)