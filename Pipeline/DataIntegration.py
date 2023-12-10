import luigi
from Pipeline.DataPreprocessing import PreprocessDataTask
import pandas as pd

class IntegrateDataTask(luigi.Task):
    def requires(self):
        return PreprocessDataTask()

    def output(self):
        return luigi.LocalTarget("Data/dataintegrate_stage_output.csv")

    def run(self):
        # Implement data integration logic here
        data = pd.read_csv("Data/dataprocess_stage_output.csv")

        data.to_csv(self.output().path, index=False)