import luigi
import pandas as pd

class LoadDataTask(luigi.Task):
    def output(self):
        # Mention path for pipeline task output if needed
        return luigi.LocalTarget("Data/raw_data.csv")

    def run(self):
        # Implement data loading logic here
        df = pd.read_csv("Data/input_sample.csv")

        df.to_csv(self.output().path, index=True)
