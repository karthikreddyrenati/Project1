import luigi

from Pipeline.DataIntegration import IntegrateDataTask
from IgniteConnect import IgniteConnector
import pandas as pd
import json
def getIgniteConnection():
    instance = IgniteConnector()
    instance.connect()
    client = instance.get_client()
    return client

class StoreInIgniteTask(luigi.Task):

    def requires(self):
        return IntegrateDataTask()
    
    def getConnection(self):
        self.instance = IgniteConnector()
        self.instance.connect()
        return self.instance.get_client()

    def run(self):
        # Implement Apache Ignite storage logic here
        client = self.getConnection()

        cache_config = {
        'cache_name': 'test_cache',
        'template_name': 'replicated'
        }

        if self.file_type == 'pdb':
            with open(self.input().path, 'r') as file:
                json_data = file.read()
        elif self.file_type == 'xyz':
            with open(self.input().path, 'r') as file:
                json_data = json.load(file)

        df = pd.read_csv("Data/integrated_data.csv")
        dfAsDict = df.to_dict(orient='list')

        cache = client.get_or_create_cache("test_cache")

        for key in dfAsDict.keys():
            cache.put(key, dfAsDict[key])

    