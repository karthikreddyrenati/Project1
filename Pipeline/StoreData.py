import luigi

from Pipeline.DataIntegration import IntegrateDataTask
from IgniteConnect import IgniteConnector
import pandas as pd
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
        
        df = pd.read_csv("Data/integrated_data.csv")
        dfAsDict = df.to_dict(orient='list')

        cache = client.get_or_create_cache("test_cache")

        for key in dfAsDict.keys():
            cache.put(key, dfAsDict[key])

    