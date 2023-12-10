import luigi

from Pipeline.DataIntegration import IntegrateDataTask
from IgniteConnect import IgniteConnector
import pandas as pd
import json
from Pipeline.DataLoading import LoadDataTask
import logging
import time


logging.basicConfig(
    filename='apache_ignite.log',  # Specify the log file
     filemode='w',
    level=logging.DEBUG,      # Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR, CRITICAL)
    format='%(asctime)s - %(levelname)s - %(message)s'  # Specify the log message format
)

# def getIgniteConnection():
#     instance = IgniteConnector()
#     instance.connect()
#     client = instance.get_client()
#     return client

class StoreInIgniteTask(luigi.Task):

    # def requires(self):
    #     # return IntegrateDataTask()
    #     return LoadDataTask()
    
    def getConnection(self):
        self.instance = IgniteConnector()
        self.instance.connect()
        return self.instance.get_client()

    def run(self):

        df = pd.read_csv('Data/1h9t_traj.xyz.gz', compression='gzip', delim_whitespace=True, names=['atoms', 'x', 'y', 'z'], skiprows=2)
        # df = pd.read_csv('Data/dataload_stage_output.csv')
        unique_atoms_frame = df.groupby('atoms').agg(lambda x: x.tolist()).reset_index()

        # Implement Apache Ignite storage logic here
        client = self.getConnection()


        # cache_config = {= 
        # 'cache_name': 'test_cache',
        # 'template_name': 'replicated'
        # }

        # if self.file_type == 'pdb':
        #     with open(self.input().path, 'r') as file:
        #         json_data = file.read()
        # elif self.file_type == 'xyz':
        #     with open(self.input().path, 'r') as file:
        #         json_data = json.load(file)

        # df = pd.read_csv("Data/dataintegrate_stage_output.csv")
        dfAsDict = unique_atoms_frame.to_dict(orient='records')
        dfAsDict = {record.pop('atoms'): record for record in dfAsDict}
        # logging.info(dfAsDict)
        cache = client.get_or_create_cache("xyz_trajectory")
        logging.info("+++++++++++++++++++++++++++  Performance Evaluation +++++++++++++++++++++++++++")
        startWrite = time.time()
        for key in dfAsDict.keys():
            logging.info(dfAsDict[key])
            cache.put(key, str(dfAsDict[key]))

        endWrite = time.time()

        logging.info("Write Speed = {} secs".format(round(endWrite - startWrite,2)))


        startRead = time.time()
        xyzData = cache.get_all(dfAsDict.keys())
        endRead = time.time()

        logging.info("Data XYZ trajectory")
        # logging.info(xyzData)

        logging.info("Read Speed = {} secs".format(round(endRead - startRead,2)))

