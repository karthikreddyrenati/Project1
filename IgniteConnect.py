from pyignite import Client

class IgniteConnector:

    def __init__(self):
        self._instance = None
        

    def connect(self, host='127.0.0.1', port=10800):
        if not self._instance:
            self._instance = Client()
            self._instance.connect(host, port)
            self._instance.get_cluster().set_state(True)
            print("Ignite Client Connected.")

        else:
            print("Ignite clinet is already established.")

    def get_client(self):
        if not self._instance:
            raise Exception("Ignite client is not connected. Call connect() first.")
        return self._instance

    def close_connection(self):
        if self._instance:
            self._instance.close()