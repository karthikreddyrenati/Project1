from IgniteConnect import IgniteConnector

instance = IgniteConnector()
instance.connect()
client = instance.get_client()


cache = client.get_or_create_cache("test_cache")


prices = cache.get("Price")

print(prices)