
import luigi
from Pipeline.StoreData import StoreInIgniteTask


def main():

    luigi.build([StoreInIgniteTask()], local_scheduler=True)


if __name__ == "__main__":
    main()
