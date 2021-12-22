from pathlib import Path
import os
import datetime


def athome(path):
    return str(Path(os.environ["HOME"], path))


def today():
    return datetime.datetime.today().strftime("%d-%b-%Y")
