import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
sys.path.append('../')

file = 'test_database.csv'
from pyrrhenious import database


def test_create_database():

    ecdatabase = database.Database(file)
    ecdatabase.load_models()
    assert 'author19_anisomodel' in ecdatabase.database['entry_id'].unique(),ecdatabase.database['entry_id'].unique()

