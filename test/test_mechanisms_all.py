import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
sys.path.append('../')

file = 'test_database.csv'
from pyrrhenious import database


def test_create_database():

    ecdatabase = database.Database(file)
    ecdatabase.create_anisotropic_models()
    assert 'isotropic_model:author19[100]+author19[010]+author19[001]' in ecdatabase.database['entry_id'].unique(), ecdatabase.database['entry_id'].unique()
