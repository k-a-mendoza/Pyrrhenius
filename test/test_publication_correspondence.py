import glob

import numpy as np
import pandas as pd
import sys
import os
file = '../database/publication_database.csv'
import pyrrhenius


def test_create_database():
    print(os.getcwd())
    ecdatabase = pyrrhenius.database.Database(pyrrhenius.DEFAULT_DATABASE_PATH)
    ecdatabase.create_isotropic_models()
    assert True

