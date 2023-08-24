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

def test_all_mechanisms():
    ecdatabase = database.Database(file)
    ecdatabase.load_models()
    T = 1400
    P = 1
    X_fe = 0.3
    Cw = 100 # ppm
    CO2 = 1
    logfO2 = -5
    models = ecdatabase.get_model_names()
    for model_name in models:
        ecmodel = ecdatabase.get_model(model_name)
        row     = ecdatabase.get_model_properties(model_name)
        try:
            result = ecmodel.get_conductivity(T=T,P=P,X_fe=X_fe,Cw=Cw,Co2=CO2,logfo2=logfO2)
        except Exception as e:
            assert False, f'Problem calculating conductivity:  {row["entry_id"]}\n{e}'
        assert np.isfinite(result), f'produced infinite result: {row["entry_id"]}'
        assert result > 1e-10, f'produced too low of a conductivity: {row["entry_id"]}'
        assert result < 1e5, f'produced too high of a conductivity: {row["entry_id"]}'

def test_all_mechanisms_for_monotonically_increasing_conductivity():
    ecdatabase = database.Database(file)
    ecdatabase.load_models()
    T = np.arange(100,2200,100)
    P = 1
    X_fe = 0.3
    Cw = 100 # ppm
    CO2 = 1
    logfO2 = -5
    models = ecdatabase.get_model_names()
    for model_name in models:
        ecmodel = ecdatabase.get_model(model_name)
        result = ecmodel.get_conductivity(T=T,P=P,X_fe=X_fe,Cw=Cw,Co2=CO2,logfo2=logfO2)
        assert np.all(np.diff(result) >= 0), f"Found a non monotonically increasing mechanism {model_name}"