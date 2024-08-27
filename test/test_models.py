import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
sys.path.append('../')
import pytest
from pyrrhenius import model

file = 'test_database.csv'


@pytest.fixture
def dataframe():
    return pd.read_csv(file, encoding_errors='replace').dropna(how='all')

@pytest.fixture
def row_1():
    return pd.read_csv(file, encoding_errors='replace').dropna(how='all').iloc[0]

@pytest.fixture
def row_2():
    return pd.read_csv(file, encoding_errors='replace').dropna(how='all').iloc[1]

def test_make_basic_model(row_1):
    new_model = model.create_model_from_row(row_1)
    assert isinstance(new_model,model.ModelInterface)
    np.testing.assert_allclose(new_model.get_conductivity(T=1500),np.asarray([0.112184]),rtol=1e-05)

def test_make_database_row(row_1):
    new_model = model.create_model_from_row(row_1).get_row()
    assert isinstance(new_model, pd.Series)

def test_add_models(row_1):
    new_model1 = model.create_model_from_row(row_1)
    new_model2 = model.create_model_from_row(row_1)
    composite_model = new_model1 + new_model2
    assert isinstance(composite_model,model.CompositeModel)
    assert str(composite_model)=='author1+author1:{author1:{10^1.0(0.1) exp( -0.67(0.07)/kT)+10^1.0(0.1) exp( '+\
                                 '-0.67(0.07)/kT)}+author1:{10^1.0(0.1) exp( -0.67(0.07)/kT)+10^1.0(0.1) exp( '+\
                                  '-0.67(0.07)/kT)}}'
    np.testing.assert_allclose(composite_model.get_conductivity(T=1500), 2*np.asarray([0.112184]), rtol=1e-05)
    assert isinstance(composite_model.get_row(), pd.Series)
    assert composite_model.get_id=='author1+author1'


def test_add_dry_model(row_1):
    new_model1 = model.create_model_from_row(row_1)
    new_model2 = model.create_model_from_row(row_1)
    composite_model = new_model1 + model.DryModel(new_model2)
    assert isinstance(composite_model,model.CompositeModel)
    np.testing.assert_allclose(composite_model.get_conductivity(T=1500), 2*np.asarray([0.112184]), rtol=1e-05)
    assert isinstance(composite_model.get_row(), pd.Series)
    assert 'dry' in composite_model.get_id


def test_anisotropic_model(dataframe):
    new_model1 = model.create_model_from_row(dataframe.loc[dataframe['Entry ID'] == 'author19[100]'])
    new_model2 = model.create_model_from_row(dataframe.loc[dataframe['Entry ID'] == 'author19[010]'])
    new_model3 = model.create_model_from_row(dataframe.loc[dataframe['Entry ID'] == 'author19[001]'])
    model_list = [new_model1,new_model2,new_model3]
    anisomodel = model.AnisotropicModel([new_model1,new_model2,new_model3])
    assert anisomodel.get_id!=None
    assert anisomodel.get_id=='isotropic_model:author19[100]+author19[010]+author19[001]'
