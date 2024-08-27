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

def test_single_metadata(row_2):
    metadata = model.PublicationMetadata.from_kwargs(**row_2)
    assert metadata.title=='Test2'
    assert metadata.author=='author2'
    assert metadata.year==2.0
    assert metadata.phase_type=='phase1'
    assert metadata.water_units=='ppm'
    assert metadata.water_max==4.0
    assert metadata.water_min==1.0
    assert metadata.iron_max==1
    assert metadata.iron_min==0.5
    assert metadata.iron_average == 0.75
    assert metadata.pressure_max==5.0
    assert metadata.pressure_min==2.0
    assert metadata.pressure_average==3.5
    assert metadata.temp_min == 1300
    assert metadata.temp_max == 1773
    assert metadata.entry_id=='author2'


def test_add_metadata(row_1,row_2):
    metadata1 = model.PublicationMetadata.from_kwargs(**row_1)
    metadata2 = model.PublicationMetadata.from_kwargs(**row_2)
    metadata_composite = metadata1+metadata2
    assert metadata_composite.entry_id=='author1+author2', 'id composite is not working'
    assert metadata_composite.publication_id == 'author1+author2', 'publication composite is not working'
    assert metadata_composite.author=='author1+author2', 'author composite is not working'
    assert metadata_composite.phase_type=='phase1', 'phase not defined'
    assert metadata_composite.water_max==4.0, 'water max not working'
    assert metadata_composite.water_min==0, 'water min not working'
    assert metadata_composite.iron_max==1, 'iron max not working'
    assert metadata_composite.iron_min==0, 'iron min not working'
    assert metadata_composite.iron_average == 0.5, 'iron avg not working'
    assert metadata_composite.pressure_max==5.0
    assert metadata_composite.pressure_min==0
    assert metadata_composite.pressure_average==2.5
    assert metadata_composite.temp_min == 1273
    assert metadata_composite.temp_max == 1773


