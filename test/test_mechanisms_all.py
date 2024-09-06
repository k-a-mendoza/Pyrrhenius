import numpy as np
import pandas as pd
import sys
import pytest
from pyrrhenius.mechanisms import SingleValue, StochasticConstant

sys.path.append('../../')
sys.path.append('../')


@pytest.fixture
def mechanism():
    return SingleValue(StochasticConstant(0,0.1,type='test'))

def test_set_water_units(mechanism):
    with pytest.raises(NotImplementedError):
        mechanism.set_water_units('invalid_units') 
        mechanism.convert_water(100)

def test_convert_water(mechanism):
    mechanism.water_units = 'wtpct'
    assert mechanism.convert_water(100) == 1e-2, 'could not convert ppm into wtpct'
    mechanism.water_units = 'ppm'
    assert mechanism.convert_water(100) == 100, 'direct ppm conversion did not work.:' + str(mechanism.convert_water(100))

    

def test_convert_co2(mechanism):
    result = mechanism.convert_co2(10)
    assert result == 10

def test_convert_pressure(mechanism):
    result = mechanism.convert_pressure(5)
    assert result == 0.051818478, 'failed to convert pressure at 5 GPa to :'+ str(result)

def test_assert_pressure(mechanism):
    with pytest.raises(AssertionError):
        mechanism._parameter_assertion(None,'P')
