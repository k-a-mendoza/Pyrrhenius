import pytest
from pytest import fixture 
import numpy as np
from sklearn import base
import pyrrhenius.database as db
import pyrrhenius.mixing as mix
import os 

DEFAULT_DATABASE_PATH = os.path.join('database', 'publication_database.csv')

@fixture 
def database():
    csv = DEFAULT_DATABASE_PATH
    return db.Database(csv_file=csv)

@fixture
def conductivity_models(database):
    model1 = database.get_model('Li_18_1%plg_brine')
    model2 = database.get_model('ks_83_basalt')
    model3 = database.get_model('xu_1999_cpx')
    return model1, model2,model3

@fixture
def temperature_array():
    return np.linspace(300, 500, num=100)


@pytest.mark.parametrize("model_class", [
    mix.EffectiveMediumTheoryMixture,
    mix.HashinStrickmanBound,
    mix.ParallelModel,
    mix.SeriesModel,
    mix.ArchiesLawGlover,
    mix.TubesModel,
    mix.GeomAverage,
    mix.ArithmeticAverage
])
def test_two_phases(model_class, conductivity_models, temperature_array):
    model1, model2 = conductivity_models[1:] 
    model = model_class([model1, model2])

    base_conductivities = [x.get_conductivity(T=temperature_array) for x in [model1,model2]]
    base_conductivities = np.stack(base_conductivities,axis=-1)
    sigma1 = np.min(base_conductivities,axis=-1)
    sigma2 = np.max(base_conductivities,axis=-1)
    mixed_conductivity = model.get_conductivity(0.3, T=temperature_array)

    assert all(sigma1<=mixed_conductivity), "float fraction: Mixed conductivity is not greater than the first phase conductivity"
    assert all(mixed_conductivity<=sigma2), "float fraction: Mixed conductivity is not less than the second phase conductivity"

    mixed_conductivity = model.get_conductivity([0.3, 0.7], T=temperature_array)
    assert all(sigma1<=mixed_conductivity), "array fraction: Mixed conductivity is not greater than the first phase conductivity"
    assert all(mixed_conductivity<=sigma2), "array fraction: Mixed conductivity is not less than the second phase conductivity"
    # Test with invalid array of floats
    with pytest.raises(AssertionError):
        model.get_conductivity([0.6, 0.6],T=temperature_array)


# The error in this test suite is in the line `model1, model2 = conductivity_models[0]`. 
# `conductivity_models[0]` returns a single model, not a tuple. 
# It should be `model1 = conductivity_models[0]`.

@pytest.mark.parametrize("model_class,args", [
    pytest.param(mix.CubeModel, 'waff74-1'),
    pytest.param(mix.ArchiesLaw, None)
])
def test_one_phase_models(model_class, args, conductivity_models, temperature_array):
    model1 = conductivity_models[1] # Fixing the error here

    if args is not None:
        model = model_class(model1, args)
    else:
        model = model_class(model1)
    base_conductivity = model1.get_conductivity(T=temperature_array)
    mixed_conductivity = model.get_conductivity(0.3, T=temperature_array)

    assert all(mixed_conductivity <= base_conductivity), "float fraction: Mixed conductivity is not greater than the first phase conductivity"
    assert all(0 < mixed_conductivity), "float fraction: Mixed conductivity is not greater than zero in all cases"

    with pytest.raises(AssertionError):
        mixed_conductivity = model.get_conductivity(1.1, T=temperature_array)

@pytest.mark.parametrize("model_class", [
    mix.EffectiveMediumTheoryMixture,
    mix.HashinStrickmanBound,
    mix.ParallelModel,
    mix.SeriesModel,
    mix.GeomAverage,
    mix.ArithmeticAverage
])
def test_three_phases(model_class, conductivity_models, temperature_array):
    model = model_class(conductivity_models)

    # Test with a single float
    with pytest.raises(AssertionError):
        model.get_conductivity(0.3,T= temperature_array)
    # Test with an array of floats

    base_conductivities = [x.get_conductivity(T=temperature_array)*np.ones(len(temperature_array)) for x in conductivity_models]
    base_conductivities = np.stack(base_conductivities,axis=-1)
    sigma1 = np.min(base_conductivities,axis=-1)
    sigma2 = np.max(base_conductivities,axis=-1)
    mixed_conductivity = model.get_conductivity([0.3, 0.5, 0.2], T=temperature_array)

    assert all(sigma1<=mixed_conductivity), "float fraction: Mixed conductivity is not greater than the first phase conductivity"
    assert all(mixed_conductivity<=sigma2), "float fraction: Mixed conductivity is not less than the second phase conductivity"

   
    with pytest.raises(AssertionError):
        model.get_conductivity([0.6, 0.6, 0.6],  T=temperature_array)
    with pytest.raises(AssertionError):
        model.get_conductivity([0.6, 0.4], T=temperature_array)


