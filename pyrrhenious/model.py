import numpy as np
import pandas as pd
from dataclasses import dataclass
from inspect import signature


@dataclass
class PublicationMetadata:
    title : str
    author: str
    year: str
    doi: str
    phase_type: str
    description: str
    sample_type: str
    equation_form: str
    publication_id: str
    entry_id: str
    complete_or_partial_fit: str
    composite_or_single: str
    pressure_average_gpa: float
    pressure_min_gpa: float
    pressure_max_gpa: float
    temp_mink: float
    temp_maxk: float
    water_min: float
    water_max: float
    water_average: float
    water_calibration: str
    water_units: str
    iron_min: float
    iron_max: float
    iron_average: float
    iron_units: str

    @classmethod
    def from_kwargs(cls, **kwargs):
        # fetch the constructor's signature
        cls_fields = {field for field in signature(cls).parameters}

        # split the kwargs into native ones and new ones
        native_args, new_args = {}, {}
        for name, val in kwargs.items():
            if name in cls_fields:
                native_args[name] = val
            else:
                new_args[name] = val

        # use the native ones to create the class ...
        ret = cls(**native_args)

        # ... and add the new ones by hand
        for new_name, new_val in new_args.items():
            setattr(ret, new_name, new_val)
        return ret



class ModelInterface:

    def get_conductivity(self, **kwargs):
        pass

    def title(self):
        pass
    @staticmethod
    def determine_starting_c_type(P=None, T=None,logfo2=None, **kwargs):
        ndarrays = list(filter(lambda x: isinstance(x, np.ndarray), [P, T, logfo2]))
        if len(ndarrays)<1:
            reference_shape=(1)
        elif len(ndarrays)==1:
            reference_shape= ndarrays[0].shape
        else:
            reference_shape = ndarrays[0].shape
            assert all(array.shape == reference_shape for array in ndarrays),  f"unsure how to construct c since dimensions of T, P, and logfo2 don't match:\nP:\t{P}\nT:\t{T}\nlogfo2:\t{logfo2}"
        return np.zeros(reference_shape)
class Model(ModelInterface):
    def __init__(self, mechanisms: list, id_rows: pd.Series):
        super().__init__()
        self.metadata = PublicationMetadata.from_kwargs(**id_rows)
        self.id = id_rows['entry_id']
        self.mechanisms = mechanisms
        for m in mechanisms:
            if m.uses_water:
                m.set_water_units(self.metadata.water_units)

    @property
    def uses_water(self):
        return any(m.uses_water for m in self.mechanisms)


    def get_conductivity(self, **kwargs):
        c = self.determine_starting_c_type(**kwargs)
        for m in self.mechanisms:
            c += m.get_conductivity(**kwargs)
        return c

    @property
    def title(self):
        return f'{self.metadata.phase_type}\n{self.metadata.publication_id}\n{self.metadata.entry_id}'

class CompositeModel(ModelInterface):

    def __init__(self, models):
        super().__init__()
        self.models = models
        self.id = '_'.join(m.id for m in models)

    def get_conductivity(self, crystal_direction=None, **kwargs):
        c = self.determine_starting_c_type(**kwargs)
        for m in self.models:
            c += m.get_conductivity(**kwargs)
        return c

    @property
    def uses_water(self):
        return any(m.uses_water for m in self.models)

    @property
    def title(self):
        return '+\n'.join([f'{m.metadata.phase_type} '+\
                          f'{m.metadata.publication_id} '+\
                          f' {m.metadata.entry_id}' for m in self.models])

class AnisotropicModel(CompositeModel):

    def __init__(self, models):
        super().__init__(models)
        self.models = {m.id[-5:]: m for m in models}
        self.id = models[0].id[:-5]

    def get_conductivity(self, crystal_direction=None, averaging='geometric', **kwargs):
        if crystal_direction in self.models.keys():
            return self._get_xstal_direction(crystal_direction,**kwargs)
        elif isinstance(crystal_direction,float):
            return self._get_anisotropic_factor(**kwargs)
        else:
            return self._get_geometric_mean(**kwargs)

    def _get_xstal_direction(self, xstal_direction,**kwargs):
        return self.models[xstal_direction].get_conductivity(**kwargs)

    def _get_geometric_mean(self,**kwargs):
        c = self.determine_starting_c_type(**kwargs)
        for m in self.models:
            c *= m.get_conductivity(**kwargs)**(1/3)
        return c

    def _get_anisotropic_factor(self,factor,**kwargs):
        conductivities = [ m.get_conductivity(**kwargs) for m in self.models.keys()]
        if isinstance(conductivities[0],np.ndarray):
            min_array = np.min(np.stack(conductivities,axis=-1),axis=-1)
            max_array = np.max(np.stack(conductivities,axis=-1), axis=-1)
        else:
            min_array = min(conductivities)
            max_array = max(conductivities)
        return min_array * (1 - factor) + factor * max_array

    @property
    def title(self):
        return 'Anisotropic+\n' +''.join([f'{m.metadata.phase_type}\n'+\
                f'{m.metadata.publication_id}' for m in self.models])


