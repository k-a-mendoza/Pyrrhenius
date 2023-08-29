import numpy as np
import pandas as pd
from dataclasses import dataclass
from inspect import signature
from . import  mechanisms as mech


def _id_row_rename(col):
    col = col.replace(' ', '_').replace('(', '').replace(')', '')
    return col.lower()

def create_model_from_row(row: pd.Series or pd.DataFrame):
    if not isinstance(row,pd.Series):
        assert len(row)==1,'Make multi model via an add operation. Function only makes single row models'
        row = row.iloc[0]
    row.rename(_id_row_rename, inplace=True)
    mechs = mech.create_mechanism_from_row(row)
    metadata_row_keys = list(set(row.keys()) - set(mech.get_const_rows(row)))
    return Model(mechs, row[metadata_row_keys])


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
    pressure_average: float
    pressure_min: float
    pressure_max: float
    temp_min: float
    temp_max: float
    water_min: float
    water_max: float
    water_average: float
    water_calibration: str
    water_units: str
    iron_min: float
    iron_max: float
    iron_average: float
    iron_units: str
    crystal_direction : str

    def recalc_average_fields(self):
        avg_fields = ['iron_average','water_average','pressure_average']
        for field in avg_fields:
            max_field_name = field.split('_')[0] + '_max'
            min_field_name = field.split('_')[0] + '_min'
            max_field = getattr(self,max_field_name)
            min_field = getattr(self,min_field_name)
            if max_field is None or min_field is None:
                continue

            max_field = float(max_field)
            min_field = float(min_field)
            if not np.isnan(max_field) and not np.isnan(min_field):
                setattr(self, field,  float((float(getattr(self,max_field_name))+
                                       float(getattr(self, min_field_name)))/2))
    @classmethod
    def get_float_field_prefixes(cls):
        return ['iron','water','pressure','temp']
    @classmethod
    def get_additive_field_strings(cls):
        return ['author','entry_id','publication_id','crystal_direction']
    @classmethod
    def get_static_field_strings(cls):
        return ['phase_type']

    @classmethod
    def from_kwargs(cls, **kwargs):
        # fetch the constructor's signature
        cls_fields = {field for field in signature(cls).parameters}

        # split the kwargs into native ones and new ones
        native_args, new_args = {}, {}
        for name, val in kwargs.items():
            name_translate = _id_row_rename(name)
            if name_translate in cls_fields:
                native_args[name_translate] = val
            else:
                new_args[name_translate] = val

        # use the native ones to create the class ...
        ret = cls(**native_args)
        ret.recalc_average_fields()
        return ret

    def __add__(self,other):
        assert self.phase_type==other.phase_type, "Trying to add different phases! You probably shouldn't combine them."
        cls_fields = {field for field in signature(PublicationMetadata).parameters}
        new_params = {}
        for name in cls_fields:
            val = None
            if any(ext in name and 'calibration' not in name and 'units' not in name for ext in self.get_float_field_prefixes()):
                value1 = getattr(self,name)
                value2 = getattr(other,name)

                if value1 is None and value2 is None:
                    pass
                elif isinstance(value1,float) and np.isnan(value1) and value2 is None:
                    pass
                elif isinstance(value2,float) and np.isnan(value2) and value1 is None:
                    pass
                elif all(isinstance(x, float) for x in [value1,value2]) and np.isnan(value1) and np.isnan(value2):
                    pass
                elif 'min' in name:
                    val = np.nanmin([float(value1), float(value2)])
                elif 'max' in name:
                    val = np.nanmax([float(value1), float(value2)])
                elif 'average' in name:
                    val = np.nanmean([float(value1), float(value2)])

            elif any(ext in name for ext in self.get_additive_field_strings()):
                val = '+'.join([getattr(self,name),getattr(other,name)])

            elif any(ext in name for ext in self.get_static_field_strings()):
                val = getattr(self, name)
            if name=='crystal_direction':
                others_xstal = getattr(other,name)
                this_xstal   = getattr(self,name)
                if others_xstal=='isotropic' and this_xstal=='isotropic':
                    val = 'isotropic'
                elif others_xstal!='isotropic' and this_xstal!='isotropic' and others_xstal!=this_xstal:
                    val = 'isotropic'
                elif any(x!='isotropic' for x in [others_xstal,this_xstal]):
                    val = others_xstal if others_xstal!='isotropic' else this_xstal
                else:
                    pass
            new_params[name]=val
        new_params['phase_type'] = self.phase_type
        return PublicationMetadata.from_kwargs(**new_params)


class ModelInterface:

    def __init__(self,id_rows: pd.Series or PublicationMetadata):
        if id_rows is not None:
            if isinstance(id_rows, pd.Series):
                self.metadata = PublicationMetadata.from_kwargs(**id_rows)
            else:
                self.metadata = id_rows

    def get_conductivity(self, **kwargs):
        pass


    @staticmethod
    def determine_starting_c_type(starting_value=0, P=None, T=None,logfo2=None,Xfe=None,Cw=None,co2=None, **kwargs):
        ndarrays = list(filter(lambda x: isinstance(x, np.ndarray), [P, T, logfo2,Xfe,Cw,co2]))
        if len(ndarrays)<1:
            reference_shape=(1)
        elif len(ndarrays)==1:
            reference_shape= ndarrays[0].shape
        else:
            reference_shape = ndarrays[0].shape
            assert all(array.shape == reference_shape for array in ndarrays),  f"unsure how to construct c since dimensions of T, P, and logfo2 don't match:\nP:\t{P}\nT:\t{T}\nlogfo2:\t{logfo2}"
        return np.ones(reference_shape)*starting_value

    def get_row(self):
        row_dict = vars(self.metadata)
        row_dict['equation']=str(self)
        row_dict['ec_model']=self
        return pd.Series(row_dict)

    def __add__(self, other):
        new_metadata = other.metadata + self.metadata
        return CompositeModel([self, other],new_metadata)

    @property
    def get_id(self):
        return self.metadata.entry_id

class Model(ModelInterface):
    def __init__(self, mechanisms: list, id_rows: pd.Series or PublicationMetadata):
        super().__init__(id_rows)
        self.mechanisms = mechanisms
        for m in mechanisms:
            if m.uses_water:
                m.set_water_units(self.metadata.water_units)
    def __repr__(self):
        if len(self.mechanisms)==1:
            return self.get_id + ':{'+str(self.mechanisms[0])+'}'
        else:
            return self.get_id + ':{' +'+'.join(str(m) for m in self.mechanisms) + '}'
    @property
    def uses_water(self):
        return any(m.uses_water for m in self.mechanisms)

    def get_conductivity(self, **kwargs):
        c = self.determine_starting_c_type(**kwargs)
        for m in self.mechanisms:
            new_conductivity = m.get_conductivity(**kwargs)
            c += new_conductivity
        return c

    def get_crystal_direction(self):
        return self.metadata.crystal_direction


class CompositeModel(ModelInterface):

    def __init__(self, models, metadata):
        super().__init__(metadata)
        self.models = models

    def get_conductivity(self, crystal_direction=None, **kwargs):
        c = self.determine_starting_c_type(**kwargs)
        for m in self.models:
            new_conductivity = m.get_conductivity(**kwargs)
            c += new_conductivity
        return c

    @property
    def uses_water(self):
        return any(m.uses_water for m in self.models)

    def __repr__(self):
        return self.get_id + ':{' +'+'.join(str(m) for m in self.models) + '}'

class AnisotropicModel(ModelInterface):

    def __init__(self, models):
        super().__init__(sum([m.metadata for m in models[1:]], models[0].metadata))
        self.metadata.entry_id='isotropic_model:' +  self.metadata.entry_id
        self.models = {m.get_id: m for m in models}
        # set up metadata for anisotropic model averages

    def get_conductivity(self, crystal_direction=None, averaging='geometric', **kwargs):
        if crystal_direction is not None and any(crystal_direction in x.get_crystal_direction() for x in self.models.values()):
            return self._get_xstal_direction(crystal_direction,**kwargs)
        elif isinstance(crystal_direction,float):
            return self._get_anisotropic_factor(crystal_direction=crystal_direction,**kwargs)
        elif averaging=='geometric':
            return self._get_geometric_mean(**kwargs)
        else:
            raise NotImplementedError("current option is not implemented")

    def _get_xstal_direction(self, xstal_direction,**kwargs):
        for id, model in self.models.items():
            if xstal_direction == model.get_crystal_direction():
                return model.get_conductivity(**kwargs)
        raise NotImplementedError("crystal direction does not exist")

    def _get_geometric_mean(self,**kwargs):
        c = self.determine_starting_c_type(starting_value=1,**kwargs)
        for m in self.models.values():
            new_conductivity = m.get_conductivity(**kwargs)
            new_conductivity = new_conductivity**(1/len(self.models))
            c *= new_conductivity
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
    def __repr__(self):
        return 'mixture of:' + '*'.join(str(k) for k,m in self.models.items())

