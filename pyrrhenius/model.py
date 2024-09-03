from functools import reduce
import numpy as np
import pandas as pd
from dataclasses import dataclass
from inspect import signature
from . import  mechanisms as mech
from . import utils


def _id_row_rename(col):
    col = col.replace(' ', '_').replace('(', '').replace(')', '')
    return col.lower()

def create_model_from_row(row: pd.Series or pd.DataFrame):
    """generates a Model object from a pandas series row

    Parameters
    ----------
    row : pd.Seriesorpd.DataFrame
        a pandas row representing model parameters

    Returns
    -------
    Model
        a Pyrrhenius model
    """
    if not isinstance(row,pd.Series):
        assert len(row)==1,'Make multi model via an add operation. Function only makes single row models'
        row = row.iloc[0]
    row.rename(_id_row_rename, inplace=True)
    mechanism = mech.create_mechanism_from_row(row)
    return Model(mechanism, row)


@dataclass
class PublicationMetadata:

    """A Metadata class holding experimental information on a pyrrhenius model

    Parameters
    ----------
    title : str
        The title of the publication.
    author : str
        The author(s) of the publication. 
    year : str
        The year the publication was published.
    doi : str
        The DOI (Digital Object Identifier) of the publication.
    phase_type : str
        The type of phase the model applies to.
    description : str
        A description of the publication or model.
    sample_type : str
        The type of sample used in the experiments.
    equation_form : str
        The form of the equation used in the model.
    publication_id : str
        A unique identifier for the publication.
    entry_id : str
        A unique identifier for the specific model entry.
    complete_or_partial_fit : str
        Indicates if the model is a complete or partial fit to the data.
    composite_or_single : str
        Indicates if the model is a composite of multiple models or a single model.
    pressure_average : float
        The average pressure condition of the experiments. Units unspecified.
    pressure_min : float
        The minimum pressure condition of the experiments. Units unspecified.
    pressure_max : float
        The maximum pressure condition of the experiments. Units unspecified.
    temp_min : float
        The minimum temperature condition of the experiments in Kelvin.
    temp_max : float 
        The maximum temperature condition of the experiments in Kelvin.
    nacl_min : float
        The minimum NaCl concentration of the experiments. Units unspecified.
    nacl_max : float
        The maximum NaCl concentration of the experiments. Units unspecified.
    nacl_average : float
        The average NaCl concentration of the experiments. Units unspecified.
    na2o_min : float
        The minimum Na2O concentration of the experiments. Units unspecified. 
    na2o_max : float
        The maximum Na2O concentration of the experiments. Units unspecified.
    na2o_average : float
        The average Na2O concentration of the experiments. Units unspecified.
    sio2_min : float
        The minimum SiO2 concentration of the experiments. Units unspecified.
    sio2_max : float
        The maximum SiO2 concentration of the experiments. Units unspecified.
    sio2_average : float
        The average SiO2 concentration of the experiments. Units unspecified.
    co2_min : float
        The minimum CO2 concentration of the experiments. Units unspecified.
    co2_max : float
        The maximum CO2 concentration of the experiments. Units unspecified.
    co2_average : float
        The average CO2 concentration of the experiments. Units unspecified.
    water_min : float
        The minimum water concentration of the experiments. Units given by water_units.
    water_max : float
        The maximum water concentration of the experiments. Units given by water_units.
    water_average : float
        The average water concentration of the experiments. Units given by water_units.
    water_calibration : str
        Information on how the water concentration was calibrated or determined.
    water_units : str
        The units of the water concentration. 
    iron_min : float
        The minimum iron concentration of the experiments. Units given by iron_units.
    iron_max : float
        The maximum iron concentration of the experiments. Units given by iron_units.
    iron_average : float
        The average iron concentration of the experiments. Units given by iron_units.
    iron_units : str
        The units of the iron concentration.
    crystal_direction : str
        The crystal direction or orientation of the sample.

    Methods
    -------
    single_pressure()
        Returns True if pressure_min equals pressure_max, False otherwise.
    single_water() 
        Returns True if water_min equals water_max, False otherwise.
    single_nacl()
        Returns True if nacl_min equals nacl_max, False otherwise.
    single_co2()
        Returns True if co2_min equals co2_max, False otherwise.
    single_na2o()
        Returns True if na2o_min equals na2o_max, False otherwise.
    single_sio2()
        Returns True if sio2_min equals sio2_max, False otherwise.
    single_iron()
        Returns True if iron_min equals iron_max, False otherwise.
    convert2ppm(water_average, water_units)
        Converts the water_average to ppm based on the water_units.
    get_water_average_ppm()
        Returns the water_average converted to ppm using the water_units.
    recalc_average_fields()
        Recalculates the average fields from the min and max fields.
    from_kwargs(**kwargs)
        Creates a PublicationMetadata instance from the provided keyword arguments.

    Raises
    ------
    NotImplementedError
        if the water unit conversion suggested is not supported. 
    """
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
    nacl_min: float
    nacl_max: float
    nacl_average: float
    na2o_min: float
    na2o_max: float
    na2o_average: float
    sio2_min: float
    sio2_max: float
    sio2_average: float
    co2_min: float
    co2_max: float
    co2_average: float
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

    def single_pressure(self):
        return self.pressure_min==self.pressure_max

    def single_water(self):
        return self.water_min==self.water_max

    def single_nacl(self):
        return self.nacl_min==self.nacl_max

    def single_co2(self):
        return self.co2_min==self.co2_max

    def single_na2o(self):
        return self.na2o_min==self.na2o_max

    def single_sio2(self):
        return self.sio2_min==self.sio2_max

    def single_iron(self):
        return self.iron_min==self.iron_max

    def convert2ppm(self,water_average,water_units):
        if water_units == 'wtpct':
            return water_average * 1e4
        elif water_units == 'wtpct10x':
            return water_average * 1e5
        elif water_units == 'ppm':
            return water_average
        elif water_units == 'total_frac':
            return water_average*1e6
        else:
            raise NotImplementedError(
                f"{water_units} conversion not implemented"
            )
    def get_water_average_ppm(self):
        if self.water_units == 'wtpct':
            return self.water_average * 1e4
        elif self.water_units == 'wtpct10x':
            return self.water_average * 1e5
        elif self.water_units == 'ppm':
            return self.water_average
        elif self.water_units == 'total_frac':
            return self.water_average*1e6
        else:
            raise NotImplementedError(
                f"{self.water_units} conversion not implemented"
            )


    def recalc_average_fields(self):
        avg_fields = ['iron_average','water_average','pressure_average','nacl_average','co2_average','na2o_average','sio2_average']
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
        return ['iron','water','pressure','temp','nacl','co2','na2o','sio2']
    @classmethod
    def get_additive_field_strings(cls):
        return ['author','entry_id','publication_id','crystal_direction']

    @classmethod
    def get_units_and_calibrations(cls):
        return ['water_calibration','water_units','iron_units']
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

    def __repr__(self):
        properties = filter(lambda x: str(x[0])!='_' and x!='entry_id',self.__dict__.keys())

        return self.entry_id + '\n:' +'\n\t'.join(k+':'+str(getattr(self,k)) for k  in properties)

    def __add__(self,other):
        """
        Adds two PublicationMetadata instances together.

        Combines the metadata fields based on the following rules:
        - Float fields are combined by taking the min of mins, max of maxs, mean of avgs
        - Additive string fields are concatenated with '+'
        - Static string fields are taken from self 
        - Differing units/calibrations print a warning and self's value is used
        - crystal_direction is set to 'isotropic' if they differ and neither is 'isotropic'

        Parameters
        ----------
        other : PublicationMetadata
            The other PublicationMetadata instance to add.

        Returns
        -------
        PublicationMetadata
            A new PublicationMetadata instance with the combined fields.

        """
        if self.phase_type!=other.phase_type:
            print("Trying to add different phases! You probably shouldn't combine them.")
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
                val = '+'.join([str(getattr(self,name)),str(getattr(other,name))])

            elif any(ext in name for ext in self.get_static_field_strings()):
                val = getattr(self, name)

            elif any(ext in name for ext in self.get_units_and_calibrations()):
                for unit_cal in self.get_units_and_calibrations():
                    value1 = getattr(self,unit_cal)
                    value2 = getattr(other,unit_cal)
                    print(f'{value1} {value2}')
                    if value1 is None and value2 is None:
                        val = np.nan
                    elif (isinstance(value1,str) and isinstance(value1,str)) and value1==value2:
                        val = value1
                    elif (isinstance(value1,str) and isinstance(value1,str)) and value1!=value2:
                        print(f'merging models with different {unit_cal}:\n\t{value1} {value2}')
                        print(f'this is not advised. However, operation will continue using {value1}')
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
    """
    A base interface class for all models in the Pyrrhenius framework.

    Parameters
    ----------
    id_row : pd.Series or PublicationMetadata
        Metadata associated with the model. This can either be a pandas Series or a PublicationMetadata object.

    Attributes
    ----------
    metadata : PublicationMetadata
        Metadata associated with the model.
    mechanism : mech.Mechanism
        The mechanism used in the model.

    Methods
    -------
    get_conductivity(**kwargs)
        Calculates the conductivity based on the model's mechanism and provided parameters.
    print_required_parameters()
        Calculates the conductivity based on the model's mechanism and provided parameters.
    uses_water
        Checks if the model uses water as a parameter.
    uses_nacl
        Checks if the model uses NaCl as a parameter.
    uses_iron
        Checks if the model uses iron as a parameter.
    uses_na2o
        Checks if the model uses Na2O as a parameter.
    uses_sio2
        Checks if the model uses SiO2 as a parameter.
    uses_co2
        Checks if the model uses CO2 as a parameter.
    uses_pressure
        Checks if the model uses pressure as a parameter.
    uses_fo2
        Checks if the model uses fo2 as a parameter.
    determine_starting_c_type(starting_value=0, **kwargs)
        generates a ``float`` or ``np.ndarray`` of a shape consistent with the provided keyword arguments

    """
    def __init__(self,id_row: pd.Series or PublicationMetadata):
        if id_row is not None:
            if isinstance(id_row, pd.Series):
                self.metadata = PublicationMetadata.from_kwargs(**id_row)
            else:
                self.metadata = id_row

    def print_required_parameters(self):
        """Prints the keywords required to use get_conductivity()
        """
        print('*'*20)
        print('Required Keywords')
        print('*'*20)
        print('\nIntrinsic Arguments')
        print('-'*20)
        for test_func, name in zip([self.uses_pressure],['P']):
            if test_func:
                print('*'+name)

        print('\nVolatile Arguments')
        print('-'*20)
        for test_func, name in zip([self.uses_water,
                                    self.uses_co2],['Cw','co2']):
            if test_func:
                print('*'+name)
        print('\nElectronic Conduction Arguments')
        print('-'*20)
        for test_func, name in zip([self.uses_iron,
                                    self.uses_fo2],
                                    ['X_fe','logfo2']):
            if test_func:
                print('*'+name)
        print('\nDiffusion Conduction Arguments')
        print('-'*20)
        for test_func, name in zip([self.uses_nacl,
                                    self.uses_na2o],
                                    ['nacl','na2o']):
            if test_func:
                print('*'+name)
        print('\nPolymerizing Agent Arguments')
        print('-'*20)
        for test_func, name in zip([self.uses_sio2],
                                    ['sio2']):
            if test_func:
                print('*'+name)
       

    def get_conductivity(self,**kwargs):
        """
        Calculates the conductivity based on the model's mechanism and provided parameters.

        Parameters
        ----------
        **kwargs : dict
            Additional parameters for the conductivity calculation. See the documentation on mech.Mechanism for a full list of possible keywords

        Returns
        -------
        float or np.ndarray
            The calculated conductivity in S/m
        """
        start= self.determine_starting_c_type(starting_value=1,**kwargs)
        start*=self.mechanism.get_conductivity(**kwargs)
        return start
    
    @property
    def get_id(self):
        """gets the model id

        Returns
        -------
        str
            the model id
        """
        return self.metadata.entry_id

    @property
    def uses_water(self):
        """
        Checks if the model uses water as a parameter.

        Returns
        -------
        bool
            True if the model uses water, False otherwise.
        """
        return self.mechanism.uses_water
    
    @property
    def uses_nacl(self):
        """
        Checks if the model uses NaCl as a parameter.

        Returns
        -------
        bool
            True if the model uses NaCl, False otherwise.
        """
        return self.mechanism.uses_nacl
    
    @property
    def uses_iron(self):
        """
        Checks if the model uses iron as a parameter.

        Returns
        -------
        bool
            True if the model uses iron, False otherwise.
        """
        return self.mechanism.uses_iron
    
    @property
    def uses_na2o(self):
        """
        Checks if the model uses Na2O as a parameter.

        Returns
        -------
        bool
            True if the model uses Na2O, False otherwise.
        """
        return self.mechanism.uses_na2o

    @property
    def uses_sio2(self):
        """
        Checks if the model uses SiO2 as a parameter.

        Returns
        -------
        bool
            True if the model uses SiO2, False otherwise.
        """
        return self.mechanism.uses_sio2
    
    @property
    def uses_co2(self):
        """
        Checks if the model uses CO2 as a parameter.

        Returns
        -------
        bool
            True if the model uses CO2, False otherwise.
        """
        return self.mechanism.uses_co2
    
    @property
    def uses_pressure(self):
        """
        Checks if the model uses pressure as a parameter.

        Returns
        -------
        bool
            True if the model uses pressure, False otherwise.
        """
        return self.mechanism.uses_pressure

    @property
    def uses_fo2(self):
        """
        Checks if the model uses fo2 as a parameter.

        Returns
        -------
        bool
            True if the model uses fo2, False otherwise.
        """
        return self.mechanism.uses_fo2
    
    @staticmethod
    def determine_starting_c_type(starting_value=0, P=None, T=None,logfo2=None,Xfe=None,Cw=None,co2=None,nacl=None,na2o=None,sio2=None,**kwargs):
        """
        Generates a ``float`` or ``np.ndarray`` of a shape consistent with the provided keyword arguments

        Currently checks keywords for P, T, logfo2, Xfe, Cw, co2, nacl, na2o, and sio2

        Parameters
        ----------
        starting_value : float, optional
            The variable value for conductivity casting, defaults to 0.

        **kwargs : dict
            Additional parameters used to calculate conductivity

        Returns
        -------
        float or np.ndarray
            The determined dummy conductivity variable used to cast the ``get_conductivity`` result into the proper shape.
        """
        ndarrays = list(filter(lambda x: isinstance(x[1], np.ndarray) and x[1] is not None, zip(['P', 'T', 'logfo2','Xfe','Cw','co2','nacl','na2o','sio2'],
                                                                                          [P, T, logfo2,Xfe,Cw,co2,nacl,na2o,sio2])))
        if len(ndarrays)<1:
            reference_shape=(1)
        elif len(ndarrays)==1:
            reference_shape= ndarrays[0][1].shape
        else:
            reference_shape = ndarrays[0][1].shape
            if not all(array[1].shape == reference_shape for array in ndarrays):
                starting_array =None
                for label, array in ndarrays:
                    if starting_array is None:
                        starting_array=np.copy(array)
                    elif array is not None:
                        try:
                            starting_array*=array
                        except Exception as e:
                            print('passed arrays will not multiply together\n')
                            for label, array in ndarrays:
                                if array is not None:
                                    print(f'\t{label}: {array.shape}')
        return np.ones(reference_shape)*starting_value

    def get_row(self):
        """ returns a pd.Series representation of the model

        Returns
        -------
        pd.Series
            a pd.Series representation of the model
        """
        row_dict = vars(self.metadata)
        row_dict['equation']=str(self)
        row_dict['ec_model']=self
        row_dict['entry_id']=self.get_id
        return pd.Series(row_dict)

    def __add__(self, other):
        """
        Adds two ModelInterface instances together.
        
        metadata is combined following the ``__add__`` methodogy of the PublicationMetadata class
        The resulting returned CompositeModel object has the same functionality of ModelInterface, but with the
        overall conductivity now representing the linear combination of two model conductivities

        Parameters
        ----------
        other : ModelInterface
            The other ModelInterface instance to add.

        Returns
        -------
        CompositeModel
            A new CompositeModel with combined properties of both objects
    """
        new_metadata = other.metadata + self.metadata
        return CompositeModel([self, other],new_metadata)

    def __repr__(self):
        """
        Returns a string representation of the ModelInterface instance.
        The representation includes the entry ID and a dictionary of the model's mechanisms.

        Returns
        -------
        str
            A string representation of the ModelInterface instance.
        """
        return self.get_id + ':{' + str(self.mechanism) + '}'

   

    def generate_representative_conditions(self,use_qfm=True):
        """generates a dictionary of conditions relevant to the experiment

        Parameters
        ----------
        use_qfm : bool, optional
            use the quartz-fayallite-magnetite buffer to calculate fo2 conditions. Currently will not work any other way

        Returns
        -------
        dict :
            dictionary consisting of keys matching the model inputs and values consisting of either 1. a range of corresponding conditions or 2. a single average condition.
            1 is returned if the min and max values are different. 2. otherwise
        """
        callouts = {}
        callouts['T']=[self.metadata.temp_min,self.metadata.temp_max]
        # if uses pressure
        if self.uses_pressure:
            if self.metadata.single_pressure():
                callouts['P']=self.metadata.pressure_average
            else:
                callouts['P']=[self.metadata.pressure_min,self.metadata.pressure_max]
        if self.uses_water:
            water_units = self.metadata.water_units
            if self.metadata.single_water():
                water_average = self.metadata.water_average
                ppm = self.metadata.convert2ppm(water_average,water_units)
                callouts['Cw']=ppm
            else:
                water_min = self.metadata.water_min
                water_max = self.metadata.water_max
                ppm_min= self.metadata.convert2ppm(water_min,water_units)
                ppm_max= self.metadata.convert2ppm(water_max,water_units)
                callouts['Cw']=[ppm_min,ppm_max]
        if self.uses_iron:
            if self.metadata.single_iron():
                callouts['X_fe']=self.metadata.iron_average
            else:
                fe_min = self.metadata.iron_min
                fe_max = self.metadata.iron_max
                callouts['X_fe']=[fe_min,fe_max]

        if self.uses_fo2:
            if 'P' not in callouts.keys():
                fo2_range = utils.calc_QFM(np.asarray(callouts['T']),1)
            else:
                fo2_range = utils.calc_QFM(np.asarray(callouts['T']),np.asarray(callouts['P']))
            callouts['logfo2']= fo2_range
        if self.uses_co2:
            if self.metadata.single_co2():
                callouts['co2']=self.metadata.co2_average
            else:
                co2_min = self.metadata.co2_min
                co2_max = self.metadata.co2_max
                callouts['co2']=[co2_min, co2_max]
        if self.uses_nacl:
            if self.metadata.single_nacl():
                callouts['nacl']=self.metadata.nacl_average
            else:
                min_nacl = self.metadata.nacl_min
                max_nacl = self.metadata.nacl_max
                callouts['nacl']=[min_nacl,max_nacl]

        if self.uses_na2o:
            if self.metadata.single_na2o():
                callouts['na2o']=self.metadata.na2o_average
            else:
                min_nacl = self.metadata.na2o_min
                max_nacl = self.metadata.na2o_max
                callouts['na2o']=[min_nacl,max_nacl]

        if self.uses_sio2:
            if self.metadata.single_sio2():
                callouts['sio2']=self.metadata.sio2_average
            else:
                min_nacl = self.metadata.sio2_min
                max_nacl = self.metadata.sio2_max
                callouts['sio2']=[min_nacl,max_nacl]

        for k in callouts.keys():
            if isinstance(callouts[k],list):
                callouts[k]=np.asarray(callouts[k])
        return callouts



class Model(ModelInterface):
    """
    Base Model class used to hold both metadata and mechanisms

    Parameters
    ----------
    mechanism : mech.Mechanism
        The mechanism associated with the model.
    id_row : pd.Series or PublicationMetadata
        The metadata or identifier row for the model.

    Attributes
    ----------
    mechanism : mech.Mechanism
        The mechanism associated with the model.
    metadata : PublicationMetadata
        Metadata associated with the model.

    Methods
    -------
    get_crystal_direction()
        Returns the crystal direction from the metadata.

    """
    def __init__(self, mechanism: mech.Mechanism, id_row: pd.Series or PublicationMetadata):
        super().__init__(id_row)
        self.mechanism = mechanism
        self.mechanism.set_water_units(self.metadata.water_units)
        
    def get_crystal_direction(self):
        """returns the crystal direction for this model

        Returns
        -------
        str
            a 4 character representation of the crystal direction ``[xxx]``, ex: [011]
        """
        return self.metadata.crystal_direction


class CompositeModel(ModelInterface):


    def __init__(self, models, metadata):
        super().__init__(metadata)
        self.models = [(model.get_id,model) for model in models]

    @property
    def uses_water(self):
        return reduce(lambda x,y: x[1].uses_water | y[1].uses_water, self.models)

    @property
    def uses_nacl(self):
        return reduce(lambda x,y: x[1].uses_nacl | y[1].uses_nacl, self.models)
    
    @property
    def uses_iron(self):
        return reduce(lambda x,y: x[1].uses_iron | y[1].uses_iron, self.models)
    
    @property
    def uses_na2o(self):
        return reduce(lambda x,y: x[1].uses_na2o | y[1].uses_na2o, self.models)

    @property
    def uses_sio2(self):
        return reduce(lambda x,y: x[1].uses_sio2 | y[1].uses_sio2, self.models)
    
    @property
    def uses_co2(self):
        return reduce(lambda x,y: x[1].uses_co2 | y[1].uses_co2, self.models)
    
    @property
    def uses_pressure(self):
        return reduce(lambda x,y: x[1].uses_pressure | y[1].uses_pressure, self.models)

    @property
    def uses_fo2(self):
        return reduce(lambda x,y: x[1].uses_fo2 | y[1].uses_fo2, self.models)
    
    def get_conductivity(self, crystal_direction=None, **kwargs):
        c = self.determine_starting_c_type(**kwargs)
        for m in self.models:
            new_conductivity = m[1].get_conductivity(**kwargs)
            c += new_conductivity
        return c

    @property
    def get_id(self):
        return '+'.join([x[0] for x in self.models])

    def __repr__(self):
        return  '{' +'+'.join(x[0] for x in self.models) + '}'

class DryModel(ModelInterface):

    def __init__(self, model):
        super().__init__(model.metadata)
        self.model = model

    def get_conductivity(self, **kwargs):
        new_kwargs = {**kwargs}
        new_kwargs['Cw']=0
        return self.model.get_conductivity(**new_kwargs)

    def __repr__(self):
        return ':dry{' + self.model.get_id+'}'

    @property
    def get_id(self):
        return 'dry:{'+self.metadata.entry_id+'}'


class WaterCorrection(ModelInterface):
    """Applies a Correction factor to Cw to any encapsulated model

    """

    def __init__(self, model,correction_factor= 1.0):
        """Creates a WaterCorrection facade around an existing model

        Parameters
        ----------
        model : ModelInterface
            a Pyrrhenius model

        correction_factor : float, optional
            the correction factor for water. Defaults to 1.0
        """
        super().__init__(model.metadata)
        self.c_factor = correction_factor
        self.main_model = model

    def get_conductivity(self, **kwargs):
        """gets the conductivity of the underlying model after locally correcting the water factor 

        """
        new_kwargs = {**kwargs}
        new_kwargs['Cw']*=self.c_factor
        return self.main_model.get_conductivity(**new_kwargs)

    def __repr__(self):
        return f':Wcorr_x_{self.c_factor}'+'{' + self.main_model.get_id + '}'

    @property
    def get_id(self):
        return f':Wcorr_x_{self.c_factor}'+'{' + self.metadata.entry_id + '}'

class CachedModel(ModelInterface):

    def __init__(self, model):
        super().__init__(model.metadata)
        self.cached_c = None
        self.main_model = model

    def reset_cache(self):
        self.cached_c = None

    def set_cache(self,**kwargs):
        self.cached_c =self.main_model.get_conductivity(**kwargs)

    def get_conductivity(self, ignore_cache=False,**kwargs):
        if self.cached_c is None:
            self.set_cache(**kwargs)

        if ignore_cache:
            return_val = self.main_model.get_conductivity(**kwargs)
        else:
            return_val = self.cached_c

        return return_val

    def __repr__(self):
        return str(self.main_model)

    @property
    def get_id(self):
        return self.main_model.get_id()

class WaterPTCorrection(ModelInterface):

    def __init__(self, model,correction_factor_call):
        super().__init__(model.metadata)
        self.c_factor = correction_factor_call
        self.main_model = model

    def get_facade_model_id(self):
        return self.main_model.get_id

    def get_conductivity(self, **kwargs):
        new_kwargs = {**kwargs}
        new_kwargs['Cw']=self.c_factor(**kwargs)
        return self.main_model.get_conductivity(**new_kwargs)

    def __repr__(self):
        return f':Wcorr(P,T)_x_{self.c_factor}'+'{' + self.main_model.get_id + '}'

    @property
    def get_id(self):
        return f':Wcorr(P,T)_x_{self.c_factor}'+'{' + self.metadata.entry_id + '}'


class IsotropicMixture(ModelInterface):

    def __init__(self, models):
        super().__init__(sum([m.metadata for m in models[1:]], models[0].metadata))
        self.metadata.entry_id='isotropic_model:' +  self.metadata.entry_id
        self.metadata.crystal_direction='isotropic'
        self.models = {m.get_id: m for m in models}
        # set up metadata for anisotropic model averages

    def get_conductivity(self, crystal_direction=None, averaging='geometric', **kwargs):
        """
        Calculate and return the conductivity based on the specified crystal direction or averaging method.
    
        Parameters:
        ----------
        crystal_direction : str or float, optional
        The crystallographic direction to query for specific conductivity.
        If a float is provided, it is used to calculate an anisotropic factor between 0 and 1.
        
        averaging : str, optional
            The method to average conductivity. Default is 'geometric'. 
            Options include:
            - 'geometric': Geometric mean of the available crystallographic directions.
            - 'max_aniso': Maximum anisotropic factor (conductivity in the most conductive direction).
            - 'min_aniso': Minimum anisotropic factor (conductivity in the least conductive direction).

        **kwargs : dict
            Additional keyword arguments passed to internal methods.
    
        Returns:
        -------
        float or np.ndarray 
            The computed conductivity based on the specified method or direction.
    
        Raises:
        -------
        NotImplementedError:
            If an unsupported averaging method is provided.
        """
        if  crystal_direction is not None:
            if isinstance(crystal_direction,float) or isinstance(crystal_direction,np.ndarray) or isinstance(crystal_direction,int):
                return self._get_anisotropic_factor(crystal_direction,**kwargs)
            elif any(crystal_direction in x.get_crystal_direction() for x in self.models.values()):
                return self._get_xstal_direction(crystal_direction,**kwargs)
            else:
                raise NotImplementedError("current crystal direction appears malformed. Check documentation")
        elif averaging=='max_aniso':
            return self._get_anisotropic_factor(1,**kwargs)
        elif averaging=='min_aniso':
            return self._get_anisotropic_factor(0,**kwargs)
        elif averaging=='geometric':
            return self._get_geometric_mean(**kwargs)
        else:
            raise NotImplementedError("current option is not implemented. Choose from \'geometric\',  \'max_aniso\', \'min_aniso\', or specify a crystal direction")

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
        conductivities = [ m.get_conductivity(**kwargs) for m in self.models.values()]
        if isinstance(conductivities[0],np.ndarray):
            min_array = np.min(np.stack(conductivities,axis=-1),axis=-1)
            max_array = np.max(np.stack(conductivities,axis=-1), axis=-1)
        else:
            min_array = min(conductivities)
            max_array = max(conductivities)
        return min_array * (1 - factor) + factor * max_array

    @property
    def uses_water(self):
        return any([model.uses_water for model in self.models.values()])
    

    @property
    def uses_nacl(self):
        return any([model.uses_nacl for model in self.models.values()])
    
    @property
    def uses_iron(self):
        return any([model.uses_iron for model in self.models.values()])
    
    @property
    def uses_na2o(self):
        return any([model.uses_na2o for model in self.models.values()])

    @property
    def uses_sio2(self):
        return any([model.uses_sio2 for model in self.models.values()])
    
    @property
    def uses_co2(self):
        return any([model.uses_co2 for model in self.models.values()])
    
    @property
    def uses_pressure(self):
        return any([model.uses_pressure for model in self.models.values()])

    @property
    def uses_fo2(self):
        return any([model.uses_fo2 for model in self.models.values()])
    
    def __repr__(self):
        return 'mixture of:' + '*'.join(str(k) for k,m in self.models.items())
