
from functools import reduce
from typing import List
import numpy as np
from dataclasses import dataclass
from pyfluids import Fluid, FluidsList, Input
from abc import ABC, abstractmethod

import inspect
from . import utils as pyutils
from . import CONSTANTS



class Mechanism(ABC):
    """Base Abstract class for all electrical conductivity MultiMechanism"""


    def __init__(self,**kwargs):
        defaults = {key:value for key, value in CONSTANTS.items()}
        
        # Override defaults with any provided keyword arguments
        defaults.update(kwargs)
        
        # Assign attributes based on the final configuration
        self.__dict__.update(defaults)
        self.repr_properties = []

    def set_water_units(self, units):
        """
        Set the units for water

        Parameters
        ----------
        units : str
            The units the model uses for water. Should be either wtpct or ppm
        """

        self.water_units = units

    def convert_water(self, Cw):
        """
        Convert water value to the equation specific units

        Parameters
        ----------
        cw : float
            The water value to be converted

        Returns
        -------
        float
            The converted water value

        Raises
        ------
        NotImplementedError
            If the specified conversion units are not implemented
        """

        assert Cw is not None, "Water value must be provided"

        if self.water_units == 'wtpct':
            return Cw * 1e-4
        elif self.water_units == 'wtpct10x':
            return Cw * 1e-5
        elif self.water_units == 'ppm':
            return Cw
        elif self.water_units == 'total_frac':
            return Cw*1e-6
        else:
            raise NotImplementedError(
                f"{self.water_units} conversion not implemented"
            )

    def convert_co2(self, co2):
        """
        Convert co2 value to the specified units

        Parameters
        ----------
        co2 : float
            The co2 value to be converted

        Returns
        -------
        float
            The converted co2 value

        Raises
        ------
        NotImplementedError
            If the specified conversion units are not implemented
        """

        assert co2 is not None, "co2 value must be provided"
        return co2

    def convert_pressure(self, P):
        """
        Convert pressure value to eV units

        Parameters
        ----------
        P : float
            The pressure value to be converted in GPa

        Returns
        -------
        float
            The converted pressure value in eV
        """

        return P * self.cm3GPamol_to_ev
    
    def _parameter_assertion(self, variable, name):
        """
        asserts the provided value provided is correct

        Parameters
        ----------
        variable : float or np.ndarray
            the value to be assessedd

        name : str 
            the name of the variable type

        Raises
        -------
        assertionError
            if variable is of a nature not digestible by the Mechanism Child Class
        """
        assert variable is not None, f"{name} not provided. Pass {name} as a keyword-argument"
        assert not np.all(np.isnan(variable)), f"{name} variable contains only NaNs. Passed parameter must contain at least one non NaN value"


    @abstractmethod
    def _get_conductivity(self, **kwargs):
        """
        Calculates the conductivity of the substance.

        Listed here are the available parameters. However, most child classes may only need two or three parameters. See documentation for more details. 

        Parameters
        ----------
        T : float or np.ndarray , optional
            The temperature value (default is None)
        P : float or np.ndarray , optional
            The pressure value in GPa (default is None)
        co2 : float or np.ndarray , optional
            The CO2 value (default is None)
        Cw : float or np.ndarray , optional
            The water value (default is None)
        nacl : float or np.ndarray , optional
            The NaCl concentration in wt% (default is None)
        sio2 : float or np.ndarray , optional
            The SiO2 concentration in wt% (default is None)
        na2o : float or np.ndarray , optional
            The na2o fraction (default is None)
        logfo2 : float or np.ndarray , optional
            The log oxygen fugacity value (default is None)
        Xfe : float or np.ndarray , optional
            The iron fraction value (default is None)
        **kwargs : dict
            Additional keyword arguments

        Returns
        -------
        np.ndarray or float
            The calculated conductivity
        """
        pass
   

    def get_conductivity(self, check_assertions=True,**kwargs):
        """
        Calculate the conductivity of the mechanism.

        Listed here are the available parameters. However, most child classes may only need two or three parameters. See documentation for more details. 

        Parameters
        ----------
        check_assertions : bool, optional
            whether to check input keywords for existence. Defaults to True
        T : float or np.ndarray , optional
            The temperature value (default is None)
        P : float or np.ndarray , optional
            The pressure value in GPa (default is None)
        co2 : float or np.ndarray , optional
            The CO2 value (default is None)
        Cw : float or np.ndarray , optional
            The water value (default is None)
        nacl : float or np.ndarray , optional
            The NaCl concentration in wt% (default is None)
        sio2 : float or np.ndarray , optional
            The SiO2 concentration in wt% (default is None)
        na2o : float or np.ndarray , optional
            The na2o fraction (default is None)
        logfo2 : float or np.ndarray , optional
            The log oxygen fugacity value (default is None)
        Xfe : float or np.ndarray , optional
            The iron fraction value (default is None)
        **kwargs : dict
            Additional keyword arguments

        Returns
        -------
        np.ndarray or float
            The calculated conductivity

        Raises
        ------
        AssertionError
            error is raised if one or more keywords are not provided which are required by the mechanism
        """
        if check_assertions:
            self.parameter_assertion(**kwargs)
        return self._get_conductivity(**kwargs)


    def parameter_assertion(self, T=None, P=None, co2=None, Cw=None, nacl=None, sio2=None, logfo2=None,Xfe=None, **kwargs):
        """convenience method for asserting passed kwargs are provided if needed by the method

        Parameters
        ----------
        T : float or np.ndarray , optional
            The temperature value (default is None)
        P : float or np.ndarray , optional
            The pressure value in GPa (default is None)
        co2 : float or np.ndarray , optional
            The CO2 value (default is None)
        Cw : float or np.ndarray , optional
            The water value (default is None)
        nacl : float or np.ndarray , optional
            The NaCl concentration in wt% (default is None)
        sio2 : float or np.ndarray , optional
            The SiO2 concentration in wt% (default is None)
        logfo2 : float or np.ndarray , optional
            The log oxygen fugacity value (default is None)
        Xfe : float or np.ndarray , optional
            The oxygen fraction value (default is None)
        **kwargs : dict
            Additional keyword arguments

        Raises
        ------
        AssertionError
            error is raised if one or more keywords are not provided which are required by the mechanism
        """
        if self.uses_temperature:
            self._parameter_assertion(T, 'T')
        if self.uses_pressure:
            self._parameter_assertion(P, 'P')
        if self.uses_water:
           self._parameter_assertion(Cw, 'Cw')
        if self.uses_nacl:
            self._parameter_assertion(nacl, 'nacl')
        if self.uses_sio2:
            self._parameter_assertion(sio2, 'sio2')
        if self.uses_co2:
            self._parameter_assertion(co2, 'co2')
        if self.uses_fo2:
            self._parameter_assertion(logfo2, 'logfo2')
        if self.uses_iron:
            self._parameter_assertion(Xfe, 'Xfe')

    @property
    def uses_temperature(self):
        """
        Returns whether the mechanism uses temperature or not

        Returns
        -------
        bool
            True if the mechanism uses temperature, False otherwise
        """
        return False

    @property
    def uses_water(self):
        """
        Returns whether the mechanism uses water or not

        Returns
        -------
        bool
            True if the mechanism uses water, False otherwise
        """
        return False

    @property
    def uses_nacl(self):
        """
        Returns whether the mechanism uses NaCl concentration or not

        Returns
        -------
        bool
            True if the mechanism uses NaCl, False otherwise
        """
        return False
    
    @property
    def uses_sio2(self):
        """
        Returns whether the mechanism uses NaCl concentration or not

        Returns
        -------
        bool
            True if the mechanism uses NaCl, False otherwise
        """
        return False
    
    @property
    def uses_na2o(self):
        """
        Returns whether the mechanism uses na2o concentration or not

        Returns
        -------
        bool
            True if the mechanism uses NaCl, False otherwise
        """
        return False

    @property
    def uses_co2(self):
        """
        Returns whether the mechanism uses co2 or not

        Returns
        -------
        bool
            True if the mechanism uses co2, False otherwise
        """
        return False

    @property
    def uses_iron(self):
        """
        Returns whether the mechanism uses iron or not

        Returns
        -------
        bool
            True if the mechanism uses iron, False otherwise
        """
        return False

    @property
    def uses_pressure(self):
        """
        Returns whether the mechanism uses pressure or not

        Returns
        -------
        bool
            True if the mechanism uses pressure, False otherwise
        """
        return False

    @property
    def uses_fo2(self):
        """
        Returns whether the mechanism uses fo2 or not

        Returns
        -------
        bool
            True if the mechanism uses fo2, False otherwise
        """
        return False
    
    @abstractmethod
    def __repr__(self):
        """String representation of a Mechanism class

        All __repr__ children must implement a routine which creates a human-readable math symbol corresponding to the represented equation. 

        Returns
        -------
        equation_str
            the string equation representation
        """
        pass

    @classmethod
    def n_args(cls):
        """
        Returns the number of parameters needed to define the class instance.

        Returns
        -------
        int
            The number of  positional parameters needed.

        """
        return _count_positional_args_required(cls.__init__)
    
    def __add__(self, other):
        return MultiMechanism(self,other)
    
class MultiMechanism(Mechanism):
    """A combination of two or more mechanisms

    class is used when the overall electric conductivity is a linear combination of multiple mechanisms

    Parameters
    ----------
    mechanism_list : List[Mechanism]
        a list of mechanisms
    """
      

    def __init__(self,mechanism_list: List[Mechanism]):
        super().__init__()
        self.multi_mechanisms = mechanism_list 

    def _get_conductivity(self, **kwargs):
        """
        calculates the conductivity

        Parameters
        ----------
        **kwargs : dict, optional
            parameters required by the object to calculate conductivity.

        Returns
        -------
        np.ndarray or float
            The calculated conductivity
        """
        return sum(mechanism._get_conductivity(**kwargs) for mechanism in self.multi_mechanisms)


   
    def set_water_units(self,*args,**kwargs):
        """Sets the units for water. 

        """
        for mechanism in self.multi_mechanisms:
            mechanism.set_water_units(*args,**kwargs)
        
    @property
    def uses_water(self):
        """
        Returns whether the mechanism uses water or not

        Returns
        -------
        bool
            True if the mechanism uses water, False otherwise
        """
        return any(mechanism.uses_water for mechanism in self.multi_mechanisms)

    @property
    def uses_nacl(self):
        """
        Returns whether the mechanism uses NaCl concentration or not

        Returns
        -------
        bool
            True if the mechanism uses NaCl, False otherwise
        """
        return any(mechanism.uses_nacl for mechanism in self.multi_mechanisms)
    
    @property
    def uses_sio2(self):
        """
        Returns whether the mechanism uses sio2 concentration or not

        Returns
        -------
        bool
            True if the mechanism uses sio2, False otherwise
        """
        return any(mechanism.uses_sio2 for mechanism in self.multi_mechanisms)
    
    
    @property
    def uses_na2o(self):
        """
        Returns whether the mechanism uses Na2O concentration or not

        Returns
        -------
        bool
            True if the mechanism uses Na2O, False otherwise
        """
        return any(mechanism.uses_na2o for mechanism in self.multi_mechanisms)

    @property
    def uses_co2(self):
        """
        Returns whether the mechanism uses co2 or not

        Returns
        -------
        bool
            True if the mechanism uses co2, False otherwise
        """
        return any(mechanism.uses_co2 for mechanism in self.multi_mechanisms)

    @property
    def uses_iron(self):
        """
        Returns whether the mechanism uses iron or not

        Returns
        -------
        bool
            True if the mechanism uses iron, False otherwise
        """
        return any(mechanism.uses_iron for mechanism in self.multi_mechanisms)

    @property
    def uses_pressure(self):
        """
        Returns whether the mechanism uses pressure or not

        Returns
        -------
        bool
            True if the mechanism uses pressure, False otherwise
        """
        return any(mechanism.uses_pressure for mechanism in self.multi_mechanisms)

    @property
    def uses_fo2(self):
        """
        Returns whether the mechanism uses fo2 or not

        Returns
        -------
        bool
            True if the mechanism uses fo2, False otherwise
        """
        return any(mechanism.uses_fo2 for mechanism in self.multi_mechanisms)
    
    def __repr__(self):
        return reduce(lambda x,y: str(x) + ' + ' + str(y), self.multi_mechanisms)

   

    
@dataclass
class StochasticConstant:
    """
    Class for storing Equation Constants. Can be used to generate a random constant value based on provided uncertainties
    
    Attributes
    ----------
    mean : float
        The mean value of the constant
    stdev : float
        The standard deviation of the constant
    type : str
        The type of the constant
    log : bool, optional
        Determines if the constant value should be calculated using a logarithmic scale,
        default is False
    isnan : bool, optional
        Determines if the constant value can be NaN,
        default is False
    
    Methods
    -------
    get_value(self, sample=False, **kwargs)
        Calculates a constant value based on the provided uncertainties
        Parameters:
            sample : bool, optional
                Determines if the value should be sampled randomly based on the mean and standard deviation,
                default is False
            **kwargs : additional keyword arguments, optional
                Additional parameters that can be used in the calculation of the constant,
                not used in the current implementation
        Returns:
            float
                The calculated constant value
        
    __repr__(self)
        Returns a string representation of the StochasticConstant object
        Returns:
            str
                A string representation of the object
        
    """
    
    mean : float
    stdev : float
    type : str
    log: bool = False
    isnan: bool = False
    
    def get_value(self, sample=False, **kwargs):
        """
        Calculates a constant value based on the provided uncertainties
        
        Parameters
        ----------
        sample : bool, optional
            Determines if the value should be sampled randomly based on the mean and standard deviation,
            default is False
        **kwargs : additional keyword arguments, optional
            Additional parameters that can be used in the calculation of the constant,
            not used in the current implementation
        
        Returns
        -------
        float
            The calculated constant value
        """
        if sample:
            val = np.random.normal(loc=float(self.mean), scale=float(self.stdev))
        else:
            val = self.mean

        if self.log:
            return 10**val
        return val
    
    def __repr__(self):
        """
        Returns a string representation of the StochasticConstant object
        
        Returns
        -------
        str
            A string representation of the object
        """
        if self.log:
            return f'10^{self.mean}({self.stdev})'
        else:
            return f'{self.mean}({self.stdev})'


class SingleValue(Mechanism):
    """Electrical conductivity model which represents single values independent of P,T,fO2, or volatiles


    Parameters
    ----------
    value : StochasticConstant
        The value represented as a StochasticConstant

    Returns
    -------
    float 
        the conductivity value in s/m
    """
    def __init__(self,value):
        super().__init__()
        self.value=value

    def _get_conductivity(self,**kwargs):
        """
        Get the conductivity.

        Returns
        -------
        float or np.ndarray
            The conductivity in units of s/m

        """
        return self.value.get_value(**kwargs)
    
    def __repr__(self):
        return str(self.value)
    

    
class ArrheniusSimple(Mechanism):
    """
    Represents a basic Arrhenius mechanism of the form:

    sigma = a * exp(-h / (k * T))

    Parameters
    ----------
    preexp : StochasticConstant
        the pre-exponential factor 

    enthalpy : StochasticConstant
        Activation enthalpy in units of eV.

    Attributes
    ----------
    preexp : StochasticConstant
        Pre-exponential factor for conductivity.

    enthalpy : StochasticConstant
        Enthalpy for conductivity.

    n_constants : int
        static value of 2 represents how many constants this class requires

    Methods
    -------
    get_preexp(**kwargs)
        Get the value of the pre-exponential factor.

    get_enthalpy(*args, **kwargs)
        Get the value of the enthalpy.

    get_conductivity(T=None, **kwargs)
        Get the conductivity.

    """

    def __init__(self, preexp: StochasticConstant, enthalpy: StochasticConstant):
        """
        Initializes the ArrheniusSimple class

        Parameters
        ----------
        preexp : StochasticConstant
            the preexponential constant
        enthalpy : StochasticConstant
            the enthalpy constant
        """
        super().__init__()
        self.preexp = preexp
        self.enthalpy = enthalpy

    def get_preexp(self, **kwargs):
        """
        Get the value of the pre-exponential factor.

        Parameters
        ----------
        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The value of the pre-exponential factor.

        """
        return self.preexp.get_value(**kwargs)

    def get_enthalpy(self, **kwargs):
        """
        Get the value of the enthalpy.

        Parameters
        ----------
        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The value of the enthalpy in units of eV.

        """
        return self.enthalpy.get_value(**kwargs)

    def _get_conductivity(self, T=None, **kwargs):
        """
        Gets the conductivity.

        Parameters
        ----------
        T : float, optional
            The temperature value in units of kelvin.

        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float or np.ndarray
            The conductivity in units of s/m

        """
        preexp = self.get_preexp(T=T, **kwargs)
        enthalpy = self.get_enthalpy(T=T, **kwargs)
        val = preexp * np.exp(-enthalpy / (self.kb * T))
        return val

    def __repr__(self):
        return f'{self.preexp} exp( -{self.enthalpy}/kT)'


class LinkedArrhenius(ArrheniusSimple):
    """
    Abstract Variation of the ArrheniusSimple class which is primarily used for linear volatile-dependent activation and preexponential factors. 
    See Sifre et al. 2014 for an example.

    Ea     = a * exp( b * volatile) + c
    sigma0 = np.exp(Ea * d + e)
    sigma  = sigma * np.exp(Ea/rT)

    Parameters
    ----------
    const1 : StochasticConstant
        Constant parameter 1.
    const2 : StochasticConstant
        Constant parameter 2.
    const3 : StochasticConstant
        Constant parameter 3.
    const4 : StochasticConstant
        Constant parameter 4.
    const5 : StochasticConstant
        Constant parameter 5.

    Attributes
    ----------
    const1 : StochasticConstant
        Constant parameter 1.
    const2 : StochasticConstant
        Constant parameter 2.
    const3 : StochasticConstant
        Constant parameter 3.
    const4 : StochasticConstant
        Constant parameter 4.
    const5 : StochasticConstant
        Constant parameter 5.

    Methods
    -------
    _get_enthalpy(volatile, **kwargs)
        Calculate the value of the enthalpy.
    get_preexp(**kwargs)
        Get the value of the pre-exponential factor.
    get_conductivity(T=None, **kwargs)
        Get the conductivity.

    """

    def __init__(self, const1: StochasticConstant, const2: StochasticConstant, const3: StochasticConstant,
                 const4: StochasticConstant, const5: StochasticConstant):
        super().__init__(None, None)
        self.const1 = const1
        self.const2 = const2
        self.const3 = const3
        self.const4 = const4
        self.const5 = const5
    def _get_volatile(self,**kwargs):
        raise NotImplementedError("need to implement this function to use the mechanism")
    
    def get_enthalpy(self, **kwargs):
        """
        Calculates enthalpy.

        Parameters
        ----------
        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The value of the enthalpy.

        """
        a = self.const1.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        c = self.const3.get_value(**kwargs)
        volatile = self._get_volatile(**kwargs)
        return a * np.exp(b * volatile) + c

    def get_preexp(self,**kwargs):
        """
        calculates the preexponential value

        Parameters
        ----------
        volatile : float or np.ndarray
            the volatile value needed for computation of the model

        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The value of the pre-exponential factor.

        """
        ea = self.get_enthalpy(**kwargs)
        d = self.const4.get_value(**kwargs)
        e = self.const5.get_value(**kwargs)
        return np.exp(ea * d + e)

    def _get_conductivity(self, T=None, **kwargs):
        """
        Gets the conductivity.

        Parameters
        ----------
        T : float, optional
            The temperature value in Kelvin
        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float or np.ndarray
            The conductivity in s/m.

        """
        preexp = self.get_preexp(**kwargs)
        h = self.get_enthalpy(**kwargs) * 1e-3
        return preexp * np.exp(-h / (self.r * T))

    @classmethod
    def n_args(cls):
        """
        Returns the number of parameters needed to define the class instance.

        Returns
        -------
        int
            The number of  positional parameters needed.

        """
        return _count_positional_args_required(LinkedArrhenius.__init__)
        

class LinkedArrheniusWet(LinkedArrhenius):
    """
    Variation of the LinkedArrhenius class that is specifically used for water-dependent conductivity
    represents the following mechanism parameterization:
    Ea     = a * exp( b * volatile) + c
    sigma0 = np.exp(Ea * d + e)
    sigma  = sigma * np.exp(Ea/rT)

    Parameters
    ----------
    *args : StochasticConstant
        Variable number of positional arguments that are StochasticConstant instances.

    Attributes
    ----------
    const1 : StochasticConstant
        Constant parameter 1.
    const2 : StochasticConstant
        Constant parameter 2.
    const3 : StochasticConstant
        Constant parameter 3.
    const4 : StochasticConstant
        Constant parameter 4.
    const5 : StochasticConstant
        Constant parameter 5.

    Properties
    ----------
    uses_water : bool
        Returns True 

    Methods
    -------
    get_enthalpy(Cw=None, **kwargs)
        calculates the enthalpy value based on water contents.

    """

    def __init__(self, *args):
        """
        Initializes the LinkedArrheniusWet class

        Parameters
        ----------
        *args : StochasticConstant
            the arguments to initialize the class
        """
        super().__init__(*args)

    @property
    def uses_water(self):
        """
        Returns True if water is used in the calculations.

        Returns
        -------
        bool
            True if water is used, False otherwise.

        """
        return True

    def _get_volatile(self,Cw=None, **kwargs):
        """

        Gets and converts the value of the volatile needed for this mechanism

        Parameters
        ----------
        Cw : float or np.ndarray

        kwargs : dict [optional]

        Returns
        -------

        """
        cw = self.convert_water(Cw)
        return cw





class LinkedArrheniusCO2(LinkedArrhenius):
    """
    Variation of the LinkedArrhenius class that is specifically used for CO2-dependent conductivity
    represents the following mechanism parameterization:
    Ea     = a * exp( b * volatile) + c
    sigma0 = np.exp(Ea * d + e)
    sigma  = sigma * np.exp(Ea/rT)

    Parameters
    ----------
    *args : StochasticConstant
        Variable number of positional arguments that are StochasticConstant instances.

    Attributes
    ----------
    const1 : StochasticConstant
        Constant parameter 1.
    const2 : StochasticConstant
        Constant parameter 2.
    const3 : StochasticConstant
        Constant parameter 3.
    const4 : StochasticConstant
        Constant parameter 4.
    const5 : StochasticConstant
        Constant parameter 5.

    Properties
    ----------
    uses_co2: bool
        Returns True 

    Methods
    -------
    get_enthalpy(Cw=None, **kwargs)
        calculates the enthalpy value based on water contents.

    """
    def __init__(self, *args):
        """
        Initializes the LinkedArrheniusCO2 class

        Parameters
        ----------
        *args : StochasticConstant
            the arguments to initialize the class
        """
        super().__init__(*args)

    @property
    def uses_co2(self):
        return True

    def _get_volatile(self, co2=None, **kwargs):
        """

        Gets and converts the value of the volatile needed for this mechanism

        Parameters
        ----------
        co2 : float or np.ndarray

        kwargs : dict [optional]

        Returns
        -------

        """
        co2 = self.convert_co2(co2)
        return co2

class SIGMELTSPommier(ArrheniusSimple):
    """
    A class for calculating SiO2, Na2O, Water, P, and T depedent
    magma conductivity. The basics of the class were taken by web scraping
    the SIGMELTS portal: https://calcul-isto.cnrs-orleans.fr/apps/sigmelts/

    Which has a corresponding publication @doi: 10.1016/j.cageo.2011.01.002

    The parameterization is based on a natural-log space regression of the form:

    delta H = a + b / (c*na2O + d) - e * na2O + f * P + g * Cw^2

    ln(sigma) = h + i*siO2 + j*T + k*P + l* natural_log(Cw) + m*na2O + n * Cw + o * delta H / rT

    Where SiO2, Na2O, and Cw are in wt%, P is in MPa, and T is in Kelvin.

    no reported error on constants was available. I modify this class to have a reported error of 0.4 log units.
    Additionally, there are three forms of the constants, piecewise implemented across different silica activities. 

    To make it compatible with Pyrrhenious, four transforms are needed:

    1. conversion of MPa into GPa

    2. conversion of the enthaly term  o * delta H / rT into  o * delta H / kT

    3. conversion of ppm Water inputs to wt% water

    4. conversion of these factors into simple conductivity. 

    Parameters
    ----------
    const1 : StochasticConstant
        standard state enthalpy

    const2 : StochasticConstant
        term1 enthalpy na2o

    const3 : StochasticConstant
        term2 enthalpy na2o

    const4 : StochasticConstant
        term3 enthalpy na2o

    const5 : StochasticConstant
        term4 enthalpy na2o

    const6 : StochasticConstant
        activation volume

    const7 : StochasticConstant
        squared water constant

    const8 : StochasticConstant
        base conductivity

    const9 : StochasticConstant
        silica activity

    const10 : StochasticConstant
        temperature species availability modifier

    const11 : StochasticConstant
        pressure species availability modifier

    const12 : StochasticConstant
        log water factor

    const13 : StochasticConstant
        na2O activity

    const14 : StochasticConstant
        linear water activity

    const15 : StochasticConstant
        static enthalpy modifier

    const16 : upper silica bound

    const17 : lower silica bound

    """
    def __init__(self, const1 : StochasticConstant,const2 : StochasticConstant,const3 : StochasticConstant,const4 : StochasticConstant,
                 const5 : StochasticConstant,const6 : StochasticConstant,const7 : StochasticConstant,const8 : StochasticConstant,
                 const9 : StochasticConstant,const10 : StochasticConstant,const11 : StochasticConstant,const12 : StochasticConstant,
                 const13 : StochasticConstant,const14 : StochasticConstant,const15 : StochasticConstant,const16 : StochasticConstant,
                 const17 : StochasticConstant):
        super().__init__(None, None)
        self.h_h0 = const1.get_value()
        self.h_na2o1 = const2.get_value()
        self.h_na2o2 = const3.get_value()
        self.h_na2o3 = const4.get_value()
        self.h_na2o4 = const5.get_value()
        self.h_dV = const6.get_value()
        self.h_Cw = const7.get_value()
        self.s_c = const8.get_value()
        self.s_sio2 = const9.get_value()
        self.s_T = const10.get_value()
        self.s_P = const11.get_value()
        self.s_logCw = const12.get_value()
        self.s_na2o = const13.get_value()
        self.s_Cw = const14.get_value()
        self.activation_modifier = const15.get_value()
        self.silica_upper = const16.get_value()
        self.silica_lower = const17.get_value()

    def get_enthalpy(self,na2o=None,Cw=None,P=None, **kwargs):
        """
        Get the value of the enthalpy.

        Parameters
        ----------

        na2o : 
            na2o equivalent concentration in wt%

        Cw: 
            water concentration, converted to wt%

        P :
            pressure in GPa

        Returns
        -------
        float or np.ndarray
            The value of the enthalpy in units of eV

        """
     
        linear_terms = self.h_h0 + self.h_na2o4*na2o + self.h_dV*P 
        nonlinear_terms = (self.h_na2o1/(self.h_na2o2*na2o + self.h_na2o3))+ self.h_Cw*self.convert_water(Cw)**2
        
        result = self.activation_modifier * (linear_terms+nonlinear_terms)*self.kb/(self.r*1e3)
        return -result

    def get_preexp(self,sio2=None,Cw=None,na2o=None,P=None,T=None, **kwargs):
        """
        Get the value of the preexponential factor.

        Parameters
        ----------
        sio2 :
            silica equivalent concentration in wt%

        na2o : 
            na2o equivalent concentration in wt%

        Cw: 
            water concentration, converted to wt%

        P :
            pressure in GPa

        T : 
            temperature in K


        Returns
        -------
        float
            The value of the enthalpy in units of eV.

        """
        # for stability, if Cw=0, set it to 0.01%

        if isinstance(Cw,float):
            cw_prime=Cw
            if Cw ==0:
                cw_prime = 1e-2
            
        else:
            cw_prime = Cw.copy()
            cw_prime[Cw==0]=1e-2
        valid_mask =  (sio2 >= self.silica_lower) &  (sio2 < self.silica_upper)
        if not np.asarray(valid_mask).any():
            return 0
        else:
            linear_terms = self.s_c+self.s_sio2*sio2 +self.s_T*T +self.s_P*P*1e3 +self.s_na2o*na2o +self.s_Cw*self.convert_water(cw_prime)

            nonlinear_terms=self.s_logCw*np.log(self.convert_water(cw_prime))
            total = linear_terms+nonlinear_terms
            terms = np.exp(total)
        terms[~valid_mask]=0
        return terms

    @property
    def uses_water(self):
        return True
    
    @property
    def uses_pressure(self):
        return True
    
    @property
    def uses_sio2(self):
        return True
    
    @property
    def uses_na2o(self):
        return True
    

class SIGMELTSGaillaird(ArrheniusSimple):
    """
    A class for calculating Water, P, and T depedent magma conductivity. The basics of the class were taken by web scraping
    the SIGMELTS portal: https://calcul-isto.cnrs-orleans.fr/apps/sigmelts/

    Which has a corresponding publication @doi: 10.1016/j.cageo.2011.01.002

    The parameterization is based on a natural-log space regression of the form:

    delta H = a*ln(water)+b + c*P

    sigma = d*ln(water)+e

    Where Cw is in wt%, P is in GPa, and T is in Kelvin.

    no reported error on constants was available. I modify this class to have a reported error of 0.4 log units.
    Additionally, there are three forms of the constants, piecewise implemented across different silica activities. 

    To make it compatible with Pyrrhenious, four transforms are needed:

    1. conversion of MPa into GPa

    2. conversion of the enthaly term  o * delta H / rT into  o * delta H / kT

    3. conversion of ppm Water inputs to wt% water

    4. conversion of these factors into simple conductivity. 

    Parameters
    ----------
    const1 : StochasticConstant
        sigma_0 offset

    const2 : StochasticConstant
        water fit constant

    const3 : StochasticConstant
        standard state enthalpy

    const4 : StochasticConstant
        water fit constant

    const5 : StochasticConstant
        activation volume


    """
    def __init__(self, const1 : StochasticConstant,const2 : StochasticConstant,const3 : StochasticConstant,const4 : StochasticConstant,
                 const5 : StochasticConstant):
        super().__init__(None, None)
        self.s_0 = const1.get_value()
        self.s_lnCw = const2.get_value()
        self.h_h0 = const3.get_value()
        self.h_lnCw = const4.get_value()
        self.h_dV = const5.get_value()
       
       
        

    def get_enthalpy(self,Cw=None,P=None, **kwargs):
        """
        Get the value of the enthalpy.

        Parameters
        ----------


        Cw: 
            water concentration, converted to wt%

        P :
            pressure in GPa

        Returns
        -------
        float or np.ndarray
            The value of the enthalpy in units of eV

        """

        linear_terms = self.h_h0 + self.h_dV*P 
        nonlinear_terms = self.h_lnCw*np.log(self.convert_water(Cw))
       
        return (linear_terms+nonlinear_terms)*self.kb/(self.r*1e3)

    def get_preexp(self,Cw=None,**kwargs):
        """
        Get the value of the preexponential factor.

        Parameters
        ----------
       
        Cw: 
            water concentration, converted to wt%


        Returns
        -------
        float
            The value of the enthalpy in units of eV.

        """
        return self.s_0 + self.s_lnCw*np.log(self.convert_water(Cw))

    @property
    def uses_water(self):
        return True
    
    @property
    def uses_pressure(self):
        return True
    


class VogelFulcherTammanWet(ArrheniusSimple):
    """
    A class that represents a Vogel-Fulcher-Tamman equation for wet substances.

    sigma = a exp( (b - c Cw^d)/k(T-e))

    Inherits from the ArrheniusSimple class.

    Parameters:
    -----------
    preexp : StochasticConstant
        The pre-exponential factor of the equation.
    enthalpy : StochasticConstant
        The enthalpy of the equation.
    const1 : StochasticConstant
        The constant that affects the water contribution.
    const2 : StochasticConstant
        The exponent of the water contribution.
    const3 : StochasticConstant
        The constant for calculating the t0 parameter.
    **kwargs
        Additional keyword arguments to be passed to the parent class constructor.

    Attributes:
    -----------
    n_constants : int
        The number of constants required for the equation.

    Methods:
    --------
    uses_water()
        Returns True, indicating that water is used in the equation.
    get_enthalpy(Cw=None, T=None, **kwargs)
        Calculates and returns the enthalpy of the equation involving water.
    get_conductivity(T=None, **kwargs)
        Calculates and returns the conductivity of the wet substance.
    """

    def __init__(self, preexp : StochasticConstant, enthalpy : StochasticConstant,const1 : StochasticConstant,
                 const2 : StochasticConstant, const3 : StochasticConstant):
        super().__init__(preexp,enthalpy)
        """
        Initializes the VogelFulcherTammanWet object.

        Parameters
        ----------
        preexp: StochasticConstant
            The pre-exponential factor of the equation.
        enthalpy: StochasticConstant
            The enthalpy of the equation.
        const1: StochasticConstant
            The constant that affects the water contribution.
        const2: StochasticConstant
            The exponent of the water contribution.
        const3: StochasticConstant
            The constant for calculating the t0 parameter.
   
        """
        super().__init__(preexp, enthalpy)
        self.const1 = const1
        self.const2 = const2
        self.const3 = const3

    @property
    def uses_water(self):
        """
        Returns True, indicating that water is used in the equation.

        Returns:
        --------
        bool
            True, indicating that water is used in the equation.
        """
        return True

    def get_enthalpy(self, Cw=None, T=None, **kwargs):
        """
        Calculates and returns the enthalpy of the equation involving water.

        Parameters:
        -----
        Cw : float, optional
            The amount of water involved in the equation. If not provided, the default value is None.
        T : float, optional
            The temperature at which the enthalpy is calculated. If not provided, the default value is None.
        **kwargs : dict
            Additional keyword arguments to be passed to the parent class method.

        Returns:
        -------
        float
            The enthalpy of the equation in eV
        """
        h = ArrheniusSimple.get_enthalpy(self,T=None,Cw=Cw,**kwargs)
        w_effect = self.const1.get_value(**kwargs)
        exponent = self.const2.get_value(**kwargs)
        return h + w_effect* self.convert_water(Cw)**exponent

    def _get_conductivity(self, T=None, **kwargs):
        """
        Calculates and returns the conductivity of the wet substance.

        Parameters:
        -----
        T : float, optional
            The temperature at which the conductivity is calculated. If not provided, the default value is None.
        **kwargs : dict
            Additional keyword arguments to be passed to the parent

        Returns:
        -------
        float or np.ndarray
            The conductivity in s/m
         """
        preexp  = self.get_preexp(T=T, **kwargs)
        enthalpy = self.get_enthalpy(T=T, **kwargs)
        t0 = self.const3.get_value(**kwargs)
        return preexp * np.exp(-enthalpy / (self.kb * (T- t0)))

class ArrheniusFugacity(ArrheniusSimple):
    """ 
    A class that represents a oxygen-fugacity dependent arrhenius preexponential constant

    sigma = a exp(b/kT) fo2^c

    Inherits from the ArrheniusSimple class.

    Parameters:
    -----------
    preexp : StochasticConstant
        The pre-exponential factor of the equation.
    exponent : StochasticConstant
        The exponential factor for the oxygen fugacity
    enthalpy : StochasticConstant
        The enthalpy value of the equation. in eV


    Methods:
    --------
    uses_fo2()
        Returns True, indicating that oxygen fugacity is used in the equation.
    get_enthalpy(Cw=None, T=None, **kwargs)
        Calculates and returns the enthalpy of the equation involving water.
    get_conductivity(T=None, **kwargs)
        Calculates and returns the conductivity of the wet substance.
    """
    def __init__(self,preexp : StochasticConstant, exponent : StochasticConstant, enthalpy : StochasticConstant):
        """
        Initializes the ArrheniusFugacity mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the preexponential constant
        exponent : StochasticConstant
            the exponent constant
        enthalpy : StochasticConstant
            the enthalpy constant
        """
        super().__init__(preexp,enthalpy)
        self.exponent = exponent

    def _get_conductivity(self, logfo2=None, **kwargs):
        """Calculates the conductivity using the Arrhenius equation with an additional factor from the oxygen fugacity.

        Parameters
        ----------
        logfo2: float
            The logarithm of the oxygen fugacity. Must not be None.
        kwargs: dict
            Additional optional keyword arguments.

        Returns
        -------
        float or np.ndarray: 
            The calculated conductivity.

        Raises
        -------
        AssertionError: 
            f logfo2 is None.

        """
        conductivity1 = super()._get_conductivity(logfo2=logfo2,**kwargs)
        return conductivity1 * (10**logfo2) ** self.exponent.get_value(logfo2=logfo2,**kwargs)

    @property
    def uses_fo2(self):
        """
        Returns True since the ArrheniusFugacity model uses oxygen fugacity.

        Returns
        --------
        bool: True.

        """
        return True

    def __repr__(self):
        return f'{self.preexp} exp( -{self.enthalpy}/kT) fO2^{self.exponent}'

@dataclass
class ArrheniusfO2(Mechanism):
    """

    sigma = exp(log10(a) - b/kT + c*log(fo2))

    """
    preexp : StochasticConstant
    enthalpy : StochasticConstant
    const : StochasticConstant

    def _get_conductivity(self, logfo2=None,T=None, **kwargs):
        """Calculates the conductivity

        Parameters
        ----------
        logfo2 : float or np.ndarray, optional
            logfo2 value, by default None
        T : float or np.ndarray, optional
            temperature value, by default None

        Returns
        -------
        conductivity: float or np.ndarray
            conductivity value in S/m
        """
        a = self.preexp.get_value(**kwargs)
        b = self.enthalpy.get_value(**kwargs)
        c = self.const.get_value(**kwargs)
        return 10**(a - b/T + c*logfo2)

    @property
    def uses_fo2(self):
        return True
    
    def __repr__(self):
        return f'exp(log10({self.preexp}) - {self.enthalpy}/kT + {self.const}*log(fo2))'

@dataclass
class ArrheniusfO22(Mechanism):
    """
    represents the following mechanism parameterization:
    sigma = log10(a + b/T + c*log10 fO2)

    """
    preexp: StochasticConstant
    const: StochasticConstant
    enthalpy: StochasticConstant

    def _get_conductivity(self, logfo2=None, T=None, **kwargs):
        """Calculates the conductivity for the ArrheniusfO22 mechanism

        Parameters
        ----------
        logfo2 : float or np.ndarray
            the logfo2 value
        T : float or np.ndarray
            the temperature value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the conductivity
        """
        a = self.preexp.get_value(**kwargs)
        b = self.enthalpy.get_value(**kwargs)
        c = self.const.get_value(**kwargs)
        return 10**(a + b / T + c*logfo2)

    @property
    def uses_fo2(self):
        return True
    
    def __repr__(self):
        return f'10**({self.preexp}) + {self.enthalpy}/T + {self.const}*log(fo2))'

class ArrheniusPreexpParam(ArrheniusSimple):
    """
    represents the following mechanism parameterization:
    sigma = A0(1 - B*P) exp( -h +P*V/kT)

    """
    def __init__(self, preexp1 : StochasticConstant, preexp2 : StochasticConstant,
                 enthalpy : StochasticConstant, volume : StochasticConstant):
        """
        Initializes the ArrheniusPreexpParam mechanism

        Parameters
        ----------
        preexp1 : StochasticConstant
            the first preexponential constant
        preexp2 : StochasticConstant
            the second preexponential constant
        enthalpy : StochasticConstant
            the enthalpy constant
        volume : StochasticConstant
            the volume constant
        """
        super().__init__(None,None)
        self.preexp1 = preexp1
        self.preexp2 = preexp2
        self.enthalpy = enthalpy
        self.volume = volume

    def get_preexp(self,P=None, **kwargs):
        """
        Calculates the preexponential constant for the ArrheniusPreexpParam mechanism

        Parameters
        ----------
        P : float or np.ndarray
            the pressure value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the preexponential constant
        """
        assert P is not None, "Pressure value must be provided"
        a0 = self.preexp1.get_value(P=P,**kwargs)
        b =  self.preexp2.get_value(P=P,**kwargs)
        return a0*(1 - b*P)

    def get_enthalpy(self, P=None,**kwargs):
        """
        Calculates the enthalpy for the ArrheniusPreexpParam mechanism

        Parameters
        ----------
        P : float or np.ndarray
            the pressure value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the enthalpy
        """
        h = self.enthalpy.get_value(**kwargs)
        v = self.volume.get_value(**kwargs)
        return h + v*self.convert_pressure(P)

    @property
    def uses_pressure(self):
        return True

    def __repr__(self):
        return f'{self.preexp1}*(1-{self.preexp2}*P)*exp( -({self.enthalpy}+P*{self.volume})/kT)'

class ArrheniusPressure(ArrheniusSimple):
    """
    represents the following mechanism parameterization:
    sigma = a exp(-(h + PV)/kt)

    """
    n_constants = 3
    def __init__(self,preexp : StochasticConstant,enthalpy : StochasticConstant,volume : StochasticConstant):
        """
        Initializes the ArrheniusPressure mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the preexponential constant
        enthalpy : StochasticConstant
            the enthalpy constant
        volume : StochasticConstant
            the volume constant
        """
        super().__init__(preexp,enthalpy)
        self.volume = volume

    def get_enthalpy(self, *args,P=None, **kwargs):
        """
        Calculates the enthalpy for the ArrheniusPressure mechanism

        Parameters
        ----------
        P : float or np.ndarray
            the pressure value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the enthalpy
        """
        enthalpy = super(ArrheniusPressure, self).get_enthalpy(**kwargs)
        return enthalpy + self.volume.get_value(**kwargs)*self.convert_pressure(P)

    @property
    def uses_pressure(self):
        return True

    def __repr__(self):
        return f'{self.preexp} exp(-({self.enthalpy} + P{self.volume})/kT)'

class IronPreexpEnthalpyArrhenius(ArrheniusSimple):
    """
    represents the following mechanism parameterization:
    sigma = a*Xfe * exp( -(b + c* Xfe^1/3 + P( d + e Xfe)/kT)

    """
    n_constants = 5

    def __init__(self, preexp: StochasticConstant, enthalpy: StochasticConstant,
                 const1: StochasticConstant, volume: StochasticConstant, const2: StochasticConstant):
        """
        Initializes the IronPreexpEnthalpyArrhenius mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the preexponential constant
        enthalpy : StochasticConstant
            the enthalpy constant
        const1 : StochasticConstant
            the first constant
        volume : StochasticConstant
            the volume constant
        const2 : StochasticConstant
            the second constant
        """
        super().__init__(None, None)
        self.preexp = preexp
        self.enthalpy = enthalpy
        self.volume = volume
        self.const1 = const1
        self.const2 = const2

    def get_preexp(self, Xfe=None, **kwargs):
        """
        Calculates the preexponential constant for the IronPreexpEnthalpyArrhenius mechanism

        Parameters
        ----------
        Xfe : float or np.ndarray
            the iron value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the preexponential constant
        """
        a = self.preexp.get_value(**kwargs)
        return a*Xfe

    def get_enthalpy(self, Xfe=None, P=None, **kwargs):
        """
        Calculates the enthalpy for the IronPreexpEnthalpyArrhenius mechanism

        Parameters
        ----------
        Xfe : float or np.ndarray
            the iron value
        P : float or np.ndarray
            the pressure value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the enthalpy
        """
        h = self.enthalpy.get_value(**kwargs)
        v = self.volume.get_value(**kwargs)
        a = self.const1.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        return h + a*Xfe**(1/3) + self.convert_pressure(P)*(v + b*Xfe)

    @property
    def uses_pressure(self):
        return True

    @property
    def uses_iron(self):
        return True
    
    def __repr__(self):
        return f'{self.preexp}*Xfe * exp( ({self.enthalpy} + {self.const1}*Xfe^1/3 + P ({self.volume}+ {self.const2}*Xfe) )/kT)'

class IronWaterArrhenius1(ArrheniusSimple):
    """
    represents the following mechanism parameterization:
    sigma = a (Cw)^b exp( -c/kt) exp(Xfe * d/kt)

    """

    def __init__(self, preexp: StochasticConstant, const: StochasticConstant,
                 enthalpy1: StochasticConstant, enthalpy2: StochasticConstant):
        """
        Initializes the IronWaterArrhenius1 mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the preexponential constant
        const : StochasticConstant
            the constant
        enthalpy1 : StochasticConstant
            the first enthalpy constant
        enthalpy2 : StochasticConstant
            the second enthalpy constant
        """
        super().__init__(preexp, enthalpy1)
        self.enthalpy2 = enthalpy2
        self.const = const


    def get_preexp(self, Xfe=None, Cw=None,T=None, **kwargs):
        """
        Calculates the preexponential constant for the IronWaterArrhenius1 mechanism

        Parameters
        ----------
        Xfe : float or np.ndarray
            the iron value
        Cw : float or np.ndarray
            the water value
        T : float or np.ndarray
            the temperature value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the preexponential constant
        """
        preexp = ArrheniusSimple.get_preexp(self,**kwargs)
        water = self.convert_water(Cw)
        b = self.const.get_value(**kwargs)
        c = self.enthalpy2.get_value(**kwargs)
        upper_enth = Xfe * c/(self.kb * T)
        return preexp * (water**b) * np.exp(upper_enth)

    @property
    def uses_iron(self):
        return True

    @property
    def uses_water(self):
        return True

    def __repr__(self):
        return f'{self.preexp} (Cw)^{self.const} * exp(-{self.enthalpy}.kT)*exp(Xfe * {self.enthalpy2}/kT)'


class IronWaterArrhenius2(ArrheniusSimple):
    """
    represents the following mechanism parameterization:
    sigma = a (Xfe+b) ^c Cw^d exp( -(e+ f(Xfe+b)+ gCw^(h))/kT)

    """


    def __init__(self, preexp: StochasticConstant, const1: StochasticConstant, const2: StochasticConstant, const3: StochasticConstant,
                 enthalpy1: StochasticConstant,const4 : StochasticConstant, const5 : StochasticConstant,const6 : StochasticConstant,
                 const7 : StochasticConstant):
        """
        Initializes the IronWaterArrhenius2 mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the preexponential constant
        const1 : StochasticConstant
            the first constant
        const2 : StochasticConstant
            the second constant
        const3 : StochasticConstant
            the third constant
        enthalpy1 : StochasticConstant
            the enthalpy constant
        const4 : StochasticConstant
            the fourth constant
        const5 : StochasticConstant
            the fifth constant
        const6 : StochasticConstant
            the sixth constant
        const7 : StochasticConstant
            the seventh constant
        """
        super().__init__(None, None)
        self.preexp = preexp
        self.enthalpy = enthalpy1
        self.const1 = const1
        self.const2 = const2
        self.const3 = const3
        self.const4 = const4
        self.const5 = const5
        self.const6 = const6
        self.const7 = const7

    def get_preexp(self, Xfe=None,Cw=None, **kwargs):
        """
        Calculates the preexponential constant for the IronWaterArrhenius2 mechanism

        Parameters
        ----------
        Xfe : float or np.ndarray
            the iron value
        Cw : float or np.ndarray
            the water value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the preexponential constant
        """
        a = self.preexp.get_value(**kwargs)
        iron_offset = self.const1.get_value(**kwargs)
        c = self.const2.get_value(**kwargs)
        d = self.const3.get_value(**kwargs)
        preexp = a * (Xfe + iron_offset) **c * self.convert_water(Cw)**d
        return preexp

    def get_enthalpy(self, Cw=None, Xfe=None, **kwargs):
        """
        Calculates the enthalpy for the IronWaterArrhenius2 mechanism

        Parameters
        ----------
        Cw : float or np.ndarray
            the water value
        Xfe : float or np.ndarray
            the iron value
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the enthalpy
        """
        h     = self.enthalpy.get_value(**kwargs)
        alpha = self.const4.get_value(**kwargs)
        iron_offset   = self.const1.get_value(**kwargs)
        iron_exponent = self.const5.get_value(**kwargs)
        beta          = self.const6.get_value(**kwargs)
        water_exp     = self.const7.get_value(**kwargs)
        water         = self.convert_water(Cw)
        total_enthalpy = h + alpha * (Xfe+iron_offset)**iron_exponent + beta * water**water_exp
        return total_enthalpy

    @property
    def uses_iron(self):
        return True

    @property
    def uses_water(self):
        return True

    def __repr__(self):
        """
        a (Xfe+b) ^c Cw^d exp( -(e+ f(Xfe+b)^g+ hCw^(i))/kT)
        Returns
        -------

        """
        return f'{self.preexp} (Xfe+{self.const1})^{self.const2} (Cw)^{self.const3} \n'+\
               f'exp( -({self.enthalpy} +{self.const4}(Xfe+{self.const1})^{self.const5}+{self.const6} Cw^{self.const7}  )/kT)'

class NerstEinstein1(Mechanism):
    """
    represents the Nernst Einstein Law. 
    sigma =  a Ch exp(b/kT)/kT

    Also used as the base class for all nernst-einstein-like mechanisms. 

    """

    def __init__(self,const1 : StochasticConstant, enthalpy : StochasticConstant):
        """
        Initializes the NerstEinstein1 mechanism

        Parameters
        ----------
        const1 : StochasticConstant
            the first constant
        enthalpy : StochasticConstant
            the enthalpy constant
        """
        super().__init__()
        self.const1 = const1
        self.enthalpy = enthalpy

    def prep_preexponential_constant(self,T=None, Cw=None,**kwargs):
        """provides the preexponential constant for the Nernst Einstein Law

        Parameters
        ----------
        T: float or np.ndarray
            the temperature value, by default None

        Cw : float or np.ndarray, optional
            the water value , by default None

        Returns
        -------
        const : float or np.ndarray 
            the preexponential constant 
        """
        a = self.const1.get_value(**kwargs)
        b = self.convert_water(Cw)
        invt = 1/(self.kb*T)
        return  a * b * invt
    
    def _get_conductivity(self, T=None,**kwargs):
        a = self.prep_preexponential_constant(T=T,**kwargs)
        h = self.enthalpy.get_value(T=T,**kwargs)

        return  a * np.exp(-h/(self.kb*T))


    @property
    def uses_water(self):
        return True
    def __repr__(self):
        return f'{self.const1} Cw * exp({self.enthalpy}/kT)/kT'


class NerstEinstein2(NerstEinstein1):
    """
    represents the Nernst Einstein Law with a different exponent for the water term
    sigma =  a Cw b q^2 exp(-h/kT)/kT

    """


    def __init__(self, const1: StochasticConstant, const2: StochasticConstant, enthalpy : StochasticConstant):
        """
        Initializes the NerstEinstein2 mechanism

        Parameters
        ----------
        const1 : StochasticConstant
            the first constant
        const2 : StochasticConstant
            the second constant
        enthalpy : StochasticConstant
            the enthalpy constant
        """
        super().__init__(const1, enthalpy)
        self.const2 = const2

    def prep_preexponential_constant(self,T=None, Cw=None,**kwargs):
        """provides the preexponential constant for the Nernst Einstein Law

        Parameters
        ----------
        T: float or np.ndarray
            the temperature value, by default None

        Cw : float or np.ndarray, optional
            the water value , by default None

        Returns
        -------
        const : float or np.ndarray 
            the preexponential constant 
        """
        a = self.const1.get_value(**kwargs)
        cw = self.convert_water(Cw)
        b = self.const2.get_value(**kwargs)
        invt = 1/(self.kb*T)
        return  a * b * self.q2 * invt * cw
    
    def __repr__(self):
        return f'{self.const1} Cw * {self.const2} * exp({self.enthalpy}/kT)/kT'


class NerstEinstein3(NerstEinstein1):
    """

    sigma =  a * b Cw^(1+c) exp( d/kT)q^2/kT

    """


    def __init__(self, const1: StochasticConstant, const2: StochasticConstant, const3: StochasticConstant, enthalpy: StochasticConstant):
        """
        Initializes the NerstEinstein3 mechanism

        Parameters
        ----------
        const1 : StochasticConstant
            the first constant
        const2 : StochasticConstant
            the second constant
        const3 : StochasticConstant
            the third constant
        enthalpy : StochasticConstant
            the enthalpy constant
        """
        super().__init__(const1, enthalpy)
        self.const2 = const2
        self.const3 = const3
    
    def prep_preexponential_constant(self,T=None, Cw=None,**kwargs):
        """provides the preexponential constant for the Nernst Einstein Law

        Parameters
        ----------
        T: float or np.ndarray
            the temperature value, by default None

        Cw : float or np.ndarray, optional
            the water value , by default None

        Returns
        -------
        const : float or np.ndarray 
            the preexponential constant 
        """
        a = self.const1.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        c = self.const3.get_value(**kwargs)

        q2 = self.q2 

        cw = self.convert_water(Cw)**(1+c)
        invt = 1/(self.kb*T)

        return  a * b * q2 * invt * cw
    


    @property
    def uses_water(self):
        return True

    def __repr__(self):
        return f'{self.const1} q^2 {self.const2} Cw^({self.const3}+1) exp(-{self.enthalpy}/kT)/kT'


@dataclass
class BrinePTDependent(Mechanism):
    """
    Represents a brine conductivity model that depends on pressure, temperature, and salinity.
    log(sigma) = a + b/T + c log(NaCl) + d log(rho) + log(A0(P,T))
    A0 = e + f*rho + g/T + h/T^2


    """
    a: StochasticConstant
    b: StochasticConstant
    c: StochasticConstant
    d: StochasticConstant
    e: StochasticConstant
    f: StochasticConstant
    g: StochasticConstant
    h: StochasticConstant

    water :  Fluid = Fluid(FluidsList.Water)


    def get_density_field(self,P=None,T=None,**kwargs):
        """
        Calculates the density field for the brine using the water object

        Parameters
        ----------
        P : float or np.ndarray
            the pressure in Pa
        T : float or np.ndarray
            the temperature in Kelvin
        **kwargs : dict
            additional keyword arguments

        Returns
        -------
        float or np.ndarray
            the density of the brine in kg/m^3
        """
        starting_var = pyutils.create_starting_variable(P=P,T=T,**kwargs)
        if isinstance(starting_var,int) or isinstance(starting_var,float):
            self.water.update(Input.pressure(P*1.0e9),Input.temperature(T-273.15))
            rho= self.water.density
        else:
            def calc_density(P,T,water):
                try:
                    water.update(Input.pressure(P * 1.0e9), Input.temperature(T-273.15))
                except:
                    return np.nan
                return water.density

            rho = np.vectorize(calc_density)(P,T,self.water)
        return rho/1e3

    def get_lambda(self,rho,T=None,**kwargs):
        """
        Calculates the lambda term for the BrinePTDependent mechanism

        Parameters
        ----------
        rho : float or np.ndarray
            the density of the brine in kg/m^3
        T : float or np.ndarray
            the temperature in Kelvin

        Returns
        -------
        float or np.ndarray
            the lambda term
        """
        val1 = self.e.get_value(**kwargs) + self.f.get_value(**kwargs)*rho
        val2 = self.g.get_value(**kwargs)/T + self.h.get_value(**kwargs)/T**2
        return val1 + val2
    
    def _get_conductivity(self,T=None,P=None,nacl=None,**kwargs):
        """
        Calculates the conductivity in s/m from this mechanism

        Parameters
        ----------
        T : np.ndarray or float
            the temperature in Kelvin

        logfo2: np.ndarray or float
            The logarithm of the oxygen fugacity. Must not be None.

        Returns
        -------
        float:
            The calculated conductivity in s/m

        Raises
        ------
        AssertionError:
            If logfo2 is None.

        """

        rho = self.get_density_field(P=P,T=T,**kwargs)
        # rho is 0.0012 for ints. 1.24 for arrays
        A0 = self.get_lambda(rho,P=P,T=T,**kwargs)
        # A0 is 1900.42 for ints. 388 for arrays
        a = self.a.get_value(**kwargs)
        # a is -0.9184 for ints
        b = self.b.get_value(**kwargs)
        # b is -872.0 for ints
        c = self.c.get_value(**kwargs)
        # c is 0.852 for ints
        d = self.d.get_value(**kwargs)
        # d is 7.6076 for ints
        # density term is -22.089 for ints
        temp_term    = b/T
        salt_term    = c*np.log10(nacl)
        density_term = d*np.log10(rho)
        return 10**(a + temp_term + salt_term + density_term + np.log10(A0))



    def __repr__(self):
        return f'{self.a} + {self.b}/T + {self.c}log(nacl) + {self.d}log(rho) + log(A0)'

    @property
    def uses_nacl(self):
        return True

    @property
    def uses_pressure(self):
        return True

@dataclass
class SEO3Term(Mechanism):
    """
    Uses a combination of ArrheniusSimple and ArrheniusFugacity models to represent a polaron hopping term for the SEO3 model of olivine
    sigma = e [a exp(b/kT) + c exp(d/kT)fO2∧e] f exp(g/kT)

    Inherits from the Mechanism class.

    Parameters:
    -----------
    preexp1 : StochasticConstant
        The pre-exponential factor for the first arrhenius law
    enth1 : StochasticConstant
        The enthalpy factor for the first arrhenius law
    preexp2 : StochasticConstant
        The pre-exponential factor for the second arrhenius law
    enth2 : StochasticConstant
        The enthalpy for the second arrhenius law
    exp : StochasticConstant
        The oxygen fugacity exponent for the second arrhenius law

    preexp3 : StochasticConstant
        The pre-exponential factor for the third arrhenius law
    enth3 : StochasticConstant
        The enthalpy for the third arrhenius law
    


    Methods:
    --------
    uses_fo2()
        Returns True, indicating that oxygen fugacity is used in the equation.

    get_conductivity(T=None, **kwargs)
        Calculates and returns the conductivity of the wet substance.

    """

    def __init__(self,preexp1 : StochasticConstant, enth1 : StochasticConstant, preexp2 : StochasticConstant,
                 enth2 : StochasticConstant, exp : StochasticConstant, preexp3 : StochasticConstant, enth3 : StochasticConstant):
        """
        Initializes the SEO3Term mechanism

        Parameters
        ----------
        preexp1 : StochasticConstant
            the pre-exponential constant for the first arrhenius law
        enth1 : StochasticConstant
            the enthalpy for the first arrhenius law
        preexp2 : StochasticConstant
            the pre-exponential constant for the second arrhenius law
        enth2 : StochasticConstant
            the enthalpy for the second arrhenius law
        exp : StochasticConstant
            the oxygen fugacity exponent for the second arrhenius law
        preexp3 : StochasticConstant
            the pre-exponential constant for the third arrhenius law
        enth3 : StochasticConstant
            the enthalpy for the third arrhenius law
        """
        super().__init__()
        self.b_species    = ArrheniusSimple(preexp1,enth1)
        self.b_species1   = ArrheniusFugacity(preexp2,exp, enth2)
        self.mu_potential = ArrheniusSimple(preexp3,enth3)

    def _get_conductivity(self, **kwargs):
        """
        Calculates the conductivity in s/m from this mechanism

        Parameters
        ----------
        T : np.ndarray or float
            the temperature in Kelvin

        logfo2: np.ndarray or float
            The logarithm of the oxygen fugacity. Must not be None.

        Returns
        -------
        float: 
            The calculated conductivity in s/m

        Raises
        ------
        AssertionError: 
            If logfo2 is None.

        """
        term1 = self.b_species.get_conductivity(**kwargs)
        term2 = self.b_species1.get_conductivity(**kwargs)
        term3 = self.mu_potential.get_conductivity(**kwargs)
        return self.q * (term1+term2)*term3

    @classmethod
    def n_args(cls):
        return ArrheniusSimple.n_args()*2 + ArrheniusFugacity.n_args()

    def __repr__(self):
        return f'e*({self.b_species} {self.b_species1})*{self.mu_potential}'

    @property
    def uses_fo2(self):
        return True


class WaterExpArrhenius1(ArrheniusSimple):
    """
    Represents an equation with the following parameterization:
    sigma = a Cw^b exp( -c /kT)

    """
    def __init__(self,preexp: StochasticConstant, const1 : StochasticConstant, enthalpy : StochasticConstant):
        """
        Initializes the WaterExpArrhenius1 mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the pre-exponential constant
        const1 : StochasticConstant
            the water exponent
        enthalpy : StochasticConstant
            the enthalpy constant
        """
        super().__init__(preexp,enthalpy)
        self.const1 = const1

    def get_preexp(self,Cw=None, **kwargs):
        """
        Calculates the preexp constant for the WaterExpArrhenius1 mechanism

        Parameters
        ----------
        Cw : float or np.ndarray
            the water concentration in mol/m^3

        Returns
        -------
        float or np.ndarray
            the preexp constant
        """
        a = self.preexp.get_value(**kwargs)
        p = self.const1.get_value(**kwargs)
        return a * self.convert_water(Cw)**p

    @property
    def uses_water(self):
        return True


class WaterExpArrhenius2(ArrheniusSimple):
    """
    Represents an equation with the following parameterization:
    sigma = a Cw exp( -(b + c C_w^(d) /kT)

    """

    def __init__(self, preexp: StochasticConstant, enthalpy: StochasticConstant,const1: StochasticConstant,
                 const2: StochasticConstant):
        """
        Initializes the WaterExpArrhenius2 mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the pre-exponential constant
        enthalpy : StochasticConstant
            the enthalpy constant
        const1 : StochasticConstant
            the water exponent
        const2 : StochasticConstant
            the water enthalpy exponent
        """
        super().__init__(preexp, enthalpy)
        self.const1 = const1
        self.const2 = const2

    def get_preexp(self, Cw=None, **kwargs):
        """
        Calculates the preexp constant for the WaterExpArrhenius2 mechanism

        Parameters
        ----------
        Cw : float or np.ndarray
            the water concentration in mol/m^3

        Returns
        -------
        float or np.ndarray
            the preexp constant
        """
        a = self.preexp.get_value(**kwargs)
        return a * self.convert_water(Cw)

    def get_enthalpy(self, Cw=None, **kwargs):
        """
        Calculates the enthalpy for the WaterExpArrhenius2 mechanism

        Parameters
        ----------
        Cw : float or np.ndarray
            the water concentration in mol/m^3

        Returns
        -------
        float or np.ndarray
            the enthalpy
        """
        a = self.enthalpy.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        c = self.const1.get_value(**kwargs)
        return a + c*self.convert_water(Cw)**(b)

    @property
    def uses_water(self):
        return True

class WaterExpArrhenius3(ArrheniusPressure):
    """
    Represents an equation with the following parameterization:
    sigma = a Cw^b exp( -(c + d C_w^(e) + P f) /kT)

    """

    def __init__(self,preexp : StochasticConstant,  waterexp_preepx : StochasticConstant,
                 enthalpy : StochasticConstant, waterenth_const : StochasticConstant,
                 waterenth_exp : StochasticConstant, volume : StochasticConstant,):
        """
        Initializes the WaterExpArrhenius3 mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the pre-exponential constant
        waterexp_preepx : StochasticConstant
            the water exponent
        enthalpy : StochasticConstant
            the enthalpy constant
        waterenth_const : StochasticConstant
            the water enthalpy constant
        waterenth_exp : StochasticConstant
            the water enthalpy exponent
        volume : StochasticConstant
            the volume constant
        """
        super().__init__(preexp,enthalpy,volume)
        self.const1 = waterexp_preepx
        self.const2 = waterenth_const
        self.const3 = waterenth_exp

    def get_preexp(self, Cw=None, **kwargs):
        """
        Calculates the preexp constant for the WaterExpArrhenius3 mechanism

        Parameters
        ----------
        Cw : float or np.ndarray
            the water concentration in mol/m^3

        Returns
        -------
        float or np.ndarray
            the preexp constant
        """
        a = ArrheniusPressure.get_preexp(self,**kwargs)
        r = self.const1.get_value(**kwargs)
        cw_converted = self.convert_water(Cw)
        return a * cw_converted**r

    def get_enthalpy(self, Cw=None, **kwargs):
        """
        Calculates the enthalpy for the WaterExpArrhenius3 mechanism

        Parameters
        ----------
        Cw : float or np.ndarray
            the water concentration in mol/m^3

        Returns
        -------
        float or np.ndarray
            the enthalpy
        """
        h = ArrheniusPressure.get_enthalpy(self, **kwargs)
        a = self.const2.get_value(**kwargs)
        r = self.const3.get_value(**kwargs)
        cw_converted = self.convert_water(Cw)
        return h + a*cw_converted**r

    @property
    def uses_water(self):
        return True

    def __repr__(self):
        return f'{self.preexp} Cw^{self.const1} exp( -({self.enthalpy} + {self.const2}Cw^{self.const3} + P{self.volume})/kT)'

class WaterExpArrheniusPressure(ArrheniusPressure):
    """
    Represents an equation with the following parameterization:
    sigma = a Cw^b exp(-(c + Pd)/kT)

    """

    def __init__(self, preexp: StochasticConstant, const1: StochasticConstant, enthalpy: StochasticConstant,
                 volume : StochasticConstant):
        """
        Initializes the WaterExpArrheniusPressure mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the pre-exponential constant
        const1 : StochasticConstant
            the water exponent
        enthalpy : StochasticConstant
            the enthalpy constant
        volume : StochasticConstant
            the volume constant
        """
        super().__init__(preexp, enthalpy, volume)
        self.const1 = const1

    def get_preexp(self, Cw=None, **kwargs):
        """
        Calculates the preexp constant for the WaterExpArrheniusPressure mechanism

        Parameters
        ----------
        Cw : float or np.ndarray
            the water concentration in mol/m^3

        Returns
        -------
        float or np.ndarray
            the preexp constant
        """
        a = self.preexp.get_value(**kwargs)
        b = self.const1.get_value(**kwargs)
        return a * self.convert_water(Cw)**b

    @property
    def uses_water(self):
        return True

    @property
    def uses_pressure(self):
        return True

    def __repr__(self):
        return f'Cw^{self.const1} {ArrheniusPressure.__repr__(self)}'

class WaterExpArrheniusInvT(WaterExpArrheniusPressure):
    """
    Represents an equation with the following parameterization:

    sigma = a Cw^b exp(-(c + Pd)/kT)/T

    """

    def __init__(self, preexp: StochasticConstant, const1: StochasticConstant, enthalpy: StochasticConstant, volume: StochasticConstant):
        """
        Initializes the WaterExpArrheniusInvT mechanism

        Parameters
        ----------
        preexp : StochasticConstant
            the pre-exponential constant
        const1 : StochasticConstant
            the water exponent
        enthalpy : StochasticConstant
            the enthalpy constant
        volume : StochasticConstant
            the volume constant
        """
        super().__init__(preexp,const1, enthalpy, volume)

    def _get_conductivity(self, T=None, **kwargs):
        """
        Calculates the conductivity in s/m from this mechanism

        Parameters
        ----------
        T : np.ndarray or float
            the temperature in Kelvin

        Returns
        -------
        c : float or np.ndarray
            the conductivity in s/m
        """
        c = super()._get_conductivity(T=T, **kwargs)
        return c/T

    @property
    def uses_water(self):
        return True

    @property
    def uses_pressure(self):
        return True

    def __repr__(self):
        return f'{WaterExpArrheniusPressure.__repr__(self)}/T'

def _count_positional_args_required(func):
    """
    counts the number of positional arguments required for a function

    Parameters
    ----------
    func : function
        the function to count the positional arguments for

    Returns
    -------
    int
        the number of positional arguments required for the function
    """
    signature = inspect.signature(func)
    empty = inspect.Parameter.empty
    total = -1
    for param in signature.parameters.values():
        if param.default is empty:
            total += 1
    return total


model_dict = {
            'arrhenius_fugacity':ArrheniusFugacity,
'arrhenius_log_fO2':ArrheniusfO2,
'arrhenius_logfO2_2':ArrheniusfO22,
'arrhenius_preexp_parameterized':ArrheniusPreexpParam,
'arrhenius_pressure':ArrheniusPressure,
'arrhenius_simple':ArrheniusSimple,
'iron_preexp_enthalpy_arrhenius':IronPreexpEnthalpyArrhenius,
'iron_water_arrhenius1':IronWaterArrhenius1,
'iron_water_arrhenius2': IronWaterArrhenius2,
'nerst_einstein1':NerstEinstein1,
'nerst_einstein_2':NerstEinstein2,
'nerst_einstien_3':NerstEinstein3,
'seo3_term':SEO3Term,
'single_value':SingleValue,
'water_exponent_arrhenius':WaterExpArrhenius1,
'water_exponent_arrhenius_pressure':WaterExpArrheniusPressure,
'water_exponent_arrhenius_pressure_invT':WaterExpArrheniusInvT,
'water_preexp_enthalpy_arrhenius':WaterExpArrhenius2,
'water_preexp_enthalpy_pressure':WaterExpArrhenius3,
'vft_wet': VogelFulcherTammanWet,
'linked_arrhenius_wet':LinkedArrheniusWet,
'linked_arrhenius_co2':LinkedArrheniusCO2,
'sigmelts_lowwater':SIGMELTSPommier,
'sigmelts_highwater':SIGMELTSGaillaird,
'brinePTdependence':BrinePTDependent}



def _count_positional_args_required(func):
    """
    counts the number of positional arguments required for a function

    Parameters
    ----------
    func : function
        the function to count the positional arguments for

    Returns
    -------
    int
        the number of positional arguments required for the function
    """
    signature = inspect.signature(func)
    empty = inspect.Parameter.empty
    total = -1
    for param in signature.parameters.values():
        if param.default is empty:
            total += 1
    return total



def _create_multiple_MultiMechanism(row):
    """ 
    creates a MultiMechanism from a row in the database

    Parameters
    ----------
    row : pd.Series
        a row in the database

    Returns
    -------
    MultiMechanism
        a MultiMechanism
    """
    potential_MultiMechanism = [x.strip() for x in row['eq_id'].split('+')]
    constants = create_constants_from_row(row)[::-1]
    mech_list = []

    for pot_mech in potential_MultiMechanism:
        n_constants = model_dict[pot_mech].n_args()
        assert n_constants <= len(constants), f"not enough constants for eqid: {row['entry_id']}\n{row}. Found {len(constants)} but should have {n_constants}"
        specific_constants = [constants.pop() for n in range(n_constants)]
        mechanism = model_dict[pot_mech](*specific_constants)
        mech_list.append(mechanism)
    
    return MultiMechanism(mech_list)


def _create_single_mechanism(row):
    """
    creates a single mechanism from a row in the database

    Parameters
    ----------
    row : pd.Series
        a row in the database

    Returns
    -------
    mechanism : Mechanism
    """
    mechanism = row['eq_id'].strip()
    target_mechanism = model_dict[mechanism]
    constants = create_constants_from_row(row)
    assert len(constants) == target_mechanism.n_args(), "Incorrect constant number defined for mechanism. " + \
                                                        "Should have: " + str(target_mechanism.n_args()) + " but found " + str(len(constants)) + " in file\n" + \
                                                        "Check database entry " + row['entry_id'] + " against " + \
                                                        str(target_mechanism)
    mechanism = target_mechanism(*constants)
    return mechanism


def create_mechanism_from_row(row) -> Mechanism:
    """
    creates a mechanism from a row in the database

    Parameters
    ----------
    row : pd.Series
        a row in the database

    Returns
    -------
    mechanism : Mechanism
        a mechanism
    """
    mechanism = row['eq_id']
    if '+' in mechanism:
        return _create_multiple_MultiMechanism(row)
    else:
        return _create_single_mechanism(row)


def get_const_rows(row):
    """
    creates a list of the row keys that are not StochasticConstants

    Parameters
    ----------
    row : pd.Series
        a row in the database

    Returns
    -------
    list of str
        a list of the row keys that are not StochasticConstants
    """
    non_constant_row_keys = ['Title','Author','Year','DOI','Phase Type','Description','Sample Type',
                             'Equation Form','Eq ID','Publication ID','Entry ID','Complete or Partial Fit',
                             'Composite or single','Pressure Average','Pressure Min','Pressure Max','Temp Min',
                             'Temp Max','Water Min','Water Max','Water Average','water calibration','Water units',
                             'Iron Min','Iron Max','Iron Average','Iron units','fO2 buffers used','Crystal Direction','Fitting Comments']
    # insert op here. r
    non_const_keys_lower_set = set(key.lower() for key in non_constant_row_keys)
    
    # Filter for keys that do not have a corresponding case-insensitive match in the non-constant row keys
    potential_row_ids = [key for key in row.keys() if key.lower() not in non_const_keys_lower_set]
    letter_rows = list(filter(lambda x: '_' not in x, potential_row_ids))
    valid_letter_rows = sorted(list(filter(lambda x: ~np.isnan(float(row[x])), letter_rows)),key=lambda x: (len(x), x))
    return valid_letter_rows


def create_constants_from_row(row):
    """
    creates a list of StochasticConstants from a row in the database

    Parameters
    ----------
    row : pd.Series
        a row in the database

    Returns
    -------
    list of StochasticConstants
        a list of StochasticConstants
    """
    letter_columns = get_const_rows(row)
    try:
        variables = [StochasticConstant(row[f'{x}'],
                                                   row[f'{x}_uncertainty'],
                                                   row[f'{x}_description'].split('_')[0],
                                                   'log' in row[f'{x}_description'])
                     for x in letter_columns]
    except Exception as e:
        print(f'problems initializing ' + row['entry_id'])
        print(e)
        variables = None
    return variables

