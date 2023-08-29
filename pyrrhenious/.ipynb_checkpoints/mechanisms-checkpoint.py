import numpy as np
from dataclasses import dataclass
import inspect
def _count_positional_args_required(func):
    signature = inspect.signature(func)
    empty = inspect.Parameter.empty
    total = -1
    for param in signature.parameters.values():
        if param.default is empty:
            total += 1
    return total
class Mechanism:
    """Base Abstract class for all electrical conductivity mechanisms"""

    k = 8.617333262e-5 # in ev/k
    q = 1.602176634e-19 # charge of electron
    q2 =(1.602176634e-19)**2
    r = 8.31446261815324e-3 # gas constant in kj/mol
    cm3GPamol_to_ev = 1/96.49 # converts pressure*cm^3/mol volume to ev

    def __init__(self):
        self.repr_properties = []

    def set_water_units(self, units):
        """
        Set the units for water

        Parameters
        ----------
        units : str
            The units the model uses for water. Should be either wtpct or ppm
        """

        self._water_units = units

    def convert_water(self, cw):
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

        assert cw is not None, "Water value must be provided"

        if self._water_units == 'wtpct':
            return cw * 1e-4
        elif self._water_units == 'wtpct10x':
            return cw * 1e-5
        elif self._water_units == 'ppm':
            return cw
        else:
            raise NotImplementedError(
                f"{self._water_units} conversion not implemented"
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

        assert P is not None, "Pressure value must be provided"
        return P * self.cm3GPamol_to_ev

    def get_conductivity(self, T=None, **kwargs):
        """
        Abstract method to calculate the conductivity

        Parameters
        ----------
        T : float, optional
            The temperature value (default is None)

        Returns
        -------
        np.ndarray or float
            The calculated conductivity
        """
        pass

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
    
    def __repr__(self):
        repr_string = f"{self.__class__.__name__} "
        for k, v in vars(self).items():

            if isinstance(v, StochasticConstant):
                repr_string += str(v) + " "

        return repr_string

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
            val = np.random.normal(loc=self.mean, scale=self.stdev)
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

    def get_conductivity(self,**kwargs):
        """
        Get the conductivity.

        Returns
        -------
        float
            The conductivity in units of s/m

        """
        return self.value.get_value(**kwargs)
    
class ArrheniousSimple(Mechanism):
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

    def get_conductivity(self, T=None, **kwargs):
        """
        Get the conductivity.

        Parameters
        ----------
        T : float, optional
            The temperature value in units of kelvin.

        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The conductivity in units of s/m

        """
        preexp = self.get_preexp(T=T, **kwargs)
        enthalpy = self.get_enthalpy(T=T, **kwargs)
        return preexp * np.exp(-enthalpy / (self.k * T))

    def __repr__(self):
        return f'{self.preexp} exp( -{self.enthalpy}/kT)'


class LinkedArrhenious(ArrheniousSimple):
    """
    Abstract Variation of the ArrheniousSimple class which is primarily used for linear volatile-dependent activation and preexponential factors. 
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

    def _get_enthalpy(self, volatile, **kwargs):
        """
        Calculates enthalpy.

        Parameters
        ----------
        volatile : float
            The volatile concentration
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
        return a * np.exp(b * volatile) + c

    def get_preexp(self, **kwargs):
        """
        calculates the preexponential value

        Parameters
        ----------
        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The value of the pre-exponential factor.

        """
        ea = self._get_enthalpy(**kwargs)
        d = self.const4.get_value(**kwargs)
        e = self.const5.get_value(**kwargs)
        return np.exp(ea * d + e)

    def get_conductivity(self, T=None, **kwargs):
        """
        Get the conductivity.

        Parameters
        ----------
        T : float, optional
            The temperature value in Kelvin
        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The conductivity in s/m.

        """
        preexp = self.get_preexp(**kwargs)
        h = self._get_enthalpy(**kwargs) * 1e-3
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
        return _count_positional_args_required(LinkedArrhenious.__init__)
        

class LinkedArrheniousWet(LinkedArrhenious):
    """
    Variation of the LinkedArrhenious class that is specifically used for water-dependent conductivity

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

    def get_enthalpy(self, Cw=None, **kwargs):
        """
        Get the value of the enthalpy considering the presence of water.

        Parameters
        ----------
        Cw : float, optional
            Concentration of water. Defaults to None.
        kwargs : dict
            Optional keyword arguments.

        Returns
        -------
        float
            The value of the enthalpy.

        """
        cw = self.convert_water(Cw)
        return self._get_enthalpy(cw)



class LinkedArrheniousCO2(LinkedArrhenious):
    """
    Variation of the LinkedArrhenious class that is specifically used for CO2-dependent conductivity

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
        super().__init__(*args)

    @property
    def uses_co2(self):
        return True

    def get_enthalpy(self, co2=None, **kwargs):
        co2 = self.convert_co2(co2)
        return self._get_enthalpy(co2)
    

class VogelFulcherTammanWet(ArrheniousSimple):
    """
    A class that represents a Vogel-Fulcher-Tamman equation for wet substances.

    sigma = a exp( (b - c Cw^d)/k(T-e))

    Inherits from the ArrheniousSimple class.

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
        h = ArrheniousSimple.get_enthalpy(self,T=None,Cw=Cw,**kwargs)
        w_effect = self.const1.get_value(**kwargs)
        exponent = self.const2.get_value(**kwargs)
        return h + w_effect* self.convert_water(Cw)**exponent

    def get_conductivity(self, T=None, **kwargs):
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
        return preexp * np.exp(-enthalpy / (self.k * (T- t0)))

class ArrheniousFugacity(ArrheniousSimple):
    """ 
    A class that represents a oxygen-fugacity dependent arrhenious preexponential constant

    sigma = a exp(b/kT) fo2^c

    Inherits from the ArrheniousSimple class.

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
        super().__init__(preexp,enthalpy)
        self.exponent = exponent

    def get_conductivity(self, logfo2=None, **kwargs):
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
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        conductivity1 = ArrheniousSimple.get_conductivity(self,**kwargs)
        return conductivity1 * (10**logfo2) ** self.exponent.get_value(**kwargs)

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
class ArrheniousfO2(Mechanism):
    """

    sigma = exp(log10(a) - b/kT + c*log(fo2))

    """
    preexp : StochasticConstant
    enthalpy : StochasticConstant
    const : StochasticConstant

    def get_conductivity(self, logfo2=None,T=None, **kwargs):
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        a = self.preexp.get_value(**kwargs)
        b = self.enthalpy.get_value(**kwargs)
        c = self.const.get_value(**kwargs)
        return 10**(a - b/T + c*logfo2)

    def uses_fo2(self):
        return True

@dataclass
class ArrheniousfO22(Mechanism):
    """
    sigma = log10(a + b/T + c*log10 fO2)

    """
    preexp: StochasticConstant
    const: StochasticConstant
    enthalpy: StochasticConstant
    def get_conductivity(self, logfo2=None, T=None, **kwargs):
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        a = self.preexp.get_value(**kwargs)
        b = self.enthalpy.get_value(**kwargs)
        c = self.const.get_value(**kwargs)
        return 10**(a + b / T + c*logfo2)

    def uses_fo2(self):
        return True

class ArrheniousPreexpParam(ArrheniousSimple):
    """

    sigma = A0(1 - B*P) exp( -h +P*V/kT)

    """
    def __init__(self, preexp1 : StochasticConstant, preexp2 : StochasticConstant,
                 enthalpy : StochasticConstant, volume : StochasticConstant):
        super().__init__(None,None)
        self.preexp1 = preexp1
        self.preexp2 = preexp2
        self.enthalpy = enthalpy
        self.volume = volume

    def get_preexp(self,P=None, **kwargs):
        assert P is not None, "Pressure value must be provided"
        a0 = self.preexp1.get_value(P=P,**kwargs)
        b =  self.preexp2.get_value(P=P,**kwargs)
        return a0*(1 - b*P)

    def get_enthalpy(self, P=None,**kwargs):
        h = self.enthalpy.get_value(**kwargs)
        v = self.volume.get_value(**kwargs)
        return h + v*self.convert_pressure(P)

    @property
    def uses_pressure(self):
        return True


class ArrheniousPressure(ArrheniousSimple):
    """

    sigma = a exp(-(h + PV)/kt)
    """
    n_constants = 3
    def __init__(self,preexp : StochasticConstant,enthalpy : StochasticConstant,volume):
        super().__init__(preexp,enthalpy)
        self.volume = volume

    def get_enthalpy(self, *args,P=None, **kwargs):
        enthalpy = super(ArrheniousPressure, self).get_enthalpy(**kwargs)
        return enthalpy + self.volume.get_value(**kwargs)*self.convert_pressure(P)

    @property
    def uses_pressure(self):
        return True

class IronPreexpEnthalpyArrhenious(ArrheniousSimple):
    """

    sigma = a*Xfe * exp( -(b + c* Xfe^1/3 + P( d + e Xfe)/kT)

    """
    n_constants = 5

    def __init__(self, preexp: StochasticConstant, enthalpy: StochasticConstant,
                 const1: StochasticConstant, volume: StochasticConstant, const2: StochasticConstant):
        super().__init__(None, None)
        self.preexp = preexp
        self.enthalpy = enthalpy
        self.volume = volume
        self.const1 = const1
        self.const2 = const2

    def get_preexp(self, X_fe=None, **kwargs):
        assert X_fe is not None, "Iron value must be provided"
        a = self.preexp.get_value(**kwargs)
        return a*X_fe

    def get_enthalpy(self, X_fe=None, P=None, **kwargs):
        assert X_fe is not None, "Iron value must be provided"
        h = self.enthalpy.get_value(**kwargs)
        v = self.volume.get_value(**kwargs)
        a = self.const1.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        return h + a*X_fe**(1/3) + self.convert_pressure(P)*(v + b*X_fe)

    @property
    def uses_pressure(self):
        return True

    @property
    def uses_iron(self):
        return True

class IronWaterArrhenious1(ArrheniousSimple):
    """
    sigma = a (Cw)^b exp( -c/kt) exp(X_fe * d/kt)

    """

    def __init__(self, preexp: StochasticConstant, const: StochasticConstant,
                 enthalpy1: StochasticConstant, enthalpy2: StochasticConstant):
        super().__init__(None, enthalpy1)
        self.preexp = preexp
        self.enthalpy2 = enthalpy2
        self.const = const


    def get_preexp(self, X_fe=None, Cw=None,T=None, **kwargs):
        assert X_fe is not None, "Iron value must be provided"
        a = self.preexp.get_value(**kwargs)
        b = self.const.get_value(**kwargs)
        c = self.enthalpy2.get_value(**kwargs)
        upper_enth = X_fe * c/(self.r * T)
        return a * ((self.convert_water(Cw))**b) * np.exp(upper_enth)

    @property
    def uses_iron(self):
        return True

    @property
    def uses_water(self):
        return True

class IronWaterArrhenious2(ArrheniousSimple):
    """
    sigma = a (X_fe+b) ^c Cw^d exp( (-e+ f(X_fe+b)+ gCw^(h))/kT)

    """
    n_constants = 8

    def __init__(self, preexp: StochasticConstant, const1: StochasticConstant,
                 const2: StochasticConstant, const3: StochasticConstant,
                 enthalpy1: StochasticConstant,const4 : StochasticConstant,
                 const5 : StochasticConstant,const6 : StochasticConstant):
        super().__init__(None, None)
        self.preexp = preexp
        self.enthalpy = enthalpy1
        self.const1 = const1
        self.const2 = const2
        self.const3 = const3
        self.const4 = const4
        self.const5 = const5
        self.const6 = const6

    def get_preexp(self, X_fe=None,Cw=None, **kwargs):
        assert X_fe is not None, "Iron value must be provided"
        a = self.preexp.get_value(**kwargs)
        b = self.const1.get_value(**kwargs)
        c = self.const2.get_value(**kwargs)
        d = self.const3.get_value(**kwargs)
        return a * (X_fe+b) **c * self.convert_water(Cw)**d

    def get_enthalpy(self, Cw=None, X_fe=None, **kwargs):
        assert X_fe is not None, "Iron value must be provided"
        b = self.const1.get_value(**kwargs)
        h = self.enthalpy.get_value(**kwargs)
        alpha = self.const4.get_value(**kwargs)
        beta = self.const5.get_value(**kwargs)
        exp = self.const6.get_value(**kwargs)
        return h + alpha * (X_fe+b)  + beta * self.convert_water(Cw)**exp

    @property
    def uses_iron(self):
        return True

    @property
    def uses_water(self):
        return True

class NerstEinstein1(Mechanism):
    """

    sigma =  a exp(b/kT) Ch q^2/kT

    """
    n_constants = 2
    def __init__(self,const1 : StochasticConstant, enthalpy : StochasticConstant):
        super().__init__()
        self.const1 = const1
        self.enthalpy = enthalpy
    def get_conductivity(self, T=None,Cw=None, **kwargs):
        a = self.const1.get_value(**kwargs)
        h = self.enthalpy.get_value(**kwargs)

        return self.convert_water(Cw) * a * np.exp(-h/(self.k*T))/(self.k*T)

    @property
    def uses_water(self):
        return True


class NerstEinstein2(Mechanism):
    """

    a Cw b q^2 exp(-h/kT)/kT

    """
    n_constants = 3

    def __init__(self, const1: StochasticConstant, const2: StochasticConstant, enthalpy : StochasticConstant):
        super().__init__()
        self.const1 = const1
        self.const2 = const2
        self.enthalpy = enthalpy

    def get_conductivity(self, T=None, Cw=None, **kwargs):
        a = self.const1.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        h = self.enthalpy.get_value(**kwargs)
        a0 = self.q2 / (self.k* T)

        enthalpy_term = np.exp(-h/(self.k*T))
        return  a * self.convert_water(Cw) * b *  enthalpy_term * a0 


    @property
    def uses_water(self):
        return True


class NerstEinstein3(Mechanism):
    """

    sigma =  a Cw^(1+b) exp( c/kT)q^2/kT

    """
    n_constants = 3

    def __init__(self, const1: StochasticConstant, const2: StochasticConstant, enthalpy: StochasticConstant):
        super().__init__()
        self.const1 = const1
        self.const2 = const2
        self.enthalpy = enthalpy

    def get_conductivity(self, T=None, Cw=None, **kwargs):
        a = self.const1.get_value(**kwargs)
        b = self.const1.get_value(**kwargs)
        h = self.enthalpy.get_value(**kwargs)
        a0 = self.q2 / self.k

        return a * self.convert_water(Cw) **(1+b) * np.exp(-h/(self.k*T)) * a0 / T


    @property
    def uses_water(self):
        return True

    def __repr__(self):
        return f'q^2 {self.const1} Cw^({self.const2}+1) exp(-{self.enthalpy}/kT)/kT'

@dataclass
class SEO3Term(Mechanism):
    """
    Uses a combination of ArrheniousSimple and ArrheniousFugacity models to represent a polaron hopping term for the SEO3 model of olivine
    sigma = e [a exp(b/kT) + c exp(d/kT)fO2∧e] f exp(g/kT)

    Inherits from the Mechanism class.

    Parameters:
    -----------
    preexp1 : StochasticConstant
        The pre-exponential factor for the first arrhenious law
    enth1 : StochasticConstant
        The enthalpy factor for the first arrhenious law
    preexp2 : StochasticConstant
        The pre-exponential factor for the second arrhenious law
    enth2 : StochasticConstant
        The enthalpy for the second arrhenious law
    exp : StochasticConstant
        The oxygen fugacity exponent for the second arrhenious law

    preexp3 : StochasticConstant
        The pre-exponential factor for the third arrhenious law
    enth3 : StochasticConstant
        The enthalpy for the third arrhenious law
    


    Methods:
    --------
    uses_fo2()
        Returns True, indicating that oxygen fugacity is used in the equation.

    get_conductivity(T=None, **kwargs)
        Calculates and returns the conductivity of the wet substance.

    """

    def __init__(self,preexp1 : StochasticConstant, enth1 : StochasticConstant, preexp2 : StochasticConstant,
                 enth2 : StochasticConstant, exp : StochasticConstant, preexp3 : StochasticConstant, enth3 : StochasticConstant):
        super().__init__()
        self.b_species    = ArrheniousSimple(preexp1,enth1)
        self.b_species1   = ArrheniousFugacity(preexp2,exp, enth2)
        self.mu_potential = ArrheniousSimple(preexp3,enth3)

    def get_conductivity(self, **kwargs):
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
        return ArrheniousSimple.n_args()*2 + ArrheniousFugacity.n_args()

    def __repr__(self):
        return f'e*({self.b_species} {self.b_species1})*{self.mu_potential}'

    @property
    def uses_fo2(self):
        return True


class WaterExpArrhenious1(ArrheniousSimple):
    """

    sigma = a Cw^b exp( -c /kT)

    """
    n_constants = 3
    def __init__(self,preexp: StochasticConstant, const1 : StochasticConstant, enthalpy : StochasticConstant):
        super().__init__(preexp,enthalpy)
        self.const1 = const1

    def get_preexp(self,Cw=None, **kwargs):
        a = self.preexp.get_value(**kwargs)
        p = self.const1.get_value(**kwargs)
        return a * self.convert_water(Cw)**p

    @property
    def uses_water(self):
        return True


class WaterExpArrhenious2(ArrheniousSimple):
    """

    sigma = a Cw exp( -(b + c C_w^(d) /kT)

    """
    n_constants = 4

    def __init__(self, preexp: StochasticConstant, enthalpy: StochasticConstant,const1: StochasticConstant,
                 const2: StochasticConstant):
        super().__init__(preexp, enthalpy)
        self.const1 = const1
        self.const2 = const2

    def get_preexp(self, Cw=None, **kwargs):
        a = self.preexp.get_value(**kwargs)
        return a * self.convert_water(Cw)

    def get_enthalpy(self, Cw=None, **kwargs):
        a = self.enthalpy.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        c = self.const1.get_value(**kwargs)
        return a + c*self.convert_water(Cw)**(b)

    @property
    def uses_water(self):
        return True

class WaterExpArrheniousPressure(ArrheniousPressure):
    """

    sigma = a Cw^b exp(-(c + Pd)/kT)

    """
    n_constants = 4

    def __init__(self, preexp: StochasticConstant, const1: StochasticConstant, enthalpy: StochasticConstant,
                 volume : StochasticConstant):
        super().__init__(preexp, enthalpy, volume)
        self.const1 = const1

    def get_preexp(self, Cw=None, **kwargs):
        a = self.preexp.get_value(**kwargs)
        b = self.const1.get_value(**kwargs)
        return a * self.convert_water(Cw)**b

    @property
    def uses_water(self):
        return True

    @property
    def uses_pressure(self):
        return True

class WaterExpArrheniousInvT(WaterExpArrheniousPressure):
    """

    sigma = a Cw^b exp(-(c + Pd)/kT)/T

    """
    n_constants = 4

    def __init__(self, preexp: StochasticConstant, const1: StochasticConstant, enthalpy: StochasticConstant,
                 volume: StochasticConstant):
        super().__init__(preexp,const1, enthalpy,volume)

    def get_conductivity(self, T=None, **kwargs):
        c = super(WaterExpArrheniousPressure).get_conductivity(T=T, **kwargs)
        return c/T

    @property
    def uses_water(self):
        return True

    @property
    def uses_pressure(self):
        return True


model_dict = {
            'arrhenious_fugacity':ArrheniousFugacity,
'arrhenious_log_fO2':ArrheniousfO2,
'arrhenious_logfO2_2':ArrheniousfO22,
'arrhenious_preexp_parameterized':ArrheniousPreexpParam,
'arrhenious_pressure':ArrheniousPressure,
'arrhenious_simple':ArrheniousSimple,
'iron_preexp_enthalpy_arrhenious':IronPreexpEnthalpyArrhenious,
'iron_water_arrhenious1':IronWaterArrhenious1,
'iron_water_arrhenious2': IronWaterArrhenious2,
'nerst_einstein1':NerstEinstein1,
'nerst_einstein_2':NerstEinstein2,
'nerst_einstien_3':NerstEinstein3,
'seo3_term':SEO3Term,
'single_value':SingleValue,
'water_exponent_arrhenious':WaterExpArrhenious1,
'water_exponent_arrhenious_pressure':WaterExpArrheniousPressure,
'water_exponent_arrhenious_pressure_invT':WaterExpArrheniousInvT,
'water_preexp_enthalpy_arrhenious':WaterExpArrhenious2,
'vft_wet': VogelFulcherTammanWet,
'linked_arrhenious_wet':LinkedArrheniousWet,
'linked_arrhenious_co2':LinkedArrheniousCO2}

