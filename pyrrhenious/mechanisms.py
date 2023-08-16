import numpy as np
from dataclasses import dataclass


class Mechanism:
    k = 8.617e-5

@dataclass
class SingleValue:
    mean : float
    stdev : float
    log : bool = False

    def get_value(self,sample=False,**kwargs):
        if sample:
            val = np.random.normal(loc=self.mean,scale=self.stdev)
        else:
            val = self.mean

        if self.log:
            return 10**val
        return val

class Val:

    def __init__(self,value):
        self.value=value

    def get_conductivity(self,**kwargs):
        return self.value.get_value(**kwargs)
class ArrheniousSimple(Mechanism):
    """

    sigma = a exp( h/kT)

    """

    def __init__(self, preexp, enthalpy, **kwargs):
        super().__init__()
        self.preexp   = preexp
        self.enthalpy = enthalpy

    def get_preexp(self, **kwargs):
        return self.preexp.get_value(**kwargs)

    def get_enthalpy(self, *args, **kwargs):
        return self.enthalpy.get_value(**kwargs)

    def get_conductivity(self, T=None, **kwargs):
        preexp = self.get_preexp(T=T, **kwargs)
        enthalpy = self.get_enthalpy(T=T, **kwargs)
        return preexp * np.exp(-enthalpy / (self.k * T))

class ArrheniousFugacity(ArrheniousSimple):
    """
    sigma = a exp(b/kT) fo2^c

    """
    def __init__(self,exponent,*args,**kwargs):
        super().__init__(*args, **kwargs)
        self.exponent = exponent

    def get_conductivity(self, logfo2=None, **kwargs):
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        conductivity1 = super(ArrheniousFugacity, self).get_conductivity(**kwargs)
        return conductivity1 * (10**logfo2) **self.exponent.get_value(**kwargs)
@dataclass
class ArrheniousfO2(Mechanism):
    """
    sigma = exp(a - b/kT + log(fo2^c))

    """
    preexp : SingleValue
    enthalpy : SingleValue
    const : SingleValue


    def get_conductivity(self, logfo2=None,T=None, **kwargs):
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        a = self.preexp.get_value(**kwargs)
        b = self.enthalpy.get_value(**kwargs)
        c = self.const.get_value(**kwargs)
        return np.exp(a - b/(self.k*T) + logfo2**c)

class ArrheniousfO22(Mechanism):
    """
    sigma = log10(a + b*log10 fO2 + c/T)

    """
    preexp: SingleValue
    const: SingleValue
    enthalpy: SingleValue

    def get_conductivity(self, logfo2=None, T=None, **kwargs):
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        a = self.preexp.get_value(**kwargs)
        b = self.enthalpy.get_value(**kwargs)
        c = self.const.get_value(**kwargs)
        return 10**(a - c / T + logfo2 ** b)
class ArrheniousPreexpParam(ArrheniousSimple):
    """

    sigma = A0(1 - B*P) exp( -h +P*V/kT)

    """

    def __init__(self, preexp1, preexp2, enthalpy, volume, **kwargs):
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
        assert P is not None, "Pressure value must be provided"
        h = self.enthalpy.get_value(**kwargs)
        v = self.volume.get_value(**kwargs)
        return h + v*P

class ArrheniousPressure(ArrheniousSimple):
    """

    sigma = a exp(-(h + PV)/kt)
    """

    def __init__(self,preexp,enthalpy,volume):
        super().__init__(preexp,enthalpy)
        self.volume = volume

    def get_enthalpy(self, *args,P=None, **kwargs):
        assert P is not None, "Pressure value must be provided"
        enthalpy = super(ArrheniousPressure, self).get_enthalpy(**kwargs)
        return enthalpy + self.volume.get_value(P=P,**kwargs)*P

class IronPreexpEnthalpyArrhenious:
    """

    sigma = a*Xfe * exp( -(b + c* Xfe^1/3 + P( d + e Xfe)/kT)

    """
    def __init__(self):
        pass
class IronWaterArrhenious1:
    """
    sigma = a (Cw)^b exp( -c/kt) exp(X_fe * d/kt)

    """
    def __init__(self):
        pass
class IronWaterArrhenious2:
    """
    sigma = a X_fe^b Cw^c exp( (-d- eX_fe- fCw^(1/3))/kT)

    """
    def __init__(self):
        pass

class NerstEinstein1:
    """

    sigma =  exp(b/kT) Ch q^2/kT

    """
    def __init__(self):
        pass


class NerstEinstein2:
    """

    C_h exp( a - b/(ln(10)RT))q^2/kT

    """
    def __init__(self):
        pass


class NerstEinstein3:
    """

    sigma =  Cw^(1+b) exp( c/kT)q^2/kT

    """
    def __init__(self):
        pass
class SEO3:

    def __init__(self):
        pass


class WaterExpArrhenious1:
    """

    sigma = a Cw^b exp( -c /kT)

    """
    def __init__(self):
        pass


class WaterExpArrhenious2:
    """

    sigma = a Cw exp( -(b + c C_w^(1/3) /kT)

    """

    def __init__(self):
        pass

class WaterExpArrheniousPressure:
    """

    sigma = a Cw^b exp(-(c + Pd)/kT)

    """

    def __init__(self):
        pass

class WaterExpArrheniousInvT:
    """

    sigma = a Cw^b exp(-(c + Pd)/kT)/T

    """

    def __init__(self):
        pass


class MechanismFactory:

    def __init__(self):
        self.models = {
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
'seo3_eq':SEO3,
'single_value':SingleValue,
'water_exponent_arrhenious':WaterExpArrhenious1,
'water_exponent_arrhenious_pressure':WaterExpArrheniousPressure,
'water_exponent_arrhenious_pressure_invT':WaterExpArrheniousInvT,
'water_preexp_enthalpy_arrhenious':WaterExpArrhenious2}

    def get_mechanism(self,str):
        return self.models[str]

