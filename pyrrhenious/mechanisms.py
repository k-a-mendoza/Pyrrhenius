import numpy as np
from dataclasses import dataclass


class Mechanism:
    k = 8.617e-5
    q2 =2*1.602176634e-19
    r = 8.314462
    cm3GPamol_to_ev = 1/96.49

    def set_water_units(self,units):
        self._water_units = units

    def convert_water(self,cw):
        assert cw is not None, "Water value must be provided"
        if self._water_units=='wtpct':
            return cw*1e-4
        elif self._water_units=='wtpct10x':
            return cw*1e-5
        elif self._water_units=='ppm':
            return cw
        else:
            raise NotImplementedError(f"{self._water_units} conversion not implemented")

    def convert_pressure(self,P):
        assert P is not None, "Pressure value must be provided"
        return P*self.cm3GPamol_to_ev
    def get_conductivity(self, T=None, **kwargs):# -> np.ndarray | float if python 3.10>
        pass

    @property
    def uses_water(self):
        return False

    @property
    def uses_iron(self):
        return False

    @property
    def uses_pressure(self):
        return False

    @property
    def uses_fo2(self):
        return False

@dataclass
class StochasticConstant:
    mean : float
    stdev : float
    type : str
    log: bool = False
    isnan: bool = False
    def get_value(self,sample=False,**kwargs):
        if sample:
            val = np.random.normal(loc=self.mean,scale=self.stdev)
        else:
            val = self.mean

        if self.log:
            return 10**val
        return val

    def __repr__(self):
        return f'{self.type}({self.mean}+-{self.stdev}: islog: {self.log})'
class SingleValue(Mechanism):
    n_constants = 1
    def __init__(self,value):
        super().__init__()
        self.value=value

    def get_conductivity(self,**kwargs):
        return self.value.get_value(**kwargs)
class ArrheniousSimple(Mechanism):
    """

    sigma = a exp( h/kT)

    """
    n_constants = 2

    def __init__(self, preexp : StochasticConstant, enthalpy : StochasticConstant, **kwargs):
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
    n_constants = 3
    def __init__(self,exponent : StochasticConstant,*args,**kwargs):
        super().__init__(*args, **kwargs)
        self.exponent = exponent

    def get_conductivity(self, logfo2=None, **kwargs):
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        conductivity1 = super(ArrheniousFugacity, self).get_conductivity(**kwargs)
        return conductivity1 * (10**logfo2) ** self.exponent.get_value(**kwargs)

    @property
    def uses_fo2(self):
        return True

@dataclass
class ArrheniousfO2(Mechanism):
    """

    sigma = exp(log10(a) - b/kT + c*log(fo2))

    """
    preexp : StochasticConstant
    enthalpy : StochasticConstant
    const : StochasticConstant

    n_constants = 3
    def get_conductivity(self, logfo2=None,T=None, **kwargs):
        assert logfo2 is not None, 'Did not provide an oxygen fugacity value!'
        a = self.preexp.get_value(**kwargs)
        b = self.enthalpy.get_value(**kwargs)
        c = self.const.get_value(**kwargs)
        print(f'{a} {b} {c}')
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
    n_constants = 3
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
    n_constants = 4
    def __init__(self, preexp1 : StochasticConstant, preexp2 : StochasticConstant,
                 enthalpy : StochasticConstant, volume : StochasticConstant, **kwargs):
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
    n_constants = 4

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
        print(self.enthalpy.get_value())
        upper_enth = X_fe * c/(self.r * T)
        return a * ((self.convert_water(Cw))**b) * np.exp(upper_enth.astype(float))

    @property
    def uses_iron(self):
        return True

    @property
    def uses_water(self):
        return True

class IronWaterArrhenious2(ArrheniousSimple):
    """
    sigma = a X_fe^b Cw^c exp( (-d- eX_fe- fCw^(1/3))/kT)

    """
    n_constants = 6

    def __init__(self, preexp: StochasticConstant, const1: StochasticConstant,
                 const2: StochasticConstant, enthalpy1: StochasticConstant,
                 const3: StochasticConstant,const4 : StochasticConstant):
        super().__init__(None, None)
        self.preexp = preexp
        self.enthalpy = enthalpy1
        self.const1 = const1
        self.const2 = const2
        self.const3 = const4
        self.const4 = const3

    def get_preexp(self, X_fe=None,Cw=None, **kwargs):
        assert X_fe is not None, "Iron value must be provided"
        a = self.preexp.get_value(**kwargs)
        b = self.const1.get_value(**kwargs)
        c = self.const2.get_value(**kwargs)
        return a * X_fe**b * self.convert_water(Cw)**c

    def get_enthalpy(self, Cw=None, X_fe=None, **kwargs):
        assert X_fe is not None, "Iron value must be provided"
        h = self.enthalpy.get_value(**kwargs)
        a = self.const3.get_value(**kwargs)
        b = self.const4.get_value(**kwargs)

        return h + a * X_fe  + b * self.convert_water(Cw)**(1/3)

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
        a0 = self.q2/self.k

        return self.convert_water(Cw) * a * np.exp(h/(self.k*T))*a0/T

    @property
    def uses_water(self):
        return True


class NerstEinstein2(Mechanism):
    """

    a Cw exp( b - c/(ln(10)RT))q^2/kT

    """
    n_constants = 3

    def __init__(self, const1: StochasticConstant, enthalpy: StochasticConstant, const2 : StochasticConstant):
        super().__init__()
        self.const1 = const1
        self.const2 = const2
        self.enthalpy = enthalpy

    def get_conductivity(self, T=None, Cw=None, **kwargs):
        a = self.const1.get_value(**kwargs)
        c = self.const1.get_value(**kwargs)
        h = self.enthalpy.get_value(**kwargs)
        a0 = self.q2 / self.k

        enthalpy_term = h - c/(self.r*np.log(10)*T)
        return  a * self.convert_water(Cw) * np.exp(enthalpy_term) * a0 / T


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

        return a * self.convert_water(Cw) **(1+b) * np.exp(h/(self.k*T)) * a0 / T


    @property
    def uses_water(self):
        return True

@dataclass
class SEO3(Mechanism):
    """

    [a∗exp(b/kT) + c∗exp(d/kT)fO2∧e]∗f∗exp(g/kT)∗h + 2∗[i∗exp(j/kT) +k∗exp(L/kT)fO2∧m]∗(n∗exp(o /kT)) p

    """
    const1: StochasticConstant
    const2: StochasticConstant
    const3: StochasticConstant
    const4: StochasticConstant
    const5: StochasticConstant
    const6: StochasticConstant
    const7: StochasticConstant
    const8: StochasticConstant
    const9: StochasticConstant
    const10: StochasticConstant
    const11: StochasticConstant
    const12: StochasticConstant
    const13: StochasticConstant
    const14: StochasticConstant
    const15: StochasticConstant
    const16: StochasticConstant
    n_constants =16

    def get_conductivity(self, T=None,logfo2=None, **kwargs):
        assert logfo2 is not None, "oxygen fugacity variable {logfo2} value must be provided"
        kT = T*self.k
        fo2 = 10**logfo2
        a = self.const1.get_value(**kwargs)
        b = self.const2.get_value(**kwargs)
        c = self.const3.get_value(**kwargs)
        d = self.const4.get_value(**kwargs)
        e = self.const5.get_value(**kwargs)

        f = self.const6.get_value(**kwargs)
        g = self.const7.get_value(**kwargs)
        h = self.const8.get_value(**kwargs)
        i = self.const9.get_value(**kwargs)
        j = self.const10.get_value(**kwargs)

        k = self.const11.get_value(**kwargs)
        l = self.const12.get_value(**kwargs)
        m = self.const13.get_value(**kwargs)
        n = self.const14.get_value(**kwargs)
        o = self.const15.get_value(**kwargs)
        p = self.const16.get_value(**kwargs)
        term1 = a*np.exp(b/kT) + c*np.exp(d/kT)*fo2**e
        term2 = f*np.exp(g/kT)*h
        term3 = i*np.exp(j/kT) +k*np.exp(l/kT)*fo2**m
        term4 = 2*(n*np.exp(o /kT))*p

        return term1*term2 + term3*term4

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

    sigma = a Cw exp( -(b + c C_w^(1/3) /kT)

    """
    n_constants = 3

    def __init__(self, preexp: StochasticConstant, enthalpy: StochasticConstant,const1: StochasticConstant):
        super().__init__(preexp, enthalpy)
        self.const1 = const1

    def get_preexp(self, Cw=None, **kwargs):
        a = self.preexp.get_value(**kwargs)
        return a * self.convert_water(Cw)

    def get_enthalpy(self, Cw=None, **kwargs):
        a = self.enthalpy.get_value(**kwargs)
        c = self.const1.get_value(**kwargs)
        return a + c*self.convert_water(Cw)**(1/3)

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
'seo3_eq':SEO3,
'single_value':SingleValue,
'water_exponent_arrhenious':WaterExpArrhenious1,
'water_exponent_arrhenious_pressure':WaterExpArrheniousPressure,
'water_exponent_arrhenious_pressure_invT':WaterExpArrheniousInvT,
'water_preexp_enthalpy_arrhenious':WaterExpArrhenious2}

