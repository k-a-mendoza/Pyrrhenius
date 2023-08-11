import numpy as np

class SifreMeltBase:
    r = 8.3145
    def __init__(self,a,b,c,d,e):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        
    def isotropic_conductivity(self,temperature=None,**kwargs):
        volatile_frac = self.correct_volatile_frac(**kwargs)
        enthalpy = self.a*np.exp(-self.b*volatile_frac)+self.c
        sigma = np.exp(self.d*enthalpy +self.e)
        
        return sigma*np.exp(-enthalpy/(self.r*temperature))
                            
    def correct_volatile_frac(self,**kwargs):
        return 1

class SilicateMelt(SifreMeltBase):

    def __init__(self):
        super().__init__(88_774,0.388,73_029,4.54e-5,5.5607)
                            
    def correct_volatile_frac(self,wtppm_H2O=0,**kwargs):
        return wtppm_H2O/1e4
    
class CarbonateMelt(SifreMeltBase):

    def __init__(self):
        super().__init__(789_166,0.1808,32_820,5.5e-5, 5.7956)
                            
    def correct_volatile_frac(self,CO2_wt=0,**kwargs):
        return CO2_wt
        
 
class SifreMelts:
    
    def __init__(self):
        self.carbonated = CarbonateMelt()
        self.hydrated = SilicateMelt()
        self.dperid = 3.3
        self.dh2O = 1.4
        self.dcarb = 2.4
        self.dbasalt = 2.8

    def isotropic_conductivity(self,**kwargs):
        carbonated =self.carbonated.isotropic_conductivity(**kwargs)
        hydrated =self.hydrated.isotropic_conductivity(**kwargs)
        return carbonated + hydrated
    
    def water_co2_wt_melt_2frac(self, wt_melt, wtppm_H2O=0, CO2_wt=0, **kwargs):
        wtH2O = self.hydrated.correct_volatile_frac(wtppm_H2O=wtppm_H2O)
        
        dmelt = wtH2O*self.dh2O/100 + self.dcarb*2*CO2_wt/100 + (1 - (wtH2O+2*CO2_wt)/100)*self.dbasalt
        x_melt = 1/(1+((100/wt_melt) -1)*dmelt/self.dperid)
        return x_melt
        