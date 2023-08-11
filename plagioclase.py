from .mineralbehavior import XStalWaterSimple, Mineral,XStalTypical
import numpy as np


def get_experiments():
    experiments = [XiaoY12Plag()]
    exps = {}
    for exp in experiments:
        exps[exp.name]=exp
    return exps
class Plagioclase(Mineral):
    
    def __init__(self,mineral_name='Plagioclase',**kwargs):
        kwargs['Plagioclase']=kwargs
        super().__init__(**kwargs)
        
class XiaoY12Plag(Plagioclase):
    
    def __init__(self):
        super().__init__(year=2012,author = "Xiaozhi Yang, Hans Keppler,"+\
                         "Catherine McCammon & Huaiwei Ni ",
                         name='XiaoY12_plag',title='Electrical conductivity of'+\
                         'orthopyroxene and plagioclase in the lower crust',
                         min_p=0.6,max_p=1.2,min_t=573,max_t=1273,hydrous=True)
    
                         
        self.mechanisms = [ XStalTypical(enthalpy=161,enthalpy_err=6,
                                preexp=4.12,preexp_err =0.34,log_sigma_err=True,
                                         convert2ev=True),
                           XStalWaterSimple(enthalpy=77,enthalpy_err=2,
                                preexp=-0.8,preexp_err =0.14,log_sigma_err=True,
                                        preexp_c_constant=0.83,
                                        preexp_c_constant_err=0.06,
                                         convert2ev=True)
                           ]
        
    def isotropic_conductivity(self,temperature=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.mechanisms:
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
    
    def max_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(x_fe=self.static_iron_frac,**kwargs)
        
    def min_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(x_fe=self.static_iron_frac,**kwargs)
    
    def polycrystalline_conductivity(self,**kwargs):
        return self.isotropic_conductivity(x_fe=self.static_iron_frac,**kwargs)