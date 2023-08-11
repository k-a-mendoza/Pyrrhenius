from .mineralbehavior import XStalWaterSimple, Mineral,XStalTypical
import numpy as np


def get_experiments():
    experiments = [Sun19Granulite_a()]
    exps = {}
    for exp in experiments:
        exps[exp.name]=exp
    return exps


class Granulite(Mineral):
    
    def __init__(self,mineral_name='Granulite',**kwargs):
        kwargs['Granulite']=kwargs
        super().__init__(**kwargs)
        
    def isotropic_conductivity(self,**kwargs):
        return self.mechanisms[0].get_conductivity(**kwargs)
    
    def polycrystalline_conductivity(self,**kwargs):
        return self.mechanisms[0].get_conductivity(**kwargs)
    
class Sun19Granulite_a(Granulite):
    
    def __init__(self):
        super().__init__(year=2019,author = "Wenquing Sun, Lidon Dai, Heping Li"+\
                         "Haiying Hu & Changcai Liu",
                         name='Sun19_2GPa_14.8\%_Bt_granulite',
                         title='Effect of temperature, pressure'+\
                        'and chemical composition on the electrical conductivity of granulite and geophysical implications',
                         min_p=2,max_p=1,min_t=573,max_t=1073,hydrous=False)
    
        self.mechanisms = [ XStalTypical(enthalpy=0.68,enthalpy_err=0.01,
                                preexp=2.19,preexp_err =0.01,log_sigma_err=True,
                                         convert2ev=False)]
        
class Sun19Granulite_b(Granulite):
    
    def __init__(self):
        super().__init__(year=2019,author = "Wenquing Sun, Lidon Dai, Heping Li"+\
                         "Haiying Hu & Changcai Liu",
                         name='Sun19_1GPa_5.49\%_Bt_granulite',
                         title='Effect of temperature, pressure'+\
                        'and chemical composition on the electrical conductivity of granulite and geophysical implications',
                         min_p=1,max_p=1,min_t=573,max_t=1073,hydrous=False)
    
        self.mechanisms = [ XStalTypical(0.61,enthalpy_err=0.01,
                                preexp=1.01,preexp_err =0.01,log_sigma_err=True,
                                         convert2ev=False)+
                            XStalTypical(0.61,enthalpy_err=0.01,
                                preexp=1.01,preexp_err =0.01,log_sigma_err=True,
                                         convert2ev=False)]
        
class Sun19Granulite_c(Granulite):
    
    def __init__(self):
        super().__init__(year=2019,author = "Wenquing Sun, Lidon Dai, Heping Li"+\
                         "Haiying Hu & Changcai Liu",
                         name='Sun19_2GPa_8.75\%_Bt_granulite',
                         title='Effect of temperature, pressure'+\
                        'and chemical composition on the electrical conductivity of granulite and geophysical implications',
                         min_p=2,max_p=2,min_t=573,max_t=1073,hydrous=False)
    
        self.mechanisms = [ XStalTypical(0.61,enthalpy_err=0.01,
                                preexp=1.01,preexp_err =0.01,log_sigma_err=True,
                                         convert2ev=False)]
        
    