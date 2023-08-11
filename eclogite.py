from .mineralbehavior import XStalWaterIronEmpirical, Mineral
import numpy as np
def get_experiments():
    experiments = [Zhang16OOmph()]
    exps = {}
    for exp in experiments:
        exps[exp.name]=exp
    return exps

class Eclogite(Mineral):
    def __init__(self,mineral_name='eclogite',**kwargs):
        kwargs['mineral_name']=mineral_name
        super().__init__(**kwargs)

class Zhang16OOmph(Eclogite):
    
    def __init__(self):
        super().__init__(year=2019, author = "Baohua Zhang, Chengcheng Zhao, Jianhua Ge, Takashi Yoshino",\
                         name='Zhang19_oomph',title='Electrical Conductivity of Omphacite'+\
                         ' as a Function of Water Content and Implications for High Conductivity '+\
                         'Anomalies in the Dabie-Sulu UHPM Belts and Tibet',
                         min_p=1.5,max_p=3,min_t=500,max_t=1300,hydrous=True,iron_frac=True)
        self.mechanisms = [ XStalWaterIronEmpirical(enthalpy=1.05,enthalpy_err=0.03,beta=-0.38, alpha_err=0.05,
                                                    alpha=-0.71,beta_err=0.07,preexp_c_constant=1.1,preexp_c_constant_err=0.1,
                                                    iron_exp=0.09,iron_exp_err=0.03,water_content_type='wt%',
                                                    enthalpy_exp=1/3,enthalpy_exp_err=0,
                                preexp=1483,preexp_err = 260,log_sigma_err=False)]
        
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