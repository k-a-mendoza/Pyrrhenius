import numpy as np
from .mineralbehavior import Mineral, XStalTypical, XStalWaterEmpirical, XStalWaterSimple,\
XStalfO2Premult, XStalWaterInvT, calc_QFM, XStalWaterSimpleIronExp,OldEmpiricalLogEqn,\
XStalXfeComplex,XStalAvg

def get_experiments():
    experiments = [XuOrtho(),XuClino(),Wang08(),Huebner88a(),Huebner88b(),Zhang16Opx(),
                  DaiKarato09(),XiaoY12(),XiaoY12CPX(),Zhang12opx()]
    exps = {}
    for exp in experiments:
        exps[exp.name]=exp
    return exps

class DryPyroxenite(Mineral):
    k = 8.617e-5
    r = 8.31446e-3 #kj/mol k
    h2e = 0.01036427
    norm_v = 1/96.49 # ev mol/Gpa cm^3

    def __init__(self,mineral_name='Pyroxenite',**kwargs):
        kwargs['mineral_name']=mineral_name
        super().__init__(**kwargs)
     
    def conductivity_001(self,temperature=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.crystal_001:
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
        
    def conductivity_010(self,temperature=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.crystal_010:
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
        
    def conductivity_100(self,temperature=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.crystal_100:
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
        
    def polycrystalline_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
    
class Wang08(DryPyroxenite):
    """
    a pyroxenite, fit to a range of mantle rocks.
    static iron fraction of x_fe = 0.545
    
    """
    def __init__(self):
        super().__init__(year=2008,author="Duojun Wang, Heping Li, Li Yi, Baoping Shi",
                         name='wang08_pyx',title ="The electrical conductivity of upper-mantle"+\
                        "rocks: water content in the upper mantle",
                         min_p=2,max_p=3,min_t=1273,max_t=1573,hydrous=True,iron_frac=True)
        # from table 3
        self.static_iron_frac = 0.545
        self.iso_mechanism = [XStalWaterSimpleIronExp(iron_frac=-79, iron_frac_err=3,
                                preexp_c_constant=0.67,preexp_c_constant_err=0.07,
                                enthalpy=183,enthalpy_err=9,water_content_type='wt%',convert2ev=True,
                                preexp=5.2,preexp_err = 0.4,log_sigma_err=True)]
        
    def isotropic_conductivity(self,temperature=None,**kwargs):
        if 'x_fe' not in kwargs.keys():
            kwargs['x_fe'] = self.static_iron_frac
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.iso_mechanism:
            
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
    
    def max_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
        
    def min_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
    def polycrystalline_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
class XuOrtho(DryPyroxenite):
    """
    a pyroxenite, fit to a range of mantle rocks.
    static iron fraction of x_fe = 0.545
    
    """
    def __init__(self):
        super().__init__(year=1999,author="Yousheng Xu, Thomas J. Shankland",
                         name='XuYousheng99_opx',title ="Electrical conductivity of orthopyroxene"+\
                         "and its high pressure phase",mineral_name='Orthopyroxene',
                         min_p=4,max_p=6,min_t=1273,max_t=1673)
        # from table 3
        self.iso_mechanism = [XStalTypical(enthalpy=1.8,enthalpy_err=0.02,
                                preexp=3.72,preexp_err = 0.06,log_sigma_err=True)]
        
    def isotropic_conductivity(self,temperature=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.iso_mechanism:
            
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
    
    def max_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
        
    def min_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
    def polycrystalline_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
class XuClino(DryPyroxenite):
    """
    a pyroxenite, fit to a range of mantle rocks.
    static iron fraction of x_fe = 0.545
    
    """
    def __init__(self):
        super().__init__(year=1999,author="Yousheng Xu, Thomas J. Shankland",
                         name='XuYousheng99_cpx',title ="Electrical conductivity of orthopyroxene"+\
                         "and its high pressure phase",mineral_name='Clinopyroxene',
                         min_p=12,max_p=14,min_t=1273,max_t=1673)
        # from table 3
        self.static_iron_frac = 0.545
        self.iso_mechanism = [XStalTypical(enthalpy=1.87,enthalpy_err=0.02,
                                preexp=3.25,preexp_err = 0.07,log_sigma_err=True)]
        
    def isotropic_conductivity(self,temperature=None,**kwargs):
        if 'x_fe' not in kwargs.keys():
            kwargs['x_fe'] = self.static_iron_frac
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.iso_mechanism:
            
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
    
    def max_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
        
    def min_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
    def polycrystalline_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
class Huebner88a(DryPyroxenite):
    """
    a nearly XFe = 0 (0.027) diopside (cpx?) collected from a vug in New York
    
    """
    def __init__(self):
        super().__init__(year=1988,author="J. Stephen Huebner, Donald E. Voigt",
                         name='Huebner88_dk7_opx',title ="Electrical conductivity of diopside"+\
                         ": Evidence for oxygen vacancies",mineral_name='Clinopyroxene',
                         min_p=0,max_p=0.1,min_t=1173,max_t=1573,fO2_dependent=True,anisotropic=True)
        self.xfe = 0.027
        self.xmg  = 1-self.xfe
        self.crystal_001 = [ OldEmpiricalLogEqn(1.43, -10850, -0.18)]
        self.crystal_010 = [ OldEmpiricalLogEqn(0.88, -10170, -0.12)]
        self.crystal_100 = [ OldEmpiricalLogEqn(1.58, -11080, -0.15)]
        
    
class Huebner88b(DryPyroxenite):
    """
    XFe = 0.12 diopside (cpx?) from Brazil
    
    """
    def __init__(self):
        super().__init__(year=1988,author="J. Stephen Huebner, Donald E. Voigt",
                         name='Huebner88_mal_opx',title ="Electrical conductivity of diopside"+\
                         ": Evidence for oxygen vacancies",mineral_name='Clinopyroxene',
                         min_p=0,max_p=0.1,min_t=1173,max_t=1573,fO2_dependent=True,anisotropic=True)
        self.xfe = 0.12
        self.xmg  = 1-self.xfe
        self.crystal_100 = [ OldEmpiricalLogEqn(-0.89, -5270, -0.03)]
        self.crystal_010 = [ OldEmpiricalLogEqn(-0.25, -4640, -0.02)]
        
    def isotropic_conductivity(self,temperature=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            conductivity_010 = np.zeros(temperature.shape)
            conductivity_100 = np.zeros(temperature.shape)
        else:
            conductivity_010 = 0
            conductivity_100 = 0
            
        for mechanism_010, mechanism_100 in zip(self.crystal_010,self.crystal_100):
            conductivity_010+=mechanism_010.get_conductivity(temperature=temperature,**kwargs)
            conductivity_100+=mechanism_100.get_conductivity(temperature=temperature,**kwargs)
            
        conductivity = np.sqrt(conductivity_010*conductivity_100)
        return conductivity
    
    def max_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
        
    def min_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
    
    def polycrystalline_conductivity(self,**kwargs):
        return self.isotropic_conductivity(**kwargs)
        
class Zhang16Opx(DryPyroxenite):
    
    def __init__(self):
        super().__init__(year=2016, author = "Baohua Zhang, Takashi Yoshino",\
                         name='ZY16_opx',title='Effect of temperature, pressure and'+\
                         'iron content on the electrical conductivity of orthopyroxene',mineral_name='Orthopyroxene',
                         min_p=1.5,max_p=5,min_t=473,max_t=1773,iron_frac=True)
        self.mechanisms = [ XStalTypical(enthalpy=2.51,enthalpy_err=0.1,
                                          activation_volume=4.15,activation_volume_err=0.71,
                                preexp=885610,preexp_err = 44561,log_sigma_err=False)
                           ,
                            XStalXfeComplex(enthalpy=2.33,enthalpy_err=0.01,
                                          activation_volume=1.06,activation_volume_err=0.44,
                                preexp=163,preexp_err = 1,log_sigma_err=False,
                                           a=1.99,a_err=0.06,b=0.12,b_err=0.076)]
        
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
    
class DaiKarato09(DryPyroxenite):
    
    def __init__(self):
        super().__init__(year=2009,author = "Lidong Dai, Shun-ichiro Karato",
                         name='DK09_opx',title='Electrical conductivity of orthopyroxene:'+\
                         'Implications for the water content of the asthenosphere',mineral_name='Orthopyroxene',
                         min_p=7,max_p=8,min_t=873,max_t=1473,hydrous=True)
    
                         #147
        self.mechanisms = [ XStalTypical(enthalpy=147,enthalpy_err=5,
                                preexp=2.4,preexp_err =0.4,log_sigma_err=True,
                                         convert2ev=True),
                           
                           XStalWaterSimple(enthalpy=82,enthalpy_err=5,
                                preexp=2.6,preexp_err =0.4,log_sigma_err=True,
                                        preexp_c_constant=0.62,water_content_type='wt%',
                                        preexp_c_constant_err=0.1,
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
    
    
class XiaoY12(DryPyroxenite):
    
    def __init__(self):
        super().__init__(year=2012,author = "Xiaozhi Yang, Hans Keppler,"+\
                         "Catherine McCammon & Huaiwei Ni ",
                         name='XiaoY12_opx',title='Electrical conductivity of'+\
                         'orthopyroxene and plagioclase in the lower crust',mineral_name='Orthopyroxene',
                         min_p=0.6,max_p=1.2,min_t=573,max_t=1273,hydrous=True)
    
                         
        self.mechanisms = [ XStalTypical(enthalpy=105,enthalpy_err=3,
                                preexp=2.39,preexp_err =0.18,log_sigma_err=True,
                                         convert2ev=True),
                           XStalWaterSimple(enthalpy=81,enthalpy_err=1,
                                preexp=3.83,preexp_err =0.1,log_sigma_err=True,
                                        preexp_c_constant=0.9,water_content_type='wt%',
                                        preexp_c_constant_err=0.04,
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
    
class XiaoY12CPX(DryPyroxenite):
    """
    averaging dry and hydrous fit results.
    """
    
    def __init__(self):
        super().__init__(year=2012,author = "Xiaozhi Yang and Catherine McCammon",
                         name='YMC12_cpx',title='Fe3+-rich augite and high '+\
                         'electrical conductivity in the deep lithosphere',mineral_name='Clinopyroxene(Augite)',
                         min_p=0.6,max_p=1.2,min_t=473,max_t=1373,hydrous=True)
    
                         
        self.mechanisms = [ XStalTypical(enthalpy=90,enthalpy_err=2,
                                preexp=2.,preexp_err =0.12,log_sigma_err=True,
                                         convert2ev=True),
                           XStalWaterSimple(enthalpy=72,enthalpy_err=1,
                                preexp=3.6,preexp_err =0.14,log_sigma_err=True,
                                        preexp_c_constant=1,water_content_type='wt%',
                                        preexp_c_constant_err=0.3,
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
    
class Zhang12opx(DryPyroxenite):
    """
    averaging dry and hydrous fit results.
    """
    
    def __init__(self):
        super().__init__(year=2012,author = "Baohua Zhang, Takashi Yoshino, Xiaoping Wu, Takuya Matsuzakia"+\
                         "Shuangming Shan, Tomoo Katsura",
                         name='Zhang12_opx',title='Electrical conductivity of enstatite as a function of water content: '+\
                         'Implications for the electrical structure in the upper mantle',mineral_name='Orthopyroxene (enstatite)',
                         min_p=2,max_p=4,min_t=1000,max_t=1723,hydrous=True)
    
                         
        self.mechanisms = [ XStalTypical(enthalpy=1.88,enthalpy_err=0.07,
                                preexp=3.99,preexp_err =0.23,log_sigma_err=True),
                           XStalWaterEmpirical(enthalpy=0.84,enthalpy_err=0.03,
                                preexp=2.58,preexp_err =0.14,log_sigma_err=True,
                                preexp_c_constant=1,water_content_type='wt%',
                                preexp_c_constant_err=0.0,enthalpy_exp=1/3,enthalpy_exp_err=0,
                                enthalpy_multiplier=-0.08,enthalpy_multiplier_err=0)
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