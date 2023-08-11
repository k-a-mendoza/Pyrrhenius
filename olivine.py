import numpy as np
from .mineralbehavior import Mineral, XStalTypical, XStalWaterEmpirical, XStalWaterSimple,\
XStalfO2Premult, XStalWaterInvT, calc_QFM, DiffusionHConduction

def get_experiments():
    experiments = [Gardes14(),DuFrane05(),Novella17(),Fei2020(),SEO3(),Yoshino09(),Yoshino16()]
    exps = {}
    for exp in experiments:
        exps[exp.name]=exp
    return exps

class DryOlivine(Mineral):
    k = 8.617e-5
    r = 8.31446e-3 #kj/mol k
    h2e = 0.01036427
    norm_v = 1/96.49 # ev mol/Gpa cm^3

    def __init__(self,mineral_name='Olivine',**kwargs):
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
    
    
class Gardes14(DryOlivine):

    def __init__(self):
        super().__init__(year=2014,author="Emmanuel Gardes, Fabrice Gaillard, Pascal", name='Gardes14',
                         title ="Towards a unified hydrous olivine conductivity law",min_p=0,max_p=6,min_t=473,max_t=1973,hydrous=True)
        # from table 3
        self.crystal_100 = [ XStalTypical(preexp=5.92,preexp_err=1.99/2,enthalpy=261*self.h2e,enthalpy_err=48*self.h2e/2,log_sigma_err=True),
                             XStalTypical(preexp=2.19,preexp_err=1.09/2,enthalpy=146*self.h2e,enthalpy_err=27*self.h2e/2,log_sigma_err=True),
                             XStalWaterEmpirical(preexp=-1.28,preexp_err=0.57/2,enthalpy=92*self.h2e,enthalpy_err=9*self.h2e/2,
                                                 enthalpy_multiplier=-2.99*self.h2e, enthalpy_multiplier_err=0.69*self.h2e/2,log_sigma_err=True)]
        self.crystal_010 = [ XStalTypical(preexp=5.45,preexp_err=2.07/2,enthalpy=268*self.h2e,enthalpy_err=55*self.h2e/2,log_sigma_err=True),
                             XStalTypical(preexp=2.38,preexp_err=1.1/2,enthalpy=141*self.h2e,enthalpy_err=21*self.h2e/2,log_sigma_err=True),
                             XStalWaterEmpirical(preexp=-0.51,preexp_err=0.99/2,enthalpy=95*self.h2e,enthalpy_err=13*self.h2e/2,
                                                 enthalpy_multiplier=-1.03*self.h2e, enthalpy_multiplier_err=0.69*self.h2e/2,log_sigma_err=True)]
        self.crystal_001 = [ XStalTypical(preexp=4.67,preexp_err=1.94/2,enthalpy=234*self.h2e,enthalpy_err=47*self.h2e/2,log_sigma_err=True),
                             XStalTypical(preexp=2.51,preexp_err=1.06/2,enthalpy=146*self.h2e,enthalpy_err=24*self.h2e/2,log_sigma_err=True),
                             XStalWaterEmpirical(preexp=-1.77,preexp_err=0.69/2,enthalpy=81*self.h2e,enthalpy_err=9*self.h2e/2,
                                                 enthalpy_multiplier=-2.25*self.h2e, enthalpy_multiplier_err=0.63*self.h2e/2,log_sigma_err=True)]
        self._preexp     = np.power(10,np.asarray((5.07, 2.34,-1.17)))
        self._preexp_err = np.power(10,np.asarray((1.32, 0.67,0.45)))
        self._cw_mult        = (1,1,5)
        self._cw_mult_err    = (0,0,2)
        self._a              = (0,0,2.8)
        self._a_err          = (0,0,0.55)
        self._enthalpy       = (239, 144,89)
        self._enthalpy_err   = (46, 16,7)
    
class DuFrane05(DryOlivine):

    def __init__(self):
        super().__init__(year=2005,author="Wyatt L. Du Frane,Jeffery J. Roberts,Daniel A. Toffelmier,James A. Tyburczy", name='duFrane05',
                         title ="Anisotropic Olivine Conductivity",min_p=0.4,max_p=6,min_t=1173,max_t=1673,anisotropic=True,fO2_dependent=True)
        # from table 3
        self.crystal_100 = [ XStalTypical(preexp=0.0119,preexp_err=0.0051,enthalpy=0.322,enthalpy_err=0.048),
                             XStalfO2Premult(preexp=0.432,preexp_err=0.147,enthalpy=0.322,enthalpy_err=0.048,fO2_exp=2/11)]
                             
        self.crystal_010 = [ XStalTypical(preexp=0.0798,preexp_err=0.0147,enthalpy=0.561,enthalpy_err=0.021),
                             XStalfO2Premult(preexp=2.39,preexp_err=0.35,enthalpy=0.561,enthalpy_err=0.021,fO2_exp=2/11)]
        
        self.crystal_001 = [ XStalTypical(preexp=0.293,preexp_err=0.072,enthalpy=0.71,enthalpy_err=0.027),
                             XStalfO2Premult(preexp=15.3,preexp_err=2.9,enthalpy=0.71,enthalpy_err=0.027,fO2_exp=2/11)]
        
        
class Novella17(DryOlivine):

    def __init__(self,density_correction_factor=1/2):
        super().__init__(year=2017,author="Davide Novella, Benjamin Jacobsen, Peter K. Weber, James A. Tyburczy, Frederick J. Ryerson & Wyatt L. Du Frane", name='Novella17',
                         title ="Hydrogen Self-diffusion in single crystal olivine and electrical conductivity of the Earth's Mantle",max_p=6,min_t=1173,anisotropic=True,hydrous=True)
        # from table 3

        duFrane = DuFrane05()
        
        self.crystal_100 = duFrane.crystal_100 + [DiffusionHConduction(preexp=-0.7,preexp_err=0.9,enthalpy=229,enthalpy_err=18,
density_correction_factor=density_correction_factor,log_sigma_err=True,density=3350)]
                             
        self.crystal_010 = duFrane.crystal_010 + [DiffusionHConduction(preexp=-5,preexp_err=0.9,enthalpy=172,enthalpy_err=19,
density_correction_factor=density_correction_factor,log_sigma_err=True)]
        
        self.crystal_001 = duFrane.crystal_001 + [DiffusionHConduction(preexp=-3.5,preexp_err=0.4,enthalpy=188,enthalpy_err=8,
density_correction_factor=density_correction_factor,log_sigma_err=True)]
        
class Sun19(DryOlivine):

    def __init__(self,density_correction_factor=1/2):
        super().__init__(year=2017,author="Davide Novella, Benjamin Jacobsen, Peter K. Weber, James A. Tyburczy, Frederick J. Ryerson & Wyatt L. Du Frane", name='Novella17',
                         title ="Hydrogen Self-diffusion in single crystal olivine and electrical conductivity of the Earth's Mantle",max_p=6,min_t=1173)
        # from table 3

        duFrane = DuFrane05()
        
        self.crystal_100 = duFrane.crystal_100 + [DiffusionHConduction(preexp=-5,preexp_err=0.5,enthalpy=140,enthalpy_err=11,
density_correction_factor=density_correction_factor,log_sigma_err=True)]
                             
        self.crystal_010 = duFrane.crystal_010 + [DiffusionHConduction(preexp=-5.1,preexp_err=1,enthalpy=153,enthalpy_err=21,
density_correction_factor=density_correction_factor,log_sigma_err=True)]
        
        self.crystal_001 = duFrane.crystal_001 + [DiffusionHConduction(preexp=-8.2,preexp_err=0.9,enthalpy=95,enthalpy_err=19,
density_correction_factor=density_correction_factor,log_sigma_err=True)]
        
   
    
class Fei2020(DryOlivine):

    def __init__(self):
        super().__init__(year=2020,author="Hongzhan Fei, Dmitry Druzhbin, Tomoo Katsura", name='FeiH20',
                         title ="Anisotropic Olivine Conductivity",min_p=2,max_p=10,min_t=1450,max_t=2180,anisotropic=True,hydrous=True)
        # from table 3
        self.crystal_100 = [ XStalTypical(preexp=1.5,preexp_err=0.2,enthalpy=117,enthalpy_err=13,convert2ev=True,
                                          activation_volume = 0.5, activation_volume_err=0.7,log_sigma_err=True),
                             XStalWaterInvT(preexp=9.9,preexp_err=0.5,enthalpy=337,enthalpy_err=15,log_sigma_err=True,convert2ev=True,
                                          activation_volume = 4.2, activation_volume_err=0.4, preexp_c_constant=1.3, preexp_c_constant_err=0.2,
                           calibrate=3)]
                             
        self.crystal_010 = [ XStalTypical(preexp=2.7,preexp_err=0.2,enthalpy=163,enthalpy_err=25,convert2ev=True,
                                         activation_volume = -0.2, activation_volume_err=0.2,log_sigma_err=True),
                             XStalWaterInvT(preexp=11.6,preexp_err=0.4,enthalpy=396,enthalpy_err=25,log_sigma_err=True,convert2ev=True,
                                         activation_volume = 3.2, activation_volume_err=0.2, preexp_c_constant=1.3, preexp_c_constant_err=0.2
                                           ,
                           calibrate=3)]
        
        self.crystal_001 = [ XStalTypical(preexp=2.7,preexp_err=0.3,enthalpy=139,enthalpy_err=13,convert2ev=True,
                             activation_volume = 0.5, activation_volume_err=0.7,log_sigma_err=True),
                             XStalWaterInvT(preexp=11.8,preexp_err=0.4,enthalpy=385,enthalpy_err=15,log_sigma_err=True,convert2ev=True,
                             activation_volume = 4.1, activation_volume_err=0.47,preexp_c_constant=1.3, preexp_c_constant_err=0.2,
                           calibrate=3)]   
    
class SEO3(DryOlivine):

    def __init__(self):
        super().__init__(year=2006,author="S. Constable", name='SEO3',
                         title ="SEO3: A new model of olivine electrical conductivity",max_p=6,min_t=1173,max_t=1873,fO2_dependent=True)
        # from table 3
        self._e = 1.602e-19
        self._preexps    = np.asarray([5.06e24, 4.58e26, 12.2e-6, 2.72e-6, 3.33e24, 6.21e30])
        self._enthalpies = np.asarray([0.357,0.752,1.05,1.09,0.02,1.83])
        self._fo2_exp = 1/6
        self._preexps_err = self._preexps*0.05
        self._enthalpies_err = self._enthalpies*0.05
        
    def isotropic_conductivity(self,**kwargs):
        return self.mechanism(**kwargs)
    
    def max_anisotropic_conductivity(self,**kwargs):
        return self.mechanism(**kwargs)
        
    def min_anisotropic_conductivity(self,**kwargs):
        return self.mechanism(**kwargs)
    
    def polycrystalline_conductivity(self,temperature=None,pressure=None,extrapolate=False,**kwargs):
        return self.mechanism(**kwargs)
    
    def get_const(self):
        return list(self._preexps) + list(self._enthalpies)
                            
    def get_err(self):
        return list(self._preexps_err) + list(self._enthalpies_err)
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None,stochastic=False,log10_fO2=None,log10_fO2_offset=0,pressure=0,**kwargs):
        if stochastic:
            sig1, sig2, sig3, sig4, sig5, sig6, h1, h2, h3, h4, h5, h6  = self.sample_const()
        else:
            sig1, sig2, sig3, sig4, sig5, sig6, h1, h2, h3, h4, h5, h6  = self.get_const()
            
        if log10_fO2 is None:
            log10_fO2 = calc_QFM(temperature,pressure) + log10_fO2_offset
        fO2 = 10**log10_fO2
        fO2_6 = fO2**self._fo2_exp
        kt  = self.k*temperature
        
        bfe = sig1*np.exp(-h1/kt)
        bmg = sig2*np.exp(-h2/kt)
        ufe = sig3*np.exp(-h3/kt)
        umg = sig4*np.exp(-h4/kt)
        fo2_1 = sig5*np.exp(-h5/kt)*fO2_6
        fo2_2 = sig6*np.exp(-h6/kt)*fO2_6
        
        return (bfe+fo2_1)*ufe*self._e+2*(bmg+fo2_2)*umg*self._e
    
    
class Yoshino09(DryOlivine):

    def __init__(self):
        super().__init__(year=2009,author="Takashi Yoshino, Takuya Matsuzaki, Anton Shatskiy, Tomoo Katsura", name='Yoshino09',
                         title ="The effect of water on the electrical conductivity of olivine aggregates and its"+\
                         "implications for the electrical structure of the upper mantle",max_p=11,min_t=500,max_t=1500,hydrous=True)
        self.polycrystal_data = [ XStalTypical(preexp=4.73,preexp_err=0.53,enthalpy=2.31,enthalpy_err=0.07,log_sigma_err=True),
                                  XStalTypical(preexp=2.98,preexp_err=0.85,enthalpy=1.71,enthalpy_err=0.04,log_sigma_err=True),
                                  XStalWaterEmpirical(preexp=1.9,preexp_err=0.44,enthalpy=0.92,enthalpy_err=0.04,log_sigma_err=True,wt_pct=True,
                                                      enthalpy_multiplier=-0.16,water_content_type='wt%', 
                                                      enthalpy_multiplier_err=0.02)]
        
        
    def isotropic_conductivity(self,temperature=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            conductivity = np.zeros(temperature.shape)
        else:
            conductivity = 0
        for mechanism in self.polycrystal_data:
            conductivity+=mechanism.get_conductivity(temperature=temperature,**kwargs)
        return conductivity
    
    def max_anisotropic_conductivity(self,**kwargs):
        return  self.isotropic_conductivity(self,**kwargs)
        
    def min_anisotropic_conductivity(self,**kwargs):
        return self.isotropic_conductivity(self,**kwargs)
    
    def polycrystalline_conductivity(self,**kwargs):
        return self.isotropic_conductivity(self,**kwargs)
        
class Yoshino16(DryOlivine):
   
    def __init__(self):
        super().__init__(year=2016,author="Takashi Yoshino, Baohua Zhang, Brandon Rhymer, Chengcheng Zhao, and Hongzhan Fei", name='Yoshino16',
                         title ="Pressure dependence of electrical conductivity in forsterite",max_p=15,max_t=2100,fO2_dependent=True,anisotropic=True)
        class SEO3Wrapper:
            
            def __init__(self):
                self.mechanism = SEO3()
                
            def get_conductivity(self,**kwargs):
                return self.mechanism.mechanism(**kwargs)
            
        self.crystal_100 = [ XStalTypical(preexp=1.1e5,preexp_err=6,enthalpy=1.75,enthalpy_err=23,
                                          activation_volume=18.6,activation_volume_err=56),
                             XStalTypical(preexp=1.3e6,preexp_err=52,enthalpy=3.05,enthalpy_err=21,
                                          activation_volume=0.3,activation_volume_err=3),
                            SEO3Wrapper()]
        
        self.crystal_010 = [ XStalTypical(preexp=3.4e5,preexp_err=35,enthalpy=2.38,enthalpy_err=11,
                                          activation_volume=6.9,activation_volume_err=5),
                             XStalTypical(preexp=1.1e4,preexp_err=4,enthalpy=2.23,enthalpy_err=6,
                                          activation_volume=-0.8,activation_volume_err=3),
                            SEO3Wrapper()]
        
        self.crystal_001 = [ XStalTypical(preexp=1e7,preexp_err=14,enthalpy=2.71,enthalpy_err=22,
                                          activation_volume=8.3,activation_volume_err=17),
                             XStalTypical(preexp=1.7e6,preexp_err=19,enthalpy=2.93,enthalpy_err=21,
                                          activation_volume=-0.3,activation_volume_err=4),
                            SEO3Wrapper()]
    
    