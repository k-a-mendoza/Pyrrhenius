from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(10)

def get_experiments():
    return [Yoshino09(),Yoshino2012(),YoshinoPolyC16(),Gardes14(),DuFrane(), Xu(), Constable92(), Poe2010(), SEO3(),DaiHu1(),OldPeridotite()] 



def calc_QFM(T,P):
    P_bars = P*10000
    if isinstance(T,np.ndarray) and isinstance(P,np.ndarray):
        log_fO2 = np.ones(T.shape)
        log_fO2[T<573]  = (-26445.3/T[T<573] ) + 10.344 + 0.092 * (P_bars[T<573]-1)/T[T<573]
        log_fO2[T>=573] = (-25096.3/T[T>=573]) + 8.735  + 0.11  * (P_bars[T>=573]-1)/T[T>=573]
    elif isinstance(T,np.ndarray) and not isinstance(P,np.ndarray):
        log_fO2 = np.ones(T.shape)
        log_fO2[T<573]  = (-26445.3/T[T<573] ) + 10.344 + 0.092 * (P_bars-1)/T[T<573]
        log_fO2[T>=573] = (-25096.3/T[T>=573]) + 8.735  + 0.11  * (P_bars-1)/T[T>=573]
    elif not isinstance(T,np.ndarray) and isinstance(P,np.ndarray):
        log_fO2 = np.ones(P.shape)
        if T > 573:
            log_fO2 = (-26445.3/T ) + 10.344 + 0.092 * (P_bars-1)/T
        else:
            log_fO2 = (-25096.3/T) + 8.735  + 0.11  * (P_bars-1)/T
    else:
        if T > 573:
            log_fO2 = (-26445.3/T ) + 10.344 + 0.092 * (P_bars-1)/T
        else:
            log_fO2 = (-25096.3/T) + 8.735  + 0.11  * (P_bars-1)/T
    return log_fO2



class DryOlivine:
    k = 8.617e-5
    r = 8.31446e-3 #kj/mol k
    norm_v = 1/96.49 # ev mol/Gpa cm^3
    
    a = -25_096.3
    b = 8.735
    c = 0.110
    
    convert2bars = 10_000
    def __init__(self,publication,ptspace):
        self.ptspace     = ptspace
        self.publication = publication
        self._mechanisms = []
    
    def conductivity(self, temperature=1000, pressure=1, stochastic=False,**kwargs):
        """
        calculates the conductivity of dry olivine
        
        Parameters
        ==========
        
        temperature : double, np.ndarray
            temperature in kelvin. Can either be a double or an nd numpy array with
            the same dimensions as pressure
            
        Pressure : double, np.ndarray
            pressure in GPa. Can either be a double or an nd numpy array with the same dimensions as temperature
            
        x_fe : double, np.ndarray
            iron fraction between 0 - 1. Can either be a double or an nd numpy array with the same dimensions as 
            temperature.
            
        log10_fO2 : double, np.ndarray
            Log10 oxygen fugacity. Can either be a  double or an nd numpy array with the same dimensions 
            as temperature
            
        stochastic : bool
            whether to estimate the conductivity using best-fit constants or use the standard errors to draw a conductivity value from a distribution.
            Defaults to False
        """
        conductivity = 0
        if 'log10_fO2' not in kwargs.keys():
            kwargs['log10_fO2'] =  calc_QFM(temperature,pressure)
        for mechanism in self._mechanisms:
            conductivity = conductivity + mechanism(temperature=temperature,stochastic=stochastic,pressure=pressure, **kwargs)
        
        return conductivity
    
@dataclass  
class Publication:
    authors: str
    year : int
    title : str
    
    def toString(self):
        print(self.authors)
        print(self.title)
        print(self.year)
    
@dataclass
class ValidPTspace:
    pressure : tuple
    temperature : tuple
    
    def temp_mask(self,temperature):
        return (temperature < self.temperature[0]) | (temperature > self.temperature[1])
    
    
class OldPeridotite(DryOlivine):
    stdev = 0.05
    def __init__(self):
        publication = Publication(authors="Z. DVorak",year=1973,
                                  title="Electrical conductivity of several samples of olivinites, peridotites, and dunites, as a function of pressure and temperature")
        ptspace = ValidPTspace(pressure=(0,np.inf),temperature=[473,973])
        super().__init__(publication,ptspace)
        self.exps    = np.asarray([0.4,1.20,2.54])
        self.preexps = 10**np.asarray([-4.78,2.86,10.06])
        self.name='dvorak73'
        self._mechanisms.append(self.mechanism)
        
    def get_constants(self):
        return self.exps[0],self.exps[1],self.exps[2], self.preexps[0],self.preexps[1],self.preexps[2]
        
        
    def sample_constants(self):
        exps = self.get_constants()
        new_constants = []
        for e in exps:
            new_constants.append(np.random.normal(e,e*self.stdev))
        return new_constants
        
    def mechanism(self,temperature=1000,stochastic=False, **kwargs):
        if stochastic:
            e1,e2,e3,sig1,sig2,sig3= self.sample_constants()
        else:
            e1,e2,e3,sig1,sig2,sig3 = self.get_constants()
      
        one = sig1*np.exp(-e1/(self.k*temperature))
        two = sig2*np.exp(-e2/(self.k*temperature))
        three = sig3*np.exp(-e3/(self.k*temperature))
        
        return  three
    
    
class Yoshino2012(DryOlivine):
    stdev = 0.05 #5% standard deviation in values
    def __init__(self):
        publication = Publication(authors="Yoshino et al",year=2012,
                                  title="Effect of temperature, pressure, and iron content on the electrical conductivity of olivine and its high pressure pseudomorphs")
        ptspace = ValidPTspace(pressure=(6,19),temperature=[500,1600])
        super().__init__(publication,ptspace)
        # from table 3
        self.preexp = 525
        self.activation = 1.96
        self.volume  = -0.01 # cm^3/mol
        self.alpha = 1.49
        self.beta = 0.22
        self._mechanisms.append(self.mechanism)
        self.name = "Yoshino12"
    
    def get_constants(self):
        return self.preexp, self.activation, self.volume, self.alpha, self.beta
        
        
    def sample_constants(self):
        constants = self.get_constants()
        new_constants = []
        for constant in constants:
            sigma = abs(self.stdev*constant)
            new_constants.append(np.random.normal(constant,sigma))
        return new_constants
        
    def mechanism(self,temperature=1000,pressure=6, x_fe=0.1,stochastic=False,stochastic_iron=False, **kwargs):
        if stochastic:
            preexp, activation, volume, alpha, beta = self.sample_constants()
        else:
            preexp, activation, volume, alpha, beta = self.get_constants()
        if stochastic_iron:
            x_fes = np.random.uniform(0.05,0.2)
        else:
            x_fes = x_fe
        enthalpy = activation - alpha*np.power(x_fes,1/3) + pressure*(volume - beta*x_fes)*self.norm_v # divided by conversion factor to electron-volts
        
        return preexp*x_fes*np.exp(-enthalpy/(self.k*temperature))
    
    
    
class Gardes14(DryOlivine):

    def __init__(self):
        publication = Publication(authors="Emmanuel Gardes, Fabrice Gaillard, Pascal",year=2014,
                                  title="Towards a unified hydrous olivine conductivity law")
        ptspace = ValidPTspace(pressure=(0,6),temperature=[200+273,1700+273])
        super().__init__(publication,ptspace)
        # from table 3
        self._preexp     = np.power(10,np.asarray((5.07, 2.34,-1.17)))
        self._preexp_err = np.power(10,np.asarray((1.32, 0.67,0.45)))
        self._cw_mult        = (1,1,5)
        self._cw_mult_err    = (0,0,2)
        self._a              = (0,0,2.8)
        self._a_err          = (0,0,0.55)
        self._enthalpy       = (239, 144,89)
        self._enthalpy_err   = (46, 16,7)
        self._mechanisms.append(self._mechanism)
        self.name = "Gardes14"
    
    def _get_const(self):
        return self._preexp[0],self._preexp[1],self._preexp[2],\
               self._cw_mult[0],self._cw_mult[1],self._cw_mult[2],\
               self._a[0],self._a[1],self._a[2],\
               self._enthalpy[0],self._enthalpy[1],self._enthalpy[2]
    
    def _get_err(self):
        return self._preexp_err[0],self._preexp_err[1],self._preexp_err[2],\
               self._cw_mult_err[0],self._cw_mult_err[1],self._cw_mult_err[2],\
               self._a_err[0],self._a_err[1],self._a_err[2],\
               self._enthalpy_err[0],self._enthalpy_err[1],self._enthalpy_err[2]
    
    def _get_sample(self):
        variables = self._get_const()
        errors    = self._get_err()
        new_vars = []
        for index in range(len(variables)):
            new_vars.append(np.random.normal(variables[index],errors[index]))
        return new_vars
                            
    def _mechanism(self,temperature=1000,stochastic=False,with_hydrous=False,**kwargs):
        if stochastic:
            variables = self._get_sample()
        else:
            variables = self._get_const()
        conductivity = 0
        
        variables = np.asarray(variables).reshape((len(variables)//4,4),order='F')
        if with_hydrous:
            n = 3
        else:
            n = 2
        for i in range(n):
            sigma, cw, alpha, enthalpy    = variables[i,:].tolist()   
            conductivity+= sigma*cw*np.exp(-(enthalpy - alpha*np.power(cw,1/3))/(self.r*temperature))
            
        return conductivity
        
        
class DuFrane(DryOlivine):

    def __init__(self):
        publication = Publication(authors="Wyatt L. Du Frane,Jeffery J. Roberts,Daniel A. Toffelmier,James A. Tyburczy",year=2005,
                                  title="Anisotropic Olivine Conductivity")
        ptspace = ValidPTspace(pressure=(0.4,6),temperature=[900+273,1400+273])
        super().__init__(publication,ptspace)
        # from table 3
        self._preexps = (2.51, 0.0653)
        self._preexp_err = (0.57,0.0187)
        self._fo2power = 2/11
        self._enthalpy = 0.531
        self._enthalpy_err = 0.032
        self._mechanisms.append(self.mechanism)
        self.name = "DuFrane05"
    
    def get_const(self):
        return self._preexps[0], self._preexps[1], self._fo2power, self._enthalpy
                            
    def get_err(self):
        return self._preexp_err[0], self._preexp_err[1], 0, self._enthalpy_err
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None, log10_fO2=1,stochastic=False,**kwargs):
        if isinstance(log10_fO2,int):
            log10_fO2 = float(log10_fO2)
        fo2 = np.power(10,log10_fO2)
        if stochastic:
            sigma_1, sigma_2, fo2_power, enthalpy = self.sample_const()
        else:
            sigma_1, sigma_2, fo2_power, enthalpy = self.get_const()
        
        return (sigma_1*np.power(fo2,fo2_power) + sigma_2)*np.exp(-enthalpy/(self.k*temperature))
        
        
class Xu(DryOlivine):

    def __init__(self):
        publication = Publication(authors="Yousheng Xu, TJ Shankland, AG Duba",year=2000,
                                  title="Pressure effect on electrical conductivity of mantle olivine")
        ptspace = ValidPTspace(pressure=(4,10),temperature=[900+273,1500+273])
        super().__init__(publication,ptspace)
        # from table 3
        self._preexp = np.asarray([2127,469,245])
        self._h      = np.asarray([1.79,1.5,1.35])
        self._v      = np.asarray([0,0.86,0.96])
        self._preexp_err = np.power(10,np.asarray([0.28,0.14,0.06]))
        self._h_err = np.asarray([0.08,0.04,0.02])
        self._v_err = np.asarray([0,0.15,0.06])
        self._mechanisms.append(self.mechanism)
        self.name = "Xu00"
    
    def get_const(self):
        return self._preexp[0],self._preexp[1],self._preexp[2],\
               self._h[0], self._h[1], self._h[2],\
               self._v[0], self._v[1], self._v[2]
                            
    def get_err(self):
        return self._preexp_err[0],self._preexp_err[1],self._preexp_err[2],\
               self._h_err[0], self._h_err[1], self._h_err[2],\
               self._v_err[0], self._v_err[1], self._v_err[2]
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None, pressure = 1,stochastic=False,**kwargs):
        if stochastic:
            s1,s2,s3,u1,u2,u3,v1,v2,v3  = self.sample_const()
        else:
            s1,s2,s3,u1,u2,u3,v1,v2,v3 = self.get_const()
            
        sig1 = s1*np.exp(-(u1 + v1*pressure*self.norm_v)/(self.k*temperature))
        sig2 = s2*np.exp(-(u2 + v2*pressure*self.norm_v)/(self.k*temperature))
        sig3 = s3*np.exp(-(u3 + v3*pressure*self.norm_v)/(self.k*temperature))
        return np.power(sig1*sig2*sig3,1/3)
                            
class Constable92(DryOlivine):

    def __init__(self):
        publication = Publication(authors="S Constable, TJ Shankland, Duba Al",year=1992,
                                  title="Electrical conductivity of an isotropic olivine mantle")
        ptspace = ValidPTspace(pressure=(0,1),temperature=[720+273,1500+273])
        super().__init__(publication,ptspace)
        # from table 3
        self._preexp = (np.power(10,2.402),np.power(10,9.17))
        self._preexp_err = (self._preexp[0]*0.03, np.power(10,0.5))
        self._h = (1.6,4.25)
        self._h_err = (0.10*1.6,0.17)
        self._mechanisms.append(self.mechanism)
        self.name = "Constable92"
    
    def get_const(self):
        return self._preexp[0], self._preexp[1], self._h[0], self._h[1]
                            
    def get_err(self):
        return self._preexp_err[0], self._preexp_err[1], self._h_err[0], self._h_err[1]
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None,stochastic=False,**kwargs):
        if stochastic:
            sigma1, sigma2, u1, u2  = self.sample_const()
        else:
             sigma1, sigma2, u1, u2  = self.get_const()
                        
        first_term = sigma1*np.exp(-u1/(self.k*temperature))
        second_term = sigma2*np.exp(-u2/(self.k*temperature))
        return first_term + second_term
                        
class Yoshino09(DryOlivine):

    def __init__(self):
        publication = Publication(authors="T Yoshino, T Matsuzaki, A Shatskiy, T Katsura",year=2009,
                                  title="The effect of water on the electrical conductivity of olivine aggregates and its implications for the electrical structure of the upper mantle")
        ptspace = ValidPTspace(pressure=(0,1),temperature=[500,2000])
        super().__init__(publication,ptspace)
        # from table 3
        self._preexp = np.asarray((np.power(10,4.73),np.power(10,2.98)))
        self._preexp_err = self._preexp *0.05
        self._h = np.asarray((2.31,1.71))
        self._h_err = self._h*0.05
        self._mechanisms.append(self.mechanism)
        self.name = "Yoshino09"
    
    def get_const(self):
        return self._preexp[0], self._preexp[1], self._h[0], self._h[1]
                            
    def get_err(self):
        return self._preexp_err[0], self._preexp_err[1], self._h_err[0], self._h_err[1]
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None,stochastic=False,**kwargs):
        if stochastic:
            sigma1, sigma2, u1, u2  = self.sample_const()
        else:
             sigma1, sigma2, u1, u2  = self.get_const()
                        
        first_term  = sigma1*np.exp(-u1/(self.k*temperature))
        second_term = sigma2*np.exp(-u2/(self.k*temperature))
        return first_term + second_term
                            
                            
class Poe2010(DryOlivine):

    def __init__(self):
        publication = Publication(authors="BT Poe, C. Romanoc, F. Nestolad, JR Smythe",year=2010,
                                  title="Electrical conductivity anisotropy of dry and hydrous olivine at 8 GPa")
        ptspace = ValidPTspace(pressure=(8,8),temperature=[500,1600])
        super().__init__(publication,ptspace)
        # from table 3
        self._preexp = np.asarray((334, 13.80, 99))
        self._preexp_err = self._preexp*0.05
        self._h = np.asarray((1.46,1.12,1.29))
        self._h_err = self._h*0.05
        self._mechanisms.append(self.mechanism)
        self.name ="Poe10"
    
    def get_const(self):
        return self._preexp[0], self._preexp[1], self._preexp[2], self._h[0], self._h[1], self._h[2]
                            
    def get_err(self):
        return self._preexp_err[0], self._preexp_err[1], self._preexp_err[2], self._h_err[0], self._h_err[1], self._h_err[2]
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None,stochastic=False,**kwargs):
        if stochastic:
            sigma1, sigma2, sigma3, u1, u2, u3  = self.sample_const()
        else:
            sigma1, sigma2, sigma3, u1, u2, u3  = self.get_const()
        
        first_term  = sigma1*np.exp(-u1/(self.k*temperature))
        second_term = sigma2*np.exp(-u2/(self.k*temperature))
        third_term  = sigma3*np.exp(-u3/(self.k*temperature))
        return np.power(first_term*second_term*third_term,1/3)
    
                            
class YoshinoPolyC16(DryOlivine):
    err=0.05
    def __init__(self):
        publication = Publication(authors="T. Yoshino, B Zhang, B Rhymer, C Zhao, H Fei",year=2016,
                                  title="Pressure dependence of electrical conductivity in forsterite")
        ptspace = ValidPTspace(pressure=(3.5,15),temperature=[1300,2100])
        super().__init__(publication,ptspace)
        self._mg = np.asarray(((1e7,3.4e5,1.1e5),
                                  (2.71,2.38,1.75),
                                       (8.3,6.9,18.6)))
        self._ox = np.asarray(((1.7e6,1.1e4,1.3e6),
                                   (2.93,2.23,3.05),
                                   (-0.3,0.8,0.3)))
        self._SEO3 = SEO3()
        self._mechanisms.append(self.mechanism)
        self._mechanisms.append(self._SEO3.mechanism)
        self.ptspace.temperature[0]=self._SEO3.ptspace.temperature[0]
        self.name ="Yoshino16"
    
    def get_const(self):
        return self._mg, self._ox
                            
    def sample_const(self):
        mg, ox = self.get_const()
        return np.random.normal(mg,np.abs(mg)*self.err),np.random.normal(ox,np.abs(ox)*self.err)
    
    def mechanism(self,temperature=None,stochastic=False,pressure=1,**kwargs):
        if stochastic:
            mg, ox  = self.sample_const()
        else:
            mg, ox = self.get_const()
            
        xstal_conductances=[]
        for crystal_idx in range(3):
            sig_mg, ev_mg, v_mg = mg[:,0].tolist()
            sig_ox, ev_ox, v_ox = ox[:,0].tolist()
        
            mg_c  = sig_mg*np.exp(-(ev_mg + pressure*v_mg*self.norm_v)/(self.k*temperature))
            ox_c = sig_ox*np.exp(-(ev_ox + pressure*v_ox*self.norm_v)/(self.k*temperature))
            xstal_conductances.append(mg_c+ox_c)
    
        if isinstance(temperature,int):
            bulk = np.power(xstal_conductances[0]*xstal_conductances[1]*xstal_conductances[2],1/3)
        else:
            original_shape = temperature.shape
            ravel_length = len(temperature.ravel())
            conductances = np.zeros((ravel_length,3))
            for i in range(3):
                conductances[:,i]=xstal_conductances[i].ravel()
            bulk = np.power(conductances[:,0]*conductances[:,1]*conductances[:,2],1/3)
            bulk = bulk.reshape(original_shape)
        return bulk
    
    
class SEO3(DryOlivine):

    def __init__(self):
        publication = Publication(authors="S. Constable",year=2006,
                                  title="SO3: A new model of olivine electrical conductivity")
        ptspace = ValidPTspace(pressure=(0,0),temperature=[900+273,1600+273])
        super().__init__(publication,ptspace)
        # from table 3
        self._e = 1.602e-19
        self._preexps    = np.asarray([5.06e24, 4.58e26, 12.2e-6, 2.72e-6, 3.33e24, 6.21e30])
        self._enthalpies = np.asarray([0.357,0.752,1.05,1.09,0.02,1.83])
        self._fo2_exp = 1/6
        self._preexps_err = self._preexps*0.05
        self._enthalpies_err = self._enthalpies*0.05
        self._mechanisms.append(self.mechanism)
        self.name ="SEO3"
    
    def get_const(self):
        return self._preexps[0], self._preexps[1], self._preexps[2], self._preexps[3], self._preexps[4], self._preexps[5], \
                self._enthalpies[0],self._enthalpies[1],self._enthalpies[2],self._enthalpies[3],self._enthalpies[4],self._enthalpies[5]
                            
    def get_err(self):
        return self._preexps_err[0], self._preexps_err[1], self._preexps_err[2], self._preexps_err[3], self._preexps_err[4], self._preexps_err[5], \
                self._enthalpies_err[0],self._enthalpies_err[1],self._enthalpies_err[2],self._enthalpies_err[3],self._enthalpies_err[4],self._enthalpies_err[5]
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None,stochastic=False,log10_fO2=1,**kwargs):
        if stochastic:
            sig1, sig2, sig3, sig4, sig5, sig6, h1, h2, h3, h4, h5, h6  = self.sample_const()
        else:
            sig1, sig2, sig3, sig4, sig5, sig6, h1, h2, h3, h4, h5, h6  = self.get_const()
        
        fo2 = np.power(10,log10_fO2)
        bfe = sig1*np.exp(-h1/(self.k*temperature))
        bmg = sig2*np.exp(-h2/(self.k*temperature))
        ufe = sig3*np.exp(-h3/(self.k*temperature))
        umg = sig4*np.exp(-h4/(self.k*temperature))
        fo2_1 = sig5*np.exp(-h5/(self.k*temperature))*fo2**self._fo2_exp
        fo2_2 = sig6*np.exp(-h6/(self.k*temperature))*fo2**self._fo2_exp
        
        
        return (bfe+fo2_1)*ufe*self._e+2*(bmg+fo2_2)*umg*self._e
    
class DaiHu1(DryOlivine):

    def __init__(self):
        publication = Publication(authors="Lidong Dai Shun-ichiro Karato",year=2014,
                                  title="Influence of FeO and H on the electrical conductivity of olivine")
        ptspace = ValidPTspace(pressure=(1,10),temperature=[873,1473])
        super().__init__(publication,ptspace)
        # from table 3
        self._preexp = 2.77
        self._preexp_err = 0.38
        self._b = -1.19
        self._b_err = 0.09
        self._c = -63
        self._c_err = 4
        self._enthalpy = 162
        self._enthalpy_err = 9
        self._mechanisms.append(self.mechanism)
        self.name = "DaiH14a"
    
    def get_const(self):
        return self._preexp, self._b, self._c, self._enthalpy
                            
    def get_err(self):
        return self._preexp_err, self._b_err, self._c_err, self._enthalpy_err
    
    def sample_const(self):
        consts = self.get_const()
        errs   = self.get_err()
        new_consts = []
        for const, err in zip(consts,errs):
            new_consts.append(np.random.normal(const,err))
        return new_consts
    
    def mechanism(self,temperature=None,stochastic=False,x_fe=0.1,**kwargs):
        if stochastic:
            sigma, b, c, enthalpy = self.sample_const()
        else:
            sigma, b, c, enthalpy= self.get_const()
        
        sigma = np.power(10,sigma+b*x_fe)
        return sigma*np.exp(-enthalpy/(self.r*temperature))
 
def plot_probabilty_graph(obj,pressure=4,log10_fO2=-2.0,x_fe=0.1):
    temperature_array  = np.linspace(200,2000,num=500) + 273.15 # temp array in kelvin
    conductivity_array = np.linspace(-9,3,num=200)
    probability_array = np.zeros((200,500))
    plt.figure(figsize=(10,10))
    for i in range(3000):
        conduct = obj.conductivity(temperature=temperature_array,pressure=pressure, stochastic=True,log10_fO2=log10_fO2,x_fe=x_fe,with_hydrous=True)
        valid_conduct = conduct[conduct>0]
        valid_temps   = temperature_array[conduct>0]
        for i in range(len(valid_conduct)):
            if np.log10(valid_conduct[i])<-9:
                continue
            closest_conductivity_index = np.argmin(np.abs(conductivity_array-np.log10(valid_conduct[i])))
            closest_temperature_index  = np.argmin(np.abs(temperature_array-valid_temps[i]))
            probability_array[closest_conductivity_index,closest_temperature_index]+=1

    probability_avg = np.max(probability_array,axis=0)
    probability_array/=probability_avg
    cax=plt.pcolormesh(conductivity_array,temperature_array,probability_array.T,cmap='nipy_spectral')
    plt.colorbar(cax,label='probabilty',ticks=np.linspace(0,1,num=11))
    conduct = obj.conductivity(temperature=temperature_array,pressure=pressure, stochastic=False,log10_fO2=log10_fO2,x_fe=x_fe,with_hydrous=True)
    plt.plot(np.log10(conduct),temperature_array,label=obj.name,color='white',linestyle='--',linewidth=2)
    plt.gca().invert_yaxis()
    plt.title(obj.name)
    plt.xlim([-9,3])
    plt.legend(loc='lower left',ncol=4)
    plt.xlabel(r'Conductivity $log_{10}$(s/m)')
    plt.ylabel('temperature (K)')
    plt.show()