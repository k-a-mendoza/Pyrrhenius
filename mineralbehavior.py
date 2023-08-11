import numpy as np
from dataclasses import dataclass

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


class Mineral:

    
    def __init__(self,name=None,author=None,title=None,year=None,min_t=-np.inf,max_t=np.inf,min_p=-np.inf,max_p=np.inf,hydrous=False,anisotropic=False,iron_frac=False,fO2_dependent=False,mineral_name=None,**kwargs):
        self.name=name
        self.author=author
        self.title=title
        self.year=year
        self.min_t=min_t
        self.max_t=max_t
        self.min_p=min_p
        self.max_p=max_p
        self.hydrous=hydrous
        self.anisotropic=anisotropic
        self.iron_frac = iron_frac
        self.fO2_dependent = fO2_dependent
        self.mineral_name = mineral_name
    
    
    def isotropic_conductivity(self,extrapolate=False,**kwargs):
        """
        calculates the isotropic conductivity of the mineral. 
        
        Uses a geometric average across the three crystallographic directions to calculate the norm. 
        
        will return np.nan values if calculation is done out of range of the data
        
        """
        xstal_001, xstal_010, xstal_100 = self.get_3_xstal_conductivities(**kwargs)
        isotropic_c = np.power(xstal_001*xstal_010*xstal_100,1/3)
        if np.all(np.isnan(isotropic_c)):
            print(f"{self.name} is all nans")
        if not extrapolate:
            isotropic_c = self.apply_nanmask(isotropic_c,**kwargs)
        return isotropic_c
    
    def conductivity_001(self,*args,**kwargs):
        raise Exception('001 conductivity not implemented')
        
    def conductivity_010(self,*args,**kwargs):
        raise Exception('010 conductivity not implemented')
        
    def conductivity_100(self,*args,**kwargs):
        raise Exception('100 conductivity not implemented')
        
    def max_anisotropic_conductivity(self, extrapolate=False,**kwargs):
        """
        calculates the maximum conductivity of the mineral if aligned.
        
        
        will return np.nan values if calculation is done out of range of the data
        
        """
        xstal_001, xstal_010, xstal_100 = self.get_3_xstal_conductivities(**kwargs)
        conductivity = np.stack([xstal_001,xstal_010,xstal_100],axis=-1)
        conductivity = np.max(conductivity,axis=-1)
        if not extrapolate:
            conductivity = self.apply_nanmask(conductivity,**kwargs)
        return conductivity
    
    def get_3_xstal_conductivities(self,**kwargs):
        xstal_001 = self.conductivity_001(**kwargs)
        xstal_010 = self.conductivity_010(**kwargs)
        xstal_100 = self.conductivity_100(**kwargs)
        return xstal_001, xstal_010, xstal_100
        
    def min_anisotropic_conductivity(self,extrapolate=False,**kwargs):
        """
        calculates the maximum conductivity of the mineral if aligned.
        
        
        will return np.nan values if calculation is done out of range of the data
        
        """
        xstal_001, xstal_010, xstal_100 = self.get_3_xstal_conductivities(**kwargs)
        
        conductivity = np.stack([xstal_001,xstal_010,xstal_100],axis=-1)
        conductivity = np.min(conductivity,axis=-1)
        if not extrapolate:
            conductivity = self.apply_nanmask(conductivity,**kwargs)
        return conductivity
    
    def polycrystalline_conductivity(self,extrapolate=False,**kwargs):
        """
        calculates the conductivity of the mineral based on polycrystalline data
        
        
        will return np.nan values if calculation is done out of range of the data
        
        """
        conductivity = self._polycrystalline_conductivity(**kwargs)
        if not extrapolate:
            conductivity = self.apply_nanmask(conductivity,**kwargs)
        return conductivity
    
    def apply_nanmask(self,conductivity,temperature=None,pressure=None,**kwargs):
        if isinstance(temperature,np.ndarray):
            mask = np.ones(temperature.shape)
        elif isinstance(pressure,np.ndarray):
            mask = np.ones(temperature.shape)
        else:
            if temperature < self.min_t or temperature < self.max_t or \
               pressure < self.min_p or pressure > self.max_p:
                return conductivity
            else:
                return np.nan
            
        mask[(temperature < self.min_t) | (temperature > self.max_t) | \
             (pressure <self.min_p) | (pressure >self.max_p)] = np.nan
        return mask*conductivity
            
    
class XStalAvg:
    """
    geometrically averages conduction mechanisms. 
    
    """
    
    def __init__(self,mechanisms,**kwargs):
        self.mechanisms=mechanisms
    
    
    def get_conductivity(self,**kwargs):
        n_mechanisms = len(self.mechanisms)
        conductivity_0 = self.mechanisms[0].get_conductivity(**kwargs)
        for index in range(1,n_mechanisms):
            conductivity*=self.mechanisms[index]
        return np.power(conductivity,1/n_mechanisms)
    
class XStalTypical:
    k = 8.617e-5
    r = 8.314e-3
    norm_v = 1/96.49
    
    def __init__(self,preexp=None,preexp_err = 0,enthalpy=None, enthalpy_err=0,activation_volume=0,activation_volume_err=0,
                 log_sigma_err=False,convert2ev=False,**kwargs):
        self.preexp=preexp
        self.preexp_err = preexp_err
        if convert2ev:
            enthalpy*=0.01036
            enthalpy_err*=0.01036
        self.enthalpy=enthalpy
        self.enthalpy_err = enthalpy_err
        self.activation_volume = activation_volume
        self.activation_volume_err = activation_volume_err
        self.log_sigma_err=log_sigma_err
    
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,scale=self.activation_volume_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
        if self.log_sigma_err:
            sigma = np.power(10,sigma)
        return sigma, enthalpy, activation_volume
    
    def get_conductivity(self,temperature=None,pressure=1,**kwargs):
        sigma, enthalpy, act_v = self.sample_constants(**kwargs)
        return sigma*np.exp(-(enthalpy + act_v*pressure*self.norm_v)/(self.k*temperature))
    
    
class XStalWaterSimple(XStalTypical):

    
    def __init__(self,preexp_c_constant=1.0,preexp_c_constant_err=0,
                 water_content_type='ppm',**kwargs):
        super().__init__(**kwargs)
        self.preexp_c_constant = preexp_c_constant
        self.preexp_c_constant_err = preexp_c_constant_err
        self.water_content_type = water_content_type
    
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,scale=self.activation_volume_err)
            exp1  = np.random.normal(loc=self.preexp_c_constant,scale=self.preexp_c_constant_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
            exp1 = self.preexp_c_constant
        if self.log_sigma_err:
            sigma = 10**sigma
        return sigma, enthalpy, activation_volume, exp1
    
    def convert_wtpcnts(self,wtppm_H2O):
        water_deliverable = wtppm_H2O
        if self.water_content_type=='frac':
            water_deliverable = wtppm_H2O/1e6
        elif self.water_content_type=='wt%':
            water_deliverable = wtppm_H2O/1e4
        return water_deliverable
    
    def get_conductivity(self,temperature=None, pressure=1, wtppm_H2O=5,**kwargs):
        sigma, enthalpy, act_v, exp1 = self.sample_constants(**kwargs)
        h2o = self.convert_wtpcnts(wtppm_H2O)
        return sigma*np.power(h2o,exp1)*np.exp(-(enthalpy + act_v*pressure*self.norm_v)/(self.k*temperature))
              
class OldEmpiricalLogEqn:
    k = 8.617e-5
    r = 8.314e-3
    norm_v = 1/96.49
    
    def __init__(self,offset,t_const,fo2_const,error=0.10):
        self.offset = offset 
        self.t_const = t_const
        self.fo2_const = fo2_const
        self.error = error
        
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            offset    = np.random.normal(loc=self.offset, scale=abs(self.offset*self.error))
            t_const   = np.random.normal(loc=self.t_const, scale=abs(self.t_const*self.error))
            fo2_const = np.random.normal(loc=self.fo2_const, scale=abs(self.fo2_const*self.error))
        else:
            offset    = self.offset
            t_const   = self.t_const
            fo2_const = self.fo2_const
        
        return offset, t_const, fo2_const
    
    def get_conductivity(self,temperature=None, pressure=1,log10_fO2=None,
                         log10_fO2_offset=0,**kwargs):
        offset, t_const, fo2_const = self.sample_constants(**kwargs)
        
        if log10_fO2 is None:
            log10_fO2 = calc_QFM(temperature,pressure) + log10_fO2_offset

        result = offset + t_const/temperature + fo2_const*log10_fO2 
        return 10**result
                            
class XStalWaterEmpirical(XStalWaterSimple):
    
    def __init__(self,enthalpy_multiplier=-1,enthalpy_multiplier_err=0,
                 enthalpy_exp=1/3,enthalpy_exp_err=0,**kwargs):
        super().__init__(**kwargs)
        self.enthalpy_multiplier = enthalpy_multiplier
        self.enthalpy_multiplier_err = enthalpy_multiplier_err
        self.enthalpy_exp=enthalpy_exp
        self.enthalpy_exp_err=enthalpy_exp_err
        
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,scale=self.activation_volume_err)
            exp1  = np.random.normal(loc=self.preexp_c_constant,scale=self.preexp_c_constant_err)
            exp2  = np.random.normal(loc=self.enthalpy_exp,scale=self.enthalpy_exp_err)
            a  = np.random.normal(loc=self.enthalpy_multiplier,scale=self.enthalpy_multiplier_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
            exp1 = self.preexp_c_constant
            exp2 = self.enthalpy_exp
            a    = self.enthalpy_multiplier
        if self.log_sigma_err:
            sigma = np.power(10,sigma)
        return sigma, enthalpy, activation_volume, exp1, exp2, a
    
    def get_conductivity(self,temperature=None, pressure=1, wtppm_H2O=5,**kwargs):
        sigma, enthalpy, act_v, exp1, exp2, a = self.sample_constants(**kwargs)
        h20 = self.convert_wtpcnts(wtppm_H2O)
        enthalpy_prime = enthalpy + a*np.power(h20,exp2) + act_v*pressure*self.norm_v
        return sigma*np.power(h20,exp1)*np.exp(-enthalpy_prime/(self.k*temperature))
    
class XStalWaterIronEmpirical(XStalWaterSimple):
    
    def __init__(self,alpha=0,alpha_err=1,beta=0,beta_err=1,enthalpy_exp=1,enthalpy_exp_err=1,
                 iron_exp=1,iron_exp_err=1,iron_offset=1e-4,**kwargs):
        super().__init__(**kwargs)
        self.alpha = alpha; self.alpha_err = alpha_err
        self.beta = beta; self.beta_err = beta_err
        self.enthalpy_exp = enthalpy_exp; self.enthalpy_exp_err=enthalpy_exp_err
        self.iron_exp = iron_exp; self.iron_exp_err = iron_exp_err
        self.iron_offset = iron_offset
        
        
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma     = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy  = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,scale=self.activation_volume_err)
            exp1         = np.random.normal(loc=self.preexp_c_constant,scale=self.preexp_c_constant_err)
            iron_exp     = np.random.normal(loc=self.iron_exp,scale=self.iron_exp_err)
            enthalpy_exp = np.random.normal(loc=self.enthalpy_exp,scale=self.enthalpy_exp_err)
            alpha = np.random.normal(loc=self.alpha,scale=self.alpha_err)
            beta  = np.random.normal(loc=self.beta,scale=self.beta_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
            exp1 = self.preexp_c_constant
            iron_exp     = self.iron_exp
            enthalpy_exp = self.enthalpy_exp
            alpha = self.alpha
            beta  = self.beta
        if self.log_sigma_err:
            sigma = np.power(10,sigma)
        return sigma, enthalpy, activation_volume, exp1, iron_exp, enthalpy_exp, alpha, beta
    
    def get_conductivity(self,temperature=None, pressure=1, wtppm_H2O=5,x_fe=0.5,**kwargs):
        sigma, enthalpy, act_v, exp1, iron_exp, enthalpy_exp, alpha, beta = self.sample_constants(**kwargs)
        h2O = self.convert_wtpcnts(wtppm_H2O)
        sigma_prime    = sigma*np.power(x_fe+self.iron_offset,iron_exp)*np.power(h2O,exp1)
        enthalpy_prime = enthalpy + beta*np.power(h2O,enthalpy_exp)+alpha*x_fe + act_v*pressure*self.norm_v
        return sigma_prime*np.exp(-enthalpy_prime/(self.k*temperature))
                                                     

    
class XStalWaterSimpleIronExp(XStalWaterSimple):

    
    def __init__(self,iron_frac=0,iron_frac_err=0,convert2ev=False,**kwargs):
        super().__init__(convert2ev=convert2ev,**kwargs)
        if convert2ev:
            iron_frac*=0.01036
            iron_frac_err*=0.01036
        self.iron_frac     = iron_frac
        self.iron_frac_err = iron_frac_err
    
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,\
                                                 scale=self.activation_volume_err)
            exp1  = np.random.normal(loc=self.preexp_c_constant,\
                                     scale=self.preexp_c_constant_err)
            alpha = np.random.normal(loc=self.iron_frac,scale=self.iron_frac_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
            exp1  = self.preexp_c_constant
            alpha = self.iron_frac
          
        if self.log_sigma_err:
            sigma = np.power(10,sigma)
        return sigma, enthalpy, activation_volume, exp1, alpha
    
    def get_conductivity(self,temperature=None, pressure=1, wtppm_H2O=5,x_fe=1.0,**kwargs):
        sigma, enthalpy, act_v, exp1, alpha = self.sample_constants(**kwargs)
        h2o = self.convert_wtpcnts(wtppm_H2O)
       
        new_enthalpy = enthalpy + act_v*pressure*self.norm_v + x_fe*alpha
        
        poww =  np.power(h2o,exp1)
        
        exp  = np.exp(-(new_enthalpy)/(self.k*temperature))
        return sigma*poww*exp
    
class XStalWaterInvT(XStalWaterSimple):
    
    def __init__(self,preexp_c_constant=1.0,preexp_c_constant_err=0,calibrate=1,**kwargs):
        super().__init__(**kwargs)
        self.calibrate=calibrate
        
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,
                                                 scale=self.activation_volume_err)
            exp1  = np.random.normal(loc=self.preexp_c_constant,scale=self.preexp_c_constant_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
            exp1 = self.preexp_c_constant
        if self.log_sigma_err:
            sigma = np.power(10,sigma)
        return sigma, enthalpy, activation_volume, exp1
    
    def get_conductivity(self,temperature=None, pressure=1, wtppm_H2O=5,**kwargs):
        sigma, enthalpy, act_v, exp1 = self.sample_constants(**kwargs)
        h2o = self.convert_wtpcnts(wtppm_H2O)
        enthalpy_prime = enthalpy + act_v*pressure*self.norm_v
        sigma_prime  = sigma/temperature
        return sigma_prime*np.power(h2o*self.calibrate,exp1)*np.exp(-enthalpy_prime/(self.k*temperature))

class XStalfO2Premult(XStalTypical):
  

    def __init__(self,fO2_exp=1/6,fO2_exp_err=0,**kwargs):
        super().__init__(**kwargs)
        self.fO2_exp = fO2_exp
        self.fO2_exp_err = fO2_exp_err
    
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,scale=self.activation_volume_err)
            fO2_exp = np.random.normal(loc=self.fO2_exp,scale=self.fO2_exp_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
            fO2_exp = self.fO2_exp
        if self.log_sigma_err:
            sigma = np.power(10,sigma)
        return sigma, enthalpy, activation_volume, fO2_exp
    
    def get_conductivity(self,temperature=None,pressure=1,log10_fO2=None,log10_fO2_offset=0,**kwargs):
        sigma, enthalpy, act_v, s = self.sample_constants(**kwargs)
        if log10_fO2 is None:
            log10_fO2 = calc_QFM(temperature,pressure) + log10_fO2_offset
        fO2 = 10**log10_fO2
        
        return np.power(fO2,s)*sigma*np.exp(-(enthalpy + act_v*pressure*self.norm_v)/(self.k*temperature))
    
class XStalXfeComplex(XStalTypical):
  

    def __init__(self,a=0,b=0,a_err=0,b_err=0,**kwargs):
        super().__init__(**kwargs)
        self.a = a
        self.a_err = a_err
        self.b = b
        self.b_err = b_err
    
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy          = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            activation_volume = np.random.normal(loc=self.activation_volume,
                                                 scale=self.activation_volume_err)
            
            a = np.random.normal(loc=self.a,scale=self.a_err)
            b = np.random.normal(loc=self.b,scale=self.b_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
            activation_volume = self.activation_volume
            a = self.a
            b = self.b
        if self.log_sigma_err:
            sigma = np.power(10,sigma)
        return sigma, enthalpy, activation_volume, a, b
    
    def get_conductivity(self,temperature=None,pressure=1,x_fe=0.3,**kwargs):
        sigma, enthalpy, act_v, a, b = self.sample_constants(**kwargs)
        
        enthalpy_modified     = enthalpy-a*x_fe**(1/3)
        activation_v_modified = pressure*(act_v - b*x_fe)*self.norm_v
        
        c = sigma*x_fe*np.exp(-( enthalpy_modified + activation_v_modified)/(self.k*temperature))
        
        return c
    
class DiffusionHConduction(XStalWaterSimple):
    f=1
    q = 1.60217663e-19
    kJ = 1.381e-23
    avos = 6.0221408e23
    def __init__(self,density_correction_factor=1,**kwargs):
        super().__init__(**kwargs)
        self.dcf = density_correction_factor
        
        
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
        if self.log_sigma_err:
            sigma =10**sigma
        return sigma, enthalpy
    
    
    def get_conductivity(self,temperature=None,wtppm_H2O=5,**kwargs):
        sigma, enthalpy = self.sample_constants(**kwargs)
        diffusivity = 10**(np.log10(sigma) - enthalpy/(np.log(10)*self.r*temperature))
        h2o = self.convert_wtpcnts(wtppm_H2O)
        nerst_einstein = self.avos*self.dcf*h2o*self.q*self.q/self.kJ
        proton_conduction = nerst_einstein*diffusivity/temperature
        return proton_conduction
    
class DiffusionHConductionRexp(DiffusionHConduction):

    def __init__(self,exponent_constant_c=1.0,exponent_constant_err=0,**kwargs):
        super().__init__(**kwargs)
        self.c = exponent_constant_c
        self.c_err = exponent_constant_c_err
        
        
    def sample_constants(self,stochastic=False,**kwargs):
        if stochastic:
            sigma = np.random.normal(loc=self.preexp, scale=self.preexp_err)
            enthalpy = np.random.normal(loc=self.enthalpy,scale=self.enthalpy_err)
            exponent = np.random.normal(loc=self.c,scale=self.c_err)
        else:
            sigma = self.preexp
            enthalpy = self.enthalpy
        if self.log_sigma_err:
            sigma =10**sigma
        return sigma, enthalpy, exponent
    
    
    def get_conductivity(self,temperature=None,wtppm_H2O=5,**kwargs):
        sigma, enthalpy, exponent = self.sample_constants(**kwargs)
        h2o = self.convert_wtpcnts(wtppm_H2O)
        diffusivity = sigma*np.power(wtppm_H2O,exponent)*np.exp(-enthalpy/(self.k*temperature))
        nerst_einstein = self.dcf*h2o*self.q*self.q/self.kJ
        proton_conduction = nerst_einstein*diffusivity/temperature
        return proton_conduction