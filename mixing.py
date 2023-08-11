import numpy as np
def sum_conductivities_matrix(estimate, sigma, frac, pol):
    summation = np.zeros(estimate.shape)
    for i in range(sigma.shape[-1]):
        f_sig = sigma[:,:,i]
        numerator = (f_sig -estimate)
        denominator = (f_sig +pol*estimate)*frac[:,:,i]
        summation+=numerator/denominator
    return summation

def sum_conductivities_array(estimate, sigma, frac, pol):
    summation = np.zeros(len(estimate))
    for i in range(sigma.shape[-1]):
        f_sig = sigma[:,i]
        summation+=(f_sig -estimate)/(f_sig +pol*estimate)*frac[:,i]
    return summation

def sum_conductivities_scalar(estimate, sigma, frac, pol):
    summation = 0
    for i in range(len(sigma)):
        f_sig = sigma[i]
        summation+=(f_sig -estimate)/(f_sig +pol*estimate)*frac[i]
    return summation


class EffectiveMedium1D:
    
    def __init__(self):
        pass
    
    
    def estimate_conductivity(self,constituent_conductivities,fractions,polarization=1):
        """
        conductivities as list, fractions as list, polarization should be =1 if sphere, 0 if rod
        
        """
        p_multiplier = ((1/((1/3)*polarization))-1)
        if isinstance(constituent_conductivities,np.ndarray):
            conductivity = binary_search_array(constituent_conductivities, fractions, polarization)
        else:
            conductivity=binary_search_value(constituent_conductivities, fractions, polarization)
                
        return conductivity
    
class EffectiveMedium2D:
    
    def __init__(self):
        pass
    
    
    def estimate_conductivity(self,constituent_conductivities,fractions,polarization=1):
        """
        conductivities as list, fractions as list, polarization should be =1 if sphere, 0 if rod
        
        """
        p_multiplier = ((1/((1/3)*polarization))-1)
        conductivity = binary_search_matrix(constituent_conductivities, fractions, polarization)

        return conductivity
    
    
class Tube:
    
    def __init__(self):
        pass
    
    def estimate_conductivity(self,fluid_conductivity,bulk_conductivity,fraction_fluid):
        return 0.3*fraction_fluid*fluid_conductivity + (1-fraction_fluid)*bulk_conductivity
    
class Film:
    
    def __init__(self):
        pass
    
    def estimate_conductivity(self,fluid_conductivity,bulk_conductivity,fraction_fluid):
        inv_melt_frac = 1-fraction_fluid
        numerator = (np.power(inv_melt_frac,2/3)*fluid_conductivity - np.power(inv_melt_frac,1/3)*bulk_conductivity)
        divisor_melt    = (inv_melt_frac - np.power(inv_melt_frac,2/3))*fluid_conductivity 
        divisor_solid   = (np.power(inv_melt_frac,2/3)-inv_melt_frac -1)*bulk_conductivity
        
        return fluid_conductivity*(numerator)/(divisor_melt+divisor_solid)
        
class MeanOfFilmTube:
    
    def __init__(self):
        self.tube = Tube()
        self.film = Film()
        
    def estimate_conductivity(self,fluid_conductivity,bulk_conductivity,frac_melt):
        tube_c= self.tube.estimate_conductivity(fluid_conductivity,bulk_conductivity,frac_melt)
        cube_c= self.tube.estimate_conductivity(fluid_conductivity,bulk_conductivity,frac_melt)
        return (tube_c + cube_c)/2

    
    
class HashinStrickland:
    """
    berryman for aggregates
    
    """
    def __init__(self):
        pass
    
    
    def max_bound(self,constituent_conductivities,fractions):
        if isinstance(constituent_conductivities,np.ndarray):
            sigma_max = np.max(constituent_conductivities,axis=1)
            conductivity = np.zeros(fractions.shape[0])
            
            for i in range(fractions.shape[1]):
                divisor = constituent_conductivities[:,i]+2*sigma_max
    
                conductivity += 1/(fractions[:,i]/divisor)
        else:
            conductivity=0
            sigma_max = np.max(constituent_conductivities)
            for i in range(len(fractions)):
                conductivity += 1/(fractions[i]/(constituent_conductivities[i]+2*sigma_max)) 
        conductivity - 2*sigma_max
        return conductivity
                
    def min_bound(self,constituent_conductivities,fractions):
        if isinstance(constituent_conductivities,np.ndarray):
            sigma_min = np.nanmin(constituent_conductivities,axis=1)
            conductivity = np.zeros(fractions.shape[0])
            for i in range(fractions.shape[1]):
                divisor = constituent_conductivities[:,i]+2*sigma_min
                
                conductivity += 1/(fractions[:,i]/divisor)
        else:
            conductivity=0
            sigma_min = np.nanmin(constituent_conductivities)
            for i in range(len(fractions)):
                conductivity += 1/(fractions[i]/(constituent_conductivities[i]+2*sigma_min)) 
        conductivity-=2*sigma_min
        return conductivity
        
        
def binary_search_array(constituent_conductivities, fractions, polarization,tolerance=1e-15):
    broadcast_vector = np.ones(constituent_conductivities.shape)
    lower_bound   = np.min(constituent_conductivities,axis=1)
    upper_bound   = np.max(constituent_conductivities,axis=1)
    mb, sign, res = _evaluate_bounds(constituent_conductivities, fractions, polarization,upper_bound,lower_bound,sum_conductivities_array)
    
    lower_prediction =  sum_conductivities(lower_bound, constituent_conductivities, fractions, polarization)
    upper_prediction = sum_conductivities(upper_bound, constituent_conductivities, fractions, polarization)
    iterations=0
    maxit=1e4
    while np.any(tolerance<abs(res)) and iterations<maxit:
        mb, sign, res = _evaluate_bounds(constituent_conductivities, fractions, polarization,upper_bound,lower_bound,sum_conductivities_array)
        prediction_sign = -np.sign(sign)
        lower_bound, upper_bound = _adjust_bounds(upper_bound, lower_bound,
                                                  mb, prediction_sign)
        iterations+=1
    return mb

def binary_search_matrix(constituent_conductivities, fractions, polarization,tolerance=1e-15):
    broadcast_vector = np.ones(constituent_conductivities.shape)
    lower_bound   = np.min(constituent_conductivities,axis=2)
    upper_bound   = np.max(constituent_conductivities,axis=2)
    mb, sign, res = _evaluate_bounds(constituent_conductivities, fractions, polarization,upper_bound,lower_bound,sum_conductivities_matrix)
    
    lower_prediction = sum_conductivities_matrix(lower_bound, constituent_conductivities, fractions, polarization)
    upper_prediction = sum_conductivities_matrix(upper_bound, constituent_conductivities, fractions, polarization)
    invalid_mask = (lower_prediction > 0) | (upper_prediction < 0) | np.isnan(lower_prediction) | np.isnan(upper_prediction)
    iterations=0
    maxit=1e4
    while np.any(tolerance<res[~invalid_mask]) and iterations<max_iterations:
        mb, sign, res = _evaluate_bounds(constituent_conductivities, fractions, polarization,upper_bound,lower_bound,sum_conductivities_matrix)
        prediction_sign = -np.sign(sign)
        lower_bound, upper_bound = _adjust_bounds(upper_bound, lower_bound,
                                                  mb, prediction_sign)
        iterations+=1
    return mb

def binary_search_value(constituent_conductivities, fractions, polarization,tolerance=1e-15):
    lower_bound   = np.min(constituent_conductivities)
    upper_bound   = np.max(constituent_conductivities)
    summation = sum_conductivities_scalar(constituent_conductivities, fractions, polarization)
    mb = (lower_bound+upper_bound)/2
    res = abs(summation)
    
    lower_prediction = sum_conductivities_scalar(lower_bound, fractions, polarization)
    upper_prediction = sum_conductivities_scalar(upper_bound, fractions, polarization)
    iterations=0
    while np.any(tolerance<res):
        summation = sum_conductivities_scalar(constituent_conductivities, fractions, polarization)
        mb = (lower_bound+upper_bound)/2
        res = abs(summation)
        print(res[0])
        prediction_sign = np.sign(summation)
        lower_bound, upper_bound = _adjust_bounds(upper_bound, lower_bound,
                                                  mb, prediction_sign)
        iterations+=1
    return mb
        
def _evaluate_bounds(constituent_conductivities, fractions, polarization,upper_bound,lower_bound,func):
    bound_average = np.mean(np.stack([lower_bound,upper_bound],axis=-1),axis=-1)
        
    middle_prediction = func(bound_average, constituent_conductivities, fractions, polarization)
    residual      = np.abs(middle_prediction)
    return bound_average, np.sign(middle_prediction), middle_prediction

def _adjust_bounds(upper_bound, lower_bound, bound_average, prediction_sign):
    if 1 in prediction_sign:
        upper_bound[prediction_sign>0] = bound_average[prediction_sign>0]
    if -1 in prediction_sign:
        lower_bound[prediction_sign<0] = bound_average[prediction_sign<0]
    return lower_bound, upper_bound