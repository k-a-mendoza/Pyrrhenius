import numpy as np
from .model import WaterPTCorrection
from scipy.optimize import minimize, root_scalar, root
from collections import OrderedDict



def binary_search(experiment, target_values, bounds, kwargs,
                  tolerance=1e-3, relative=True,debug=True,max_iterations = 50):
    broadcast_vector = np.ones(target_values.shape)
    lower_bound = bounds[0] * broadcast_vector
    upper_bound = bounds[1] * broadcast_vector
    bound_average = lower_bound+upper_bound
    bound_average/=2

    upper_prediction = experiment(upper_bound, **kwargs) - target_values
    lower_prediction = experiment(lower_bound, **kwargs) - target_values
    middle_prediction = experiment(bound_average, **kwargs) - target_values
    new_nans = np.isnan(lower_prediction*upper_prediction*middle_prediction)
    original_bounded_regions = lower_prediction*upper_prediction<0
    valid_mask = original_bounded_regions & ~new_nans
    valid_points = 100 * np.count_nonzero(valid_mask) / len(target_values.ravel())
    if debug:
        print(f'Starting Solver:\nIteration NaNs: {100 * np.count_nonzero(new_nans) / len(target_values.ravel())} %'+\
        f' In-bounds: {100 * np.count_nonzero(original_bounded_regions) / len(target_values.ravel())} %'+\
        f' Valid: {valid_points} %')

    iterations = 0

    def relative_cost(middle_residual,original_field):
        return abs(middle_residual) / original_field

    def absolute_cost(middle_residual,original_field):
        return abs(middle_residual)

    if relative:
        cost_function = relative_cost
    else:
        cost_function = absolute_cost

    while np.any(tolerance < cost_function(middle_prediction,target_values)[valid_mask]) and iterations < max_iterations:
        bound_average = lower_bound+upper_bound
        bound_average/=2
        upper_prediction  = experiment(upper_bound, **kwargs)   - target_values
        lower_prediction  = experiment(lower_bound, **kwargs)   - target_values
        middle_prediction = experiment(bound_average, **kwargs) - target_values
        
       
        change_mask = upper_prediction * middle_prediction < 0
 
        lower_bound[change_mask] = bound_average[change_mask]
        upper_bound[~change_mask] = bound_average[~change_mask]
        

        new_nans = np.isnan(lower_prediction*upper_prediction*middle_prediction)
        bounded_regions = lower_prediction*upper_prediction<0
        valid_mask   = bounded_regions & ~new_nans
        valid_points = 100 * np.count_nonzero(valid_mask) / len(target_values.ravel())
        if debug:
            print(f'At {iterations}.  Max Relative Residual: {np.nanmax(relative_cost(middle_prediction,target_values)[valid_mask]):2.2e}')
            print(f' Max Absolute Residual: {np.nanmax(absolute_cost(middle_prediction,target_values)[valid_mask]):2.2e}')
            print(f' Iteration NaNs: {100 * np.count_nonzero(new_nans) / len(target_values.ravel())} %')
            print(f' In-bounds: {100 * np.count_nonzero(bounded_regions) / len(target_values.ravel())} %')
            print(f' Valid: {valid_points} %')
        iterations += 1
    if debug:
        print(f"{iterations} iterations to reach {tolerance}")
    bound_average[~original_bounded_regions]=np.nan
    return bound_average


def _find_new_bounds(bound_average, lower_bound, lower_prediction, middle_prediction, upper_bound, upper_prediction):
    product_upper = upper_prediction * middle_prediction < 0
    product_lower = lower_prediction * middle_prediction < 0
    if product_upper.any():
        lower_bound[product_upper] = bound_average[product_upper]
    if product_lower.any():
        upper_bound[product_lower] = bound_average[product_lower]
    bound_average = (lower_bound + upper_bound) / 2
    return bound_average, lower_bound, upper_bound


def create_starting_variable(P=None, T=None, Cw=None, X_fe=None,
                             logfo2=None, co2=None,start_value =0,**kwargs):
    kwargs = [P, T, Cw, X_fe, logfo2, co2]
    if any(isinstance(arg, np.ndarray) for arg in kwargs):
        arrays = filter(lambda x : isinstance(x,np.ndarray),kwargs)
        shapes = [x.shape for x in arrays]
        if len(set(shapes)) > 0:
            assert all(x == shapes[0] for x in shapes), f"Some passed shapes are different. {shapes}"
            return np.ones(shapes[0])*start_value
    else:
        start_value = np.ones((1))*start_value

    return start_value

def access_array_indices(P=None, T=None, Cw=None, X_fe=None, logfo2=None, co2=None,**kwargs):
    kwargs = [P, T, Cw, X_fe, logfo2, co2]
    if any(isinstance(arg, np.ndarray) for arg in kwargs):
        arrays = filter(lambda x : isinstance(x,np.ndarray),kwargs)
        shapes = [x.shape for x in arrays]
        if len(set(shapes)) > 0:
            assert all(x == shapes[0] for x in shapes), f"Some passed shapes are different. {shapes}"
            return np.zeros(shapes[0])

    return 0


class BruggemanSolver:
    def __init__(self, other, **kwargs):
        self.experiment = other
        self.constituent_conductivities = [x.get_conductivity(**kwargs) for x in other.phase_list]
        if isinstance(self.constituent_conductivities[0],float):
            self.baseline_img = 0.0
        else:
            self.baseline_img = np.zeros(self.constituent_conductivities[0].shape,dtype=float)

    def __call__(self, bulk):
        baseline_img =  np.copy(self.baseline_img)
        for phase_index in range(len(self.constituent_conductivities)):

            fraction     = self.experiment.phase_fractions[phase_index]
            polarization = self.experiment.polarization_factors[phase_index]
            conductivity_slice = self.constituent_conductivities[phase_index]

            numerator = conductivity_slice - bulk
            denominator = conductivity_slice + (1 / polarization - 1) * bulk
            new_value = fraction * numerator / denominator
            baseline_img += new_value

        return baseline_img

    def __call__single(self, bulk,**kwargs):
        self.baseline_img = create_starting_variable(**kwargs)
        for phase_index in range(len(self.constituent_conductivities)):
            fraction = self.experiment.phase_fractions[phase_index]

            polarization = self.experiment.polarization_factors[phase_index]
            conductivity_slice = self.constituent_conductivities[phase_index]
            numerator = conductivity_slice - bulk
            denominator = conductivity_slice + (1 / polarization - 1) * bulk

            self.baseline_img += fraction * numerator / denominator

        return self.baseline_img


class BruggemanSymmetrical:

    mixing_model='EMT1'
    def __init__(self,phase_list=None,phase_fractions=None,polarization_factors=None,**kwargs):
        if polarization_factors is None:
            pol_fac = [1/3 for x in phase_list]
        else:
            pol_fac = []
            for x in polarization_factors:
                if x is None or x > 1/3 or x <0:
                    pol_fac.append(1/3)
                else:
                    pol_fac.append(x)

        if phase_fractions is None:
            phase_fractions = [1/len(phase_list)]*len(phase_list)
        self.phase_list           = phase_list
        self.phase_fractions      = phase_fractions
        self.polarization_factors = pol_fac

    def get_conductivity(self,phase_relevant_kwargs=None,**kwargs):
        if phase_relevant_kwargs is None:
            constituent_conductivities = np.stack([x.get_conductivity(**kwargs) \
                                                   for x in self.phase_list], axis=-1)
        else:
            constituent_conductivities = np.stack([x.get_conductivity(**{**kwargs,**phase_kwargs}) \
                                                   for x,phase_kwargs in zip(self.phase_list,\
                                                                        phase_relevant_kwargs)], axis=-1)

        bounds_min = np.min(constituent_conductivities,axis=-1)
        bounds_max = np.max(constituent_conductivities,axis=-1)
        minibruggeman = BruggemanSolver(self, **kwargs)
        if constituent_conductivities.ndim==1:
            result = root_scalar(minibruggeman,x0 = np.zeros(bounds_max.shape),
                                 bracket=(bounds_min,bounds_max)).root
        else:
            result = binary_search(minibruggeman,np.zeros(bounds_max.shape),
                      bounds=(bounds_min, bounds_max),
                      kwargs={},
                      tolerance=1e-8, relative=True,debug=False)

        return result




class HashinShtrikmanUpper:

    mixing_model = 'HS+'
    def __init__(self,matrix,inclusion,**kwargs):
        self.matrix = matrix
        self.inclusion=inclusion


    def get_conductivity(self,phi,**kwargs):

        c_matrix = self.matrix.get_conductivity(**kwargs)
        c_inclusion = self.inclusion.get_conductivity(**kwargs)

        denominator = 1/(c_matrix - c_inclusion) + phi/(3*c_inclusion)
        return c_inclusion + (1-phi)/denominator

    @property
    def metadata(self):
        return self.inclusion.metadata+self.matrix.metadata


class HashinShtrikmanLower:

    mixing_model = 'HS-'
    def __init__(self,matrix,inclusion,**kwargs):
        self.matrix = matrix
        self.inclusion=inclusion


    def get_conductivity(self,phi,**kwargs):

        c_matrix = self.matrix.get_conductivity(**kwargs)
        c_inclusion = self.inclusion.get_conductivity(**kwargs)

        denominator = 1/(c_inclusion - c_matrix) + (1-phi)/(3*c_matrix)
        return c_matrix + phi/denominator

class ParallelModel:
    """
        https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GC006530
    """
    def __init__(self, matrix, inclusion):
        self.matrix = matrix
        self.inclusion = inclusion

    def get_conductivity(self, phi, **kwargs):
        phi_conjugate = 1- phi
        c_inclusion = self.inclusion.get_conductivity(**kwargs)
        c_matrix = self.matrix.get_conductivity(**kwargs)

        return phi*c_inclusion + c_matrix * phi_conjugate

class SeriesModel:
    """
    https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GC006530
    """
    def __init__(self, matrix, inclusion):
        self.matrix = matrix
        self.inclusion = inclusion

    def get_conductivity(self, phi, **kwargs):
        phi_conjugate = 1- phi
        c_inclusion = self.inclusion.get_conductivity(**kwargs)
        c_matrix = self.matrix.get_conductivity(**kwargs)
        numerator = c_inclusion*c_matrix
        denominator = phi_conjugate*c_inclusion + c_matrix * phi
        return numerator / denominator


class CubeModel:
    mixing_model = 'Cube'
    def __init__(self, matrix_conductivity_model, fluid_conduction_model,\
                 enforce_fluid_is_always_larger=False,variant='waff74-1'):
        self.variant = variant
        self.matrix = matrix_conductivity_model
        self.fluid = fluid_conduction_model
        self.enforce_fluid_is_always_larger=enforce_fluid_is_always_larger

    def get_conductivity(self, phi,**kwargs):
        """

        if variant: waff74-1:

                 bulk = (1 -(1-phi)^(2/3) ) sigma_melt

        elif variant waff74-2:
            numerator = solid * melt * ( 1 - phi)^(2/3)
            denominator = melt * (1 - phi)^(1/3) + solid * (1 - (1-phi)^(1/3))

        (1+2*phi)solid

        if variant partzsh00:
        (look up the paper)

        does waff have 2 models for cubes?
        Parameters
        ----------
        phi
        variant
        kwargs

        Returns
        -------

        """

        c_inclusion = self.fluid.get_conductivity(**kwargs)
        if self.variant=='waff74-1':
            value = (1 -(1-phi)**(2/3) )*c_inclusion
        elif self.variant=='waff74-2':
            c_matrix  = self.matrix.get_conductivity(**kwargs)
            numerator   = c_matrix * c_inclusion * ( 1 - phi)**(2/3)
            denominator = c_inclusion * (1 - phi)**(1/3) + c_matrix * (1 - (1-phi)**(1/3))
            value = numerator/denominator
            value += c_inclusion* (1 - (1-phi)**(2/3))
        elif self.variant=='P00':
            c_matrix  = self.matrix.get_conductivity(**kwargs)
            cube_side_length = (1 - phi)**(1/3)
            value1 = (1 - cube_side_length)/c_inclusion
            value2 = cube_side_length/(c_inclusion*(1-cube_side_length**2) + \
                                       c_matrix*cube_side_length**2)
            value = 1/(value1+value2)

        return value


class ArchiesLaw:
    mixing_model = 'archiesLaw'
    def __init__(self, matrix_conductivity_model, fluid_conduction_model, c, n):
        self.matrix = matrix_conductivity_model
        self.fluid = fluid_conduction_model
        self.c = c
        self.n = n

    def get_conductivity(self, phi, **kwargs):
        c_inclusion = self.fluid.get_conductivity(**kwargs)

        return (self.c*phi**self.n)*c_inclusion


class ArchiesLawGlover:
    mixing_model = 'archiesLaw'
    """
    suggested values for m are:
    m: 1.9 for brine in peridotite: Huang et al., (2021)
    m: 1.05 for silicic melts Gaillard et al., (2005)
    m: 1.2 to match basalt vs peridotite Yoshino et al., (2010)


    """
    def __init__(self, matrix_conductivity_model, fluid_conduction_model, m=1.9):
        self.matrix = matrix_conductivity_model
        self.fluid = fluid_conduction_model
        self.m = m

    def get_conductivity(self, phi=None, **kwargs):
        c_inclusion = self.fluid.get_conductivity(**kwargs)
        c_matrix    = self.matrix.get_conductivity(**kwargs)
        p = np.log(1 - phi**self.m)/np.log(1 - phi)
        term1 = c_matrix*(1 - phi)**p
        term2 = c_inclusion*phi**self.m
        return term1 + term2


class TubesModel:

    mixing_model = 'Tube'
    def __init__(self, matrix_conductivity_model, fluid_conduction_model):
        self.matrix = matrix_conductivity_model
        self.fluid = fluid_conduction_model

    def get_conductivity(self, phi, **kwargs):
        c_inclusion = self.fluid.get_conductivity(**kwargs)
        c_matrix = self.matrix.get_conductivity(**kwargs)

        first = phi*c_inclusion/3
        second = (1-phi)*c_matrix

        return first+second

class GeomAverage:

    mixing_model = 'GM'
    def __init__(self, phases = None, phase_fractions = None):
        if phase_fractions is None:
            self.phase_fractions = [1/len(phases)]*len(phases)
        else:
            self.phase_fractions = phase_fractions
        self.phases = phases

    def get_conductivity(self,phase_relevant_kwargs=None, **kwargs):
        starting_var = create_starting_variable(start_value=1,**kwargs)
        if phase_relevant_kwargs is None:
            for phase, fraction in zip(self.phases,self.phase_fractions):
                starting_var*=phase.get_conductivity(**kwargs)**(fraction)
        else:
            for phase, fraction,phase_kwargs in zip(self.phases,self.phase_fractions,phase_relevant_kwargs):
                starting_var*=phase.get_conductivity(**{**kwargs,**phase_kwargs})**(fraction)


        return starting_var

class ArithmeticAverage:

    mixing_model = 'GM'
    def __init__(self, phases = None, phase_fractions = None):
        if phase_fractions is None:
            self.phase_fractions = [1/len(phases)]*len(phases)
        else:
            self.phase_fractions = phase_fractions
        self.phases = phases

    def get_conductivity(self, **kwargs):
        starting_var = create_starting_variable(**kwargs)
        for phase, fraction in zip(self.phases,self.phase_fractions):
            starting_var+=fraction*phase.get_conductivity(**kwargs)**(fraction)

        return starting_var
