import numpy as np
from .model import WaterPTCorrection
from scipy.optimize import minimize, root_scalar, root
from collections import OrderedDict
from . import model 
from typing import Union, List
from abc import ABC, abstractmethod
from . import utils as pyutils



class AbstractMixingModel(ABC):
    """
    Abstract class for all mixing models
    """

    def __init__(self,phases : List[model.ModelInterface]):
        """
        Initializes the AbstractMixingModel class

        Parameters
        ----------
        phases : List[model.ModelInterface]
            the phases to mix
        """
        if isinstance(phases,list) or isinstance(phases,tuple):
            self.phases = phases
        else:
            self.phases = [phases]

    def prep_phase_fractions(self,fractions):
        """
        Prepares the phase fractions for the mixture

        Parameters
        ----------
        phase_fractions : Union[np.ndarray,float, List[float]]
            the volume fractions of the conductive phase

        Raises
        ------
        AssertionError
            if the phase fractions do not sum to 1

        Returns
        -------
        np.ndarray
            the volume fractions of each phase 
        """
        if isinstance(fractions,float):
            assert len(self.phases)<3, "not enough phase_fractions provided"
            phase_fractions = np.ones(2)
            phase_fractions[0]=fractions
            phase_fractions[1]=1-fractions
        elif isinstance(fractions,list) or isinstance(fractions):
            assert len(self.phases)==len(fractions),"passed fraction and phase list length mismatch."
            phase_fractions = np.array(fractions)

        assert sum(phase_fractions)==1.0, "phase_fractions don't add up to one"
        return phase_fractions
    
    @abstractmethod 
    def get_conductivity(self, **kwargs):
        pass

class BruggemanSolver:
    """
    Bruggeman solver for the effective medium theory

    The Bruggeman symmetrical equation is the most widely used EMT equation. It is sometimes known as the Bruggemane Landauer equation. 
    The system is considered to be an aggregate of unconnected units

    The bulk conductivity of the system is parameterized as:

    sum_i f_i * (sigma_i - sigma_bulk) / (sigma_i + (1/p_i - 1) * sigma_bulk) = 0 

    and thus, sigma_bulk must be solved for. To generalize over N-phases, sigma_bulk is solved for iteratively using the binary search algorithm. 

    """

    def __init__(self, models,polarizations,phase_fractions,**kwargs):
        """
        Initializes the BruggemanSolver class

        Parameters
        ----------
        models : List[model.ModelInterface]
            the models to mix
        polarizations : List[float]
            the polarizations of the models
        phase_fractions : List[float]
            the volume fractions of the models
        kwargs : dict
            the keyword arguments to pass to the conductivity function
        """
        self.phase_fractions = phase_fractions
        starting_var    = pyutils.create_starting_variable(starting_value=1,**kwargs)
        self.constituent_conductivities = [x.get_conductivity(**kwargs)*starting_var for x in models]
        self.constituent_polarizations  = polarizations
        self.baseline_img =  pyutils.create_starting_variable(**kwargs)

    def __call__(self, bulk):
        """
        Evaluates the residual of the Bruggeman equation given an estimated bulk conductivity. 

        Parameters
        ----------
        bulk : float or np.ndarray
            the bulk conductivity to evaluate the residual at

        Returns
        -------
        float or np.ndarray
            the residual of the Bruggeman equation
        Raises
        ------
        AssertionError
            if the bulk conductivity is less than the minimum constituent conductivity or greater than the maximum constituent conductivity
        """
        assert all(bulk >= np.min(np.stack(self.constituent_conductivities,axis=-1),axis=-1)), "Bulk conductivity is less than the minimum constituent conductivity. Cannot solve for bulk conductivity."
        assert all(bulk <= np.max(np.stack(self.constituent_conductivities,axis=-1),axis=-1)), "Bulk conductivity is greater than the maximum constituent conductivity. Cannot solve for bulk conductivity."
        baseline_img =  np.copy(self.baseline_img)
        for conductivity, P, volume_fraction in zip(self.constituent_conductivities,
                                                                self.constituent_polarizations,
                                                                self.phase_fractions):

            numerator = conductivity - bulk
            denominator = conductivity + (1 / P - 1) * bulk
            new_value = volume_fraction * numerator / denominator
            baseline_img += new_value

        return baseline_img



class EffectiveMediumTheoryMixture(AbstractMixingModel):
    """
    Mixing model based on the effective medium theory

    The effective medium theory is a method to predict the macroscopic conductivity of a mixture of materials. 
    It is based on the assumption that the mixture is a random distribution of the individual phases.

    """

    def __init__(self,phases=None,type='Bruggeman'):
        """
        Initializes the EffectiveMediumTheoryMixture class

        Parameters
        ----------
        phase_list : list
            the list of phases to mix
        type : str
            the type of the effective medium theory to use. Options are currently limited to 'Bruggeman', but may extend to others once
            the relevant literature is evaluated. 
        kwargs : dict
            the keyword arguments to pass to the conductivity function

        Raises
        ------
        AssertionError
            if the type is not supported
        """
        super().__init__(phases)
        if type=='Bruggeman':
            self.solver = BruggemanSolver
        else:
            raise ValueError(f"Type {type} is not supported for the effective medium theory. Please choose 'Bruggeman'.")
        

    def set_polarization_factors(self,polarization_factors):
        """
        Sets the polarization factors for the phases

        Parameters
        ----------
        phase_list : list
            the list of phases
        polarization_factors : list
            the polarization factors of the phases

        Returns
        -------
        list
            the polarization factors of the phases
        Raises
        ------
        AssertionError
            if the polarization factors are not between 0 and 1/3
        """
        if polarization_factors is None:
            pol_fac = [1/3]*len(self.phases)
        else:
            pol_fac = []
            for x in polarization_factors:
                assert x > 0 and x < 1/3, "Polarization factors must be between 0 and 1/3. Current factors: {x}"
                if x is None:
                    pol_fac.append(1/3)
                else:
                    pol_fac.append(x)
        return pol_fac

    def get_conductivity(self,phase_fractions=None,polarization_factors=None,**kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        phase_fractions : list
            the fractions of the phases
        polarization_factors : list
            the polarization factors of the phases. Defaults to 1/3 for all phases. 
            For reference, 1/3 is the polarization factor for a sphere and 0 is the polarization factor for a rod. 
        kwargs : dict
            the keyword arguments to pass to the conductivity function for the mixture

        Returns
        -------
        float or np.ndarray
            the conductivity of the mixture
        """
        phase_fractions            = self.prep_phase_fractions(phase_fractions)
        polarization_factors       = self.set_polarization_factors(polarization_factors)
        starting_var    = pyutils.create_starting_variable(starting_value=1,**kwargs)
        constituent_conductivities = np.stack([x.get_conductivity(**kwargs)*starting_var \
                                                   for x in self.phases], axis=-1)

        bounds_min = np.min(constituent_conductivities,axis=-1)
        bounds_max = np.max(constituent_conductivities,axis=-1)
        solver = self.solver(self.phases,polarization_factors,phase_fractions,**kwargs)
        if constituent_conductivities.ndim==1:
            result = root_scalar(solver,x0 = np.zeros(bounds_max.shape),
                                 bracket=(bounds_min,bounds_max)).root
        else:
            result = pyutils.binary_search(solver,np.zeros(bounds_max.shape),
                      bounds=(bounds_min, bounds_max),
                      kwargs={},
                      tolerance=1e-8, relative=True,debug=False)

        return result



class HashinStrickmanBound(AbstractMixingModel):
    """
    Represents the Hashin-Shtrikman bound for the conductivity of a mixture

    The Hashin-Shtrikman bound simulates the conductivity of spheres embedded in a matrix. 


    summation_term = sum phase_fractions / (sigma_i +2*sigma_{min/max})
    sigma_total    =  summation_term^-1 - 2*sigma_{min/max}

    """
    def __init__(self,phases : List[model.ModelInterface]):
        """
        Initializes the HashinShtrikmanLower class

        Parameters
        ----------
        sigma_0 : model.ModelInterface
            the conductivity model of the matrix
        sigma_1 : model.ModelInterface
            the conductivity model of the inclusion
        """
        super().__init__(phases)

    def get_conductivity(self,fractions: Union[np.ndarray,float, List[float]],type='max',**kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fractions : Union[np.ndarray,float, List[float]]
            the volume fractions of the conductive phase
        type : str
            the type of the Hashin-Shtrikman bound to use. Options are 'max' or 'min'.
        kwargs : dict
            the keyword arguments to pass to the conductivity function

        Returns
        -------
        float or np.ndarray
            the conductivity of the mixture
        """
        phase_fractions = self.prep_phase_fractions(fractions)
        starting_var    = pyutils.create_starting_variable(starting_value=1,**kwargs)
        conductivities = [x.get_conductivity(**kwargs)*starting_var for x in self.phases]
        if type=='max':
            sigma_1 = np.max(np.stack(conductivities,axis=-1),axis=-1)
        elif type=='min':
            sigma_1 = np.min(np.stack(conductivities,axis=-1),axis=-1)

        summation_term = pyutils.create_starting_variable(**kwargs)
        for phase, c in zip(phase_fractions,conductivities):
            summation_term+= phase / (c + 2*sigma_1)
        sigma_total =  1/summation_term - 2*sigma_1
        return sigma_total
    

class ParallelModel(AbstractMixingModel):
    """
    Represents the conductivity of two phases connected in parallel.
    The parameterization looks like:
    
    sigma_total = sum  phase_fraction_i * conductivity_i

    """
    def __init__(self,phases : List[model.ModelInterface]):
       """
       Initializes the ParallelModel class

       Parameters
       ----------
       phases : List[model.ModelInterface]
           the phases to connect in parallel
       """
       super().__init__(phases)

    def get_conductivity(self, fractions, **kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fractions : Union[np.ndarray,float, List[float]]
            the volume fractions of the conductive phase
      
        kwargs : dict
            the keyword arguments to pass to the conductivity function

        Returns
        -------
        float or np.ndarray
            the conductivity of the mixture
        """
        phase_fractions = self.prep_phase_fractions(fractions)
        summation_term = pyutils.create_starting_variable(**kwargs)
        for fraction, phase in zip(phase_fractions,self.phases):
            summation_term+=fraction*phase.get_conductivity(**kwargs)
        return summation_term
    

class SeriesModel(AbstractMixingModel):
    """
    Represents the conductivity of several phases connected in series. 
    The parameterization looks like:
    1/sigma_total = sum phase_fraction_i / conductivity_i

    """
    def __init__(self, phases):
        """
        Initializes the SeriesModel class

        Parameters
        ----------
        phases : List[model.ModelInterface]
            the phases to connect in series
        """
        super().__init__(phases)

    def get_conductivity(self, fractions, **kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fractions : float or np.ndarray
            the porosity of the mixture/ the volume fraction of the inclusion
        kwargs : dict
            the keyword arguments to pass to the conductivity function
        """
        phase_fractions = self.prep_phase_fractions(fractions)
        summation_term = pyutils.create_starting_variable(**kwargs)
        for fraction, phase in zip(phase_fractions,self.phases):
            summation_term+=fraction/phase.get_conductivity(**kwargs)
        return 1/summation_term


class CubeModel(AbstractMixingModel):
    """
    Represents the conductivity of a lattice of cubes separated by thin sheets. Also known as the thin film model.

    Two variants are implemented:
    waff74-1:
    sigma_total = (1 -(1-phi)^(2/3) ) * sigma_melt

    waff74-2:
    sigma_total = solid * melt * ( 1 - phi)^(2/3) / ( melt * (1 - phi)^(1/3) + solid * (1 - (1-phi)^(1/3)) )
    """
    mixing_model = 'Cube'
    def __init__(self, phases,variant='waff74-1'):
        """
        Initializes the CubeModel class

        Parameters
        ----------
        phases : List[model.ModelInterface]
            a list of models corresponding to the consitutent phases. The model assuems the first phase 
            is the matrix phase, and the second is the inclusion. Alternatively a single phase passed will be interpreted to 
            represent that of the fluid, but will fail if the variant selected requires two phases
       
        variant : str
            the variant of the cube model to use
            Two variants are implemented:
            -waff74-1:
             sigma_total = (1 -(1-phi)^(2/3) ) * sigma_melt

            -waff74-2:
             sigma_total = solid * melt * ( 1 - phi)^(2/3) / ( melt * (1 - phi)^(1/3) + solid * (1 - (1-phi)^(1/3)) )

            -P00:
             sigma_total = 1/( (1 - phi) / sigma_solid + phi / sigma_fluid )

        """
        super().__init__(phases)
        self.variant = variant

    def get_conductivity(self, fractions,**kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fractions : float or np.ndarray
            the volume fractions of the conductive phase
        kwargs : dict
            the keyword arguments to pass to the conductivity function

        Returns
        -------
        float or np.ndarray
            the conductivity of the mixture
        """
        phase_fractions = self.prep_phase_fractions(fractions)
        phi = phase_fractions[0]
        c_inclusion = self.phases[0].get_conductivity(**kwargs)
        if self.variant=='waff74-1':
            value = (1 -(1-phi)**(2/3) )*c_inclusion
        elif self.variant=='waff74-2':
            c_matrix  = self.phases[0].get_conductivity(**kwargs)
            numerator   = c_matrix * c_inclusion * ( 1 - phi)**(2/3)
            denominator = c_inclusion * (1 - phi)**(1/3) + c_matrix * (1 - (1-phi)**(1/3))
            value = numerator/denominator
            value += c_inclusion* (1 - (1-phi)**(2/3))
        elif self.variant=='P00':
            c_matrix  = self.phases[0].get_conductivity(**kwargs)
            cube_side_length = (1 - phi)**(1/3)
            value1 = (1 - cube_side_length)/c_inclusion
            value2 = cube_side_length/(c_inclusion*(1-cube_side_length**2) + \
                                       c_matrix*cube_side_length**2)
            value = 1/(value1+value2)

        return value


class ArchiesLaw(AbstractMixingModel):
    """
    Represents the conductivity of a mixture using Archie's law

    The parameterization looks like:
    sigma_total = c*phi^n * sigma_fluid
    """
    def __init__(self, phases, c=1, n=1):
        """
        Initializes the ArchiesLaw class

        Parameters
        ----------
        matrix_conductivity_model : model.ModelInterface
            the conductivity model of the matrix
        fluid_conduction_model : model.ModelInterface
            the conductivity model of the fluid
        c : float
            the coefficient of the parameterization. Defaults to 1
        n : float
            the exponent of the parameterization. Defaults to 1
        """
        super().__init__(phases)
        self.c = c
        self.n = n

    def get_conductivity(self, fraction, **kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fraction : float or np.ndarray
            the porosity of the mixture
        kwargs : dict
            the keyword arguments to pass to the conductivity function

        Returns
        -------
        float or np.ndarray
            the conductivity of the mixture
        """
        phase_fraction = self.prep_phase_fractions(fraction)[0]
        c_inclusion = self.phases[0].get_conductivity(**kwargs)

        return (self.c*phase_fraction**self.n)*c_inclusion


class ArchiesLawGlover(AbstractMixingModel):
    """
    a reparameterization of archies law to reduce the number of unknown coefficients first described by Glover et al., (2000)

    the parameterization looks like:
    sigma_total = C_m * (1 - phi)^p + C_f * phi^m

    where p: 
    p = log(1 - phi^m)/log(1 - phi)

    suggested values for m are:
    m: 1.9 for brine in peridotite: Huang et al., (2021)
    m: 1.05 for silicic melts Gaillard et al., (2005)
    m: 1.2 to match basalt vs peridotite Yoshino et al., (2010)


    """
    def __init__(self, phases, m=1.9):
        """
        Initializes the ArchiesLawGlover class

        Parameters
        ----------
        phases: model.ModelInterface
            the list of phases for the model. Assumes the first phase is the inclusion and the second is the matrix
        m : float
            the exponent of the parameterization
        """
        super().__init__(phases)
        self.m = m

    def get_conductivity(self, fractions, **kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fractions : float or np.ndarray
            the volume fractions of the mixture
        kwargs : dict
            the keyword arguments to pass to the conductivity function

        Returns
        -------
        float or np.ndarray
            the conductivity of the mixture
        """
        phase_fraction = self.prep_phase_fractions(fractions)
        c_inclusion = self.phases[0].get_conductivity(**kwargs)
        c_matrix    = self.phases[1].get_conductivity(**kwargs)
        p = np.log(1 - phase_fraction[0]**self.m)/np.log(phase_fraction[1])
        term1 = c_matrix*phase_fraction[1]**p
        term2 = c_inclusion*phase_fraction[0]**self.m
        return term1 + term2


class TubesModel(AbstractMixingModel):
    """
    Represents the conductivity of a tube-like lattice embedded in a matrix

    The parameterization looks like:
    sigma_total = C_m * (1 - phi) + C_f * phi
    """

    def __init__(self, phases : model.ModelInterface):
        """
        Initializes the TubesModel class

        Parameters
        ----------
        phases: model.ModelInterface
            the list of phases for the model. Assumes the first phase is the inclusion and the second is the matrix

        """
        super().__init__(phases)
 

    def get_conductivity(self, fractions, **kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fractions : float or np.ndarray
            the phase fractions of the mixture
        kwargs : dict
            the keyword arguments to pass to the conductivity function
        """
        phase_fraction = self.prep_phase_fractions(fractions)
        c_inclusion = self.phases[1].get_conductivity(**kwargs)
        c_matrix = self.phases[0].get_conductivity(**kwargs)

        first  = phase_fraction[0]*c_inclusion/3
        second = phase_fraction[1]*c_matrix

        return first+second

class GeomAverage(AbstractMixingModel):
    """
    Represents the geometric average of the conductivities of the phases

    The geometric average is defined as:
    sigma_total = exp(1/N * sum(log(C_i)))
    or alternatively as
    sigma_total = (prod(C_i^(f_i)))
    """

    def __init__(self, phases : List[model.ModelInterface]):
        """
        Initializes the GeometricAverage class

        Parameters
        ----------
        phases : List[model.ModelInterface]
            the phases to average
        """
        super().__init__(phases)

    def get_conductivity(self,fractions=None, **kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        fractions : dict
            the volume fractions of each phase
        kwargs : dict
            the keyword arguments to pass to the conductivity function
        """
        starting_var = pyutils.create_starting_variable(starting_value=1,**kwargs)
        phase_fractions = self.prep_phase_fractions(fractions)
        conductivities = [starting_var*phase.get_conductivity(**kwargs) for phase in self.phases]
        for c, fraction in zip(conductivities,phase_fractions):
            starting_var*= c**fraction

        return starting_var

class ArithmeticAverage(AbstractMixingModel):
    """
    Represents the arithmetic average of the conductivities of the phases

    The arithmetic average is defined as:
    sigma_total = sum(f_i * C_i)
    """

    def __init__(self, phases):
        """
        Initializes the ArithmeticAverage class

        Parameters
        ----------
        phases : list
            the phases to average
        phase_fractions : list
            the fractions of the phases
        """
        super().__init__(phases)
        

    def get_conductivity(self,fractions, **kwargs):
        """
        Calculates the conductivity of the mixture

        Parameters
        ----------
        kwargs : dict
            the keyword arguments to pass to the conductivity function

        Returns
        -------
        float or np.ndarray
            the conductivity of the mixture
        """
        starting_var = pyutils.create_starting_variable(**kwargs)
        phase_fractions = self.prep_phase_fractions(fractions)
        for phase, fraction in zip(self.phases,phase_fractions):
            starting_var+=fraction*phase.get_conductivity(**kwargs)

        return starting_var
