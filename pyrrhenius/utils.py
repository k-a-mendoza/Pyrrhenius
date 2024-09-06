import numpy as np




def _find_new_bounds(bound_average, lower_bound, lower_prediction, middle_prediction, upper_bound, upper_prediction):
    product_upper = upper_prediction * middle_prediction < 0
    product_lower = lower_prediction * middle_prediction < 0
    if product_upper.any():
        lower_bound[product_upper] = bound_average[product_upper]
    if product_lower.any():
        upper_bound[product_lower] = bound_average[product_lower]
    bound_average = (lower_bound + upper_bound) / 2
    return bound_average, lower_bound, upper_bound


def create_starting_variable(starting_value=0, P=None, T=None,logfo2=None,Xfe=None,Cw=None,co2=None,nacl=None,na2o=None,sio2=None,**kwargs):
    """
    Generates a ``float`` or ``np.ndarray`` of a shape consistent with the provided keyword arguments

    Currently checks keywords for P, T, logfo2, Xfe, Cw, co2, nacl, na2o, and sio2

    Parameters
    ----------
    starting_value : float, optional
        The variable value for conductivity casting, defaults to 0.

    **kwargs : dict
        Additional parameters used to calculate conductivity

    Returns
    -------
    float or np.ndarray
            The determined dummy conductivity variable used to cast the ``get_conductivity`` result into the proper shape.
    """
    ndarrays = list(filter(lambda x: isinstance(x[1], np.ndarray) and x[1] is not None, zip(['P', 'T', 'logfo2','Xfe','Cw','co2','nacl','na2o','sio2'],
                                                                                          [P, T, logfo2,Xfe,Cw,co2,nacl,na2o,sio2])))
    if len(ndarrays)<1:
        reference_shape=(1)
    elif len(ndarrays)==1:
        reference_shape= ndarrays[0][1].shape
    else:
        reference_shape = ndarrays[0][1].shape
        if not all(array[1].shape == reference_shape for array in ndarrays):
            starting_array =None
            for label, array in ndarrays:
                if starting_array is None:
                    starting_array=np.copy(array)
                elif array is not None:
                    try:
                        starting_array*=array
                    except Exception as e:
                        print('passed arrays will not multiply together\n')
                        for label, array in ndarrays:
                            if array is not None:
                                print(f'\t{label}: {array.shape}')
    return np.ones(reference_shape)*starting_value

def access_array_indices(P=None, T=None, Cw=None, X_fe=None, logfo2=None, co2=None,**kwargs):
    kwargs = [P, T, Cw, X_fe, logfo2, co2]
    if any(isinstance(arg, np.ndarray) for arg in kwargs):
        arrays = filter(lambda x : isinstance(x,np.ndarray),kwargs)
        shapes = [x.shape for x in arrays]
        if len(set(shapes)) > 0:
            assert all(x == shapes[0] for x in shapes), f"Some passed shapes are different. {shapes}"
            return np.zeros(shapes[0])

    return 0

def binary_search(experiment, target_values, bounds, kwargs,
                  tolerance=1e-3, relative=True,debug=True,max_iterations = 50):
    """
    Binary search algorithm used to 

    Parameters
    ----------
    experiment : function
        the conductivity model
    target_values : np.ndarray
        the target values of the conductivity model
    bounds : tuple
        the bounds of the conductivity model
    kwargs : dict
        the keyword arguments to pass to the conductivity model
    tolerance : float
        the tolerance of the conductivity model
    relative : bool
        if True, the tolerance is relative
    debug : bool
        if True, the debug information is printed
    max_iterations : int
        the maximum number of iterations

    Returns
    -------
    np.ndarray
        the root of the conductivity model
    """
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



def calc_QFM(T,P):
    """
    Calculate the Quartz-Fayallite-Magnetite buffer across Pressure and Temperature Conditons
    
    Adapted from the plotting scripts of
    Iacovino, Kayla. (2022). Calculate fO2 Buffer (1.5). Zenodo. https://doi.org/10.5281/zenodo.7387196
    
    Parameters
    ----------
    T : float or np.ndarray
        Temperature in Kelvin
    P : float or np.ndarray
        Pressure in bars
    
    Returns
    -------
    float or np.ndarray
        The calculated QFM buffer value in log10(fo2) units
    """
    a = -26445.3
    aa = -25096.3
    b = 10.344
    bb = 8.735
    c = 0.092
    cc = 0.11

    P_bars = P*10000
    if not isinstance(T,np.ndarray):
        if T < 573:
            val = a/T + b + c *(P_bars-1)/T
        else:
            val = aa/T + bb + cc*(P_bars-1)/T
    elif isinstance(T,np.ndarray) and isinstance(P,np.ndarray):

        val = np.ones(T.shape)
        low_t_mask = T < 573
        high_t_mask = T >=573
        val_low  = a / T[low_t_mask] + b + c * (P_bars[low_t_mask] - 1) / T[low_t_mask]
        val_high = aa / T[high_t_mask] + bb + cc * (P_bars[high_t_mask] - 1) / T[high_t_mask]
        val[low_t_mask]  = val_low
        val[high_t_mask] = val_high

    elif isinstance(T,np.ndarray) and not isinstance(P,np.ndarray):

        val = np.ones(T.shape)
        low_t_mask = T < 573
        high_t_mask = T >=573
        val_low  = a / T[low_t_mask] + b + c * (P_bars - 1) / T[low_t_mask]
        val_high = aa / T[high_t_mask] + bb + cc * (P_bars - 1) / T[high_t_mask]
        val[low_t_mask]  = val_low
        val[high_t_mask] = val_high


    return val


def calc_IW(T,P):
    """
    Calculate the Iron-Wustite buffer based on Pressure and Temperature 
    
    Adapted from the plotting scripts of
    Iacovino, Kayla. (2022). Calculate fO2 Buffer (1.5). Zenodo. https://doi.org/10.5281/zenodo.7387196

    Parameters
    ----------
    T : float or np.ndarray
        Temperature in Kelvin
    P : float or np.ndarray
        Pressure in bars
    
    Returns
    -------
    float or np.ndarray
        The calculated IW buffer value in log10(fo2) units
    
    """
    a0= 6.54106 
    a1= 0.0012324
    b0= -28163.6
    b1= 546.32
    b2= -1.13412
    b3= 0.0019274   

    return a0+a1*P + (b0+b1*P-b2*P**2+b3*P**3)/T


def format_ax_arrhenian_space(ax,xlim=[5,20],ylim=[1e-4,1],linear_x_units='K',inverse_x_units='K',
                yunits=r'$\sigma$ (S/m)',xscale=1e4,linear_major_ticks=None,linear_minor_ticks=None,
                inverse_ticks=None):
    """
    Format an axis for an Arrhenius plot with a top linear scale and a bottom inverse temperature scale.
    The y axis is shared with an assumed logarithmic scale.

    For use with a matplotlib axis. 
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib Axes object to format.
    xlim : list of float, optional
        The limits for the inverse temperature x-axis in 1e4/T (K) [t_min, t_max]. Default is [5, 20].
    ylim : list of float, optional
        The limits for the y-axis. Default is [1e-4, 1].
    linear_x_units : str, optional
        The units for the linear temperature scale. Default is 'K'.
    inverse_x_units : str, optional
        The units for the inverse temperature scale. Default is 'K'.
    yunits : str, optional
        The units for the y-axis. Default is r'$\sigma$'.
    xscale : float, optional
        The scaling factor for the inverse temperature axis. Default is 1e4.
    linear_major_ticks : array-like, optional
        Major tick positions for the linear temperature axis. Default is None, which auto-generates ticks.
    linear_minor_ticks : array-like, optional
        Minor tick positions for the linear temperature axis. Default is None, which auto-generates ticks.
    inverse_ticks : array-like, optional
        Tick positions for the inverse temperature axis. Default is None, which auto-generates ticks.
    
    Returns
    -------
    ax2
        a matplotlib.axes object sharing the same y value
    
    Example
    -------
    >>> import pyrrhenius.utils as pyrhutil
    >>> fig, ax = plt.subplots()
    >>> pyrhutil.format_ax_arrhenian_space(ax)
    >>> plt.show()
    """
    ax2 = ax.twiny()
    circ = r'$^\circ$'
    xscale_label = f'{xscale:.0e}'.replace('+0', '')
    xlim = sorted(xlim)
    t_min = min(xlim)
    t_max = max(xlim)
    ax.set_ylim(ylim)
    ax.set_yscale('log')
    ax.set_xlabel(f'{xscale_label}/T ({circ}{inverse_x_units})', loc='left')
    ax.set_ylabel(yunits)
    ax2.set_xlabel(f'T ({circ}{linear_x_units})', loc='left')

    if linear_x_units == 'K':
        min_k_location = np.round(xscale / t_max,-2)
        max_k_location = np.round(xscale / t_min,-2)
        if linear_major_ticks is None:
            linear_major_ticks = np.arange(min_k_location,max_k_location + 200, 200)
            linear_x_tick_major_locations = xscale / linear_major_ticks
            linear_x_tick_major_labels = linear_major_ticks
        else:
            linear_x_tick_major_locations = xscale / linear_major_ticks
            linear_x_tick_major_labels = linear_major_ticks

        if linear_minor_ticks is None:
            linear_minor_ticks = np.arange(min_k_location,max_k_location + 100, 100)
            linear_x_tick_minor_locations = xscale / linear_minor_ticks
        else:
            linear_x_tick_minor_locations = xscale / linear_minor_ticks

    elif linear_x_units == 'C':
        min_c_location = np.round(xscale / t_max -273.15,-2)
        max_c_location = np.round(xscale / t_min -273.15,-2)
        if linear_major_ticks is None:
           
            linear_major_ticks_c = np.arange(min_c_location, max_c_location + 200, 200)
            linear_major_ticks_k_location = xscale/(linear_major_ticks_c_location+273.15)
            
            linear_x_tick_major_locations = linear_major_ticks_k_location
            linear_x_tick_major_labels    = linear_major_ticks_c
        else:
            linear_x_tick_major_locations = xscale / (linear_major_ticks + 273.15)
            linear_x_tick_major_labels    = linear_major_ticks

        if linear_minor_ticks is None:
            linear_major_ticks_c = np.arange(min_c_location, max_c_location + 100, 100)
            linear_x_tick_minor_locations = xscale/(linear_major_ticks_c+273.15)
        else:
            linear_x_tick_minor_locations = xscale / (linear_minor_ticks+273.15)

    if inverse_ticks is None:
        inverse_ticks = np.arange(xlim[0],xlim[1]+1,1)

    ax.set_xticks(inverse_ticks)
    ax.set_xlim([t_min,t_max])
    
    # Configure the linear x-axis ticks (top axis)
    ax2.set_xticks(linear_x_tick_major_locations)
    ax2.set_xticklabels(linear_x_tick_major_labels,minor=False)
    ax2.set_xlim([t_min, t_max])
    # Add minor ticks to the top axis if provided
    ax2.set_xticks(linear_x_tick_minor_locations, minor=True)
    ax2.tick_params(axis='x', which='minor')
    return ax2
