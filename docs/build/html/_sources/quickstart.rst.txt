=====================
Pyrrhenius Quickstart
=====================

Basic Usage
-----------

The basic workflow of using pyrrhenius is structured like so:

1. Import necessary modules

.. jupyter-execute::
    :hide-code:

    import numpy as np
    np.set_printoptions(edgeitems=3, infstr='inf', linewidth=75, nanstr='nan', 
                       precision=8, suppress=False, threshold=3, formatter=None)

.. jupyter-execute::

   import pyrrhenius.database as phsd

2. Create a pyrrhenius database object

.. jupyter-execute::

   ecdatabase = phsd.Database()

3. Create a :py:class:`pyrrhenius.model.Model` object corresponding to the desired ``model_id``. In this case the
model id corresponding to the SEO3 model of olivine (insert citation) is loaded

.. jupyter-execute::

   model = ecdatabase.get_model('SEO3_ol')

4. Call the ``get_conductivity(*args,**kwargs)`` method on the model with the relevant parameters. 

.. jupyter-execute::

   conductivity = model.get_conductivity(T=1000, P=1.0, logfo2=10**-11)
   print('*'*20)
   print('Calculated conductivity (S/m) at T=1000 K, P=1 GPa,and fO2=10^-11 bars:')
   print(conductivity)

Keywords can be of type ``float`` or :py:class:`numpy.ndarray`'s. The latter functionality is especially useful for computing conductivities across a vectorized grid.
In this case, using meshgrid on a dummy variable (``zz``) and a temperature array `T` creates a temperature array of size  (10x10). 

 .. jupyter-execute::

   import numpy as np

   T = np.linspace(500,1500,num=11)
   z = np.arange(0,11,1) # dummy variable for demonstration purposes only 
   tt, xx = np.meshgrid(T,z)
   print(tt.shape)
   tt

As long as the input keyword-arguments are directly broadcastable, pyrrhenius can use mixed ``float`` and ``numpy.ndarray`` objects, vectorizing them as needed.

.. jupyter-execute::

   model.get_conductivity(T=tt, P=1.0, logfo2=10**-11)


Accessing Database Options
--------------------------

Pyrrhenius currently ships with a ``.csv`` database which is loaded by default.

.. jupyter-execute::

   import pyrrhenius.database as phsd

   ecdatabase = phsd.Database()

Once the database object has been created, you can use the ``get_phases()``, ``get_model_list_for_phase()``, and the ``get_model()`` methods
to specify which model to load. 

.. jupyter-execute::

   ecdatabase.get_phases()


.. jupyter-execute::

   ecdatabase.get_model_list_for_phase('granite')


.. jupyter-execute::

   ecmodel = ecdatabase.get_model('han_23_HD_granite')
   ecmodel


Isotropic Models
----------------
The default database comes with a number of anisotropic models, visible as ``model_id``'s with "[xxx]" strings appended to the end. To get an isotropic model, first tell 
the database to generate isotropic models via ``create_isotropic_models()``, then examine the available models

.. jupyter-execute::

   before_isotropic_calculation = ecdatabase.get_model_list_for_phase('plagioclase')
   ecdatabase.create_isotropic_models()
   after_isotropic_calculation = ecdatabase.get_model_list_for_phase('plagioclase')
   print('*'*20)
   print('Before Isotropic Calculation')
   print('*'*20)
   print(*before_isotropic_calculation,sep='\n')
   print('*'*20)
   print('After Isotropic Calculation')
   print('*'*20)
   print(*after_isotropic_calculation,sep='\n')

You should see that calling ``create_isotropic_models()`` on the database procedurally creates new ``model_id``'s where multiple crystal directions are present for the same base id. These procedurally generated new models are identified by a prepended ``isotropic:`` string. They can now be accessed in the same way as default models

.. jupyter-execute::

    ecmodel = ecdatabase.get_model('isotropic_model:yang_12b_plag[100]+yang_12b_plag[010]+yang_12b_plag[001]')
    conductivity = ecmodel.get_conductivity(T=1000, P=1.0)
    conductivity

Mixing Models
-------------

Pyrrhenius provides several N phase mixing models which are accessed via the ``mixing`` module. Since the interfaces for these mixing models 
can be different, consult the documentation prior to using them.

.. jupyter-execute::

    import pyrrhenius.mixing as pyhmix

    brine_id = 'Li_18_1%plg_brine'
    plag_id = 'isotropic_model:yang_12b_plag[100]+yang_12b_plag[010]+yang_12b_plag[001]'

    brine_model = ecdatabase.get_model(brine_id)
    plag_model  = ecdatabase.get_model(plag_id)

    # The HashinStrikman mixing model needs to be initialized with a matrix and inclusion ecmodel
    hashinshtrikman_matrix = pyhmix.HashinShtrikmanUpper(plag_model,brine_model)
    # The Geometric Average model requires intitialization with a phase and phase fraction list. 
    geometric_mixed_matrix = pyhmix.GeomAverage(phases=[brine_model,plag_model],
                                                phase_fractions=[0.05,0.95])

    # Only the HS model in this example requires a provided phase fraction (0.05), positional argument. 
    hs_conductivity = hashinshtrikman_matrix.get_conductivity(0.05,T=1000)
    gm_conductivity = geometric_mixed_matrix.get_conductivity(T=1000)

    # Also calculate endmember phase conductivities for comparison 
    plagioclase_conductivity = plag_model.get_conductivity(T=1000)
    brine_conductivity        = brine_model.get_conductivity(T=1000)

    print(f'HS: {hs_conductivity} GM:{gm_conductivity}')
    print(f'Plag: {plagioclase_conductivity} Brine:{brine_conductivity}')

Metadata Access
---------------
Most pyrrhenius objects come equipped with a ``metadata`` object which describes the source publication, experimental conditions, and calibration settings used to create the model

.. jupyter-execute::

    plag_model.metadata

``metadata`` objects can be used by the parent model to produce input data representative of the experimental conditions

.. jupyter-execute::

    plag_model.generate_representative_conditions()

you can use the output from ``generate_representative_conditions()`` to construct your own input arrays, or directly evaluate the condition dictionary within the model itself

.. jupyter-execute::

    condition_dict = plag_model.generate_representative_conditions()
    plag_model.get_conductivity(**condition_dict)

Plotting Utilities
------------------

Since most experimental petrologists conduct their parameter fitting in :math:`\log_{10}(\sigma), \frac{1}{T}` space, pyrrhenius provides a convenience plotting method to format a `matplotlib.Axis` for a similar plotting space 

.. jupyter-execute::

    import matplotlib.pyplot as plt
    import numpy as np
    import pyrrhenius.mixing as pyhmix
    import pyrrhenius.database as phsd
    import pyrrhenius.utils as pyhutils

    ecdatabase = phsd.Database()
    
    # endmember models
    brine_id = 'Li_18_1%plg_brine'
    plag_id = 'Li_18_wet_plag'

    brine_model = ecdatabase.get_model(brine_id)
    plag_model  = ecdatabase.get_model(plag_id)

    # The HashinStrikman mixing model needs to be initialized with a matrix and inclusion ecmodel
    hashinshtrikman_matrix = pyhmix.HashinShtrikmanUpper(plag_model,brine_model)

    # provide a range of temperature conditions at which to evaluate the models 
    T = np.linspace(400,1200,num=120)

    # Only the HS model in this example requires a provided phase fraction (0.05), positional argument. 
    hs_5pct = hashinshtrikman_matrix.get_conductivity(0.05,T=T)
    hs_1pct = hashinshtrikman_matrix.get_conductivity(0.01,T=T)

    # Also calculate endmember phase conductivities for comparison 
    plagioclase_conductivity = plag_model.get_conductivity(T=T)
    brine_conductivity        = brine_model.get_conductivity(T=T)

    # set up matplotlib plotting 

    fig, ax = plt.subplots()
    linear_major_ticks = np.asarray([2000,1400,1100,900,800,700,600,500,400])
    pyhutils.format_ax_arrhenian_space(ax,linear_major_ticks=linear_major_ticks)

    ax.plot(1e4/T,plagioclase_conductivity,color='purple',label='Li et al., (2018) Plagioclase')
    ax.plot(1e4/T,brine_conductivity,color='orange',label='Li et al., (2018) 1% NaCl Plagioclase-Equilibrated Brine')
    ax.plot(1e4/T,hs_1pct,color='blue',label='Li et al., (2018) 1% NaCl, 1% Vol Plagioclase-Equilibrated Brine',linestyle=':')
    ax.plot(1e4/T,hs_5pct,color='blue',label='Li et al., (2018) 1% NaCl, 5% Vol Plagioclase-Equilibrated Brine',linestyle='--')
    ax.set_title('A Plagioclase, Brine, and two HS mixed Plag-Brine Models')
    ax.legend()
    fig.show()




