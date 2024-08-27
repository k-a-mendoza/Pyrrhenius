=========================
Overview: Why Pyrrhenius?
=========================

Pyrrhenius provides a powerful, flexible, and intuitive framework for working with mineral electrical conductivity models and data in Python. It introduces a modular, object-oriented approach that allows scientists to efficiently compute conductivities while seamlessly keeping track of metadata about the minerals, experimental conditions, and model parameters.

Motivation
----------


Over the past two decades, electromagnetic geophysical techniques have greatly improved in resolution and spatial extent. At the same time, mineral physicists have produced hundreds of models relating the physiochemical state of common minerals to their electrical conductivity.

Pyrrhenius was designed to provide a means to tractably merge electrical modeling with laboratory observations, allowing geoscientists to flexibly use, compare, validate, and extend what is currently a heterogeneous corpus of mineral physics research. By abstracting away the implementation details and unit conversions, Pyrrhenius enables researchers to focus on science.

The choice of Python and integration with the scientific Python ecosystem makes Pyrrhenius accessible to a wide audience and interoperable with existing electromagnetic geophysics libraries like MtPy and SimPEG. The open-source, object-oriented design also facilitates community contributions to expand the database and functionality over time.


The Arrhenius Equation
----------------------

.. math:: 
   :label: arrhenius

   K = A \cdot \exp\left(\frac{-\Delta H}{k_b T}\right)   



The Arrhenius equation :eq:`arrhenius` is a fundamental formula in chemistry that describes the temperature dependence of reaction rates :math:`(K)` based on an activation energy :math:`\Delta H`. As it turns out, the movement of electrons is a type of reaction, and thus can be described in a similar manner

.. math::   
   :label: arrhenius_electric

   \sigma = \sigma_0 \cdot \exp\left(\frac{-\Delta H}{k_b T}\right) 


Where: 

.. list-table:: 
   :widths: 33 33 33
   :header-rows: 1

   * - Symbol
     - Description
     - Units
   * - :math:`\sigma`
     - Electric Conductivity
     - Siements per meter (S/m)
   * - :math:`\sigma_0`
     - The Preexponental Factor
     - Siements per meter (S/m)
   * - :math:`\Delta H`
     - The Reaction Enthalpy
     - Electron-Volts (eV)
   * - :math:`k_b`
     - Boltzmann's Constant
     - Electron-Volts per Kelvin (eV/K)
   * - :math:`T`
     - Absolute Temperature
     - Kelvin (K)
   
An alternate formulation of the Arrhenius equation uses the gas constant :math:`R` instead of :math:`k_b`. This substitution results in a unit change for :math:`\Delta H`, but otherwise the relationship is the same. 

An Arrhenian Catastrophe
------------------------

The values of both :math:`\Delta H` and :math:`\sigma_0` may themselves be dependent on :math:`\mathbf{PX}` *Physiochemical state* (redox state, electron carrier density, pressure, and other thermodynamic parameters), such that it is more appropriate to rewrite :eq:`arrhenius_electric` as

.. math:: 
    
    \sigma = \sigma_0(\mathbf{PX}) \cdot \exp\left(\frac{-\Delta H(\mathbf{PX})}{k_b T}\right)

Where :math:`X_j` indicates a dependence on the material chemistry and :math:`fO_2` (oxygen fugacity) is used as a proxy for the redox state.
Furthermore, within a single geologic material multiple reactions may result in the movement of electrons. The linear combination of several of these mechanisms is then needed to describe the mineral's total conductivity

.. math::

    \sigma_{total} =\sum_i \sigma_i(\mathbf{PX})  \cdot \exp\left(\frac{-\Delta H_i(\mathbf{PX})}{k_b T}\right)


While an arrhenius equation can be used to describe the conductivity of materials, there are likely numerous unexpected complexities arising from unknown dependencies on numerous variables. This presents several design challenges, as there is a large hereogeneity in both the published parameter values and mathematical form of arrhenian-type equations. Lets look at a few type examples of these
where numerical constants have been replaced with subscripted :math:`\alpha_i`'s. Also,the following symbols are used:

* :math:`X_{fe}` for Iron Fraction 

* :math:`C_w` for Water Content 

* :math:`P` for Pressure


**Omphacite**

Liu et al., (2021) "Electrical conductivity of omphacite and garnet indicates limited deep water recycling by crust subduction".10.1016/j.epsl.2021.116784

.. math::

   \sigma = \alpha_0 (X_{fe} + \alpha_1)^{\alpha_2} C_w^{\alpha_3} \exp\left( \frac{-\alpha_4 - \alpha_5\cdot(X_{fe} + \alpha_6) - \alpha_7\cdot C_w^{\alpha_8}}{k_bT} \right)

**Orthopyroxene**

Zhang & Yoshino (2016) "Effect of temperature, pressure and iron content on the electrical conductivity of orthopyroxene". 10.1007/s00410-016-1315-z

.. math::

   \sigma = \alpha_0 \exp\left( \frac{-(\alpha_1 + \alpha_2 \cdot P)}{k_bT} \right) 
   + \alpha_3 \cdot X_{fe} \cdot \exp\left( \frac{-(\alpha_4 + \alpha_5 \cdot X_{fe}^{\alpha_6} + P \cdot (\alpha_7 + \alpha_8 \cdot X_{fe}))}{k_b T} \right)

**Silicate Melt**

Pommier & Le-Trong (2011). "SIGMELTS: A web portal for electrical conductivity calculations in geosciences". 10.1016/j.cageo.2011.01.002

.. math::

   \sigma = \left( \exp(\alpha_0 \cdot \ln(C_w) + \alpha_1) \right) \cdot \exp\left( \frac{-(\alpha_2 \cdot \ln(C_w) + \alpha_3 + \alpha_4 \cdot P)}{k_bT} \right)

Clearly, experimentalists suggest more or less complicated relationships are needed to parameterize electric conductivity in terms of geologically relevant chemistries and thermodynamic conditions. The upside of this approach is that the sensitivity of electric conductivity to mineral-relevant factors can be modeled. The downside is that it becomes hard to create a non hard coded database of parameterizations.

Since the :math:`\Delta H` and :math:`\sigma_0`'s of the original arrhenian equation might have oft repeated paramterizations, one way to simplify the matter is to encapsulate all Enthalpy and preexponential-like constants into a object oriented hierarchy. Lets apply this concept to the previous equations and see what happens. Each bolded constant is assumed to depend on Physiochemical states. Their subscripts indicate a unique category of *meta*-parameter 

**Omphacite**


.. math::

   \sigma = \mathbf{\alpha_0}\exp\left( \frac{-\mathbf{\Delta H_0}}{k_bT} \right)

**Orthopyroxene**

.. math::

   \sigma = \alpha_0 \exp\left( \frac{\mathbf{\Delta H_0}}{k_bT} \right) + \alpha_1 \cdot \exp\left( \frac{ \mathbf{\Delta H_1}}{k_b T} \right)

**Silicate Melt**

.. math::

   \sigma = \mathbf{\alpha_0}\cdot \exp\left( \frac{\mathbf{\Delta H_0}}{k_bT} \right)

So as long as each *meta*-parameter is handled appropriately, each mineral can be represented by a simplified internal state essentially corresponding to an arrhenian equation. 

We can do one better though. Each arrhenian-type equation can be substituted by an Arrhenian object:

**Omphacite**


.. math::

   \sigma = \mathbf{Arr}(\sigma_0,\Delta H_0)

**Orthopyroxene**


.. math::

   \sigma = \mathbf{Arr}(\sigma_0,\Delta H_0) + \mathbf{Arr}(\sigma_1,\Delta H_1)

**Silicate Melt**


.. math::

   \sigma = \mathbf{Arr}(\sigma_0,\Delta H_0) 

Use of the arrhenian object seems to simplify the relationships to either linear superpositions of equations, or single equation objects, provided the input parameters are implemented correctly. Thus we can see that by adopting an *object oriented* approach, representing a heterogeneous inventory of possible conductivity parameterizations can be dramatically simplified. 


Simplifying Design with Objects
---------------------------------

Object-oriented programming (OOP) is a programming paradigm that organizes software design around data, or objects, rather than functions and logic. The main principles of OOP include encapsulation (aggregating behaviors in isolation), abstraction (providing high-level interfaces to more complex operations), inheritance (providing variants of objects with small changes in behavior), and polymorphism (enforcement of a single outward-facing structure to several heterogeneous internal ones). Use of OOP can allow for simplified, maintainable, extensible, and debuggable code. 

Pyrrhenius utilizes OOP design across four levels of abstraction:

- **Level 1** Empirically-fitted constants (:py:class:`pyrrhenius.mechanisms.StochasticConstants`) are objects, allowing for both mean value sampling and random sampling via bootstraping using a single keyword-argument.  

- **Level 2** Variants of the Arrhenius equation (:py:class:`pyrrhenius.mechanisms.Mechanism` and children) are coded in an inheritance heirarchy, using :py:class:`pyrrhenius.mechanisms.StochasticConstants` to calculate arrhenian equation parameters. Most child class differences revolve around providing a more complex enthalpy or preexponential constant.

- **Level 3** Mineral physics models are represented by :py:class:`pyrrhenius.model.Model` objects, which contain a single :py:class:`pyrrhenius.mechanisms.Mechanism` or other :py:class:`pyrrhenius.model.Model` objects. :py:class:`pyrrhenius.model.Model`'s also use a metadata property object :py:class:`pyrrhenius.model.PublicationMetadata`, taking care of metadata management so you don't have to.

- **Level 4** A database  object :py:class:`pyrrhenius.database.Database` which provides a flexible interface to access specific :py:class:`pyrrhenius.model.Model`'s reported by the literature. 

Because each level of abstraction utilizes a near-identical interface, more complicated mineral models involving compound mechanisms or models are easily created using the `+` operator. These unified interfaces also allow for easy application of mineral mixing models.  Thus, adapting pyrrhenius to your needs is usually easy to express in-code and does not require modifications to the codebase.

Key Features 
------------

Pyrrhenius offers several key features that enable expressive and efficient computations of mineral electric conductivity:

- A curated collection of over 100 models spanning a wide range of minerals and conditions

- An extensible spreadsheet-based database framework for accessing models and associated metadata, so substituting your own work-groups database is as easy as providing an alternate ``.csv`` file.

- Calculate electric conductivity across any mechanism combinations of mechanisms, model or combinations of models using the ``get_conductivity(*args,**kwargs)`` method.

- Linearly combine metadata, mechanisms, and models using the ``+`` operator

- Stochastic sample model parameters by using the optional keyword-argument ``sample=True``

- Easily apply Cubes, Tubes, hashinStrikman, Effective Medium Theory, and Archie's Law phase mixing Variants

- Access human readable metadata via the ``.metadata`` property of all objects

- Integration with the broader scientific Python ecosystem including NumPy, SciPy, and Pandas

With Pyrrhenius, you can go from a scattered collection of published conductivity models to science-driven insights in just a few lines of code
The immediate payoff of using Pyrrhenius is that you'll write less code while worrying less about errors. The long-term payoff is that you'll be able to entertain more advanced geophysical analyses while easily recognizing what you did weeks to months ago.