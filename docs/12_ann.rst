.. _ann:

========================================
Artificial Neural Network Potentials
========================================

FLAME implements a range of artificial neural network (ANN) potentials.
In the original Behler-Parrinello approach,
the output values of ANNs are considered as atomic energies
which sum up to the total energy of a system.
In the charge equilibration via neural network
technique (CENT),
the ANN serve as an intemediate step 
to include long-range electrostatic interactions.
The parameters related to an ANN must be set
in the block ``ann``, as described below.

ann options
=================

general ``ann`` parameters
------------------------------------------

**subtask**: (string) Sets the task to perform within the ANN schemes.

   default: ``No default value.``

   options:

      ``train``: routines to perform an ANN training.

.. _ref-ann-approach:

**approach**: (string) Determines the ANN technique to be employed.

   default: ``atombased``
    
   options: 

      ``atombased``: the original Behler-Parrinello method, where the total energy is the sum of atomic energies :cite:`Behler2011`
       
      ``cent1``: the charge equilibration via neural network technique :cite:`ghasemi_interatomic_2015`

``train`` parameters
--------------------------------

**optimizer**: (string) Method for the optimization during an ANN training process.


   default: ``No default value.``

   options: 
   
      ``rivals``: Currently the only available method, based on a modification of the extended Kalman filter.

**nstep_opt**: (integer) Number of extended Kalman filter steps.

    default: ``100``


**nconf_rmse**: (integer) Number of configurations for which energy
and force RMSE values will be calculated and reported in the
output.

    default: ``0``

**ampl_rand**: (real) The amplitude of random numbers used
as the initial values of the ANN weights. Arbitrary units.

    default: ``1.d0``

**symfunc**: (string) Determines how the evaluation of the symmetry functions
is treated.

   default: ``only_calculate``

   options: 
   
      ``only_calculate``: values are stored in memory only

      ``write``: values are written to disk

      ``read``:  values are read from disk from a previous run


``cent1`` parameters
--------------------------------
**syslinsolver**: (string) Determines what method to use
to solve the system of linear equations for the charge equilibration process (CEP).

   default: ``direct``

   options: 
   
      ``direct``: direct approach, e.g., Gaussian elemination

      ``operator``: iterative approach by convex optimization

**nstep_cep**: (integer) Maximum number of steps allowed in the self-consistent CEP
required by the CENT approach.

    default: ``200``

**alphax_q**: (real) Stepsize in the CEP optimization.

    default: ``1.d0``

**qgnrmtol**: (real) Convergence criterion for the norm of the gradient
in the CEP when **syslinsolver**  is set to  ``operator``.
Tight values, e.g. ``1.d-5``, may be required for accurate forces.
Units in Ha/Bohr.

    default: ``5.d-4``

**freeBC_direct**: (logical) Selectively employs the
``direct`` method for the **syslinsolver** 
for cluster structures only 
during the CEP. This is important when the reference (training
or validation) data set contains both
periodic and molecular configurations.

   default: ``False``
      

Additional files required for **ann**: ``train``
================================================

**list_posinp_train.yaml**: Contains a single keyword, **files**, with
a list of configurations for training the ANN.

   **files**: (list of strings) List of filenames of the atomic structure files in ``yaml`` format.
   For more information on YAML structure file format, see :ref:`yamlstructure`.

**list_posinp_valid.yaml**: Contains a single keyword, **files**, with
a list of configurations for validating the ANN.

   **files**: (list of strings) List of filenames of the atomic structure files in ``yaml`` format.
   For more information on YAML structure file format, see :ref:`yamlstructure`.

**SE.ann.input.yaml**: ANN parameter files for every chemical element ``SE``
in the system.
They contain element-based parameters of the ANN potential as well as
parameters of the symmetry functions.
They contain two keys, ``main`` and ``symfunc`` (see below).
All paramters are in atomic units.

``main`` parameters
--------------------------------

**nodes**: (list of integers) Determines the architecture of the ANN.
For example, ``[3, 5]`` means
that the ANN has two hidden layers (number of elements) with
three and five nodes, respectively.
Currently, architectures with only one or two hidden
layers are implemented, where the latter is well tested
and has been employed in several applications.

   default: ``No default value.``

**rcut**: (real) The cutoff radius used for the symmetry functions. Units in Bohr.

   default: ``No default value.``

**ampl_chi**: (real) Determines the amplitude of the
hyperbolic tangent function used to map the value of the ANN output nodes
to the atomic electronegativity.
We recommend to employ ``1.d0``, and using smaller values is strongly discouraged.

   default: ``No default value.``

**prefactor_chi**: (real) Determines the prefactor of the argument of the
hyperbolic tangent function used to map the value of the ANN output nodes
to the atomic electronegativity.
We recommend using ``1.d0``.

   default: ``No default value.``

**ener_ref**: (real) Determines the reference energy value.
We recommend to set it to the energy of an isolated atom so that
ANN trains indeed the formation energies.

   default: ``No default value.``

**gausswidth**: (real) Determines the width of the Gaussian function
representing the atomic charge density in CENT.
We recommend to try different values in the range
between ``1.d0`` and ``3.d0``, which correspond to 
most atomic radii. Units in Bohr.

   default: ``No default value.``

**hardness**: (real) Determines the atomic hardness by which
to control how much charge, approximately, to expect to be
collected by this type of atom.
Similar to **gausswidth**, we recommend using
a value in a physically meaningful range given
in textbooks. Units in :math:`{\textrm{Ha}}/{\textrm{Bohr}^{3}}`.

   default: ``No default value.``

**chi0**: (real) Determines the offset of atomic electronegativity. Arbitrary units.

   default: ``No default value.``

**method**: (string) Determines the type of symmetry functions
used as the atomic environment descriptor.
Currently, only symmetry functions of type ``behler`` are implemented.

   default: ``No default value.``

   options: 
   
      ``behler``: For more details see Ref. :cite:`Behler2011`

``symfunc`` parameters
--------------------------------
Two types of symmetry functions, radial and angular, are implemented in FLAME.
A complete description and comparison between these symmetry functions is
given in Ref.:cite:`Behler2011`.
In the parameter file, this information is provided
line by line, in a specific format:

The radial symmetry function in FLAME is called `g02` and has two parameters,
the exponent to control the broadness of the Gaussian function, and
the offset that determines the center of the Gaussian function.
Using a finite offset is not well tested and we recommend to set it to zero.
A `g02` symmetry function is defined by the key `g02_` appended by
a zero-padded enumeration, e.g. `g02_001`.
The value of the key is exponent, offset, 0.0, 0.0, and atom type,
all separated by spaces.
The two zeros at position 2 and 3 are the lower and upper bounds of the symmetry function
for all training data sets.
The zeros cannot be omitted.
The last item is the atomic species of the surrounding atom.

The angular symmetry function is of type `g05` and contains three parameters.
The first parameter, the exponent of the Gaussian function, is similar to that of `g02`.
The other two parameters are the prefactor of the cosine function
and the value of the power.
Similar to `g02`, the parameters must be followed by two zero.
The last two entries per line indicate the two 
atomic species surrounding the center atom.
An example of an ANN parameter file for type Na in sodium chloride
system is given below:

  .. literalinclude:: Na.ann.input.yaml

