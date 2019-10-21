.. _ann:

========================================
Artificial Neural Network Potentials
========================================

FLAME implements a range of artificial neural network (ANN) potentials.
In high dimensional ANN potentials, the total
energy is decomposed into atomic energies.
The output values of ANNs are considered as atomic energies
in the original Behler-Parrinello approach, while
the charge equilibration via neural network
technique (CENT) 
includes long-range electrostatic interactions.
The parameters related to ANN must be set
in the block ``ann`` described below.

ann options
=================

general ``ann`` parameters
------------------------------------------

**subtask**: (string) Sets the task to perform within the ANN schemes.

   default: ``No default value.``

   options:

      ``train``: routines to perform an ANN training.

..  warning:: Other methods?

.. _ref-ann-approach:

**approach**: (string) Determines the ANN technique to be employed.

   default: ``atombased``
    
   options: 

      ``atombased``: total energy is the sum of atomic energies, the original Behler-Parrinello method
       
      ``cent1``: charge equilibration via neural network :cite:`ghasemi_interatomic_2015`

``train`` parameters
--------------------------------

**optimizer**: (string) Method for the optimization during the ANN training process.


   default: ``No default value.``

   options: 
   
      ``rivals``: Well tested and stable optimizer based on a modification of extended Kalman filter.

..  warning:: Other methods?

**nstep_opt**: (integer) Number of extended Kalman filter steps

    default: ``100``


**nconf_rmse**: (integer) Number of configurations for which energy
and force RMSE values will be calculated and reported in the
output.

    default: ``0``

**ampl_rand**: (real) The amplitude of random numbers used
as the initial values of the ANN weights.

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
to solve the system of linear equations for the CEP.

   default: ``direct``

   options: 
   
      ``direct``: direct approach, e.g., Gaussian elemination

      ``operator``: iterative approach by convex optimization

**nstep_cep**: (integer) Maximum number of steps in self-consistent
charge equilibration process (CEP) required by the
CENT approach.

    default: ``200``

**alphax_q**: (real) Stepsize in CEP optimization.

    default: ``1.d0``

**qgnrmtol**: (real) Convergence criterion for the norm of the gradient
in CEP when **syslinsolver**  is set top  ``operator``.
Tight values, e.g. ``1.d-5``, may be required for accurate forces.

    default: ``5.d-4``

**freeBC_direct**: (logical) Selectively selects the
``direct`` method for the **syslinsolver** 
for cluster structures only 
during the CEP. This is important when reference (training
or validation) data set contains both
periodic and molecular configurations.

   default: ``False``
      

extra files required by ann train
======================================

**list_posinp_train.yaml**: Contains a list of files in which
configurations of the ANN train task are stored.
The file must contain the key ``files`` whose value
is a list of YAML structure files containing the structures.
For more information on YAML structure file format, see :ref:`yamlstructure`.

**list_posinp_valid.yaml**: Contains a list of files in which
configurations for the validation of the ANN train task are stored.
The file must contain the key ``files`` whose value
is a list of YAML structure files containing the structures.
For more information on YAML structure file format, see :ref:`yamlstructure`.

**SE.ann.input.yaml**: SE is the symbol of the element.
This file must exist for each type of atom in this system.
It contains element based parameters of ANN potential as well as
parameters of the symmetry functions.
It has two keys, ``main`` and ``symfunc``.
All paramters are in atomic units.

``main`` parameters
--------------------------------

**nodes**: Determines the architecture of the ANN.
It is a list of integers, for example [3,5] means
the ANN will have two hidden layers with
three and five nodes, respectively.
Currently, architectures with only one or two hidden
layers are implemented where the latter is well tested
and employed in several different applications.

   default: ``No default value.``

**rcut**: (real) The cutoff radius used for the symmetry functions.

   default: ``No default value.``

**ampl_chi**: (real) Determines the amplitude of the
hyperbolic tangent function used to map the value of ANN output node
to the atomic electronegativity.
We recommend to use 1.0 and smaller values are strongly discouraged.

   default: ``No default value.``

**prefactor_chi**: (real) Determines the prefactor of the argument of the
hyperbolic tangent function used to map the value of ANN output node
to the atomic electronegativity.
We recommend to use 1.0.

   default: ``No default value.``

**ener_ref**: (real) Determines the reference values for energy.
We recommend to set it to the energy of an isolated atom so that
ANN trains indeed the formation energies.

   default: ``No default value.``

**gausswidth**: (real) Determines the width of the Gaussian function
representing the atomic charge density in CENT.
We recommend you to try different values in the range
between 1.0 and 3.0 which most atomic radii are.

   default: ``No default value.``

**hardness**: (real) Determines the atomic hardness by which
you can control how much charge, approximately, you expect to be
collected by this type of atom.
Similar to ``gausswidth``, we recommend you to fine tune for
an optimal value around physically meaningful values given
in textbooks.

   default: ``No default value.``

**chi0**: (real) Determines the offset for atomic electronegativity.

   default: ``No default value.``

**method**: (string) Determines the type of symmetry function
to be used as the atomic environment descriptor.
Currently, only symmetry functions of type ``behler`` are functional.

   default: ``No default value.``

   options: 
   
      ``behler``: For more details see Ref. :cite:`Behler2011`

``symfunc`` parameters
--------------------------------
Two types of symmetry functions, one radial and one angular, are implemented FLAME.
A complete description and comparison between these symmetry functions is
given in :cite:`Behler2011`.
The radial symmetry function in FLAME is `g02` which has two parameters,
the exponent to control the broadness of the Gaussian function and
the offset that determines the center of the Gaussian function.
The offset is not well tested and we recommend to set it to zero.
A `g02` symmetry functions is defined by key `g02_` appended by
a zero-padded enumeration, e.g. `g02_001`.
The value of the key is exponent, offset, 0.0, 0.0, and atom type
all separated by spaces.
The two zeros are the lower and upper bounds of the symmetry function
for all training data set.
The zeros cannot be ignored.
The last item is the atom type of the surrounding atom.
The angular symmetry function is of type `g05` which contains three parameters.
The exponent of the Gaussian function is similar to that of `g02`.
The next two parameters are the prefactor of the cosine function
and the value of the power.
Similar to `g02`, the parameters must be followed by two zero.
In the last, two atom types indicating the types of the two atoms
surrounding the center atom.
An example of an ANN parameter file for type Na in sodium chloride
system is given below:

  .. literalinclude:: Na.ann.input.yaml

