.. _ann:

========================================
Artificial Neural Network Potentials
========================================

Artificial neural network (ANN) potentials are
a class of machine learning interatomic potentials.
In high dimensional ANN potentials, the total
energy is decomposed into atomic energies.
The output values of ANNs are considered as atomic energies.
In order to use such an approach in FLAME,
see below :ref:`approach <ref-ann-approach>`.
The equilibration via neural network
technique (CENT) is another method implemented
in the FLAME code.
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

.. _ref-ann-approach:


**approach**: (string) Determines the ANN technuque to be employed.

   default: ``atombased``
    
   options: 

      ``atombased``: total energy is the sum of atomic energies
       
      ``cent1``: charge equilibration via neural network :cite:`ghasemi_interatomic_2015`

``train`` parameters
--------------------------------

**optimizer**: (string) Method for the optimization during the ANN training process.


   default: ``No default value.``

   options: 
   
      ``rivals`` Well tested and stable optimizer based on a modification of extended Kalman filter.

**nstep_opt**: (integer) Number of extended Kalman filter steps

    default: ``100``


**nconf_rmse**: (integer) Number of configurations for which energy
and force RMSE values will be calculated and reported in the
output.

    default: ``0``

**ampl_rand**: (real) The amplitude of random numbers used
as the initial values for ANN weights.

    default: ``1.d0``

**symfunc**: (string) Determines how the evaluation of the symmetry functions
is treated.

   default: ``only_calculate``

      options: ``only_calculate`` values are stored in memory only

      options: ``write`` values are written to disk

      options: ``read``  values are read from disk from a previous run



``cent1`` parameters
--------------------------------
**syslinsolver**: (string) Determines what method to use
to solve the system of linear equations for the CEP.

   default: ``direct``

   options: 
   
      ``direct`` direct approach, e.g., Gaussian elemination

      ``operator`` iterative approach

**nstep_cep**: (integer) Maximum number of steps in self-consistent
charge equilibration process (CEP) required by the
CENT approach.

    default: ``200``

**alphax_q**: (real) Stepsize in CEP optimization.

    default: ``1.d0``

**qgnrmtol**: (real) Convergence criterion for the norm of the gradient
in CEP when **syslinsolver**  is ``operator``.
Tight values, e.g. ``1.d-5``, may be required for accurate forces.

    default: ``5.d-4``

**freeBC_direct**: (logical) Selectively selects the
``direct`` method for the **syslinsolver** 
for cluster structures only 
during the CEP. This is important when reference (training
or validation) data set contains both
periodic and molecular configurations.

   default: ``False``
      
