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

**subtask**: If task=ann, then different subtasks are available,
e.g. ``subtask: train`` invokes routines to perform an ANN training.
There is no default for this key.

.. _ref-ann-approach:

**approach**: If ``atombased`` then, the total energy is the summation of
atomic energies. If ``cent1`` then, the CENT approach as explained
in :cite:`ghasemi_interatomic_2015`.

    default: ``atombased``

**optimizer**: The method to be used for the optimization during a training process.
There are several methods, however, only
``optimizer: rivals`` is well tested.
It is a modification of extended Kalman filter.
There is no default for this key.

**nstep_opt**: Number of extended Kalman filter steps

    default: ``100``

**nstep_cep**: Maximum number of steps in self-consistent
charge equilibration process (CEP) required by the
CENT approach.

    default: ``200``

**alphax_q**: Stepsize in CEP optimization.

    default: ``1.0``

**nconf_rmse**: Number of configurations for which energy
and force RMSE values will be calculated and reported in the
output.

    default: ``0``

**ampl_rand**: The amplitude of random numbers used
as the initial values for ANN weights.

    default: ``1.0``

**symfunc**: It determines whether symmetry functions
values to be saved in memory (``only_calculate``),
to be written into files (``write``) to be uesd later,
to be read from files (``read``) that are written by a previous run.

    default: ``only_calculate``

**syslinsolver**: It determines whether a ``direct``
method should be employed for CEP or an iterative
(``operator``).

    default: ``direct``

**qgnrmtol**: The convergence criterion be applied to the norm of gradient
in CEP when an iterative method is invoked.
Tight values, e.g. 1.E-5, may be required for accurate forces.

    default: ``5.E-4``

**freeBC_direct**: ``True`` or ``False``.
If ``True``, then for task equal to train,
CEP is performed using a ``direct`` method
for cluster structures.
This is important when reference (training
or validation) data set contains both
periodic and molecular configurations.

    default: ``False``

