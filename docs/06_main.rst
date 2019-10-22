.. _main:

=============
Main Block
=============

The options available in the **main** block of the *flame_in.yaml* file determines the 
overall task to be performed in FLAME and sets up the 
atomistic simulation environment.

**task**: (string) Determines the main FLAME task. 

   default: ``No default value.``

   options:


        ``geopt``: Local geometry optimization. The parameters of this
        task can be set within the block :ref:`geopt <geopt>`.

        ``saddle``: Saddle point search.
        The parameters of this task can be set within the
        block :ref:`saddle <saddle>`.

        ``dynamics``: Molecular dynamics simulation.
        The parameters of this task can be set within the
        block :ref:`dynamics <dynamics>`.

        ``ann``: All tasks related to the training and evaluation  
        of artificial neural network potentials. Subtasks can be selected in
        the block :ref:`ann <ann>`.

        ``single_point``: Perform a single point calculation
        to evaluate the energy, forces, and stresses for one or more configurations.
        Allows linking to external sampling codes. The parameters for this taks
        can be specified in the block :ref:`single_point <single_point>`.

        ``minhocao``:  Perform a global optimization calculation based on the
        minima hopping method. The parameters for this task can be set within
        the block :ref:`minhocao <minhocao>`.
        
..        ``conf_comp``: Compare atomic structures to determine their similarities and
..        dissimilarities. The parameters of this task can be set within the
..        block :ref:`conf_comp <conf_comp>`.
..        ``minhopp``: Perform a global optimization calculation based on the
..        minima hopping method. The parameters of this task can be set within
..        the block :ref:`minhopp <minhopp>`.
..        ``genconf``: Generate configurations based on
..        subtask chosen in the block :ref:`genconf <genconf>`.


**two_level_geopt**: (logical) Determines if geometry optimizations
are performed with two accuracy levels.
If ``True``, the block **geopt_prec** must
be present in *flame_in.yaml*.

    default: ``False``

**verbosity**: (integer) Verbosity of output data.

    default: ``0``

    options: ``0, 1, 2, 3`` Increasing number for higher verbosity.


**verbosity_mh**: (integer) Verbosity of the output related to the task ``minhocao``.

    default: ``3``

    options: ``0, 1, 2, 3`` Increasing number for higher verbosity.


**nat**: (integer) Number of atoms involved in the simulation.

   default: ``0``

**types**: (list of strings of length [number of atomic types]) 
Character of the atoms involved in the simulation.
Internally, the atomic types are enumerated, starting from 1.

   default: ``No default value.``

**znucl**: (list of integers of length [number of atomic types]) 
The atomic number of the atoms involved in the simulation,
charge of the nuclei.
Corresponds to an equivalent input of the **types** keyword.

   default: ``No default value.``

**amass**: (list of reals of length [number of atomic types]) 
Overrides the physical atomic masses in dynamics simulations.
Particularly useful to reduce the spectrum of the
vibrational eigenfrequencies in MD escape trials during
minima hopping runs.

   default: ``No default value.``

**typat**: (list of integers of length **nat**) Indexes
of the atomic types from 1 to number of atomic types.

   default: ``No default value.``

**findsym**: (logical) Activates symmetry detection of crystalline solids.
FLAME must be compiled and linked with SPGLIB.

    default: ``False``

**rng_type**:  (string) Not used for production runs, only used for development
and regression tests.

**seed**: (integer) Seed value to initialize the random number generator.

   default: ``-2``

**pressure**: (real) External pressure. In units of GPa

   default: ``0.d0``

**params_new**: (logical) Enables parsing the ``params_new.in`` for
the **task** ``minhocao`` runs. This functionality is available to allow backwards
compatibility with earlier versions of ``minhocao``.
This option will be removed in the future.

    default: ``False``

..   nrun_lammps                         : 0
..   nat                                 : 0
..   pressure                            : 0.0
..   findsym                             : False
..   finddos                             : False
..   params_new                          : False
..   verbosity                           : 0
..   verbose                             : 1
..   rng_type                            : only_for_tests
..   seed                                : -2
..   finddos                             : False
..   params_new                          : False
