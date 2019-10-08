.. _main:

=============
Main Block
=============

Options of the **main** block in the *flame_in.yaml* file determines the 
overall task to be performed in FLAME and sets up the 
atomistic simulation environment.

**task**: (task) Determines the main FLAME task. 

   default: ``No default value.``

   options:

        ``minhopp``: Perform a global optimization calculation based on the
        minima hopping method. The parameters of this task can be set within
        the block :ref:`minhopp <minhopp>`.

        ``geopt``: Local geometry optimization. The parameters of this
        task can be set within the block :ref:`geopt <geopt>`.

        ``saddle``: Saddle point search.
        The parameters of this task can be set within the
        block :ref:`saddle <saddle>`.

        ``dynamics``: Molecular dynamics simulation.
        The parameters of this task can be set within the
        block :ref:`dynamics <dynamics>`.

        ``conf_comp``: Compare atomic structures to determine their similarities and
        dissimilarities. The parameters of this task can be set within the
        block :ref:`conf_comp <conf_comp>`.

        ``ann``: All tasks related to the training and evaluation of 
        relevant to artificial neural networks. Subtasks can be selected in
        the block :ref:`ann <ann>`.

        ``genconf``: Generate configurations based on
        subtask chosen in the block :ref:`genconf <genconf>`.

        ``single_point``: Perform a single point calculation
        to evaluate the energy, forces, and stresses for one or more configurations.

        ``minhocao``:  Perform a global optimization calculation based on the
        minima hopping method. The parameters of this task can be set within
        the block :ref:`minhocao <minhocao>`.
        
**two_level_geopt**: (logical) Determines if geometry optimizations
are performed with two accoracy levels.
If ``True``, the block **geopt_prec** must
be present in *flame_in.yaml*.

    default: ``False``

**verbosity**: (integer) Verbosity of output data.

    default: ``0``

    options: ``0, 1, 2, 3`` Increasing number for higher verbosity.

**verbose**: (integer) Verbosity of output data.

    default: ``0``

    options: ``0, 1, 2, 3`` Increasing number for higher verbosity.

**types**: (list of string) Character of the atoms involved in the simulation.

   default: ``No default value.``

**nat**: (integer) Number of atoms involved in the simulation.

   default: ``0``


**findsym**: (logical) Activates symmetry detection of crystalline solids.
FLAME must be compiled and linked with SpgLib.

    default: ``False``

**rng_type**: (string)

   default: ``intrinsic``

**seed**: (integer) Seed value to initialze the random number generator.

   default: ``-2``

**pressure**: (real) External pressure. In units of GPa

   default: ``0.d0``

..   nrun_lammps                         : 0
..   nat                                 : 0
..   pressure                            : 0.0
..   findsym                             : False
..   finddos                             : False
..   params_new                          : False
