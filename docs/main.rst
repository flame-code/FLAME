.. _main:

====
Main
====

Options of section **main** in flame_in.yaml file.

**task**: The task which is supposed to be executed by Alborz.

    default: There is no default value. One must select one from the list below.

    options:

        ``minhopp``: It performs a global optimization calculation based on the
        minima hopping method. The parameters of this task can be set within
        the block :ref:`minhopp <minhopp>`.

        ``geopt``: It performs a geometry optimization run. The parameters of this
        task can be set within the block :ref:`geopt <geopt>`.

        ``saddle_1s``: It performs a saddle point search calculation.
        The parameters of this task can be set within the
        block :ref:`saddle_1s <saddle_1s>`.

        ``dynamics``: It performs a Molecular dynamics simulation.
        The parameters of this task can be set within the
        block :ref:`dynamics <dynamics>`.

        ``conf_comp``: It compares configurations for similarities and
        dissimilarities. The parameters of this task can be set within the
        block :ref:`conf_comp <conf_comp>`.

        ``ann``: This task can perform several type of calculations
        relevant to Artificial Neural Network based on subtask chosen in
        the block :ref:`ann <ann>`.

        ``genconf``: This task can generate configurations based on
        subtask chosen in the block :ref:`genconf <genconf>`.

        ``testforces``: It examines atomic forces based on selected method
        which can be set :ref:`testforces <testforces>`.

        ``single_point``: It performs single point calculation
        (energy and forces) for one or more configurations.

        ``misc``: It solves poisson's equation for a charge density given
        in cube format.

**two_level_geopt**: It specifies the geometry optimization should
be in level or two. It it is .true., the block geopt_prec must be
provided in addition to block geopt.

    default: ``.false.``

**verbosity**: The verbosity of information printed in screen.
It is an integer value which the smaller value means more information
on screen. The smallest possible is 0 and the largest value may vary
in future.

    default: ``0``


