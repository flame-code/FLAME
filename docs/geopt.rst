.. _geopt:

==================================
Optimizer parameters
==================================

The optimizer parameters can be set in this block.
Some of these parameters are used by all optimizers and
some are only relevant to some of optimizers.

[geopt] options
=================

**method**: The optimization method for geometry optimization.

    default: ``There is no default value.``

    options:

        ``sd``: steepest descent method with energy feedback.

        ``sdcg``: steepest descent method followed by conjugate gradient method.

        ``sddiis``: steepest descent method followed by DIIS method.

        ``nlbfgs``: Nocedal's implementation of LBFGS.

        ``bfgs``: Ghasemi's implementation of BFGS in which the moves along soft
        modes are controlled.

        ``fire``: The FIRE method.

**fmaxtol**: Maximum absolute value component of force vector.
This parameter must be specified.

    default: ``There is no default value.``

**alphax**: The standard step size used in most optimizers. It is system dependent.

    default: ``There is no default value.``

**condnum**: The condition number of the system as expected by the user. This value is
used in Ghasemi's controlled BFGS.

    default: ``10.0``

**precaution**: In Ghasemi's controlled BFGS, moves are gradually affected by
the quasi-Hessian rather than the beginning.
The parameter *precaution* controls the rate of this gradual process.

    default: ``normal``

    options:

        ``high``: A high precaution is employed and the quasi-Hessian affects are
        applied with a high rate.

        ``normal``: A  normal precaution is employed and the quasi-Hessian affects are
        applied with a medium rate.

        ``low``: A  normal precaution is employed and the quasi-Hessian affects are
        applied with a low rate.

**lprint**: A logical variable. If it is set to true, information at each iteration is printed.

    default: ``.false.``

**dt_start**: The starting time step used in the FIRE method.

    default: ``0.005``

**nit**: The maximum number of iterations.

    default: ``1000``

**dxmax**: Maximum displacement for each component. It is not implemented in all methods yet.

    default: ``0.1``

**nsatur**: Number of saturation steps.

    default: ``5``

**cellrelax**: A logical variable. If it is set to true, a variable cell shape optimization is performed.

    default: ``false``

**funits**: The factor which scales some parameters relevant in energy.
Currently only useful for lennard-jones, since it is the only potential with different units.

    default: ``1.0``

