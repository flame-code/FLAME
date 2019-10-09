.. _geopt:

==================================
Optimizer
==================================

The geometry optimizer parameters can be set in this block.
Some parameters are used by all methods, while others
apply only to specific optimizers.


geopt parameters
======================

general ``geopt`` parameters
-----------------------------

**method**: (string) Optimization method.

    default: ``No default value.``

    options:

        ``sd``: steepest descent method with energy feedback.

        ``sdcg``: steepest descent method followed by conjugate gradient method.

        ``sddiis``: steepest descent method followed by DIIS method.

        ``nlbfgs``: Nocedal's implementation of LBFGS :cite:`Nocedal`.

        ``qbfgs``: Quantum Espresso's version of BFGS for variable cell shapes :cite:`Espresso`.

        ``bfgs``: Ghasemi's implementation of BFGS with preferential moves along soft
        modes.

        ``fire``: The Fast Inertial Relaxation Engine method :cite:`bitzek_structurel_2006`.

**fmaxtol**: (real) Convergence parameter. The opimization will terminate 
as soon as the maximum absolute component of force vector falls below this value. In units of Ha/Bohr.
Must be specified.

    default: ``No default value.``

**alphax**: (real) Sandard step size used in most optimizers. Optimal value depends on the system and method.

    default: ``No default value.``


**lprint**: (logical) Verbosity setting. If true, detailed information at each iteration is printed.

    default: ``False``

**nit**: (integer) Maximum number of iterations.

    default: ``1000``

**dxmax**: (real) Maximum displacement for each atomic component. Units in Bohr. Not implemented in all methods yet.

    default: ``1.d-1``

**funits**: (real) Factor to scale energy and force units. 
Currently only useful for Lennard-Jones potential, 
where the units can be scaled to be comparable to
other potentials.

    default: ``1.d0``

**cellrelax**: (logical) Activates variable cell shape relaxations.
 If ``True``, a variable cell shape optimization is performed.

    default: ``False``


**strfact**: (real) The stress tensor scaled by this factor and treated like forces in the process of optimization.
Only relevant if **cellrelax** is ``True``.

    default: ``1.d2``

**geoext**: (logical) Some atomic simulation packages (e.g., LAMMPS, GULP, etc.) come with their
own implementations of geometry optimizers. If  ``True``, these  extrenal optimizers 
will be used. Not implemented for all extrenal codes.

    default: ``False``

``bfgs`` parameters
---------------------

**condnum**: (real) Predicted condition number of the system.

    default: ``1.d1``

**precaution**: (string) Weighting of the quasi-Hessian. 
In Ghasemi's controlled BFGS, moves are gradually affected by
the quasi-Hessian instead of the initial, diagonal matrix.
The parameter *precaution* controls the rate of this gradual transition.

    default: ``normal``

    options:

        ``high``: A high precaution is employed and the quasi-Hessian is
        applied with a high rate.

        ``normal``: A  normal precaution is employed and the quasi-Hessian is
        applied with a medium rate.

        ``low``: A  normal precaution is employed and the quasi-Hessian is
        applied with a low rate.


``fire`` parameters
---------------------

**dt_start**: (real) Initial time step. Arbitrary units.

    default: ``No default value.``

**dt_min**: (real) Minimal time step. Arbitrary units. 

    default: ``1.d0``

**dt_max**: (real) Maximal time step. Arbitrary units. 

    default: ``8.d1``


``qbfgs`` parameters
---------------------

**qbfgsndim**: (integer) Number of old forces and displacements vectors used in the
PULAY mixing of the residual vectors obtained on the basis
of the inverse hessian matrix given by the BFGS algorithm.
When bfgs_ndim = 1, the standard quasi-Newton BFGS method is
used.

    default: ``1``

**qbfgstri**: (real) Initial ionic displacement in the structural relaxation.

    default: ``5.d-1``

**qbfgstrmin**: (real) Minimum ionic displacement in the structural relaxation
BFGS is reset when trust_radius < trust_radius_min.

    default: ``1.d-3``

**qbfgstrmax**: (real) Maximum ionic displacement in the structural relaxation.

    default: ``8.d-1``

**qbfgsw1**: (real) Parameters used in line search based on the Wolfe conditions.

    default: ``1.d-2``

**qbfgsw2**: (real) Parameters used in line search based on the Wolfe conditions.

    default: ``5.d-1``

