.. _saddle:

===========
Saddle
===========

One and two-sided saddle point search methods 
are implemented in FLAME. Note that not all 
parameters apply to 
all available saddle point serach methods
when **task** is set to ``saddle``.


saddle parameters
=================================

general ``saddle`` parameters
------------------------------------------

**method**: (string) Method for the saddle point search.

    default: ``No default value.``

    options:
        
         ``bar_saddle``:  The bar-saddle method based from Ref.:cite:`schaefer_minima_2014`.
    
         ``dimer```:
    
..  warning:: Options and their descriptions are missing


..  warning:: Sections for other methods need to be added


``bar_saddle`` parameters
-----------------------------------

**dbar**: (real) The legth of the bar/dimer in the bar-saddle method. In units of Bohr in configurational space.

   default: ``1.d-1``

**nstep**: (integer) Maximal number of steps in for the saddle search.

   default: ``1000``

**stepsize**: (real) The initial step size for the saddle point search iterations. Will be adjusted based on a
a feedback mechanism. Arbitratry units.

   default: ``20.d0``

**bar_contract**: (logical) Activates the contraction of the bar size during
the optimization. The final bar size can be set with **dbar_contracted**.
The contraction occurs gradually over **nstep_contract** iterations.

   default: ``True``

**dbar_contracted**: (real) The final bar size after contraction. In units of Bohr in configurational space.

   default: ``1.d-3``

**nstep_contract**: (integer) Number of iterations to contract the bar from the
initial value **dbar** to the final length **dbar_contracted**. 

   default: ``40``

**fnrmtol_coarse**: (real)  The force norm tolerance that has to be
reached to initiallize the contraction of the bar. Units in Ha/Bohr.

   default: ``1.d-2``

**fnrmtol_contracted**: (real) The force norm tolerance 
for the convergence of the final, contracted bar. Units in Ha/Bohr.

   default: ``5.d-4``


..  warning:: Following parameters not yet assigned

Not yet assigned
--------------------

**list_random_displace**: List of atoms to be displaced randomly before initiating
the saddle optimization in dimer method.

   default: ``No default value.``

**dimsep**: The dimer separation in the dimer method. It has no default value
and it must be set if dimer method is invoked.

   default: ``No default value.``

**ampl**: The amplitude of random displacement applied to atoms listed by
the key **list_random_displace**.

   default: ``No default value.``

**np_splsad**: np_splsad-1 is the number of moving anchor points in the splined saddle method.
Having the default value set, the splined saddle will run with only one moving achor point.

   default: ``2``

**np_neb**: np_neb-1 is the number of moving images in the nudged elastic band (NEB) method.
Having the default value set, NEB will run with only one moving image.

   default: ``2``

**ns2**: Number of extra points to be added to the number of anchor points in order to
assist the maximization process, the inner loop in the splined saddle method.

   default: ``0``

**vdtol**: The convergence criterion in the maximization process in
the inner loop in the splined saddle method.

   default: ``1.d-1``

**htol**: The smallest value of the normalized pathway parameter between
two neighboring points in the maximization process of the splined saddle method.

   default: ``2.d-2``

**alphax**: The step size of the optimizer in the NEB and splined saddle methods.

   default: ``5.d-1``

**docineb**: if ``yes``, it does climbing image NEB.

   default: ``no``

**doneb**: if ``yes``, it performs an NEB calculations. No default value so
it must be set by ``yes`` or ``no``.

   default: ``No default value.``

**pickbestanchorpoints**: If ``no``, anchor points are distributed uniformly
in the beginning of simulation. If ``yes``, anchor points are initially
selected to favor higher energy points based on estimates obtained
by an interpolation. This is not well tested and we recommend you to set it
to ``no``.

   default: ``No default value.``

**runstat**: It determines whether it is a new run or a restart of a previous run.

   default: new
   options:

         ``new``: A new run so NEB images or splined saddle anchor points set at the beginning of the run.
         ``restart``: A restart run so NEB images or splined saddle anchor points to be read from a file,
         not tested yet, so we do not recommend it for now.

**typintpol**: The type of interpolation in the maximization process in the splined saddle method.

   default: ``cubic``
   options:

         ``cubic``: Natural cubic splines
         ``quintic``: A spline using fifth-order polynomial. This is unstable except for simple pathways.

**fcalls_max**: The maximum number of calls to force evaluation.

   default: ``100``

**fmaxtol_splsad**: The convergence criterion for the saddle optimization
in the splined saddle method.

   default: ``2.d-4``

**fmaxtol_neb**: The convergence criterion for the saddle optimization
in the NEB method.

   default: ``2.d-2``

**opt_method**: The optimization method used in the saddle point search
when using NEB or the splined saddle method.

   default: ``SD``
   options:

         ``SD``: The steepest descent method.
         ``BFGS``: The Broyden–Fletcher–Goldfarb–Shanno (BFGS) method.

