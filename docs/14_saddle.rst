.. _saddle:

===========
Saddle
===========

One and two-sided saddle point search methods 
are implemented in FLAME. Note that not all 
parameters apply to 
all available saddle point search methods
when **task** is set to ``saddle``.


saddle parameters
=================================

general ``saddle`` parameters
------------------------------------------

**method**: (string) Method for the saddle point search.

    default: ``No default value.``

    options:
        
         ``bar_saddle``:  The bar-saddle method from Ref. :cite:`schaefer_minima_2014`.
    
         ``dimer``: The dimer method from Ref. :cite:`henkelman_dimer_1999`.
    
         ``splined_saddle``: The splined saddle method from Ref. :cite:`ghasemi_an_2011`.
    
``bar_saddle`` parameters
-----------------------------------

**dbar**: (real) The length of the bar/dimer in the bar-saddle method. In units of Bohr in configurational space.

   default: ``1.d-1``

**nstep**: (integer) Maximal number of steps for the saddle search.

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
reached to commence the contraction of the bar. Units in Ha/Bohr.

   default: ``1.d-2``

**fnrmtol_contracted**: (real) The force norm tolerance 
for the convergence of the final, contracted bar. Units in Ha/Bohr.

   default: ``5.d-4``


``dimer`` parameters
-----------------------------------

**list_random_displace**: (list of integers) List of atomic indices that are 
displaced randomly before initiating
the saddle optimization in the dimer method.

   default: ``No default value.``

**ampl**: (real) The amplitude of the random displacement applied to 
the atoms listed by the key **list_random_displace**. In units of Bohr in configurational space.

   default: ``No default value.``

**dimsep**: (real) The dimer separation in the dimer method. 
No default value given, must be set when the dimer method is invoked.
In units of Bohr in configurational space.

   default: ``No default value.``


``splined_saddle`` parameters
-----------------------------------

**np_splsad**: (integer) Number of anchor points.
**np_splsad**-1 is the number of moving anchor points in the splined saddle method.
Note: with the default value of ``2``, the splined saddle will only have one moving achor point.

   default: ``2``

**np_neb**: (integer) Number of NEB images. 
**np_neb**-1 is the number of moving images in the nudged elastic band (NEB) method.
Note: with the default value of ``2``, the NEB will run with only one moving image.

   default: ``2``

**ns2**: (integer) Number of extra points added to the number of anchor points in order to
assist the maximization process, the inner loop of the splined saddle method.

   default: ``0``

**vdtol**: (real) The convergence criterion in the maximization process in
the inner loop of the splined saddle method. Units in Ha.

   default: ``1.d-1``

**htol**: (real) The smallest value of the normalized pathway parameter between
two neighboring points in the maximization process of the splined saddle method.
Arbitrary units.


   default: ``2.d-2``

**alphax**: (real) The step size of the optimizer in the NEB and splined saddle methods.
Arbitrary units.

   default: ``5.d-1``

**docineb**: (string) Activates the climbing image NEB.

   default: ``no``

   options: 

         ``no``: NEB without climbing image

         ``yes``: NEB with climbing image

**doneb**: (string) Activates NEB calculations.

   default: ``No default value.``

   options: 

         ``no``: No NEB calcualtion is performed

         ``yes``: NEB calculation is performed

**pickbestanchorpoints**: (string) Activates an automated selection of favored anchor points. 
This feature is not well tested and we recommend setting this parameter to ``no``.

   default: ``No default value.``

   options: 

         ``no``: anchor points are distributed uniformly in the beginning of simulation

         ``yes``: anchor points are initially selected to favor higher energy points based on estimates obtained by an interpolation.


**runstat**: (string) Determines whether a new or a restart run is performed.

   default: ``new``

   options:

         ``new``: New run, the NEB images or splined saddle anchor points are initialized at the beginning of the run.

         ``restart``: Restart run, the NEB images or splined saddle anchor points are read from a file. Restart runs have not yet been well tested,
         so we currently do not recommend using it.

**typintpol**: (string) The type of interpolation in the maximization process of the splined saddle method.

   default: ``cubic``

   options:

         ``cubic``: Natural cubic splines.

         ``quintic``: A spline using fifth-order polynomials. Rather unstable except for simple pathways.

**fcalls_max**: (integer) The maximum number of calls to force evaluation.

   default: ``100``

**fmaxtol_splsad**: (real) The convergence criterion for the saddle optimization
in the splined saddle method. Units in Ha/Bohr.

   default: ``2.d-4``

**fmaxtol_neb**: (real) The convergence criterion for the saddle optimization
in the NEB method. Units in Ha/Bohr.

   default: ``2.d-2``

**opt_method**: The optimization method used in the saddle point search
when using NEB or the splined saddle method.

   default: ``SD``

   options:

         ``SD``: The steepest descent method.

         ``BFGS``: The Broyden–Fletcher–Goldfarb–Shanno (BFGS) method.

