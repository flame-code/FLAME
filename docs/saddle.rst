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

**list_random_displace**: (...)

   default: ``No default value.``

   options: 

**dimsep**: (real)                           

   default: ``No default value.``

**ampl**: (real)

   default: ``No default value.``

**np_splsad**: (integer)

   default: ``2``

**np_neb**: (integer)

   default: ``2``

**ns2**: (integer)

   default: ``0``

**vdtol**: (real)

   default: ``1.d-1``

**dt**: (real) 

   default: ``3.d-2``

**htol**:  (real)

   default: ``2.d-2``

**alphax**: (real)

   default: ``5.d-1``

**hybrid**: ()

   default: no

**docineb**: ()

   default: no

**doneb**: ()

   default: unknown

**pickbestanchorpoints**: ()

   default: unknown

**runstat**: ()

   default: new

**typintpol**: (string)

   default: ``cubic``

**fcalls_max**: (integer)

   default: ``100``

**fmaxtol_splsad**: (real)

   default: ``2.d-4``

**fmaxtol_neb**: (real)

   default: ``2.d-2``

**opt_method**: (string) 

   default: ``SD``

   options:

      ``SD``
