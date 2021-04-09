indeControl of the precision :f:mod:`f_precision`
=================================================

Work in progress: here the automodule documentation of f_precision should appear

.. f:automodule:: f_precisions


Now I write some customized information here: 

.. f:function:: f_loc(x)

      Inquiry Function which identifies the address of the scalar object
      associated to a unknown quantity. 
      Portable version of the :f:func:`loc()` function of GNU fortran standard

      :p x[in]: 
 
	 Variable of any type; it may be a subroutine or a fortran intrinsic or derived type. It cannot be a nullified pointer 
	 
      :r integer(f_address) f_loc: Address of the object

      :Example:
	 .. code-block:: fortran

	    program test_f_loc
	       use f_precisions, only: f_address,f_loc
	       implicit none
               integer(f_address) :: i
	       !note that altough this variable is uninitialized it 
               !already has an address in the stack
               real :: r 
      	    
               i = f_loc(r)
               print *, i
            end program test_f_loc

       .. todo::

	  Some text 2


.. f:currentmodule:: 
