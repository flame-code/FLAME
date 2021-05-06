Mapping and Lists: the :f:mod:`dictionaries` module
===================================================

.. f:currentmodule::

.. f:currentmodule:: dictionaries

.. f:type:: dictionary
   
    Opaque type of the :f:mod:`dictionaries` module, that is used to serialize and store data information.
    When using it as a entity in a given scoping unit, it has to be declared with the `pointer` attribute, as it 
    refers to internal storage data, not directly accessible.
    For this reason its status is undefined after declaration. Each dictionary has to be therefore iniitalized via 
    calls to :f:subr:`dict_init` function, of by association to the :f:func:`dict_new` function.


.. list-table:: Public procedures
   :widths: auto
   :header-rows: 1

   * - Constructors/Destructors:  
     - Accessors:	       
   * - :f:subr:`dict_init`
     - :f:func:`dict_len`     
   * - :f:func:`dict_new`
     - :f:func:`dict_size`    
   * - :f:func:`dict_free`	      
     - :f:func:`dict_key`     
   * - :f:func:`list_new`
     - :f:func:`dict_item` 

.. f:function:: dict_free(dict [,dictA, dictB, ...])

    Free and nullify the dictionary. Optionally accepts up to ten different dictionaries as arguments.

    :p dictionary dict [inout,pointer]: 

       The dictionary to be freed. Nullified on output. The same applied to the optional arguments :f:var:`dict1`, :f:var:`dict2`, etc.
	 
    :Example:
	 .. code-block:: fortran
         
  	    program freedict
	      use dictionaries
	      implicit none
              type(dictionary), pointer :: dict1, dict2
	      [...]
	      call dict_free(dict1)
	      call dict_free(dict2)
	      !or instead
	      !call dict_free(dict1,dict2)
            end program freedict

.. f:autosubroutine:: dict_init

.. f:function:: list_new([[entity-list] [entity,  [,entity] ...]])

    Prepare storage data in the dictionries database ans provides a pointer to it. The function f:func:`list_new` conceives the database as a list.
    Optionally a comma-separated list of entities can be provided, or, alternatively an array of them. 
    If more thant one argument is present, or if an array is provided as a single argument, 
    the stored dictionary should be interpreted as a list, in arguments or array element order, respectively.

    Each `entity` is described as a 

    `.item.`  value

    The variable `value` may be any scalar object of intrinsic type, or a :f:type:`dictionary` `pointer` variable. The variable might also be a fortran
    array of intrinsic type and rank one, in which case :f:func:`list_new` returns a pointer to a list whose items correspond to `value` elements in array element order.

    :r dictionary list_new [inout,pointer]: reference to the stored data

    :Example:

       .. literalinclude:: /../tests/flib/dicts.f90
           :language: fortran
           :start-after: example of a list
	   :end-before: end of example of a list
       .. code-block:: fortran

            use dictionaries
	    real(f_double) :: var
  	    real(f_double), dimension(3) :: arr
	    type(dictionary), pointer :: list
	    ![...]
	    list => list_new(.item. arr)
	    !retrieval
	    val = list // 0 !this corresponds to val=arr(1)
	    val = list // 1 !this corresponds to val=arr(2)
	    val = list // 2 !this corresponds to val=arr(3)
	    call dict_free(list)

.. f:function:: dict_new([[entity-list] [entity,  [,entity] ...]])

    Prepare storage data in the dictionries database ans provides a pointer to it.
    When called without arguments the association to the :f:func:`dict_new` function is identical to a :f:subr:`dict_init` routine.
    Optionally a comma-separated list of entities can be provided, or, alternatively an array of them. 
    If more thant one argument is present, or if an array is provided as a single argument, 
    the stored dictionary should be interpreted as a ordered mapping, in arguments or array element order, respectively.

    Each `entity` is described as a 

    key `.is.` value

    where key is a string. 
    The variable value may be any scalar object of intrinsic type, or a :f:type:`dictionary` `pointer` variable.

    :r dictionary dict_new [inout,pointer]: reference to the stored data

    :Example:

       .. literalinclude:: /../tests/flib/dicts.f90
           :language: fortran
           :start-after: example of dict_new usage1
	   :end-before: end of example of dict_new usage1
       .. literalinclude:: /../tests/flib/dicts.f90
           :language: fortran
           :start-after: perform an iterator on a scalar
	   :end-before: end of perform an iterator
       .. literalinclude:: /../tests/flib/dicts.f90
           :language: fortran
           :start-after: new test, build dictionary
	   :end-before: end of new test, build dictionary
       .. literalinclude:: /../tests/flib/dicts.f90
           :language: fortran
           :start-after: another comprehensive test
	   :end-before: end of another comprehensive test


.. f:autofunction:: dict_len

.. f:autofunction:: dict_size

.. f:autofunction:: dict_key

.. f:autofunction:: dict_item

.. f:autofunction:: dict_value
 

and again

.. f:currentmodule:: 


and yet again

.. f:autofunction:: f_err_check
