.. _structures:

===========
Structures
===========

FLAME can handle various forms of atomic structure
file formats. The preferred format is ``yaml``,
but other formats are supported as well.
The subdirectory ``utils/python/`` contains
scripts to convert structural data from wealth input formats
to a ``yaml`` file.


``yaml`` structure format
=================================

Keyword-based structural information format.
Must start with the keyword **conf**, followed by
the following parameters.

**bc**: (string) Boundary condition.

   default: ``No default value.``

   options:

       ``free``: Free boundary conditions for molecules and clusters.

       ``slab``: Slab boundary conditions for surfaces or interfaces.
        
       ``bulk``: Periodic boundary conditions for crystals and solids.
    
**cell**: (list of 3 lists) Three cell vectors.

   default: ``No default value.``

**nat**: (int)  Number of atoms.

   default: ``No default value.``

**units_length**: (string) Length units used in the structure.

   default: ``angstrom``

   options:

      ``angstrom```: units in Angstrom

      ``atomic```: units in Bohr

**coord**: (list of lists)  For every atom in the system,
a list of the form
``[x, y, z, AT, CCC]`` must be given. ``x, y, z`` are the cartesian 
coordinates of  the atom, while ``AT`` is the symbol of the particular atomic type,
and ``CCC`` determines if its ``x, y, z`` degree of freedom are allowed to
change during the simulation.


E.g., ``[0., 0., 0., Si, TTT]`` is a silicon atom at the origin, free to move in all dimensions,
while ``[1., 1., 1., C, TFF]`` is a carbon atomot at ``(1, 1, 1)`` only allowed to move along the ``x``
direcion.

   default: ``No default value.``

Example: bulk calcium fluoride
----------------------------------
.. code-block:: yaml

    conf:
      nat                                 : 12
      bc                                  : bulk
      units_length                        : angstrom
      cell:
      -  [5.462, 0.0, 0.0]
      -  [0.0, 5.462, 0.0]
      -  [0.0, 0.0, 5.462]
      coord:
      -  [0.0, 0.0, 0.0, Ca, TTT]
      -  [2.731, 2.731, 0.0, Ca, TTT]
      -  [2.731, 0.0, 2.731, Ca, TTT]
      -  [0.0, 2.731, 2.731, Ca, TTT]
      -  [1.3655, 1.3655, 1.3655, F, TTT]
      -  [4.0965, 1.3655, 1.3655, F, TTT]
      -  [1.3655, 4.0965, 1.3655, F, TTT]
      -  [4.0965, 4.0965, 1.3655, F, TTT]
      -  [1.3655, 1.3655, 4.0965, F, TTT]
      -  [4.0965, 1.3655, 4.0965, F, TTT]
      -  [1.3655, 4.0965, 4.0965, F, TTT]
      -  [4.0965, 4.0965, 4.0965, F, TTT]



``ascii`` structure format
=================================



