
=================================
Atomic Structure Files
=================================

FLAME can handle various atomic structure
file formats. The preferred format is ``yaml``,
but other formats are supported as well.
The subdirectory ``utils/python/`` contains
scripts to convert structural data from a wealth of input formats
to a ``yaml`` file.


.. _yamlstructure:

``yaml`` Format
-----------------------------

Keyword-based structural information format.
Must start with a block with keyword **conf**, followed by
the following parameters.

**bc**: (string) Boundary condition.

   default: ``No default value.``

   options:

       ``free``: Free boundary conditions for molecules and clusters

       ``slab``: Slab boundary conditions for surfaces, interfaces, and 2D materials
        
       ``bulk``: Periodic boundary conditions for crystals and solids
    
**cell**: (list of three lists) Three cell vectors.

   default: ``No default value.``

**nat**: (int)  Number of atoms.

   default: ``No default value.``

**units_length**: (string) Length units used in the structure.

   default: ``atomic``

   options:

      ``angstrom``: units in Angstrom

      ``atomic``: units in Bohr

**coord**: (list of ``nat`` lists)  For every atom in the system,
a list of the form
``[x, y, z, AT, CCC]`` must be given. ``x, y, z`` are the cartesian 
coordinates of  the atom, while ``AT`` is the symbol of the particular atomic type,
and ``CCC`` determines if its ``x, y, z`` degree of freedom are allowed to
change during the simulation.


E.g., ``[0., 0., 0., Si, TTT]`` is a silicon atom at the origin, free to move in all dimensions,
while ``[1., 1., 1., C, TFF]`` is a carbon atom at ``(1, 1, 1)``, only allowed to move along the ``x``
direction.

Example
**********************************

Calcium fluoride in a crystal lattice.

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



``ascii``  Format
-----------------------------

The ascii format is described here in detail: http://inac.cea.fr/sp2m/L_Sim/V_Sim/sample.html#sample_ascii
It is the native format of the **task** ``minhocao`` and the only format 
that currently supports the constraint of individual lattice parameters.
Periodic cell vectors have to be provided to describe the simulation cell.
The cell consisting of the cell vectors ``a, b, c`` 
has to be rotated such that ``a`` points along the ``x`` direction, 
and ``b`` lies within the ``x-y`` plane.
See below figure for the orientation of the cell
and the projections ``dxx``, ``dyx``, ``dyy``, ``dzx``, ``dzy``, and ``dzz``.


.. image:: repere-ascii.png
   :scale: 50 %
   :align: center

**line 1**: (integer) Number of atoms in the system, ``Nat``.

**line 2**: (real, real, real) ``dxx`` ``dyx`` ``dyy`` values

**line 3**: (real, real, real) ``dzx`` ``dzy`` ``dzz`` values

**lines 4 -- n**: Only lines starting with ``#keyword:`` are allowed 
and read/interpreted accordingly.
The ``#keyword:`` tag must be followed by the following options.

   ``reduced``: the atomic coordinates are treated as reduced coordinates, i.e., with respect to the cell vectors.

   ``fixlat a b c alpha beta gamma vol``: specific components of the lattice vectors can be fixed. 
   A ``True`` (or ``T``) value
   fixes the degree of freedom, while ``False`` (or ``F``) does not fix it.
   The last parameter (``vol``) is used for fixed-cell-shape-variable-volume simulations.
   I.e., ``#keyword: fixlat  F F T T T F F`` allows the length of the ``a`` and ``b`` vectors to change, while keeping the 
   length of the ``c`` vector constant. At the same time, only the angle between ``a`` and ``b`` is allowed to change.
   This setting is particularly useful for simulations of 2D materials or surfaces.

**lines n+1 -- Nat**: ``Nat`` lines with the coordinates of every atom of the form ``x, y, z, AT, c``.
The ``x, y, z`` coordintes, followed by the chemical symbol ``AT``, and optionally ``f`` for fixed atom.
Note that the reduced coordinates will be fixed instead of the Cartesian one if the
system is periodic

Example
**********************************

Calcium fluoride in a crystal lattice, with selectively fixed lattice parameters.
The Ca atoms are not allowed to move.

.. code-block:: none
      
   12
   5.4620E+00   0.0000E+00   5.4620E+00
   0.0000E+00   0.0000E+00   5.4620E+00
   #keywords: fixlat F F T T T F F 
   0.0000E+00   0.0000E+00   0.0000E+00   Ca f
   2.7310E+00   2.7310E+00   0.0000E+00   Ca f
   2.7310E+00   0.0000E+00   2.7310E+00   Ca f
   0.0000E+00   2.7310E+00   2.7310E+00   Ca f
   1.3655E+00   1.3655E+00   1.3655E+00    F 
   4.0965E+00   1.3655E+00   1.3655E+00    F 
   1.3655E+00   4.0965E+00   1.3655E+00    F 
   4.0965E+00   4.0965E+00   1.3655E+00    F 
   1.3655E+00   1.3655E+00   4.0965E+00    F 
   4.0965E+00   1.3655E+00   4.0965E+00    F 
   1.3655E+00   4.0965E+00   4.0965E+00    F 
   4.0965E+00   4.0965E+00   4.0965E+00    F 



