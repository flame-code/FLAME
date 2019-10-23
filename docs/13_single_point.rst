.. _single_point:

===========
Singlepoint
===========

The single point task performs a sole energy and force evaluation
given a ``posinp.yaml`` atomic structure file without changing the atomic positions.
This task is especially useful to establish  
a simple interface between an external 
sampling code and FLAME,
which require energy/force evaluation given a atomic 
structure file. 
A prominent example for such an application is
the computation of force sets of displaced 
supercell structures in phonon calculations, e.g., using Phonopy
or ShengBTE.
Further, this task allows to use a potential implemented in 
the FLAME code as a black-box engine, either by directly 
linking to the flame library,
or by using the socket i-Pi interface :cite:`ceriotti_i-pi_2014`.

``single_point`` parameters
=================================



**print_force**: (logical) The forces will be written into a separate output file if set to ``True``.

    default: ``False``

**format**: (string) Fortran format string used to print the forces.

    default: ``No default value.``


**usesocket**: (logical) Activates the communication scheme over unix or TCP sockets. The i-Pi protocol is used.
Note: the server side of the socket communication has to be initialized first before 
FLAME can connect to the socket. A *posinp.yaml* must be provided that coincides
in terms of composition and ordering of the elements with the system run at the i-Pi host.

    default: ``False``

**sockinet**: (integer) Selects Unix socket or internet (TCP) socket.

    default: ``0``

    options:
        
        ``0``: Unix socket

        ``1``: internet (TCP) socket

**sockport**: (integer) Socket port number.

   default: ``0``


**sockhost**: (string) Socket address. If **sockinet** is ``0``, a string with the **sockhost** name will be
created in a temporary directory. Otherwise, a valid IP address must be provided (`127.0.0.1` for localhost).

    default: ``0``

