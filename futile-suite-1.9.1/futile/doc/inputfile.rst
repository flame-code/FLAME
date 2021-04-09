Automatic definition of the input dictionary the :f:mod:`f_input_file` module
=============================================================================

In a Fortran program it is often useful to have a high-level handling of the input
variables. The Fortran specification provides the ``namelist`` approach, which is
relatively easy to use. However such an approach presents few drawbacks of portability
and code-intrusivity, which limit its usage in a multi-language context.
With this module we provide a set of rules to automatically parse, inspect, verify and
convert a input file from the ``yaml`` format into a ``futile`` dictionary.
With the usage of the :f:mod:`f_input_file` module developers might easily
write _specifications_ for their program input file.

Such module is in tight connection with its python counterpart, :py:mod:`futile.Inputvars`, which
employs the same conventions and make possible a full interplay between a ``yaml`` input file and
a python dictionary.

.. f:automodule:: f_input_file

