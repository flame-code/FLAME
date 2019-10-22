.. _flame_in:

============================
The ``flame_in.yaml`` file
============================

The contents of *flame_in.yaml* is arranged in labeled blocks.
Python-style indentation is used to indicate nesting.
Blocks contain keyword--parameter pairs
as well as subblocks.

Following top-level blocks are available.

.. hlist::
   :columns: 1

   * :ref:`main <main>`
   * :ref:`potential <potential>`
   * :ref:`fingerprint <fingerprint>`
   * :ref:`geopt <geopt>`
   * :ref:`geopt_prec <geopt_prec>`
   * :ref:`dynamics <dynamics>`
   * :ref:`ann <ann>`
   * :ref:`single_point <single_point>`
   * :ref:`saddle <saddle>`
   * :ref:`saddle_opt <saddle_opt>`
   * :ref:`minhocao <minhocao>`

Only the *main* block must be present in every FLAME run.

