.. _flame_in:

=========================
Contents of flame_in.yaml
=========================

Contents of *flame_in.yaml* are arranged in blocks and in each
block a set of parameters can be assigned the desired value. 
The label of each block appears in square bracket and the
block continues to the begining of next block.
Each parameter must be followed by an equal sign.
Notice that flame_in.yaml is case sensitive.

The list of blocks is as following:

.. hlist::
   :columns: 1

   * :ref:`main <main>`
   * :ref:`minhopp <minhopp>`
   * :ref:`geopt <geopt>`
   * :ref:`potential <potential>`
   * :ref:`dynamics <dynamics>`
   * :ref:`genconf <genconf>`
   * :ref:`conf_comp <conf_comp>`
   * :ref:`testforces <testforces>`
   * :ref:`ann <ann>`
   * :ref:`saddle_1s <saddle_1s>`
   * :ref:`geopt_prec <geopt_prec>`
   * :ref:`saddle_1s_opt <saddle_1s_opt>`
   * :ref:`minhocao <minhocao>`

Among all these block only block main must be present and the
rest are optional.

