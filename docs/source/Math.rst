Math
====

Here we describe the mathematical framework.

Inner product
-------------

All inner products are defined as

.. math::
    
   \langle a | b \rangle = 4 \Delta f \Re \sum_{i} \frac{\tilde{a}(f_i) \tilde{b}^*(f_i)}{S_n(f_i)},

where :math:`\tilde{a}(f)` and :math:`\tilde{b}(f)` are the Fourier transforms of the time-domain signals,
and :math:`\Delta f` is the frequency resolution. The sum is over frequency bins.
There is always a minimum frequency, maximum frequency and frequency resolution. Which can 
be adjusted to the user's preferences in the options of the functions ``SNR`` and ``DALITensors``.

``DALI`` tensors
----------------

The ``DALI`` Likelihood expansion is defined as

.. math::

   \ln L(\theta) = - \frac{1}{2} \langle h_i | h_j \rangle \Delta \theta^{i j} \\
   -\left(
    \frac{1}{2! 1!} \langle h_{ij} | h_k \rangle \Delta \theta^{i j k} - 
    \frac{1}{2 \, 2!^2} \langle h_{ij} | h_{kl} \rangle \Delta \theta^{i j k l} 
   \right) \\ 
   -\left(
    \frac{1}{3! 1!} \langle h_{ijk} | h_l \rangle \Delta \theta^{i j k l} - 
    \frac{1}{3! 2!} \langle h_{ijk} | h_{lm} \rangle \Delta \theta^{i j k l m} + 
    \frac{1}{2 \, 3!^2} \langle h_{ijk} | h_{lmn} \rangle \Delta \theta^{i j k l m n}
   \right) + ...,

where :math:`\Delta^{i_1 ... i_N } = \Delta \theta^{i_1} ... \Delta \theta^{i_N}`. It can be 
proved that the order N contribution to the ``DALI`` expansion is given by

.. math::

   \frac{1}{N! 1!} \langle h_{i_1 ... i_N} | h_{j_1} \rangle \Delta \theta^{i_1 ... i_N j_1} + \dots + 
   \frac{1}{2 \, N!^2} \langle h_{i_1 ... i_N} | h_{j_1 ... j_N} \rangle \Delta \theta^{i_1 ... i_N j_1 ... j_N}.

Presently, the framework can calculate the ``DALI`` tensors up to derivative order N=4 
with the function ``DALITensors``.

Taylor form
-----------

Taylor form refers to the grouping of the terms in the ``DALI`` expansion in powers of 
\Delta \theta^i. This form is useful for sampling the likelihood, as it allows for a more
efficient evaluation algorithm. For instance, to derivative order N=3 the Taylor Form would be 

