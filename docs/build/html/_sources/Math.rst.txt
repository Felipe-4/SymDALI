Math
====

Here we describe the mathematical framework of the ``DALI`` formalism.

Inner product
-------------

All inner products are defined as

.. math::
    
   \langle a | b \rangle = 4 \, \Delta f \, \sum_{i} \Re \left[ \frac{\tilde{a}(f_i) \tilde{b}^*(f_i)}{S_n(f_i)} \right],

where :math:`\tilde{a}(f)` and :math:`\tilde{b}(f)` are the Fourier transforms of the time-domain signals
and :math:`\Delta f` is the size of the frequency bin. There is always a minimum 
frequency, maximum frequency and a grid resolution. These can be adjusted to 
the user's preferences in the options of the functions ``SNR`` and ``DALITensors``.

``DALI`` tensors
----------------

The ``DALI`` Likelihood expansion, up to derivative order 3, is defined as

.. math::

   \ln \mathcal{L}(\theta) = &- \frac{1}{2} \langle h_i | h_j \rangle \Delta \theta^{i j} \\
   &-\left(
    \frac{1}{2! 1!} \langle h_{ij} | h_k \rangle \Delta \theta^{i j k} + 
    \frac{1}{2 \, 2!^2} \langle h_{ij} | h_{kl} \rangle \Delta \theta^{i j k l} 
   \right) \\ 
   &-\left(
    \frac{1}{3! 1!} \langle h_{ijk} | h_l \rangle \Delta \theta^{i j k l} + 
    \frac{1}{3! 2!} \langle h_{ijk} | h_{lm} \rangle \Delta \theta^{i j k l m} + 
    \frac{1}{2 \, 3!^2} \langle h_{ijk} | h_{lmn} \rangle \Delta \theta^{i j k l m n}
   \right),

where :math:`\Delta \theta^{i_1 ... i_N } = \Delta \theta^{i_1} ... \Delta \theta^{i_N}`. It can be 
proved that the order N contribution to the ``DALI`` expansion is given by

.. math::

   \frac{1}{N! 1!} \langle h_{i_1 ... i_N} | h_{j_1} \rangle \Delta \theta^{i_1 ... i_N j_1} + \dots + 
   \frac{1}{2 \, N!^2} \langle h_{i_1 ... i_N} | h_{j_1 ... j_N} \rangle \Delta \theta^{i_1 ... i_N j_1 ... j_N}.

``SymDALI`` can calculate the ``DALI`` tensors up to derivative order N=4 
with the function ``DALITensors``. 


In the case of a network of detectors all ``DALI`` tensors are summed over the detectors: 

.. math::

   &\langle h_i | h_j \rangle = \sum_a \langle h_i | h_j \rangle_a \\ 
   &\langle h_i | h_{j k} \rangle  = \sum_a \langle h_i | h_{j k} \rangle_a \\
   & \qquad \qquad \dots, 
where the inner product with subscript :math:`a` is the usual inner product calculated
with the noise power spectral density of the :math:`a`-th detector.

.. _taylor_section:

Taylor form
-----------

Taylor form refers to the grouping of the terms in the ``DALI`` expansion in powers of 
:math:`\Delta \theta^i`. This form is useful for sampling the likelihood, as it allows for a more
efficient evaluation algorithm. For instance, to derivative order N=3 the Taylor Form would be 

.. math::

   \ln \mathcal{L}(\theta) = &- \frac{1}{2} \langle h_i | h_j \rangle \Delta \theta^{i j} \\
   &-\frac{1}{2! 1!} \langle h_{ij} | h_k \rangle \Delta \theta^{i j k} \\ 
   &-\left(\frac{1}{2 \, 2!^2} \langle h_{ij} | h_{kl} \rangle  + 
   \frac{1}{3! 1!} \langle h_{ijk} | h_l \rangle \right) \, \Delta \theta^{i j k l} \\ 
   &-\frac{1}{3! 2!} \langle h_{ijk} | h_{lm} \rangle \Delta \theta^{i j k l m} \\ 
   &-\frac{1}{2 \, 3!^2} \langle h_{ijk} | h_{lmn} \rangle \Delta \theta^{i j k l m n}
   + ...

Once the ``DALI`` tensors are calculated, the Taylor form can be obtained with the 
function ``TaylorForm``. 