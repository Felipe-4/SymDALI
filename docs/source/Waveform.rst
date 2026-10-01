Waveform
========

The Wolfram Language code has support for the following waveform models:

* ``IMRPhenomD``: non-precessing, aligned-spin, inspiral-merger-ringdown
  waveform model for binary black holes (BBH), see `1508.07250 <https://arxiv.org/abs/1508.07250>`_ and `1508.07253 <https://arxiv.org/abs/1508.07253>`_.
* ``IMRPhenomHM``: non-precessing, aligned-spin, inspiral-merger-ringdown
  waveform model for BBH with higher modes, see `1708.00404 <https://arxiv.org/abs/1708.00404>`_.

These waveform models are implemented as functions that return the plus and
cross polarizations of the GW signal, given a set of parameters.

``hphcIMRPhenomD``
------------------

.. code-block:: mathematica

   hphcIMRPhenomD[f, Mc, q, s1z, s2z, iota, tc, phiRef, invdL]

Calculates the polarizations :math:`\{h_+, h_\times\}` of the ``IMRPhenomD`` approximant.

**Parameters**:

   ``f``: Frequency vector, where ``f[[i]]`` is in [Hz].

   ``Mc``: Chirp mass [solar mass].

   ``q``: Mass ratio in (0,1].

   ``s1z``: Aligned spin, in (-1,1).

   ``s2z``: Aligned spin, in (-1,1).

   ``iota``: Inclination angle, in [0, :math:`\pi`].

   ``tc``: Coalescence time [s].

   ``phiRef``: Reference phase, in [0, 2 :math:`\pi`].

   ``invdL``: Inverse luminosity distance [1/Gpc].

**Returns**: 
   The polarizations :math:`\{h_+, h_\times\}` of the ``IMRPhenomD`` approximant. 
   The output has dimensions ``{2, N}``, where ``N = Length[f]``. 



``hphcIMRPhenomHM``
------------------

.. code-block:: mathematica

   hphcIMRPhenomHM[f, Mc, q, s1z, s2z, iota, tc, phiRef, invdL]

Calculates the polarizations :math:`\{h_+, h_\times\}` of the ``IMRPhenomHM`` approximant.

**Parameters**:

   ``f``: Frequency vector, where ``f[[i]]`` is in [Hz].

   ``Mc``: Chirp mass [solar mass].

   ``q``: Mass ratio in (0,1].

   ``s1z``: Aligned spin, in (-1,1).

   ``s2z``: Aligned spin, in (-1,1).

   ``iota``: Inclination angle, in [0, :math:`\pi`].

   ``tc``: Coalescence time [s].

   ``phiRef``: Reference phase, in [0, 2 :math:`\pi`].

   ``invdL``: Inverse luminosity distance [1/Gpc].

**Returns**
   The polarizations :math:`\{h_+, h_\times\}` of the ``IMRPhenomHM`` approximant. 
   The output has dimensions ``{2, N}``, where ``N = Length[f]``. 

