.. SymDALI documentation master file, created by
   sphinx-quickstart on Fri Sep 25 10:51:30 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============
``SymDALI`` is a hybrid framework, using Wolfram Language and Python, to perform approximations of the 
gravitational wave (GW) likelihood through the Fisher matrix and its higher order extension: Derivative 
Approximation for LIkelihoods (``DALI``) [1]_, [2]_. It can calculate the Fisher matrix and ``DALI`` tensors for any
terrestrial network of GW detectors. However, there is built-in support for the LIGO-Virgo-KAGRA (LVK) 
network, Einstein Telescope (ET) and Cosmic Explorer (CE) detectors. 

Wolfram Language
-----------------
Most of the framework is implemented in Wolfram Language. Its functionalities involve the calculation
of the signal-to-noise ratio (SNR) of arbitrary networks; the plus and cross polarizations 
of the GW signal; pattern functions; the Fisher matrix and the higher order ``DALI`` tensors for arbitrary
networks. The waveform derivatives are calculated with Wolfram Language's symbolic 
differentiation, providing accurate numerical results, and compiled to the Wolfram Virtual Machine (WVM), 
for fast evaluation. The calculation of the Fisher matrix for the ``IMRPhenomHM`` waveform model takes less than a second on a Macbook Air M1.

Python
------
The Python part of the framework contains sampling utilities to perform parameter estimation with the Fisher
matrix and ``DALI`` tensors. The sampling stage is always necessary when working 
with higher order ``DALI`` tensors. Although the Fisher formalism allows for parameter estimation with matrix inversion,
it has been shown that sampling the Fisher likelihood with exact priors produces better results 
[3]_, [4]_, [5]_, [6]_, and therefore it is recommended.


.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   Waveform
   Detectors
   quickstart
   user_guide
   examples
   api


References
----------

.. [1] `1401.6892 <https://arxiv.org/abs/1401.6892>`_
.. [2] `2203.02670 <https://arxiv.org/abs/2203.02670>`_
.. [3] `2205.02499 <https://arxiv.org/abs/2205.02499>`_
.. [4] `2307.10154 <https://arxiv.org/abs/2307.10154>`_
.. [5] `2404.16103 <https://arxiv.org/abs/2404.16103>`_
.. [6] `2510.16955 <https://arxiv.org/abs/2510.16955>`_ 