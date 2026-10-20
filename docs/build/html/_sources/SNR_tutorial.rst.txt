SNR tutorial
============

The signal-to-noise ratio (SNR) can be calculated using the function
``SNR``. Its inputs are the waveform model, the source parameters, and
the detector network. The output is the total network SNR.

.. code-block:: mathematica

   SNR[Appr, fp, {det1, det2, ...}]

Calculates the SNR for the waveform model ``Appr``, with source parameters
``fp``, and the detector network ``{det1, det2, ...}``.

**Parameters**

* ``Appr``: Waveform model. Supported values are ``"IMRPhenomD"`` and
  ``"IMRPhenomHM"``.
* ``fp``: Association containing the source parameters.
* ``{det1, det2, ...}``: List of associations representing the detectors
  in the network.

**Returns**

The total SNR of the detector network.


The convention for ``fp`` follows section :ref:`fiducial` and the convention for ``deti`` 
follows section :ref:`detector_asso`.

Example
-------

The following example calculates the SNR using the ``IMRPhenomD``
waveform model:

.. code-block:: mathematica

   fp = <|
       "dec" -> 3.14/2, "\[Phi]" -> 3.4, "\[Psi]" -> 0.1,
       "\[ScriptCapitalM]c" -> 0.0116, "q" -> 0.2,
       "s1z" -> 0.1, "s2z" -> 0.1, "\[Iota]" -> 2.5,
       "1/dL" -> 1, "tc" -> 10, "\[Phi]ref" -> 0.
   |>;

   H1 = <|
       "Position" -> DetectorVertex["H1"],
       "DetectorTensor" -> DetectorTensor["H1"],
       "ASD" -> ASD["L1H1-O5"]
   |>;

   L1 = <|
       "Position" -> DetectorVertex["L1"],
       "DetectorTensor" -> DetectorTensor["L1"],
       "ASD" -> ASD["L1H1-O5"]
   |>;

   snr = SNR["IMRPhenomD", fp, {H1, L1}];


Options
-------

The function also accepts options to specify the frequency range, resolution 
and whether to return the SNR of each detector separately.
The default settings are shown below:

.. code-block:: mathematica

   snr = SNR[
       "IMRPhenomD", fp, {H1, L1},
       "fmin" -> 10, "fmax" -> 1024,
       "res" -> 1000, "AllSNRs" -> False
   ]

When ``"AllSNRs" -> True`` the function returns a list containing the
SNR of each detector in the network. The order of SNRs is the same of the detector 
network, in this case ``{H1, L1}``.