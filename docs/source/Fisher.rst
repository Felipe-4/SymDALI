.. _fisher_tutorial:

Fisher tutorial
===============
The Fisher matrix can be calculated using the function
``DALITensors``. Its inputs are the waveform model, the source parameters, the 
detector network and the order for ``DALI`` tensors. Since we want the 
Fisher matrix, we will set the order to 1. The output is the Fisher matrix of the network.

.. code-block:: mathematica

   DALITensors[Appr, fp, {det1, det2, ...}, 1]

Calculates the Fisher matrix for the waveform model ``Appr``, with source parameters
``fp``, and the detector network ``{det1, det2, ...}``.

**Parameters**

   ``Appr``: Waveform model. Supported values are ``"IMRPhenomD"`` and ``"IMRPhenomHM"``.

   ``fp``: Association containing the source parameters.

   ``{det1, det2, ...}``: List of associations representing the detectors in the network.

**Returns**

The Fisher matrix of the detector network.

The convention for ``fp`` follows section :ref:`fiducial` and the convention for ``deti`` 
follows section :ref:`detector_asso`.

Example
-------

The following example calculates the Fisher matrix using the ``IMRPhenomD``
waveform model:

.. code-block:: mathematica

   fp = <|
       "dec" -> 3.14/2, "phi" -> 3.4, "psi" -> 0.1,
       "Mc" -> 0.0116, "q" -> 0.2,
       "s1z" -> 0.1, "s2z" -> 0.1, 
       "iota" -> 2.5, "invdL" -> 1, "tc" -> 10, "phiRef" -> 0.
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

   FM = DALITensors["IMRPhenomD", fp, {H1, L1}, 1];

We note that ``FM`` has dimensions ``{1,1,121}``:

.. code-block:: mathematica

   Dimensions[FM]
    {1, 1, 121}

If we want the familiar matrix dimensions we can use 

.. code-block:: mathematica

   FM2 = ArrayReshape[FM, {11, 11}];

To invert the Fisher matrix we can use the ``CholeskyInverse`` function:

.. code-block:: mathematica

   {Cov, epsilon} = CholeskyInverse[FM2];

``Cov`` is the Fisher inverse calculated through a Cholesky decomposition. 
``epsilon`` is the error on the inversion. We defined ``epsilon``, following `2207.02771 <https://arxiv.org/pdf/2207.02771>`_, 
to be the largest error between the product of the Fisher matrix and its inverse with respect to the identity:

.. code-block:: mathematica

   diff1 = FM2.Cov - IdentityMatrix[11]//Abs;
   diff2 = Cov.FM2- IdentityMatrix[11]//Abs;
   epsilon = Max[{diff1, diff2}];

The recommendation is that ``epsilon`` should be less than ``0.05``. 

The derivatives that are used to calculate the Fisher matrix are with respect to
the variables in ``fp``. Their order, that is which variable corresponds to which
dimension, can be accessed through the function ``EmptyAssociation``. For example, 
the following code gives us the exact correspondence between waveform variables and 
matrix dimensions: 

.. code-block:: mathematica:

    EmptyAssociation["Aligned"]
    <|
        "dec" -> 1, "phi" -> 2, "psi" -> 3, "Mc" -> 4,
        "q" -> 5, "s1z" -> 6, "s2z" -> 7, "iota" -> 8, "invdL" -> 9, 
        "tc" -> 10, "phiRef" -> 11
    |>


.. _options_fisher:

Options
-------

The function also accepts options to specify the frequency range, 
resolution and whether to return the Fisher matrix of each detector separately.
The default settings are shown below:

.. code-block:: mathematica

   FM = DALITensors[
       "IMRPhenomD", fp, {H1, L1}, 1,
       "fmin" -> 20, "fmax" -> 1024, "res" -> 1000,
       "AllFisherMatrices" -> False
   ];

When ``AllFisherMatrices`` is set to ``True``, the output will be a 
list of the previous output form:

.. code-block:: mathematica

   FM = DALITensors[
       "IMRPhenomD", fp, {H1, L1}, 1,
       "fmin" -> 20, "fmax" -> 1024, "res" -> 1000,
       "AllFisherMatrices" -> True
   ];
   Dimensions[FM]
    {2, 1, 1, 121}

The order of Fisher matrices is the same of the detector network, in this case ``{H1, L1}``. 

For the ``IMRPhenomD`` model is possible to calculate the Fisher matrix with derivatives
of the waveform on ``TIGER`` parameters:

.. math::
    \left(\delta \varphi_{-2}, \, \delta \varphi_{0}, \delta \varphi_{1}, \delta \varphi_{2},
    \delta \varphi_{3}, \delta \varphi_{4}, \delta \varphi_{5l}, \delta \varphi_{6}, \delta \varphi_{6l}, \delta \varphi_{7} \right)

The way of including these extra parameters is passing them as additional keys for the 
fiducial point association. For instance

.. code-block:: mathematica

    fp = <|
       "dec" -> 3.14/2, "phi" -> 3.4, "psi" -> 0.1,
       "Mc" -> 0.0116, "q" -> 0.2,
       "s1z" -> 0.1, "s2z" -> 0.1, 
       "iota" -> 2.5, "invdL" -> 1, "tc" -> 10, "phiRef" -> 0.,
       "\[Delta]\[CurlyPhi]0" -> 0
   |>;

when passed to ``DALITensors`` will include derivatives of :math:`\delta \varphi_0`. The
variable :math:`\delta \varphi_0` will correspond to the dimension 12 of the Fisher matrix,
keeping the relation of GR variables and dimension intact. You can pass more than one 
parameter:

.. code-block:: mathematica

    fp = <|
       "dec" -> 3.14/2, "phi" -> 3.4, "psi" -> 0.1,
       "Mc" -> 0.0116, "q" -> 0.2,
       "s1z" -> 0.1, "s2z" -> 0.1, 
       "iota" -> 2.5, "invdL" -> 1, "tc" -> 10, "phiRef" -> 0.,
        "\[Delta]\[CurlyPhi]2" -> 0, "\[Delta]\[CurlyPhi]0" -> 0
   |>;

In which case :math:`\delta \varphi_0` will correspond to dimension 12 and :math:`\delta \varphi_2`
to dimension 13. The ordering is defined by the sequence 

.. math::
    \left(\delta \varphi_{-2}, \, \delta \varphi_{0}, \delta \varphi_{1}, \delta \varphi_{2},
    \delta \varphi_{3}, \delta \varphi_{4}, \delta \varphi_{5l}, \delta \varphi_{6}, \delta \varphi_{6l}, \delta \varphi_{7} \right),

where lower numbers will always correspond to lower dimensions. The syntax for all 
``TIGER`` parameters is the same of ``Options[hphcIMRPhenomD]``.
