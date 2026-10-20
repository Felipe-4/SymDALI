.. _dali_tensors:

DALI tensors tutorial
=====================

Higher order ``DALI`` tensors are calculated with the function ``DALITensors``. This 
function takes as input the waveform model, the source parameters, the detector network 
and the derivative order for ``DALI``.

.. code-block:: mathematica

   DALITensors[Appr, fp, {det1, det2, ...}, n]

Calculates the ``DALI`` tensors for the waveform model ``Appr``, with source parameters
``fp`` and the detector network ``{det1, det2, ...}`` up to derivative order ``n``.

**Parameters**

   ``Appr``: Waveform model. Supported values are ``"IMRPhenomD"`` and ``"IMRPhenomHM"``.

   ``fp``: Association containing the source parameters.

   ``{det1, det2, ...}``: List of associations representing the detectors in the network.

   ``n``: Derivative order for ``DALI`` tensors. ``SymDALI`` has support up to ``n=4``.

**Returns**

The ``DALI`` tensors for the specified derivative order.

The convention for ``fp`` follows section :ref:`fiducial` and the convention for ``deti`` 
follows section :ref:`detector_asso`.

.. _dali_tensors_ex:

Example
-------

The following example calculates the third order ``DALI`` tensors using ``IMRPhenomHM``
at the same fiducial point and detector network of Figure 5 in `2609.38399 <https://arxiv.org/abs/2609.38399>`_

.. code-block:: mathematica

   fp = <|
       "dec" -> -0.883787, "\[Phi]" -> 6.18412, "\[Psi]" -> 0.471846,
       "\[ScriptCapitalM]c" -> 190.258, "q" -> 0.935815,
       "s1z" -> 0.107562, "s2z" -> -0.0254365, "\[Iota]" -> 2.60021,
       "1/dL" -> 0.0284797, "tc" -> 61093.87906420612, "\[Phi]ref" -> 5.93578
   |>;

   ET1 = <|
	   "Position"-> DetectorVertex["ET-S"], 
	   "DetectorTensor"->DetectorTensor["ET-S-L"], 
	   "ASD"->ASD["ET-15-lfhf"]
   |>;

   ET2 = <|
	   "Position"->DetectorVertex["ET-MR"], 
	   "DetectorTensor"->DetectorTensor["ET-MR-L-45"], 
	   "ASD"->ASD["ET-15-lfhf"]
   |>;

   CE = <|
	   "Position"-> DetectorVertex["CE-NM"], 
	   "DetectorTensor"->DetectorTensor["CE-NM"], 
	   "ASD"->ASD["CE-40"]
   |>;

   (*Given the total mass of the injection after 93 Hz the waveform is 0*)
   fcut = 93.;

   (*We take Delta_f to be 1/Delta_t == 1/2 for this signal*)
   res = Length@Range[5., 93., 1/2.] + 1;

   tensors = DALITensors["IMRPhenomHM", fp, {ET1, ET2, CE}, 3, "fmin"->5, "fmax"->93, "res"->res];

We note that the structure of the output is a depth-2 list of vectors: 

.. code-block:: mathematica

   Table[
        tensors[[i, j]] // Length,
        {i, 1, 2},
        {j, 1, i}
    ]
    {{121}, {726, 4356}, {3146, 18876, 81796}}

Each tensor in the ``DALI`` formalism is represented as a one dimensional vector. We recognize the 
first element as the Fisher matrix in 11 dimensions, with 121 elements. The second element has 
the representation of :math:`\langle h_{ij} | h_k  \rangle` and :math:`\langle h_{ij} | h_{kl}  \rangle` with 726 
elements with 4356 elements respectvely. In summary, the output of ``DALITensors`` is 
schematically equivalent to the following structure:

.. math::

   \Big\{& \\
        &\quad \{ \langle h_{i} | h_j  \rangle  \}, \\
        &\quad \{ \langle h_{ij} | h_k  \rangle  , \,  \langle h_{ij} | h_{kl}  \rangle \} , \\
        &\quad \dots \\
   \Big\}&,

where each quantity :math:`\langle h_{i_1 \dots i_n} | h_{j_1 \dots j_m}  \rangle` 
is represented as a one dimensional vector. The relation between waveform variables and 
dimensions of the tensors follows that of the Fisher matrix case. As explained in :ref:`fisher_tutorial`
this relation can be obtained by calling ``EmptyAssociation["Aligned"]``.


The function ``DALITensors`` does eliminate
some redundant elements in the tensors, but not all of them. For instance, 
the Fisher matrix representation does contain 121 elements, even though it is symmetric and 
only 66 of them are independent. The same applies to the higher order tensors. 
As pointed out in section :ref:`taylor_section`, for likelihood sampling it is more 
efficient to use the Taylor form of the ``DALI`` expansion, which can be obtained with 
the function ``TaylorForm``

.. code-block:: mathematica

   taylor = TaylorForm[tensors];

Now, we observe that this quantity carries only the independent elements of all tensors:

.. code-block:: mathematica

   Table[
        taylor[[i]] // Length,
        {i, 1, 3}
    ]
    {66, 286, 1001, 3003, 8008}

In fact, the output ``taylor`` defined above is a representation of the following structure:

.. math::

   \Big\{& \\
        &\quad -\frac{1}{2} \{ \langle h_{i} | h_j  \rangle  \}, \\
        &\quad - \left(\frac{1}{2 \, 2!^2} \langle h_{ij} | h_{kl} \rangle  + \frac{1}{3! 1!} \langle h_{ijk} | h_l \rangle \right) \\ 
        &\quad -  \frac{1}{3! 2!} \langle h_{ijk} | h_{lm} \rangle  \\
        &\quad - \frac{1}{2 \, 3!^2} \langle h_{ijk} | h_{lmn} \rangle \\ 
   \Big\}&

where the labels were all symmetrized and only the LI components are kept. This explains why 
we got ``{66, 286, 1001, 3003, 8008}`` as the lengths of the elements of ``taylor``. These 
are the numbers of LI components in a totally symmetric tensor of rank (2, 3, 4, 5, 6) in 11 dimensions.

After calculating the Taylor form, the recommendation is to save it in a text file:

.. code-block:: mathematica

   Export["taylor.txt", Flatten[taylor], "Table"];

This file is used in the ``Python`` part of the framework to sample the likelihood, which is described 
in section :ref:`sampling`.

Options
-------

The options of ``DALITensors`` are already described in the Fisher tutorial in section
:ref:`options_fisher`. They extend naturally to the higher order case.