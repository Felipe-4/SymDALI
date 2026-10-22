.. _fiducial:

Fiducial point
==============

In ``SymDALI`` in order to calculate the SNR of a network, the Fisher matrix or 
higher order ``DALI`` tensors an injection value for the source is required. Here we
described the convention adopted by the functions ``SNR`` and ``DALITensors``.

An injection is specified by an association in the following form:

.. code-block:: mathematica

   fp = <|
        "dec"->3.14/2,"phi"->3.4,"psi"->0.1,
        "Mc"->0.0116,"q"->0.2, "s1z"->0.1,"s2z"->0.1,
        "iota"->2.5,"invdL"->1,"tc"->10, "phiRef"->0.
    |>;

where the parameters correspond to 

**Parameters**:

    ``"dec"``: declination, in [:math:`-\pi/2`,  :math:`\pi/2`].

    ``"phi"``: azimuthal angle, in [0,  2 :math:`\pi`].

    ``"psi"``: polarization angle, in [0,  :math:`\pi`].

    ``"Mc"``: Chirp mass [:math:`M_\odot`].

    ``"q"``: Mass ratio, in (0,1].

    ``"s1z"``: Aligned spin of the heaviest BH, in (-1,1).

    ``"s2z"``: Aligned spin of the lightest BH, in (-1,1).

    ``"iota"``: Inclination angle, in [0, :math:`\pi`].

    ``"tc"``: Coalescence time [s].

    ``"phiRef"``: Reference phase, in [0, 2 :math:`\pi`].

    ``"invdL"``: Inverse luminosity distance [1/Gpc].

For practical reasons, an empty association with this convention is available through
the function ``EmptyAssociation``. For example:

.. code-block:: mathematica:

    EmptyAssociation["Aligned"]
    <|
        "dec" -> 1, "phi" -> 2, "psi" -> 3, "Mc" -> 4,
        "q" -> 5, "s1z" -> 6, "s2z" -> 7, "iota" -> 8, "invdL" -> 9, 
        "tc" -> 10, "phiRef" -> 11
    |>

