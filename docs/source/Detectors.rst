Detectors 
=========================
Here we describe how the Wolfram Language code characterizes a detector, how to access
the built-in detectors and how to define a new detector.

Characterization
----------------


A detector is specified by three properties: its location, orientation, and
noise power spectral density :math:`S_n(f)`. Our convention
is to use the WGS-84 model to define the detector's position and orientation.

The position is given by the Cartesian coordinates of the detector's vertex
in meters. For example, the LIGO Livingston detector position vector is:

.. code-block:: mathematica

   PosL1 = {-74276., -5.49628*10^6, 3.22426*10^6};

The orientation is given by the detector tensor:

.. math::

   D_{ij} = \frac{1}{2} \left(nx_i nx_j - ny_i ny_j\right),

where :math:`nx` and :math:`ny` are the unit vectors along the two detector
arms. For the LIGO Livingston detector, the orientation is given by the
matrix:

.. code-block:: mathematica

   DL1 = {
       {0.411318, 0.14021, 0.247279},
       {0.14021, -0.108998, -0.181597},
       {0.247279, -0.181597, -0.302236}
   };

The noise properties are represented by an
``InterpolatingFunction`` of the amplitude spectral density (ASD),
:math:`\sqrt{S_n(f)}`, which can be generated using ``Interpolation``:

.. code-block:: mathematica

   data = Import["AplusDesign.txt", "Data"];
   ASDL1 = Interpolation[data];

The functions in ``SymDALI`` require detectors to be specified by quantities
such as ``PosL1``, ``DL1``, and ``ASDL1``. These quantities can be defined
for any detector. 

Buit-in detectors
-----------------

``SymDALI`` has built-in support for the
LIGO-Virgo-KAGRA (LVK) network, Einstein Telescope (ET), and Cosmic Explorer
(CE). Predefined detector positions are provided by ``DetectorVertex``. Predefined
orientations are provided by ``DetectorTensor``, and built-in ASDs are
provided by ``ASD``. The usage messages of these three functions (for instance ``?DetectorVertex``) provide 
a list of the available detectors and their corresponding names.

For example, the previous definitions can instead be generated with:

.. code-block:: mathematica

   PosL1 = DetectorVertex["L1"];
   DL1 = DetectorTensor["L1"];
   ASDL1 = ASD["L1H1-O5"];


Defining a new detector
-----------------------

To define the position and orientation of a new detector, we can use the functions
``FromGeocentricCoordinates`` and ``ArmDirection``.

.. code-block:: mathematica

   FromGeocentricCoordinates[{h, lat, lon}]

Returns the Cartesian coordinates ``{x, y, z}`` of a point.

**Parameters**:

    ``h`` :Height of the vertex above the reference ellipsoid, in meters. A value of zero can be used for simplicity.

    ``lat``: Vertex latitude in radians.

    ``lon``: Vertex longitude in radians.

**Returns**:
    The Cartesian coordinates ``{x, y, z}`` of the point, in meters.

.. code-block:: mathematica

   ArmDirection[
       {lat, lon},
       {{altitudex, azimuthx}, {altitudey, azimuthy}}
   ]

Returns the unit vectors along the two detector arms.

**Parameters**:

    ``lat``: Vertex latitude in radians.

    ``lon``: Vertex longitude in radians.

    ``altitudex``: Altitude angle of the first arm, in radians.

    ``azimuthx``: Azimuth angle of the first arm, in radians.

    ``altitudey``: Altitude angle of the second arm, in radians.

    ``azimuthy``: Azimuth angle of the second arm, in radians.

**Returns**:
    The unit vectors ``{nx, ny}`` along the two detector arms.

The altitude and azimuth angles follow the convention used in
`LALDetectors.h <https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lal/lib/tools/LALDetectors.h>`_:

**Azimuth**: Angle measured clockwise from north of the direction vector.

**Altitude**: The angle the direction vector makes with the local horizontal plane. Positive values indicate directions above the horizontal.

The following example shows how to define the position and orientation of
LIGO Livingston using these functions. The reference values are taken from
`LALDetectors.h <https://git.ligo.org/lscsoft/lalsuite/-/blob/master/lal/lib/tools/LALDetectors.h>`_:

.. code-block:: mathematica

   PosL1 = FromGeocentricCoordinates[{-6.574, 0.53342313506, -1.58430937078}];

   {nx, ny} = ArmDirection[
       {0.53342313506, -1.58430937078},
       {
           {-0.00031210000, 4.40317772346},
           {-0.00061070000, 2.83238139666}
       }
   ];

   DL1 = (TensorProduct[nx, nx] - TensorProduct[ny, ny])/2;

A new ASD can be defined by importing the corresponding data and applying
``Interpolation``.


Detector association
--------------------

The functions ``SNR`` and ``DALITensors`` take a list of one or more
detectors as one of their inputs. Each detector is specified by an
association with the following keys:

.. code-block:: mathematica

   L1 = <|
       "Position" -> PosL1,
       "Orientation" -> DL1,
       "ASD" -> ASDL1
   |>;


Pattern functions
-----------------


The ...