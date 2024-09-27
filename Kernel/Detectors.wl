(* ::Package:: *)

BeginPackage["FelipeBarbosa`SymDALI`Detectors`"]


FpFc::usage =
"
FpFc[H1, \[Alpha], \[Delta], \[Psi], GMST] gives the antena pattern functions {F+, Fx}, for:
- right ascencion: \[Alpha]
- declination: \[Delta]
- polarization angle: \[Psi] (geocentric frame)
- Greenwich mean sidereal time: GMST
(angles must be given in radians)
Obs.: H1 is a string,
H1: Ligo Hanford
L1: Ligo Livingston
V1: Virgo
";


Begin["`Private`"];


LSO::usage="LSO[M] gives the GW frequency at the last stable circular orbit for
a test mass in a Schwarzschild spacetime of mass M
M[SolarMass]";


LSO[M_] = With[{G = (Quantity[("GravitationalConstant")/("SpeedOfLight")^3]//UnitConvert[#, ("Seconds")/("SolarMass")]&)[[1]]}, (6^(3/2) M G \[Pi])^-1];


DetectorTensor::usage = "Normal vectors taken from https://arxiv.org/pdf/gr-qc/0008066 which seem to agree with the ones in
https://www.ligo.org/scientists/GW100916/detectors.txt";


DetectorTensor["H1"] = With[
    {nx= {-0.2239, 0.7998, 0.5569}, ny = {-0.9140, 0.0261, -0.4049}},
    (nx\[TensorProduct]nx - ny\[TensorProduct]ny)/2
];

DetectorTensor["L1"] = With[
   {nx = {\[Minus]0.9546,\[Minus]0.1416,\[Minus]0.2622}, ny = {+0.2977,\[Minus]0.4879,\[Minus]0.8205} },
   0.5 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];

DetectorTensor["V1"] = With[
    {nx = {\[Minus]0.7005,+0.2085,+0.6826}, ny = {\[Minus]0.0538,\[Minus]0.9691,+0.2408} },
	0.5 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];


(* ::Section:: *)
(*Antenna Pattern Functions*)


FpFc[detector_String, \[Alpha]_, \[Delta]_, \[Psi]_, GMST_] := Module[
    {\[Theta] = \[Pi]/2 - \[Delta], \[Phi] = \[Alpha] - GMST, eplus, ecross, R1, R2, T, D, eplusGeoFrame, ecrossGeoFrame},
	(*Polarization frame tensors*)
    eplus = {{1,0,0}, {0,-1,0}, {0,0,0}};
    ecross = {{0,1,0}, {1,0,0}, {0,0,0}};
    
    (*R1: Geocentric frame -> Wave-Frame (apart of a parity transf.)
	  R2: Wave-frame -> polarization frame*)
    R1 = List[
        {-Sin[\[Phi]], Cos[\[Phi]], 0}, (*\hat{\[Phi]} decomposed in \hat{e}_i*)
        {Cos[\[Theta]] Cos[\[Phi]], Cos[\[Theta]] Sin[\[Phi]], -Sin[\[Theta]]}, (*\hat{\[Theta]} decomposed in \hat{e}_i*)
        {Sin[\[Theta]] Cos[\[Phi]], Sin[\[Theta]] Sin[\[Phi]], Cos[\[Theta]]} (*\hat{r}decomposed in \hat{e}_i*)
    ];
    R2 = RotationMatrix[-\[Psi], {0,0,1}];
    T = R2 . R1//Simplify;
    eplusGeoFrame = (T\[Transpose] . eplus . T//Simplify);
    ecrossGeoFrame = (T\[Transpose] . ecross . T//Simplify);
    
    D = DetectorTensor[detector];
    (* D_{ij} e_{ij} = Tr[D.e\[Transpose]]  *)
    Tr/@{D . eplusGeoFrame\[Transpose], D . ecrossGeoFrame\[Transpose]}//Simplify
]


End[];


EndPackage[]
