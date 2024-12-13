(* ::Package:: *)

BeginPackage["FelipeBarbosa`SymDALI`Detectors`"]


epec::usage = "epec[\[Alpha], sin\[Delta], \[Psi], GMST] returns the {\!\(\*SubscriptBox[\(e\), \(+\)]\), \!\(\*SubscriptBox[\(e\), \(x\)]\)} polarization tensors in the geocentric frame, assuming
a GW coming from a sky direction with right ascencion \[Alpha], declination with Sin sin\[Delta] a polarization angle \[Psi] at 
the GMST time GMST.";

\[CapitalDelta]t::usage = "\[CapitalDelta]t[pos, \[Alpha], sin\[Delta], GMST] returns the time interval that it takes a GW (with propagation speed c)
to go from the detector to the Earth's center, assuming the wave is coming from the direction \[Alpha], sin\[Delta] at GMST time
GMST and the detector is in the position pos. The position should be given in metters following the convention of 
https://www.ligo.org/scientists/GW100916/detectors.txt 
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


epec[\[Alpha]_, sin\[Delta]_, \[Psi]_, GMST_] := Module[
    {\[Theta] = \[Pi]/2 - \[Delta], \[Phi] = \[Alpha] - GMST, eplus, ecross, R1, R2, T, D, eplusGeoFrame, ecrossGeoFrame, sinvar},
	(*Polarization frame tensors*)
    eplus = {{1,0,0}, {0,-1,0}, {0,0,0}};
    ecross = {{0,1,0}, {1,0,0}, {0,0,0}};
    
    (*R1: Geocentric frame -> Wave-Frame (apart of a parity transf.)
	  R2: Wave-frame -> polarization frame*)
    R1 = List[ (*Writing it already as a function of sin\[Delta]: *)
        {-Sin[\[Phi]], Cos[\[Phi]], 0}, (*\hat{\[Phi]} decomposed in \hat{e}_i*)
        {sin\[Delta] Cos[\[Phi]],  sin\[Delta] Sin[\[Phi]],  -Sqrt[1- sin\[Delta]^2]}, (*\hat{\[Theta]} decomposed in \hat{e}_i*)
        {Sqrt[1- sin\[Delta]^2] Cos[\[Phi]], Sqrt[1- sin\[Delta]^2] Sin[\[Phi]], sin\[Delta]} (*\hat{r}decomposed in \hat{e}_i*)
    ];
    
    R2 = RotationMatrix[-\[Psi], {0,0,1}];
    T = R2 . R1//Simplify;
    eplusGeoFrame = (T\[Transpose] . eplus . T);
    ecrossGeoFrame = (T\[Transpose] . ecross . T);
    
    {eplusGeoFrame,ecrossGeoFrame}//FullSimplify
    
    
   (* D = 0.5 (nx\[TensorProduct]nx - ny\[TensorProduct]ny);
    (*D_{ij} e_{ij} = Tr[D.e\[Transpose]]*)
    Tr/@{D . eplusGeoFrame\[Transpose], D . ecrossGeoFrame\[Transpose]}//Simplify
   *)
]


\[CapitalDelta]t[pos_, \[Alpha]_, sin\[Delta]_, GMST_] := Module[
	{\[Phi] = \[Alpha]-GMST, r, c = UnitConvert["SpeedOfLight"][[1]]},
	
	(*radial unit vector: {Cos[\[Phi]] Sin[\[Theta]],Sin[\[Theta]] Sin[\[Phi]],Cos[\[Theta]]}, with \[Theta] = \[Pi]/2 - \[Delta] *)
	r = {Sqrt[1- sin\[Delta]^2] Cos[\[Phi]], Sqrt[1- sin\[Delta]^2] Sin[\[Phi]], sin\[Delta]};
	
	-(pos . r/c)
	
]//Simplify;


End[];


EndPackage[]
