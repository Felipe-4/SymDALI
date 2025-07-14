(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`Detectors`"]


Unprotect[DetectorTensor, DetectorVertex, ASD, ArmDirection, DetectorPosition];


ClearAll[DetectorTensor, DetectorVertex, ASD, ArmDirection, DetectorPosition]


DetectorTensor::usage = "```DetectorTensor[\"det\"]``` returns the detector tensor \!\(\*SubscriptBox[\(D\), \(ij\)]\) = 0.5 (\!\(\*SubscriptBox[\(n\), \(i\)]\) \!\(\*SubscriptBox[\(n\), \(j\)]\) - \!\(\*SubscriptBox[\(n\), \(j\)]\) \!\(\*SubscriptBox[\(n\), \(i\)]\)),
in the form of a 3x3 matrix.
availabe detectors: {\"H1\", \"L1\", \"V1\", \"K\", \"ET1\", \"ET2\", \"ET3\"}
Data taken from lal suite 
(https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html)";


DetectorVertex::usage="```DetectorVertex[\"det\"]``` returns the detector position {x,y,z} [meters].
availabe detectors: {\"H1\", \"L1\", \"V1\", \"K\", \"ET1\", \"ET2\", \"ET3\"}
Data taken from lal suite 
(https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html)";


ASD::usage="```ASD[\"name\"]``` returns \!\(\*SqrtBox[\(\*SubscriptBox[\(S\), \(n\)] \((f)\)\)]\) in the form of an ```InterpolationFunction```.
available \"names\":
	\"ET-D\":  public ET-D (https://www.et-gw.eu/index.php/etsensitivities).
	\"CE-20\":  baseline 20 km detector, \"compact binary tuned\" (https://dcc.cosmicexplorer.org/CE-T2000017/public).
	\"CE-20-pm\":  20 km detector tuned for post-merger signals (https://dcc.cosmicexplorer.org/CE-T2000017/public).
	\"CE-40\":  baseline 40 km detector (https://dcc.cosmicexplorer.org/CE-T2000017/public).
	\"CE-40-lf\": 40 km detector tuned for low-freqency signals (https://dcc.cosmicexplorer.org/CE-T2000017/public).";


ArmDirection::usage="```ArmDirection[{\[CurlyPhi], \[Lambda]}, { {\[Omega]1, \[Psi]1},  {\[Omega]2, \[Psi]2} }]``` returns {n1, n2}, where \"ni\": {nix,niy,niz} unit vector in the direction of arm \"i\"
\[CurlyPhi]: vertex latitude [radians]
\[Lambda]: vertex longitude [radians]
\[Omega]i: tilt of arm \"i\" above the horizontal plane [radians] 
\[Psi]i: Azimuth of arm \"i\" [radians] (North of East)

coordinates follow WGS-84 model used in lal suite (https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html)";


DetectorPosition::usage="```DetectorPosition[h, \[Phi], \[Lambda]]``` returns {x,y,z} position of the detector.

h: displacement along the local vertical [meters]
\[CurlyPhi]: latitude [radians]
\[Lambda]: longitude [radians]

coordinates follow WGS-84 model used in lal suite (https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html)";


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*DetectorTensor*)


DetectorTensor = <||>;


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


DetectorTensor["ET1"] = Module[
	{nx =  {-0.70045821479, 0.20848948619,0.68256166277}, ny = {-0.39681482542,-0.73500471881, 0.54982366052},d},
	
	d = 1/2 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];

DetectorTensor["ET2"] = Module[
	{nx = {0.30364338937, -0.94349420500, -0.13273800225 }, ny = {0.70045821479, -0.20848948619, -0.68256166277},d},
	
	d = 1/2 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];

DetectorTensor["ET3"] = Module[
	{nx = {0.39681482542, 0.73500471881, -0.54982366052 }, ny = {-0.30364338937, 0.94349420500, 0.13273800225},d},
	
	d = 1/2 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];


DetectorTensor["K"] = Module[
	{nx = {-0.3759040,-0.8361583, 0.3994189}, ny = {0.7164378,0.01114076, 0.6975620},d},
	
	d = 1/2 (nx\[TensorProduct]nx-ny\[TensorProduct]ny)
];


Protect[DetectorTensor];


(* ::Subsection::Closed:: *)
(*DetectorVertex*)


DetectorVertex = <||>;


DetectorVertex["H1"] = {-2.16141492636 10^6,  -3.83469517889 10^6   , 4.60035022664 10^6};
DetectorVertex["L1"] = {-7.42760447238 10^4, -5.49628371971 10^6 ,  3.22425701744 10^6 };
DetectorVertex["V"] = {4.54637409900 10^6, 8.42989697626 10^5, 4.37857696241 10^6};
DetectorVertex["K"] = {-3777336.024,  3484898.411, 3765313.697};
DetectorVertex["ET1"] = {4.54637409900 10^6, 8.42989697626 10^5,4.37857696241 10^6};
DetectorVertex["ET2"] = {4.53936951685 10^6, 8.45074592488 10^5, 4.38540257904 10^6};
DetectorVertex["ET3"] = { 4.54240595075 10^6, 8.35639650438 10^5, 4.38407519902 10^6};


Protect[DetectorVertex];


(* ::Subsection::Closed:: *)
(*ASDs*)


ASDdir = With[
	{pacletDir= FindFile["FelipeBarbosa`SymDALI`"]//FileNameDrop[#,-2]&},
	FileNameJoin[{pacletDir, "/Data/ASD"}]
];


ASD =<||>;


Module[
	{data},
	
	(*ET-D*)
	data = Import[FileNameJoin[{ASDdir,"ET-0000A-18_ETDSensitivityCurveTxtFile.txt"}], "Data"];
	data = data[[All,{1,4}]];
	ASD["ET-D"] = Interpolation[data];
	
	(*CE-20*)
	data = Import[FileNameJoin[{ASDdir,"cosmic_explorer_20km_strain.txt"}], "Data"]//N;
	ASD["CE-20"] = Interpolation[data];
	 
	(*CE-40-lf*)
	data = Import[FileNameJoin[{ASDdir,"cosmic_explorer_40km_lf_strain.txt"}], "Data"]//N;
	ASD["CE-40-lf"] = Interpolation[data];
	 
	(*CE-20-pm*)
	data = Import[FileNameJoin[{ASDdir, "cosmic_explorer_20km_pm_strain.txt"}], "Data"]//N;
	ASD["CE-20-pm"] =Interpolation[data];
	
	data = Import[FileNameJoin[{ASDdir, "cosmic_explorer_strain.txt"}], "Data"]//N;
	ASD["CE-40"] =Interpolation[data];
]


Protect[ASD];


(* ::Subsection::Closed:: *)
(*ArmDirection and DetectorPosition*)


ArmDirection[{\[CurlyPhi]_, \[Lambda]_}, {{\[Omega]1_, \[Psi]1_}, {\[Omega]2_, \[Psi]2_}}] := Module[
	{comps1, comps2, n1, n2, \[Theta] = \[Pi]/2 - \[CurlyPhi], \[Phi] = \[Lambda], rhat, \[Phi]hat, \[Theta]hat},
	
	(* compsi: {r, \[Theta], \[Phi]} components in the base {rhat, \[Theta]hat, \[Phi]hat} of armi*)
	comps1 = {Sin[\[Omega]1], - Cos[\[Omega]1] Sin[\[Psi]1], Cos[\[Omega]1] Cos[\[Psi]1]}; 
	comps2 = {Sin[\[Omega]2], - Cos[\[Omega]2] Sin[\[Psi]2], Cos[\[Omega]2] Cos[\[Psi]2]};
	
	(*rhat, \[Theta]hat, \[Phi]hat components in {x,y,z}*)
	rhat = {Sin[\[Theta]] Cos[\[Phi]], Sin[\[Theta]] Sin[\[Phi]], Cos[\[Theta]]};
	\[Theta]hat = {Cos[\[Theta]] Cos[\[Phi]], Cos[\[Theta]] Sin[\[Phi]], - Sin[\[Theta]]};
	\[Phi]hat = {-Sin[\[Phi]], Cos[\[Phi]], 0};
	
	
	n1 = comps1[[1]]rhat + comps1[[2]]\[Theta]hat + comps1[[3]]*\[Phi]hat;
	n2 = comps2[[1]]rhat + comps2[[2]]\[Theta]hat + comps2[[3]]*\[Phi]hat;
	
	{n1, n2}//N
]


DetectorPosition[h_, \[CurlyPhi]_, \[Lambda]_] := Module[
	{R, a =6378137 , b = 6356752.314, X,Y,Z},
	
	R = a^2/Sqrt[a^2 Cos[\[CurlyPhi]]^2 + b^2 (Sin[\[CurlyPhi]]^2) ];
	
	X = (R+h) Cos[\[CurlyPhi]] Cos[\[Lambda]];
	Y = (R+h) Cos[\[CurlyPhi]] Sin[\[Lambda]];
	Z =  (0.993306 R +h) Sin[\[CurlyPhi]];
	
	{X,Y,Z}
]


Protect[ArmDirection]; Protect[DetectorPosition];


(* ::Section::Closed:: *)
(*Package Footer*)


End[];


EndPackage[]
