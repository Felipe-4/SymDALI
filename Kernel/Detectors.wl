(* ::Package:: *)

(* ::Section::Closed:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`Detectors`"]


Unprotect[DetectorTensor, DetectorVertex, ASD, ArmDirection, FromGeocentricCoordinates];


ClearAll[DetectorTensor, DetectorVertex, ASD, ArmDirection, FromGeocentricCoordinates];


DetectorTensor::usage = "```DetectorTensor[\"det\"]``` returns the detector tensor \!\(\*SubscriptBox[\(D\), \(ij\)]\) = 0.5 (\!\(\*SubscriptBox[\(nx\), \(i\)]\) \!\(\*SubscriptBox[\(nx\), \(j\)]\) - \!\(\*SubscriptBox[\(ny\), \(j\)]\) \!\(\*SubscriptBox[\(ny\), \(i\)]\)),
in the form of a 3x3 matrix.
availabe \"det\": 
	{\"H1\", \"L1\", \"V1\", \"K\", \"ET1\", \"ET2\", \"ET3\", \"I1\"}   (standards from \!\(\*TemplateBox[{\"\\\"LALDetectors.h\\\"\", \"https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html\"},\n\"HyperlinkURL\"]\)).

	\"ET-S-L\":  L-Shaped Detector in Sardinia (lat  = 40\[Degree] 31',  lon  = 9\[Degree] 25').
	\"ETi-S\":   Arms of the triangle configuration in Sardinia, with i=1,2,3 (assuming colocated vertices).
	\"ET-MR-L-0\":  L-Shaped Detector in Meuse-Rhine with \[Alpha]=0\[Degree] regarding \"ET-S-L\" , see \!\(\*TemplateBox[{\"\\\"2303.15923\\\"\", \"https://arxiv.org/pdf/2303.15923\"},\n\"HyperlinkURL\"]\) (lat = 50\[Degree] 43' 23'', lon = 5\[Degree] 55' 14'').
	\"ET-MR-L-45\":  L-Shaped Detector in Meuse-Rhine with \[Alpha]=45\[Degree] regarding \"ET-S-L\" , see \!\(\*TemplateBox[{\"\\\"2303.15923\\\"\", \"https://arxiv.org/pdf/2303.15923\"},\n\"HyperlinkURL\"]\).
	\"CE-I\": CE L-shaped Detector in Idaho (see Table III \!\(\*TemplateBox[{\"\\\"2010.15202\\\"\", \"https://arxiv.org/pdf/2010.15202\"},\n\"HyperlinkURL\"]\)).
	\"CE-NM\": CE L-shaped Detector in New Mexico (see Table III \!\(\*TemplateBox[{\"\\\"2010.15202\\\"\", \"https://arxiv.org/pdf/2010.15202\"},\n\"HyperlinkURL\"]\)).

Check  ```ArmDirection``` if you need to define a new Dij.";


DetectorVertex::usage="```DetectorVertex[\"det\"]``` returns the vertex position {x,y,z} [meters].
availabe \"det\": 
	{\"H1\", \"L1\", \"V1\", \"K\", \"ET1\", \"ET2\", \"ET3\", \"I1\"}  (standards from \!\(\*TemplateBox[{\"\\\"LALDetectors.h\\\"\", \"https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html\"},\n\"HyperlinkURL\"]\)).

	\"ET-S\": ET Sardinia location  (lat  = 40\[Degree] 31',  lon  = 9\[Degree] 25').
	\"ET-MR\": ET Meuse-Rhine location  (lat = 50\[Degree] 43' 23'', lon = 5\[Degree] 55' 14'').
	\"CE-I\": CE Idaho location (see Table III \!\(\*TemplateBox[{\"\\\"2010.15202\\\"\", \"https://arxiv.org/pdf/2010.15202\"},\n\"HyperlinkURL\"]\)).
	\"CE-NM\": CE New Mexico location (see Table III \!\(\*TemplateBox[{\"\\\"2010.15202\\\"\", \"https://arxiv.org/pdf/2010.15202\"},\n\"HyperlinkURL\"]\)).

Check  ```FromGeocentricCoordinates```  if you need to define a new vertex.";


(*This is usefull so that ```?ASD``` doesn't take too long to evaluate*)
Attributes[ASD] = {ReadProtected};


ASD::usage="```ASD[\"name\"]``` returns \!\(\*SqrtBox[\(\*SubscriptBox[\(S\), \(n\)] \((f)\)\)]\) in the form of an ```InterpolationFunction```.
available \"names\":
	\"ET-D\":  public ET-D (https://www.et-gw.eu/index.php/etsensitivities).
	\"CE-20\":  baseline 20 km detector, \"compact binary tuned\" (https://dcc.cosmicexplorer.org/CE-T2000017/public).
	\"CE-20-pm\":  20 km detector tuned for post-merger signals (https://dcc.cosmicexplorer.org/CE-T2000017/public).
	\"CE-40\":  baseline 40 km detector (https://dcc.cosmicexplorer.org/CE-T2000017/public).
	\"CE-40-lf\": 40 km detector tuned for low-freqency signals (https://dcc.cosmicexplorer.org/CE-T2000017/public).
	\"ET-n-hf\": 'n' can be 10,15 or 20 (arm length kms). 'hf' means high-frequency (https://apps.et-gw.eu/tds/?r=18213).
	\"ET-n-lf\": 'n' can be 10,15 or 20 (arm length kms). 'lf' means low-frequency (https://apps.et-gw.eu/tds/?r=18213).
	\"ET-n-lfhf\": 'n' can be 10, 15 or 20 (arm length kms). 'lfhf means combined, xylophone, configuration (https://apps.et-gw.eu/tds/?r=18213).
	\"L1H1-05\":  LIGO A+ Design target for O5 (https://dcc.ligo.org/LIGO-T2000012/public).
	\"K-80Mpc\": Kagra for O5 simullations 80 Mpc (https://dcc.ligo.org/LIGO-T2000012/public).
	\"V1-O5\": Virgo target sensitivity O5 low noise (https://dcc.ligo.org/LIGO-T2000012/public).
";


ArmDirection::usage="```ArmDirection[{\[CurlyPhi], \[Lambda]}, {{\[Omega]1, \[Psi]1}, {\[Omega]2, \[Psi]2}}]``` returns ```{n1,n2}``` where ```ni = {nix,niy,niz}```, the unit vector in the direction of the arm i.
\[CurlyPhi]: vertex latitude [radians]
\[Lambda]: vertex longitude [radians]
\[Omega]i: Altitude of the arm i [radians]
\[Psi]i: Azimuth of arm i [radians] 

coordinates follow WGS-84 model used in \!\(\*TemplateBox[{\"\\\"LALDetectors.h\\\"\", \"https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html\"},\n\"HyperlinkURL\"]\)";


FromGeocentricCoordinates::usage="```FromGeocentricCoordinates[{h, \[CurlyPhi], \[Lambda]}]``` returns the {x,y,z} position of the detector.

h: displacement along the local vertical [meters]
\[CurlyPhi]: latitude [radians]
\[Lambda]: longitude [radians]

coordinates follow WGS-84 model used in \!\(\*TemplateBox[{\"\\\"LALDetectors.h\\\"\", \"https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_detectors_8h_source.html\"},\n\"HyperlinkURL\"]\)";


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Subsection:: *)
(*ArmDirection and DetectorPosition*)


ArmDirection//Clear

ArmDirection[{\[CurlyPhi]_, \[Lambda]_}, {{\[Omega]1_, \[Psi]1_}, {\[Omega]2_, \[Psi]2_}}] := Module[
	{comps1,comps2, n1, n2, \[Theta] = \[Pi]/2 - \[CurlyPhi], \[Phi] = \[Lambda], rhat, \[Phi]hat, \[Theta]hat},
	
	(* compsi: {r, \[Theta], \[Phi]} components in the base {rhat, \[Theta]hat, \[Phi]hat} of armi*)
	comps1 = {Sin[\[Omega]1], - Cos[\[Omega]1] Cos[\[Psi]1], Cos[\[Omega]1] Sin[\[Psi]1]}; 
	comps2 = {Sin[\[Omega]2], - Cos[\[Omega]2] Cos[\[Psi]2], Cos[\[Omega]2] Sin[\[Psi]2]}; 
	
	(*rhat, \[Theta]hat, \[Phi]hat components in {x,y,z}*)
	rhat = {Sin[\[Theta]] Cos[\[Phi]], Sin[\[Theta]] Sin[\[Phi]], Cos[\[Theta]]};
	\[Theta]hat = {Cos[\[Theta]] Cos[\[Phi]], Cos[\[Theta]] Sin[\[Phi]], - Sin[\[Theta]]};
	\[Phi]hat = {-Sin[\[Phi]], Cos[\[Phi]], 0};
	
	
	n1 = comps1[[1]] rhat + comps1[[2]] \[Theta]hat + comps1[[3]]*\[Phi]hat;
	n2 = comps2[[1]] rhat + comps2[[2]] \[Theta]hat + comps2[[3]]*\[Phi]hat;
	
	
	{n1, n2}//N
]


FromGeocentricCoordinates//Clear

FromGeocentricCoordinates[{h_, \[CurlyPhi]_, \[Lambda]_}] := Module[
	{R, a =6378137 , b = 6356752.314, X,Y,Z},
	
	R = a^2/Sqrt[a^2 Cos[\[CurlyPhi]]^2 + b^2 (Sin[\[CurlyPhi]]^2) ];
	
	X = (R+h) Cos[\[CurlyPhi]] Cos[\[Lambda]];
	Y = (R+h) Cos[\[CurlyPhi]] Sin[\[Lambda]];
	Z =  (0.993306 R +h) Sin[\[CurlyPhi]];
	
	{X,Y,Z}
]


Protect[ArmDirection]; Protect[FromGeocentricCoordinates];


(* ::Subsection:: *)
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


DetectorTensor["I1"] = Module[
	{nx={0.38496278183, -0.39387275094, 0.83466634811}, ny={0.89838844906, -0.04722636126, -0.43665531647}},
	1/2 (nx\[TensorProduct]nx-ny\[TensorProduct]ny)
];


Module[
	{
		lat = (40 + 31/60) Degree, lon = (9 + 25/60) Degree, n1, n2
	}, 
	
	(*sardinia L shape detector*)
	{n1, n2} = ArmDirection[
		{lat, lon}, 
		{{0, 0}, {0, \[Pi]/2}}
	];
	
	DetectorTensor["ET-S-L"] = (n1\[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
]


?ArmDirection


(* ::Text:: *)
(*triangle + L-shape detectors below in Sardinia:*)


Module[
	{
		lat = (40 +31/60) Degree, lon = (9 + 25/60) Degree, n1, n2
	}, 
	
	(*ET1 triangle detector*)
	{n1, n2} = ArmDirection[
		{lat, lon}, 
		{{0, 30 Degree}, {0, \[Pi]/2}}
	];
	
	DetectorTensor["ET1-S"] = (n1\[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
	
	(*ET2 triangle detector*)
	{n1, n2} = ArmDirection[
		{lat, lon}, 
		{{0, 150 Degree}, {0, 210 Degree}}
	];
	
	DetectorTensor["ET2-S"] = (n1\[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
	
	(*ET3 triangle detector*)
	{n1, n2} = ArmDirection[
		{lat, lon}, 
		{{0, 270 Degree}, {0, 330 Degree}}
	];
	
	DetectorTensor["ET3-S"] = (n1\[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
]


(* ::Text:: *)
(*Below you have the Meuse-Rhine L-shaped detectors, \[Alpha]=0  \[And] \[Alpha]=45*)


Module[
	{lat, lon, n1, n2},
	
	lat = (50 + 43/60 + 23/3600) Degree;
	lon = (5 + 55/60 + 14/3600) Degree;
	
	(*Meuse-Rhine L-shaped aligned*)
	{n1, n2} = ArmDirection[
		{lat, lon},
		{{0,0}, {0, \[Pi]/2}}
	];
	
	DetectorTensor["ET-MR-L-0"] = (n1 \[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
	
	
	(*Meuse-Rhine L-shaped \[Alpha]=45*)
	{n1, n2} = ArmDirection[
		{lat, lon},
		{{0, 45 Degree}, {0, 135 Degree}}
	];
	
	DetectorTensor["ET-MR-L-45"] = (n1 \[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
]


(* ::Text:: *)
(*Below we give CE Idaho and CE New Mexico: *)


Module[
	{lat, lon, n1,n2, \[Gamma]dx, \[Gamma]dy, \[Psi]y, \[Psi]x},
	(*CE-IDAHO*)
	
	lat =  0.764918; lon = \[Minus]1.969170;
	\[Gamma]dy = 0; \[Gamma]dx = \[Gamma]dy-\[Pi]/2;
	\[Psi]y = \[Pi]/2-\[Gamma]dy; \[Psi]x = \[Pi]/2 - \[Gamma]dx;
	
	{n1, n2} = ArmDirection[{lat, lon}, {{0, \[Psi]x}, {0, \[Psi]y}}];
	
	DetectorTensor["CE-I"] = (n1\[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
	
	
	
	lat = 0.578751; lon =\[Minus]1.858430;
	\[Gamma]dy = \[Minus]1.047200; \[Gamma]dx = \[Gamma]dy-\[Pi]/2;
	\[Psi]y = \[Pi]/2-\[Gamma]dy; \[Psi]x = \[Pi]/2 - \[Gamma]dx;
	
	{n1, n2} = ArmDirection[{lat, lon}, {{0, \[Psi]x}, {0, \[Psi]y}}];
	
	DetectorTensor["CE-NM"] = (n1\[TensorProduct]n1 - n2\[TensorProduct]n2)/2;
]


Protect[DetectorTensor];


(* ::Subsection::Closed:: *)
(*DetectorVertex*)


DetectorVertex = <||>;


DetectorVertex["H1"] = {-2.16141492636 10^6,  -3.83469517889 10^6   , 4.60035022664 10^6};
DetectorVertex["L1"] = {-7.42760447238 10^4, -5.49628371971 10^6 ,  3.22425701744 10^6 };
DetectorVertex["V1"] = {4.54637409900 10^6, 8.42989697626 10^5, 4.37857696241 10^6};
DetectorVertex["K"] = {-3777336.024,  3484898.411, 3765313.697};
DetectorVertex["ET1"] = {4.54637409900 10^6, 8.42989697626 10^5,4.37857696241 10^6};
DetectorVertex["ET2"] = {4.53936951685 10^6, 8.45074592488 10^5, 4.38540257904 10^6};
DetectorVertex["ET3"] = { 4.54240595075 10^6, 8.35639650438 10^5, 4.38407519902 10^6};
DetectorVertex["I1"] = {1.34897115479 10^6, 5.85742826577 10^6, 2.12756925209 10^6};
DetectorVertex["ET-S"] = {4.790201377886417`*^6,794444.3662803674`,4.121768605202005`*^6};
DetectorVertex["ET-MR"] = {4.0243452115988736`*^6,417334.87759974424`,4.914100068789029`*^6};


Module[
	{lat, lon},
	(*CE-IDAHO*)
	
	lat =  0.764918; lon = \[Minus]1.969170;
	
	DetectorVertex["CE-I"] = FromGeocentricCoordinates[{0, lat, lon}];
	
	lat = 0.578751; lon =\[Minus]1.858430;
	
	DetectorVertex["CE-NM"] = FromGeocentricCoordinates[{0, lat, lon}];
]


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
	data = Import[FileNameJoin[{ASDdir,"ET-0000A-18_ETDSensitivityCurveTxtFile.txt"}], "Data"]//N;
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
	
	(*ET 10 km arms*)
	data = Import[FileNameJoin[{ASDdir, "18213_ET10kmcolumns.txt"}], "Data"]//N;
	data[[All,{2,3,4}]] = data[[All,{2,3,4}]]//Sqrt; (*These files carry the PSD, not ASD*)
	
		(*HF only*)
		ASD["ET-10-hf"] = Interpolation[data[[All,{1,2}]], InterpolationOrder->1];
		(*LF only*)
		ASD["ET-10-lf"] = Interpolation[data[[All,{1,3}]], InterpolationOrder->1];
		(*HFLF*)
		ASD["ET-10-lfhf"] = Interpolation[data[[All,{1,4}]], InterpolationOrder->1];
		
	
	(*ET 15 km arms*)
	data = Import[FileNameJoin[{ASDdir, "18213_ET15kmcolumns.txt"}], "Data"]//N;
	data[[All,{2,3,4}]] = data[[All,{2,3,4}]]//Sqrt; (*These files carry the PSD, not ASD*)
		(*HF only*)
		ASD["ET-15-hf"] = Interpolation[data[[All,{1,2}]], InterpolationOrder->1];
		(*LF only*)
		ASD["ET-15-lf"] = Interpolation[data[[All,{1,3}]], InterpolationOrder->1];
		(*HFLF*)
		ASD["ET-15-lfhf"] = Interpolation[data[[All,{1,4}]], InterpolationOrder->1];
		
	(*ET 20 km arms*)
	data = Import[FileNameJoin[{ASDdir, "18213_ET20kmcolumns.txt"}], "Data"]//N;
	data[[All,{2,3,4}]] = data[[All,{2,3,4}]]//Sqrt; (*These files carry the PSD, not ASD*)
		(*HF only*)
		ASD["ET-20-hf"] = Interpolation[data[[All,{1,2}]],InterpolationOrder->1];
		(*LF only*)
		ASD["ET-20-lf"] = Interpolation[data[[All,{1,3}]], InterpolationOrder->1];
		(*HFLF*)
		ASD["ET-20-lfhf"] = Interpolation[data[[All,{1,4}]], InterpolationOrder->1];
		
	(*A+Design LIGO 05*)
	data = Import[FileNameJoin[{ASDdir, "AplusDesign.txt"}], "Data"]//N;
	ASD["L1H1-O5"] = Interpolation[data];
	
	(*kagra 80 mpc*)
	data = Import[FileNameJoin[{ASDdir, "kagra_80Mpc.txt"}], "Data"]//N;
	ASD["K-80Mpc"] = Interpolation[data];
	
	(*Virgo low noise*)
	data = Import[FileNameJoin[{ASDdir, "avirgo_O5low_NEW.txt"}], "Data"]//N;
	ASD["V1-O5"] = Interpolation[data];
	
	
		
]


Import[FileNameJoin[{ASDdir, "avirgo_O5low_NEW.txt"}], "Data"]


Protect[ASD];


(* ::Section::Closed:: *)
(*Package Footer*)


End[];


EndPackage[]
