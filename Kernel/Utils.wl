(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`Utils`", {"FelipeBarbosa`SymDALI`DerivativeTools`"}];


Unprotect[
	PatternFunctions, hphcIMRPhenomD, EmptyAssociation, SNR, hphcIMRPhenomHM
];

ClearAll[
	PatternFunctions, hphcIMRPhenomD, EmptyAssociation, SNR, hphcIMRPhenomHM
];


PatternFunctions::usage="PatternFunctions[f, \[Delta], \[Phi], \[Psi], pi, Dij] returns the pattern functions 
{\!\(\*SubscriptBox[\(F\), \(+\)]\), \!\(\*SubscriptBox[\(F\), \(x\)]\)}.

f: frequency vector, f[[i]] is in [Hz];
\[Delta]: declination \[Element] [-\[Pi]/2, \[Pi]/2];
\[Phi]: azimuthal angle \[Element] [0, 2 \[Pi]];
\[Psi]: polarization angle \[Element] [0, \[Pi]];
pi: position vector of the detector, see ```DetectorVertex```;
Dij: 3x3 detector tensor, see ```DetectorTensor```.";


hphcIMRPhenomD::usage="hphcIMRPhenomD[f, \[ScriptCapitalM]c, q, s1z, s2z, \[Iota], tc, \[Phi]ref, invdL] returns the polarizations
{\!\(\*SubscriptBox[\(h\), \(+\)]\), \!\(\*SubscriptBox[\(h\), \(x\)]\)} of the \"IMRPhenomD\" approximant.

f: frequency vector, f[[i]] is in [Hz];
\[ScriptCapitalM]c: chirp mass, [solar mass];
q: mass ratio, (0,1];
s1z: aligned spin \[Element] (-1,1);
s2z: aligned spin \[Element] (-1,1);
\[Iota]: inclination angle \[Element] [0, \[Pi]];
tc: coalescence time [s];
\[Phi]ref: reference phase \[Element] [0, 2 \[Pi]];
invdL: inverse of luminosity distance [1/Gpc].";


hphcIMRPhenomHM::usage="hphcIMRPhenomHM[f, \[ScriptCapitalM]c, q, s1z, s2z, \[Iota], tc, \[Phi]ref, invdL] returns the polarizations
{\!\(\*SubscriptBox[\(h\), \(+\)]\), \!\(\*SubscriptBox[\(h\), \(x\)]\)} of the \"IMRPhenomHM\" approximant.

f: frequency vector, f[[i]] is in [Hz];
\[ScriptCapitalM]c: chirp mass, [solar mass];
q: mass ratio, (0,1];
s1z: aligned spin \[Element] (-1,1);
s2z: aligned spin \[Element] (-1,1);
\[Iota]: inclination angle \[Element] [0, \[Pi]];
tc: coalescence time [s];
\[Phi]ref: reference phase \[Element] [0, 2 \[Pi]];
invdL: inverse of luminosity distance [1/Gpc].";


EmptyAssociation::usage="EmptyAssociation[name] returns an association with the correct syntax for the 
```DALITensors``` and ````SNR``` functions.

name: \"Aligned\" or \"Detector\"";


SNR::usage="SNR[Appr, fp, {det1, det2,...}] returns the network SNR of detectors \"{det1, det2,...}\"
for the fiducial \"fp\" using the waveform approximant \"Appr\".

Appr: \"IMRPhenomD\" or \"IMRPhenomHM\";
fp: fiducial in the form of an Association see ```EmptyAssociation[\"Aligned\"]```;
det_i: detector information in the form of an Association, see ```EmptyAssociation[\"Detector\"]```.

Options: 

\"res\": resolution of the frequency grid, defaults to 1000;
\"fmin\": smallest frequency in the grid [Hz], defaults to 10;
\"fmax\": largest frequency in the grid [Hz], defaults to 1024;
\"AllSNRs\": whether to return the SNR of different detectors separatelly (boolean), defaults to False.";


Begin["Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*SymRules, NRules and iFunctions for Fisher calculation*)


NRules = <||>;


(*Module[
	{symPv2, symD, symFpFc, nPv2, nD, nFpFc, m},
	m = Quiet[
		DerivativeRulesLoad/@{(*"IMRPhenomPv2",*) "IMRPhenomD", "Detectors"},
		{Part::partw, Part::take}
	];
	
	
	{(*NRules["IMRPhenomPv2"],*) NRules["IMRPhenomD"], NRules["FpFc"]} = m[[All,2]];
	
];*)


FpFcHead = FelipeBarbosa`SymDALI`DALICoefficients`Private`name;

D\[ScriptA]IMR = FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[ScriptCapitalA]IMR;
D\[CapitalPhi]IMR = FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[CapitalPsi];

HMhphc := FelipeBarbosa`SymDALI`DALICoefficients`Private`HMhphc;


ClearAll@NRules;


(* ::Subsection::Closed:: *)
(*(Public) PatternFunctions*)


PatternFunctions[$f_, $\[Delta]_, $\[Phi]_, $\[Psi]_, $pi_, $Dij_]/;(
	VectorQ[$f, NumberQ] && MatrixQ[$Dij, NumberQ] && VectorQ[$pi, NumberQ] && 
	-\[Pi]/2 <= $\[Delta] <= \[Pi]/2 && 0<=$\[Phi]<=2 \[Pi] && 0<= $\[Psi] <= \[Pi] 
):= Module[
	{LI = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]], \[Theta] = \[Pi]/2-$\[Delta], dij},
	
	dij = Extract[$Dij, LI];
	
	FpFcHead[$f, \[Theta], $\[Phi], $\[Psi], $pi, dij]	
]


Protect[PatternFunctions];


(* ::Subsection::Closed:: *)
(*(Public) hphc's*)


(* ::Subsubsection::Closed:: *)
(*PhenomD*)


Options[hphcIMRPhenomD] = {
	"\[Delta]\[CurlyPhi]-2"->0, "\[Delta]\[CurlyPhi]0"->0, "\[Delta]\[CurlyPhi]1"->0,"\[Delta]\[CurlyPhi]2"->0,"\[Delta]\[CurlyPhi]3"->0,"\[Delta]\[CurlyPhi]4"->0,"\[Delta]\[CurlyPhi]5l"->0,"\[Delta]\[CurlyPhi]6"->0,"\[Delta]\[CurlyPhi]6l"->0,"\[Delta]\[CurlyPhi]7"->0
};


hphcIMRPhenomD[$f_, $\[ScriptCapitalM]c_, $q_, $s1z_, $s2z_, $\[Iota]_, $tc_, $\[Phi]ref_, $invdL_, OptionsPattern[]]/;(
	VectorQ[$f, NumberQ] && Min[$f] > 0 && $\[ScriptCapitalM]c>0 && -1<=$s1z<=1 && -1<=$s2z<=1 &&
	0<=$\[Iota]<=\[Pi] && 0<=$\[Phi]ref<= 2 \[Pi] && $invdL>0
) := Module[
	{
		\[Omega], \[Omega]ref, G = 4.9254664969309`3.6105383994801805*^-6,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7, auxh,
		Amp, \[Phi], $M, $\[Chi]s, $\[Chi]a, $\[Eta]
	},
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7} = OptionValue[hphcIMRPhenomD, Options[hphcIMRPhenomD][[All,1]]];
	
	$\[Eta] = $q/(1+$q)^2;
	$M = $\[ScriptCapitalM]c ($q/(1+$q)^2)^(-3/5);
	
	\[Omega] = $f G $M;
	\[Omega]ref = Min[$f] G $M;
	$\[Chi]s = ($s1z+$s2z)/2;
	$\[Chi]a = ($s1z-$s2z)/2;
	
	
	auxh = Exp[-I (2 \[Pi] $f $tc - 2 $\[Phi]ref)]*$invdL*$M^2;
	
	Amp = D\[ScriptA]IMR[\[Omega], $\[Eta], $\[Chi]s, $\[Chi]a, $\[Iota]];
	
	\[Phi] = D\[CapitalPhi]IMR[
		\[Omega], \[Omega]ref, $\[Eta], $\[Chi]s, $\[Chi]a,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7
	];
	
	{
		Amp[[All, 1]] Exp[-I \[Phi]]*auxh, 
		Amp[[All, 2]] Exp[-I \[Phi]]*auxh
	}
	
]


Protect[hphcIMRPhenomD];


(* ::Subsubsection::Closed:: *)
(*PhenomHM*)


hphcIMRPhenomHM[
	$f_, $\[ScriptCapitalM]c_, $q_, $s1z_, $s2z_, $\[Iota]_, $tc_, $\[Phi]ref_, $invdL_
] := Module[
	{G = 4.925490947641267`*^-6, \[Omega], \[Omega]ref, M, hphc, $\[Chi]s, $\[Chi]a, $\[Eta]},
	
	$\[Eta] = $q/($q+1)^2;
	
	M = $\[ScriptCapitalM]c $\[Eta]^(-3/5);
	\[Omega] = $f M G; 
	\[Omega]ref = Min[\[Omega]];
	$\[Chi]s = ($s1z+ $s2z)/2;
	$\[Chi]a = ($s1z - $s2z)/2;
	
	hphc = M^2 HMhphc[\[Omega], \[Omega]ref, $\[Eta], $\[Chi]s, $\[Chi]a, $\[Iota], $\[Phi]ref]*Exp[-I 2 \[Pi] $f $tc]*$invdL;
	
	{
		hphc[[All,1]],
		hphc[[All,2]]
	}
];


Protect[hphcIMRPhenomHM];


(* ::Subsection::Closed:: *)
(*(Public) EmptyAssociation*)


EmptyAssociation["Aligned"] = <|
	"dec"->1,"\[Phi]"->2,"\[Psi]"->3,"\[ScriptCapitalM]c"->4,"q"->5,"s1z"->6,"s2z"->7,"\[Iota]"->8,"1/dL"->9,"tc"->10,"\[Phi]ref"->11
|>;

(*EmptyAssociation["Precessing"] = <|
	"\[Theta]"->1,"\[Phi]"->2,"\[Psi]"->3,"m1"->4,"m2"->5,"s1x"->6,"s1y"->7,"s1z"->8,"s2x"->9,"s2y"->10,"s2z"->11,"\[Iota]"->12,"dL"->13,"tc"->14,"\[Phi]ref"->15
|>;*)

EmptyAssociation["Detector"] = <|
	"Position"->{0,0,0},
	"DetectorTensor"-> ConstantArray[0, {3,3}],
	"ASD"->Interpolation[Table[{i,1}, {i, 1, 5000}]]
|>;


Protect[EmptyAssociation]


(* ::Subsection:: *)
(*SNR*)


(* ::Subsubsection::Closed:: *)
(*RetrieveFiducial*)


vars["Aligned"] = {
	"\[ScriptCapitalM]c", "q", "s1z", "s2z",
	"\[Iota]", "dec", "\[Phi]", "\[Psi]", "1/dL", 
	"tc", "\[Phi]ref"
};


GenErrorMessage["Aligned", True] := Null;

GenErrorMessage["Aligned", False] := Throw[
"\n
PhenomD variables should satisfy the following conditions: \n
0 < \[ScriptCapitalM]c < \[Infinity] && 0 q <= 1 && -1 \[LessEqual] s1z\[LessEqual]1 && -1\[LessEqual]s2z\[LessEqual]1 && 0\[LessEqual]\[Iota]\[LessEqual]\[Pi] && -\[Pi]/2\[LessEqual]dec\[LessEqual]\[Pi]/2 && \n
0\[LessEqual]\[Phi]\[LessEqual]2\[Pi] && 0\[LessEqual]\[Psi]\[LessEqual]\[Pi] && 0 < 1/dL < \[Infinity] && -\[Infinity]\[LessEqual]tc\[LessEqual]\[Infinity] && 0\[LessEqual]\[Phi]ref\[LessEqual]2 \[Pi]
\n"
];


RetrieveFiducial["Aligned", fp_Association]/;(
	(*there should not be more than 12 variables*)
	Length@Keys[fp] <= 26 && 
	(*the keys must be contained in vars["IMRPhenomD"]*)
	ContainsAll[vars["Aligned"],Keys[fp]]
) := Module[
	{\[ScriptCapitalM]c, q, s1z, s2z, \[Iota], dec, \[Phi], \[Psi], invdL, tc, \[Phi]ref, \[Delta]p, test, res},
	
	{\[ScriptCapitalM]c, q, s1z, s2z, \[Iota], dec, \[Phi], \[Psi], invdL, tc, \[Phi]ref} = fp/@{
		"\[ScriptCapitalM]c", "q", "s1z", "s2z",
		"\[Iota]", "dec", "\[Phi]", "\[Psi]", "1/dL", 
		"tc", "\[Phi]ref"
	};
	
	
	(*Test variables:*)
	test = Thread@LessEqual[
		{10^-6, 0, -1, -1, 0, -\[Pi]/2, 0,0, 10^-11, -Infinity, 0},
		{\[ScriptCapitalM]c, q, s1z, s2z, \[Iota], dec, \[Phi], \[Psi], invdL, tc, \[Phi]ref},
		{Infinity, 1, 1,1, \[Pi], \[Pi]/2, 2 \[Pi], \[Pi], Infinity,  Infinity, 2 \[Pi]}
	];
	
	(*this collapses to True or False*)
	test = (And@@test)//TrueQ;
	
	GenErrorMessage["Aligned", test];
	
	If[
		Length[Keys[fp]] === 11,
		res = {{\[Pi]/2-dec, \[Phi], \[Psi]}, {\[ScriptCapitalM]c, q, s1z, s2z, \[Iota], tc, \[Phi]ref, invdL}},
		
		res = {{\[Pi]/2- dec, \[Phi], \[Psi]}, {\[ScriptCapitalM]c, q, s1z, s2z, \[Iota], tc, \[Phi]ref, invdL}}
	]
	(*$\[ScriptCapitalM]c_, $q_, $s1z_, $s2z_, $\[Iota]_, $tc_, $\[Phi]ref_, $invdL*)
]


LIDij = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]];


(* ::Subsubsection:: *)
(*SNR*)


ClearAll@idot

idot[m1_, m2_] := m1[[1]] m2[[1]] + m1[[2]] m2[[2]]


convert["IMRPhenomD"]  = "Aligned"
convert["IMRPhenomHM"]  = "Aligned"


testAligned["IMRPhenomD"] = True;
testAligned["IMRPhenomHM"] = True;


iSNR[
	approximant_String, fiducialPoint_Association, 
	{detectors__Association},
	fmin_, fmax_, res_, AllFisherMatrices_]/;(
	approximant === "IMRPhenomD" ||approximant ===  "IMRPhenomHM" 
) := Module[
	{
		FiducialFpFc,  Fiducialhphc, iFpFcHead = FpFcHead, detecs, fpDetectors, 
		FpFcDetectors, IFMAX, \[CapitalDelta]f, fvec, totalM, hphc, M, \[Eta], hs, PSDs, SNRs, AlignedOrP = convert[approximant],
		isAligned = testAligned[approximant]
	},
	
	
	(*get fiducial points*)
	{FiducialFpFc,  Fiducialhphc} = With[
	
		{l = RetrieveFiducial[AlignedOrP, fiducialPoint]},
		If[
			Cases[l, Missing, Infinity, Heads->True]==={},
			l,
			Throw["Wrong variables for the fiducial point"]
		]
	];
	
	totalM = If[
		isAligned, 
		Fiducialhphc[[1]] (Fiducialhphc[[2]]/(Fiducialhphc[[2]]+1)^2)^(-3/5)
	];
	
	IFMAX = Min@{fmax, 0.2/(4.9254664969309`3.6105383994801805*^-6 totalM)}//Round;
	
	\[CapitalDelta]f = (IFMAX-fmin)/(res-1);
	
	fvec = Range[fmin, IFMAX, \[CapitalDelta]f];
	
	(*non-zero values for \[Delta]p do not affect the SNR:*)
	
	Which[
		approximant==="IMRPhenomD", 
		Fiducialhphc = Fiducialhphc[[1;;8]];
		hphc = hphcIMRPhenomD@@Join[{fvec}, Fiducialhphc],
		
		approximant==="IMRPhenomHM",
		Fiducialhphc = Fiducialhphc[[1;;8]];
		hphc = hphcIMRPhenomHM@@Join[{fvec}, Fiducialhphc]
	];
	(*D: {\[ScriptCapitalM]c, q, s1z, s2z, \[Iota], tc, \[Phi]ref, 1/dL}*)
	
	
	
	(*Check detector keys*)
	If[
		Cases[{#["Position"], #["DetectorTensor"], #["ASD"]}&/@{detectors}, Missing, Infinity, Heads->True]==={},
		0,
		Throw["Wrong Detector Keys"]
	];
	
	fpDetectors = Join[
		FiducialFpFc, 
		{#["Position"], Extract[#["DetectorTensor"], LIDij]}
	]&/@{detectors};
	
	(*this should be an array of real numbers*)
	If[
		VectorQ[Flatten[fpDetectors],RealValuedNumberQ]===False,
		Throw["DetectorTensor and Position should contain only real numbers"]
	];
	
	FpFcDetectors = iFpFcHead[fvec, Sequence@@#]&/@fpDetectors;
	
	(*all hs*)
	hs = idot[hphc, #]&/@FpFcDetectors;
	
	PSDs = (#["ASD"][fvec])&/@{detectors}; 
	PSDs = PSDs^2;
	
	SNRs = MapThread[
		4 \[CapitalDelta]f Re[(Conjugate[#1]*#1)/#2//Total]&,
		{hs, PSDs}
	];
	
	
	If[
		AllFisherMatrices===True,
		Sqrt[SNRs],
		Sqrt[Total[SNRs]]
	]
		
]


Options[SNR] = {
	"res" -> 1000,
	"fmin" -> 10,
	"fmax" -> 1024,
	"AllSNRs" -> False
};


SNR[x__, OptionsPattern[]] := Module[
	{fmin, fmax, res, AllFisherMatrices},
	{fmin, fmax, res, AllFisherMatrices} = OptionValue[SNR, #]&/@{"fmin", "fmax", "res", "AllSNRs"};
	Catch[iSNR[x, fmin, fmax, res, AllFisherMatrices]]
]


Protect[SNR]


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
