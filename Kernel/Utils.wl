(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`Utils`", {"FelipeBarbosa`SymDALI`DerivativeTools`"}];


Unprotect[
	PatternFunctions, hphcIMRPhenomD, hphcIMRPhenomPv2,EmptyAssociation, SNR
];

ClearAll[
	PatternFunctions, hphcIMRPhenomD, hphcIMRPhenomPv2,EmptyAssociation, SNR
];


PatternFunctions::usage="```PatternFunctions[f, \[Delta], \[Alpha], \[Psi], GMST, pi, Dij]``` returns the pattern functions 
{\!\(\*SubscriptBox[\(F\), \(+\)]\), \!\(\*SubscriptBox[\(F\), \(x\)]\)} in the form of a matrix with dimensions {2, Length[f]}. 
The definition includes the time delay from the detector to the Earth's center.

f: frequency vector, {f1,f2,f3,...} [Hz].
\[Delta]: declination in the range [-\[Pi]/2, \[Pi]/2]  (related to polar angle \[Theta] by \[Theta] = \[Pi]/2 - \[Delta]).
\[Alpha]: right ascension in the range [0, 2 \[Pi]]  (related to the azimuthal angle \[Phi] by \[Phi] = \[Alpha] - GMST).
\[Psi]: polarization angle in the range [0, \[Pi]].
GMST: Greenwich mean sideral time in the range [0, 2 \[Pi]]. 
pi: position vector, {x,y,z} [meters], of the detector vertex  (same coordinate system of  ```DetectorVertex```).
Dij: detector tensor, \!\(\*SubscriptBox[\(D\), \(ij\)]\) = 0.5 (\!\(\*SubscriptBox[\(n\), \(i\)]\) \!\(\*SubscriptBox[\(n\), \(j\)]\) - \!\(\*SubscriptBox[\(n\), \(j\)]\) \!\(\*SubscriptBox[\(n\), \(i\)]\)), in the form of a 3x3 matrix   (same coordinate system of  ```DetectorTensor```).";


hphcIMRPhenomD::usage="```hphcIMRPhenomD[f, M, \[Eta], s1z, s2z, \[Iota], tc, \[Phi]ref, dL]``` returns the plus and cross 
polarizations {\!\(\*SubscriptBox[\(h\), \(+\)]\), \!\(\*SubscriptBox[\(h\), \(x\)]\)} in the form of a matrix with dimensions {2, Length[f]} for the IMRPhenomD approximant.

f: frequency vector, {f1,f2,f3,...} [Hz].
M: total mass [solar masses] (M>0).
\[Eta]: symmetric mass ratio, in the range (0, 0.25].
s1z: dimensionless spin component (direction of orbital \!\(\*OverscriptBox[\(L\), \(->\)]\)) of the heaviest BH, in the range [-1,1].
s2z: dimensionless spin component (direction of orbital \!\(\*OverscriptBox[\(L\), \(->\)]\)) of the lightest BH, in the range [-1,1].
\[Iota]: inclination angle, in the range [0, \[Pi]].
tc: coalescence time [seconds].
\[Phi]ref: reference frequency, in the range [0, 2 \[Pi]].
dL: luminosity distance [Gpc].

Obs.: GR deviation parameters from the ```TIGER``` framework can be passed in the form of Options: 
```hphcIMRPhenomD[f, M, \[Eta], s1z, s2z, \[Iota], tc, \[Phi]ref, dL, \"\[Delta]\[CurlyPhi]-2\"->1, \"\[Delta]\[Beta]3\"->1]```.
By default all deviation parameters are 0.";





(*hphcIMRPhenomPv2::usage="```hphcIMRPhenomPv2[f, M, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], tc, \[Phi]ref, dL]``` returns 
the plus and cross polarizations {\!\(\*SubscriptBox[\(h\), \(+\)]\), \!\(\*SubscriptBox[\(h\), \(x\)]\)} in the form of a matrix with dimensions {2, Length[f]} for the 
IMRPhenomD approximant.

f: frequency vector, {f1,f2,f3,...} [Hz].
M: total mass [solar masses] (M>0).
\[Eta]: symmetric mass ratio, in the range (0, 0.25].
s1x: dimensionless spin component in the x diretion of the heaviest BH, in the range [-1,1].
s1y: dimensionless spin component in the y diretion of the heaviest BH, in the range [-1,1].
s1z: dimensionless spin component (direction of orbital \!\(\*OverscriptBox[\(L\), \(->\)]\)) of the heaviest BH, in the range [-1,1].
s2x: dimensionless spin component in the x direction of the lightest BH, in the range [-1,1].
s2y: dimensionless spin component in the y direction of the lightest BH, in the range [-1,1].
s2z: dimensionless spin component (direction of orbital \!\(\*OverscriptBox[\(L\), \(->\)]\)) of the lightest BH, in the range [-1,1].
\[Iota]: inclination angle, in the range [0, \[Pi]].
tc: coalescence time [seconds].
\[Phi]ref: reference frequency, in the range [0, 2 \[Pi]].
dL: luminosity distance [Gpc].

Obs.: GR deviation parameters from the ```TIGER``` framework can be passed in the form of Options: 
```hphcIMRPhenomPv2[f, M, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], tc, \[Phi]ref, dL, \"\[Delta]\[CurlyPhi]-2\"->1, \"\[Delta]\[Beta]3\"->1]```.
By default all deviation parameters are 0.";*)


EmptyAssociation::usage="```EmptyAssociation[\"name\"]``` returns an association with the correct variable syntax to be 
filled and feed to ```FisherMatrix```.

Available \"name\": 
	\"IMRPhenomPv2\".
	\"IMRPhenomD\".
	\"Detector\"
";


Begin["Private`"];


(* ::Section:: *)
(*Definitions*)


(* ::Subsection:: *)
(*SymRules, NRules and iFunctions for Fisher calculation*)


NRules = <||>;


Module[
	{symPv2, symD, symFpFc, nPv2, nD, nFpFc, m},
	m = Quiet[
		DerivativeRulesLoad/@{(*"IMRPhenomPv2",*) "IMRPhenomD", "Detectors"},
		{Part::partw, Part::take}
	];
	
	
	{(*NRules["IMRPhenomPv2"],*) NRules["IMRPhenomD"], NRules["FpFc"]} = m[[All,2]];
	
];


FpFcHead = NRules["FpFc"][FpFc][[1,2,0]];
D\[ScriptA]IMR = NRules["IMRPhenomD"][\[ScriptA]IMR][[1,2,0]];
D\[CapitalPhi]IMR = NRules["IMRPhenomD"][\[CapitalPhi]IMR][[1,2,0]];
(*Pv2\[ScriptA]IMR =  NRules["IMRPhenomPv2"][\[ScriptA]IMR][[1,2,0]];
Pv2\[CapitalPhi]IMR =  NRules["IMRPhenomPv2"][\[CapitalPhi]IMR][[1,2,0]];*)


ClearAll@NRules;


(* ::Subsection:: *)
(*(Public) PatternFunctions*)


PatternFunctions[$f_, $\[Delta]_, $\[Alpha]_, $\[Psi]_, $GMST_, $pi_, $Dij_]/;(
	VectorQ[$f, NumberQ] && MatrixQ[$Dij, NumberQ] && VectorQ[$pi, NumberQ] && 
	-\[Pi]/2 <= $\[Delta] <= \[Pi]/2 && 0<=$\[Alpha]<=2 \[Pi] && 0<= $\[Psi] <= \[Pi] 
):= Module[
	{LI = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]], \[Theta] = \[Pi]/2-$\[Delta], \[Phi] = $\[Alpha]-$GMST, dij},
	
	dij = Extract[$Dij, LI];
	
	FpFcHead[$f, \[Theta], \[Phi], $\[Psi], $pi, dij]	
]


Protect[PatternFunctions];


(* ::Subsection::Closed:: *)
(*(Public) hphc's*)


(* ::Subsubsection::Closed:: *)
(*PhenomD*)


Options[hphcIMRPhenomD] = {
	"\[Delta]\[CurlyPhi]-2"->0, "\[Delta]\[CurlyPhi]0"->0, "\[Delta]\[CurlyPhi]1"->0,"\[Delta]\[CurlyPhi]2"->0,"\[Delta]\[CurlyPhi]3"->0,"\[Delta]\[CurlyPhi]4"->0,"\[Delta]\[CurlyPhi]5l"->0,"\[Delta]\[CurlyPhi]6"->0,"\[Delta]\[CurlyPhi]6l"->0,"\[Delta]\[CurlyPhi]7"->0,
	"\[Delta]\[Beta]2"->0,"\[Delta]\[Beta]3"->0,
	"\[Delta]\[Alpha]2"->0,"\[Delta]\[Alpha]3"->0,"\[Delta]\[Alpha]4"->0
};


hphcIMRPhenomD[$f_, $\[ScriptCapitalM]c_, $\[Delta]_, $\[Chi]s_, $\[Chi]a_, $\[Iota]_, $tc_, $\[Phi]ref_, $invdL_, OptionsPattern[]]/;(
	VectorQ[$f, NumberQ] && Min[$f]>0 && $\[ScriptCapitalM]c>0 && 0 <= $\[Delta] < 1 && -1<=$\[Chi]s<=1 && -1<=$\[Chi]a<=1 &&
	0<=$\[Iota]<=\[Pi] && 0<=$\[Phi]ref<= 2 \[Pi] && $invdL>0
) := Module[
	{
		\[Omega], \[Omega]ref, G = 4.9254664969309`3.6105383994801805*^-6(*G/c^3 [s/solarMass]*), 
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4, auxh,
		Amp, \[Phi], $M = $\[ScriptCapitalM]c ((1-$\[Delta]^2)/4)^(-3/5)
	},
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = OptionValue[hphcIMRPhenomD, Options[hphcIMRPhenomD][[All,1]]];
	
	\[Omega] = $f G $M;
	\[Omega]ref = Min[$f] G $M;
	
	auxh = Exp[-I (2 \[Pi] $f $tc - 2 $\[Phi]ref)]*$invdL;
	
	Amp = D\[ScriptA]IMR[$f, $\[ScriptCapitalM]c, $\[Delta], $\[Chi]s, $\[Chi]a, $\[Iota]];
	
	\[Phi] = D\[CapitalPhi]IMR[
		\[Omega], \[Omega]ref, $\[Delta], $\[Chi]s, $\[Chi]a,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
	];
	
	{
		Amp[[1]] Exp[-I \[Phi]]*auxh, 
		Amp[[2]] Exp[-I \[Phi]]*auxh
	}
	
]


Protect[hphcIMRPhenomD];


(* ::Subsubsection::Closed:: *)
(*PhenomPv2*)


(*Options[hphcIMRPhenomPv2] = {
	"\[Delta]\[CurlyPhi]-2"->0, "\[Delta]\[CurlyPhi]0"->0, "\[Delta]\[CurlyPhi]1"->0,"\[Delta]\[CurlyPhi]2"->0,"\[Delta]\[CurlyPhi]3"->0,"\[Delta]\[CurlyPhi]4"->0,"\[Delta]\[CurlyPhi]5l"->0,"\[Delta]\[CurlyPhi]6"->0,"\[Delta]\[CurlyPhi]6l"->0,"\[Delta]\[CurlyPhi]7"->0,
	"\[Delta]\[Beta]2"->0,"\[Delta]\[Beta]3"->0,
	"\[Delta]\[Alpha]2"->0,"\[Delta]\[Alpha]3"->0,"\[Delta]\[Alpha]4"->0
};*)


(*hphcIMRPhenomPv2[$f_, $M_, $\[Eta]_, $s1x_, $s1y_, $s1z_, $s2x_, $s2y_, $s2z_, $\[Iota]_, $tc_, $\[Phi]ref_, $dL_, OptionsPattern[]]/;(
	VectorQ[$f, NumberQ] && Min[$f]>0 && $M>0 && 0 < $\[Eta] <= 0.25 && -1<=$s1z<=1 && -1<=$s2z<=1 &&
	0<=$\[Iota]<=\[Pi] && 0<=$\[Phi]ref<= 2 \[Pi] && $dL>0
) := Module[
	{
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4, 
		auxh,  Amp, \[Phi],hphc, m1 ,m2, fref
	},
	m1 = ($M/2) (1 + Sqrt[1 - 4 $\[Eta]]);
	m2 = ($M/2) (1 - Sqrt[1 - 4 $\[Eta]]);
	fref = Min[$f];
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = OptionValue[
		hphcIMRPhenomPv2, Options[hphcIMRPhenomPv2][[All,1]]
	];
	
	auxh = Exp[-I 2 \[Pi] $f $tc]/$dL;
	
	
	Amp = Pv2\[ScriptA]IMR[$f, fref, m1, m2, $s1x, $s1y, $s1z, $s2x, $s2y, $s2z, $\[Phi]ref, $\[Iota]];
	
	\[Phi] = Pv2\[CapitalPhi]IMR[
		$f, fref, m1, m2, $s1x, $s1y, $s1z, $s2x, $s2y, $s2z,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
	];
	
	hphc = Transpose[Amp]*Exp[-I \[Phi]]*auxh;
	
	Transpose[hphc]
]*)


(*Protect[hphcIMRPhenomPv2]*)


(* ::Subsection:: *)
(*(Public) EmptyAssociation*)


EmptyAssociation["Aligned"] = <|
	"\[Theta]"->1,"\[Phi]"->2,"\[Psi]"->3,"\[ScriptCapitalM]c"->4,"\[Delta]"->5,"\[Chi]s"->6,"\[Chi]a"->7,"\[Iota]"->8,"1/dL"->9,"tc"->10,"\[Phi]ref"->11
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


(* ::Subsubsection:: *)
(*RetrieveFiducial*)


vars["Aligned"] = {
	"\[ScriptCapitalM]c", "\[Delta]", "\[Chi]s", "\[Chi]a",
	"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "1/dL", 
	"tc", "\[Phi]ref",
	"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"
};


(*vars["IMRPhenomPv2"] = {
	"m1", "m2", "s1x","s1y","s1z","s2x", "s2y", "s2z",
	"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "dL",
	"tc", "\[Phi]ref",
	"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"
};*)


GenErrorMessage["Aligned", True] := Null;

GenErrorMessage["Aligned", False] := Throw[
"\n
PhenomD variables should satisfy the following conditions: \n
1\[LessEqual]M\[LessEqual]\[Infinity] && \!\(\*SuperscriptBox[\(10\), \(-4\)]\)\[LessEqual]\[Eta]\[LessEqual]0.25`&& -1\[LessEqual]s1z\[LessEqual]1 && -1\[LessEqual]s2z\[LessEqual]1 && 0\[LessEqual]\[Iota]\[LessEqual]\[Pi] && 0\[LessEqual]\[Theta]\[LessEqual]\[Pi] && \n
0\[LessEqual]\[Phi]\[LessEqual]2\[Pi] && 0\[LessEqual]\[Psi]\[LessEqual]\[Pi] && \!\(\*SuperscriptBox[\(10\), \(-11\)]\)\[LessEqual]dL\[LessEqual]\[Infinity] && -\[Infinity]\[LessEqual]tc\[LessEqual]\[Infinity] && 0\[LessEqual]\[Phi]ref\[LessEqual]2 \[Pi]
\n"
];

GenErrorMessage["Precessing", True] := Null;

GenErrorMessage["Precessing", False] := Throw[
"
\n
PhenomPv2 variables should satisfy the following conditions: \n
m1 >= m2 && Norm[{s1x, s1y,s1z}]<=1  && Norm[{s2x,s2y, s2z}]<=1 \n
1\[LessEqual] M \[LessEqual]\[Infinity] && \!\(\*SuperscriptBox[\(10\), \(-4\)]\)\[LessEqual]\[Eta]\[LessEqual] 0.25`&& -1\[LessEqual]s1x\[LessEqual]1 && -1\[LessEqual]s1y\[LessEqual]1 && -1\[LessEqual]s1z\[LessEqual]1 &&\n
-1\[LessEqual]s2x\[LessEqual]1 && -1\[LessEqual]s2y\[LessEqual]1 && -1\[LessEqual]s2z\[LessEqual]1 && 0\[LessEqual]\[Iota]\[LessEqual]\[Pi] && 0\[LessEqual]\[Theta]\[LessEqual]\[Pi] &&\n
0\[LessEqual]\[Phi]\[LessEqual]2\[Pi] && 0\[LessEqual]\[Psi]\[LessEqual]\[Pi] && \!\(\*SuperscriptBox[\(10\), \(-11\)]\)\[LessEqual]dL\[LessEqual]\[Infinity] && -\[Infinity]\[LessEqual]tc\[LessEqual]\[Infinity] && 0\[LessEqual]\[Phi]ref\[LessEqual]2 \[Pi] \n"
];


RetrieveFiducial["Aligned", fp_Association]/;(
	(*there should not be more than 12 variables*)
	Length@Keys[fp] <= 26 && 
	(*the keys must be contained in vars["IMRPhenomD"]*)
	ContainsAll[vars["Aligned"],Keys[fp]]
) := Module[
	{\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], \[Theta], \[Phi], \[Psi], invdL, tc, \[Phi]ref, \[Delta]p, test, res},
	
	{\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], \[Theta], \[Phi], \[Psi], invdL, tc, \[Phi]ref} = fp/@{
		"\[ScriptCapitalM]c", "\[Delta]", "\[Chi]s", "\[Chi]a",
		"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "1/dL", 
		"tc", "\[Phi]ref"
	};
	
	
	(*Test variables:*)
	test = Thread@LessEqual[
		{10^-3, 0, -1, -1, 0, 0, 0,0, 10^-11, -Infinity, 0},
		{\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], \[Theta], \[Phi], \[Psi], invdL, tc, \[Phi]ref},
		{Infinity, 0.9999, 1,1, \[Pi], \[Pi], 2 \[Pi], \[Pi], 10^20,  Infinity, 2 \[Pi]}
	];
	
	(*this collapses to True or False*)
	test = (And@@test)//TrueQ;
	
	GenErrorMessage["Aligned", test];
	
	If[
		Length[Keys[fp]] === 11,
		res = {{\[Theta], \[Phi], \[Psi]}, {\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], invdL, tc, \[Phi]ref}},
		
		\[Delta]p = DeleteElements[Keys[fp], vars["Aligned"][[1;;11]]];(*check for real value*)
		\[Delta]p = fp/@SortBy[\[Delta]p, order\[Delta]];
		
		If[VectorQ[\[Delta]p, RealValuedNumberQ] === False, Throw["\[Delta]pi values should be in the range -\[Infinity]< \[Delta]pi <\[Infinity]"]];
		
		res = {{\[Theta], \[Phi], \[Psi]}, {\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], invdL, tc, \[Phi]ref, Sequence@@\[Delta]p}}
	]
	
]


RetrieveFiducial["IMRPhenomPv2", fp_Association]/;(
	Length@Keys[fp] <= 30 && 
	
	ContainsAll[vars["IMRPhenomPv2"],Keys[fp]]
) := Module[
	{m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], \[Theta], \[Phi], \[Psi], dL, tc, \[Phi]ref, \[Delta]p, test, test2},
	
	{
	m1,m2, 
	s1x, s1y, s1z, 
	s2x, s2y, s2z, 
	\[Iota], \[Theta], \[Phi], \[Psi], dL, tc, \[Phi]ref} = fp/@{
		"m1", "m2", 
		"s1x", "s1y","s1z", 
		"s2x", "s2y", "s2z", 
		"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "dL", 
		"tc", "\[Phi]ref"
	};
	
	If[m1==m2, m2 = (0.9999) m1];
	
	test = Thread@LessEqual[
		{1, 1, Sequence@@ConstantArray[-1,6], 0,0,0,0, 10^-11, -Infinity, 0},
		{m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], \[Theta], \[Phi], \[Psi], dL, tc, \[Phi]ref},
		{Infinity, Infinity, Sequence@@ConstantArray[1, 6], \[Pi],\[Pi], 2 \[Pi], \[Pi], Infinity, Infinity, 2 \[Pi]}
	];
	
	test = (And@@test)//TrueQ;
	
	test2 = (m1 >= m2 && Norm[{s1x, s1y, s1z}] <= 1 && Norm[{s2x, s2y, s2z}] <= 1)//TrueQ;
	
	GenErrorMessage["IMRPhenomPv2", TrueQ[test&&test2]];
	
	
	If[
		Length[Keys[fp]] === 15,
		{{\[Theta], \[Phi], \[Psi]}, {m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota],dL, tc, \[Phi]ref}},
		\[Delta]p = DeleteElements[Keys[fp], vars["IMRPhenomPv2"][[1;;15]]]; 
		\[Delta]p = fp/@SortBy[\[Delta]p, order\[Delta]];
		If[VectorQ[\[Delta]p, RealValuedNumberQ]===False, Throw["\[Delta]pi value should be in the range -\[Infinity] < \[Delta]pi <\[Infinity]"]];
		{{\[Theta], \[Phi], \[Psi]}, {m1,m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota],dL, tc, \[Phi]ref, Sequence@@\[Delta]p}}
	]
]

RetrieveFiducial[x___] := Throw[$Failed, failTag[RetrieveFiducial]]


LIDij = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]];


(* ::Subsubsection:: *)
(*SNR*)


ClearAll@idot
idot[m1_, m2_] := m1[[1]] m2[[1]] + m1[[2]] m2[[2]]


convert["IMRPhenomD"]  = "Aligned"
convert["IMRPhenomPv2"] = "Precessing"


testAligned["IMRPhenomD"] = True;
testAligned["IMRPhenomPv2"] = False;


iSNR[
approximant_String, fiducialPoint_Association, 
{detectors__Association},
fmin_, fmax_, res_, AllFisherMatrices_]/;(
	approximant === "IMRPhenomD" ||approximant ===  "IMRPhenomPv2" 
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
	
	totalM = If[isAligned, Fiducialhphc[[1]] ((1-Fiducialhphc[[2]]^2)/4)^(-3/5), Fiducialhphc[[1]]+Fiducialhphc[[2]]];
	
	IFMAX = Min@{fmax, 0.2/(4.9254664969309`3.6105383994801805*^-6 totalM)}//Round;
	
	\[CapitalDelta]f = (IFMAX-fmin)/(res-1);
	
	fvec = Range[fmin, IFMAX, \[CapitalDelta]f];
	
	(*non-zero values for \[Delta]p do not affect the SNR:*)
	
	If[
		approximant==="IMRPhenomD", 
		Fiducialhphc = Fiducialhphc[[1;;8]];
		hphc = hphcIMRPhenomD@@Join[{fvec}, Fiducialhphc[[{1,2,3,4,5,7,8}]], Fiducialhphc[[{6}]]],
		Fiducialhphc = Fiducialhphc[[1;;12]];
		M = Fiducialhphc[[1]]+Fiducialhphc[[2]];
		\[Eta] = Fiducialhphc[[1]]*Fiducialhphc[[2]]/M^2;
		hphc = hphcIMRPhenomPv2@@Join[{fvec}, {M, \[Eta]}, Fiducialhphc[[{3,4,5,6,7,8,9,11,12}]], Fiducialhphc[[{10}]]]
	];
	(*PV2: {m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota],dL, tc, \[Phi]ref}*)
	(*D: {M, \[Eta], s1z, s2z, \[Iota], dL, tc, \[Phi]ref}*)
	
	
	
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
