(* ::Package:: *)

System`$D::usage="$D[{n__}, symbol][y__] is equivalent to Derivative[n__][symbol][y__].";
System`TagRule::usage="TagRule[f, g[f[x_]], 3] is equivalent in spirit to f/: g[f[x_]] = 3. At some point I will fix the formating so that you 
can type f/: g[f[x_]] -> 3. For now the syntax is TagRule[tag, lhs, rhs]";

BeginPackage["FelipeBarbosa`SymDALI`"];
ClearAll[us\[Theta]];
us\[Theta];
Begin["`Private`"]
Attributes@us\[Theta] = {Listable};
us\[Theta][x_]/; x<0 = 0; us\[Theta][x_]/; x==0 = 1/2; us\[Theta][x_]/; x>0 = 1;
Unprotect[Derivative, $D]; Clear[Derivative, $D];

Derivative[n_][us\[Theta]][x_]/;n>=1 := 0; $D[n_, us\[Theta]][x_] := 0

Protect[Derivative, $D];
End[];


EndPackage[];


<<FelipeBarbosa`SymDALI`DALICoefficients`;
<<FelipeBarbosa`SymDALI`DerivativeTools`;
<<FelipeBarbosa`SymDALI`Detectors`;
<<FelipeBarbosa`SymDALI`DALIPolynomial`;
<<FelipeBarbosa`SymDALI`Likelihood`;


(* ::Section::Closed:: *)
(*PackageHeader*)


BeginPackage["FelipeBarbosa`SymDALI`", {"FelipeBarbosa`SymDALI`DerivativeTools`"}];


Unprotect[
	SymRules, NRules,iFpFc,ihphcIMRPhenomD, ihphcIMRPhenomPv2, \[ScriptA]IMR, \[CapitalPhi]IMR, FpFc,(*1st subsection*)
	PatternFunctions, hphcIMRPhenomD, hphcIMRPhenomPv2 (*2o subsection*)
];


ClearAll[
	SymRules, NRules,iFpFc,ihphcIMRPhenomD, ihphcIMRPhenomPv2, \[ScriptA]IMR, \[CapitalPhi]IMR, FpFc,(*1st subsection*)
	PatternFunctions, hphcIMRPhenomD, hphcIMRPhenomPv2 (*2o subsection*)
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
dL: luminosity distance [Mpc].

Obs.: GR deviation parameters from the ```TIGER``` framework can be passed in the form of Options: 
```hphcIMRPhenomD[f, M, \[Eta], s1z, s2z, \[Iota], tc, \[Phi]ref, dL, \"\[Delta]\[CurlyPhi]-2\"->1, \"\[Delta]\[Beta]3\"->1]```.
By default all deviation parameters are 0.";


hphcIMRPhenomPv2::usage="```hphcIMRPhenomPv2[f, M, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], tc, \[Phi]ref, dL]``` returns 
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
dL: luminosity distance [Mpc].

Obs.: GR deviation parameters from the ```TIGER``` framework can be passed in the form of Options: 
```hphcIMRPhenomPv2[f, M, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], tc, \[Phi]ref, dL, \"\[Delta]\[CurlyPhi]-2\"->1, \"\[Delta]\[Beta]3\"->1]```.
By default all deviation parameters are 0.";


(* ::Section:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*SymRules, NRules and iFunctions for Fisher calculation*)


SymRules = <||>;
NRules = <||>;


Module[
	{symPv2, symD, symFpFc, nPv2, nD, nFpFc, m},
	m = FelipeBarbosa`SymDALI`DerivativeTools`DerivativeRulesLoad/@{"IMRPhenomPv2", "IMRPhenomD", "Detectors"};
	
	{SymRules["IMRPhenomPv2"], SymRules["IMRPhenomD"], SymRules["FpFc"]} = m[[All,1]];
	{NRules["IMRPhenomPv2"], NRules["IMRPhenomD"], NRules["FpFc"]} = m[[All,2]];
];


Protect[SymRules, NRules];


ihphcIMRPhenomD[
	M_, \[Eta]_, \[Chi]1_, \[Chi]2_,
	\[Iota]_, 
	\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_,
	fref_, f_ 
] = Block[
	{G = 4.9254664969309`3.6105383994801805*^-6 (*G/c^3 [s/solarMass]*), \[Omega], \[Omega]ref},
	\[Omega] = f M G; 
	\[Omega]ref = fref M G;
	{(1+Cos[\[Iota]]^2)/2, -I Cos[\[Iota]]}*\[ScriptA]IMR[f,M,\[Eta],\[Chi]1,\[Chi]2] Exp[
		-I \[CapitalPhi]IMR[\[Omega],\[Omega]ref,\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4]
	]
];


ihphcIMRPhenomPv2[
	m1_, m2_, 
	s1x_, s1y_, s1z_,
	s2x_, s2y_, s2z_,
	\[Phi]ref_, \[Iota]_,
	\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_,
	fref_, f_
] = \[ScriptA]IMR[f, fref, m1, m2, s1x,s1y,s1z,s2x,s2y,s2z,\[Phi]ref,\[Iota]]*Exp[
	-I \[CapitalPhi]IMR[f,fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4]
];


iFpFc[
	\[Theta]_, \[Phi]_, \[Psi]_,
	pi_, Dij_,
	f_
] = FpFc[f,\[Theta],\[Phi],\[Psi],pi,Dij];


Protect[ihphcIMRPhenomD, ihphcIMRPhenomPv2, iFpFc, \[CapitalPhi]IMR, \[ScriptA]IMR, FpFc];


(* ::Subsection::Closed:: *)
(*public FpFc and hphc*)


PatternFunctions[f_, \[Delta]_, \[Alpha]_, \[Psi]_, GMST_, pi_, Dij_]/;(
	VectorQ[f, NumberQ] && MatrixQ[Dij, NumberQ] && VectorQ[pi, NumberQ] && 
	-\[Pi]/2 <= \[Delta] <= \[Pi]/2 && 0<=\[Alpha]<=2 \[Pi] && 0<= \[Psi] <= \[Pi] && 0 <= GMST <= 2 \[Pi]
):= Module[
	{LI = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]], \[Theta] = \[Pi]/2-\[Delta], \[Phi] = \[Alpha]-GMST, dij},
	
	dij = Extract[Dij, LI];
	
	NRules["FpFc"]["FpFc"][[1,2,0]][f, \[Theta], \[Phi], \[Psi], pi, dij]	
]


Protect[PatternFunctions];


Options[hphcIMRPhenomD] = {
	"\[Delta]\[CurlyPhi]-2"->0, "\[Delta]\[CurlyPhi]0"->0, "\[Delta]\[CurlyPhi]1"->0,"\[Delta]\[CurlyPhi]2"->0,"\[Delta]\[CurlyPhi]3"->0,"\[Delta]\[CurlyPhi]4"->0,"\[Delta]\[CurlyPhi]5l"->0,"\[Delta]\[CurlyPhi]6"->0,"\[Delta]\[CurlyPhi]6l"->0,"\[Delta]\[CurlyPhi]7"->0,
	"\[Delta]\[Beta]2"->0,"\[Delta]\[Beta]3"->0,
	"\[Delta]\[Alpha]2"->0,"\[Delta]\[Alpha]3"->0,"\[Delta]\[Alpha]4"->0
};


hphcIMRPhenomD[f_, M_, \[Eta]_, s1z_, s2z_, \[Iota]_, tc_, \[Phi]ref_, dL_, OptionsPattern[]]/;(
	VectorQ[f, NumberQ] && Min[f]>0 && M>0 && 0 < \[Eta] <= 0.25 && -1<=s1z<=1 && -1<=s2z<=1 &&
	0<=\[Iota]<=\[Pi] && 0<=\[Phi]ref<= 2 \[Pi] && dL>0
) := Module[
	{
		\[Omega], \[Omega]ref, G = 4.9254664969309`3.6105383994801805*^-6(*G/c^3 [s/solarMass]*), 
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4, auxh, idL = dL 10^-3 (*Mpc->Gpc*),
		Amp, \[Phi],h
	},
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = OptionValue[hphcIMRPhenomD, Options[hphcIMRPhenomD][[All,1]]];
	
	\[Omega] = f G M;
	\[Omega]ref = Min[f] G M;
	
	auxh = Exp[-I (2 \[Pi] f tc - 2 \[Phi]ref)]/idL;
	
	Amp = NRules["IMRPhenomD"]["\[ScriptA]IMR"][[1,2,0]][f, M, \[Eta], s1z, s2z];
	
	\[Phi] = NRules["IMRPhenomD"]["\[CapitalPhi]IMR"][[1,2,0]][
		\[Omega], \[Omega]ref, \[Eta], s1z, s2z,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
	];
	
	h = Amp*Exp[-I \[Phi]]*auxh;
	
	{
		(1+Cos[\[Iota]]^2)/2 h, 
		-I Cos[\[Iota]] h
	}
	
]


Protect[hphcIMRPhenomD];


Options[hphcIMRPhenomPv2] = {
	"\[Delta]\[CurlyPhi]-2"->0, "\[Delta]\[CurlyPhi]0"->0, "\[Delta]\[CurlyPhi]1"->0,"\[Delta]\[CurlyPhi]2"->0,"\[Delta]\[CurlyPhi]3"->0,"\[Delta]\[CurlyPhi]4"->0,"\[Delta]\[CurlyPhi]5l"->0,"\[Delta]\[CurlyPhi]6"->0,"\[Delta]\[CurlyPhi]6l"->0,"\[Delta]\[CurlyPhi]7"->0,
	"\[Delta]\[Beta]2"->0,"\[Delta]\[Beta]3"->0,
	"\[Delta]\[Alpha]2"->0,"\[Delta]\[Alpha]3"->0,"\[Delta]\[Alpha]4"->0
};


hphcIMRPhenomPv2[f_, M_, \[Eta]_, s1x_, s1y_, s1z_, s2x_, s2y_, s2z_, \[Iota]_, tc_, \[Phi]ref_, dL_, OptionsPattern[]]/;(
	VectorQ[f, NumberQ] && Min[f]>0 && M>0 && 0 < \[Eta] <= 0.25 && -1<=s1z<=1 && -1<=s2z<=1 &&
	0<=\[Iota]<=\[Pi] && 0<=\[Phi]ref<= 2 \[Pi] && dL>0
) := Module[
	{
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4, 
		auxh, idL = dL 10^-3 (*Mpc->Gpc*), Amp, \[Phi],hphc, m1 ,m2, fref
	},
	m1 = (M/2) (1 + Sqrt[1 - 4 \[Eta]]);
	m2 = (M/2) (1 - Sqrt[1 - 4 \[Eta]]);
	fref = Min[f];
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = OptionValue[hphcIMRPhenomD, Options[hphcIMRPhenomD][[All,1]]];
	
	auxh = Exp[-I 2 \[Pi] f tc]/idL;
	
	
	Amp = NRules["IMRPhenomPv2"]["\[ScriptA]IMR"][[1,2,0]][f, fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]ref, \[Iota]];
	
	\[Phi] = NRules["IMRPhenomPv2"]["\[CapitalPhi]IMR"][[1,2,0]][
		f, fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
	];
	
	hphc = Transpose[Amp]*Exp[-I \[Phi]]*auxh;
	
	Transpose[hphc]
]


Protect[hphcIMRPhenomPv2]


(* ::Subsection::Closed:: *)
(*MakeSNR*)


iApAcIMRPhenomD[f_, M_, \[Eta]_, s1z_, s2z_, \[Iota]_, dL_]/;(
	VectorQ[f, NumberQ] && Min[f]>0 && M>0 && 0 < \[Eta] <= 0.25 && -1<=s1z<=1 && -1<=s2z<=1 &&
	0<=\[Iota]<=\[Pi] && 0<=\[Phi]ref<= 2 \[Pi] && dL>0
) := Module[
	{
		G = 4.9254664969309`3.6105383994801805*^-6(*G/c^3 [s/solarMass]*), 
		Amp, \[Phi],h, idL = 10^-3 dL
	},
	
	
	Amp = NRules["IMRPhenomD"]["\[ScriptA]IMR"][[1,2,0]][f, M, \[Eta], s1z, s2z];
	
	h = Amp/idL;
	
	{
		(1+Cos[\[Iota]]^2)/2 h, 
		-I Cos[\[Iota]] h
	}
]


iApAcIMRPhenomPv2[f_, M_, \[Eta]_, s1x_, s1y_, s1z_, s2x_, s2y_, s2z_, \[Iota]_, \[Phi]ref_, dL_]/;(
	VectorQ[f, NumberQ] && Min[f]>0 && M>0 && 0 < \[Eta] <= 0.25 && -1<=s1z<=1 && -1<=s2z<=1 &&
	0<=\[Iota]<=\[Pi] && 0<=\[Phi]ref<= 2 \[Pi] && dL>0
) := Module[
	{
		idL = dL 10^-3 (*Mpc->Gpc*), Amp, m1 ,m2, fref
	},
	m1 = (M/2) (1 + Sqrt[1 - 4 \[Eta]]);
	m2 = (M/2) (1 - Sqrt[1 - 4 \[Eta]]);
	fref = Min[f];
	
	Amp = NRules["IMRPhenomPv2"]["\[ScriptA]IMR"][[1,2,0]][f, fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]ref, \[Iota]]/idL
]


(*Options[MakeSNR] = {"fmax"->1024, "fmin"-> 10, "\[CapitalDelta]f"-> 1};*)


WFVars["IMRPhenomD"] = {f, M, \[Eta], s1z, s2z, \[Iota], dL};
WFVars["IMRPhenomPv2"] ={f,M,\[Eta],s1x,s1y,s1z,s2x,s2y,s2z,\[Iota],\[Phi]ref,dL};


(*MakeSNR[WF_String, dets_Association, OptionsPattern[]]/;(
	WF === "IMRPhenomD" || WF === "IMRPhenomPv2"
) := Module[
	{AmpHead, vars = WFVars[WF]},
	
	AmpHead = ToExpression["iApAc"<>WF];
	


]*)


(*Block[
	{det1, det2, det3},
	amp@@vars;
]*)





(* ::Section:: *)
(*PackageFooter*)


EndPackage[];
