(* ::Package:: *)

(* ::Section::Closed:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`Population`"];


Unprotect[
	CholeskyInverse, 
	ComovingDistance, DifferentialComovingVolume, 
	MadauDickinsonProfile, PowerLawPlusPeak, SpinDistribution
];


ClearAll[
	CholeskyInverse,
	ComovingDistance, DifferentialComovingVolume, 
	MadauDickinsonProfile, PowerLawPlusPeak, SpinDistribution
];


CholeskyInverse::usage="```CholeskyInverse[m]``` returns ```{im, \[Epsilon]}```, where ```im``` is the inverse of m 
calculated through ```CholeskyDecomposition``` and ```\[Epsilon]``` is the largest non zero number in 
```{im.m - IdentityMatrix, m.im - IdentityMatrix}```. 
```m``` has to be positive semidefinite, symmetric and purelly numeric.";


ComovingDistance::usage="```ComovingDistance[]``` returns \!\(\*SubscriptBox[\(D\), \(c\)]\)[z], where

\!\(\*SubscriptBox[\(D\), \(c\)]\)[z] \[Congruent] \!\(\*FractionBox[\(c\), \(H0\)]\) \!\(\*SuperscriptBox[SubscriptBox[\(\[Integral]\), \(0\)], \(z\)]\) \!\(\*FractionBox[\(\[DifferentialD]zp\), SqrtBox[\(\[CapitalOmega]m\\\ \*SuperscriptBox[\((1 + zp)\), \(3\)]\\\  + \\\ \(\[CapitalOmega]\[Kappa]\\\ \((1 + \*SuperscriptBox[\(zp\), \(2\)])\)\)\\\  + \\\ \((1\\\  - \\\ \[CapitalOmega]m\\\  - \\\ \[CapitalOmega]\[Kappa])\)\)]]\) [Mpc],

in the form of an ```InterpolatingFunction```. The default values are (\!\(\*TemplateBox[{\"\\\"Planck18\\\"\", \"https://arxiv.org/pdf/1807.06209\"},\n\"HyperlinkURL\"]\) Table II):
				{\[CapitalOmega]m, \[CapitalOmega]\[Kappa]} = {0.3111, 0},
				H0 = 67.66 \!\(\*FractionBox[\(km/s\), \(Mpc\)]\).
These values can be altered with the options {\"H0\", \"\[CapitalOmega]m\", \"\[CapitalOmega]\[Kappa]\"}.
The maximum value for redshift (the default is zmax=20) can be altered with the option \"zmax\". ";


DifferentialComovingVolume::usage="```DifferentialComovingVolume[]``` returns \!\(\*FractionBox[\(\[DifferentialD]Vc[z]\), \(\[DifferentialD]z\)]\), where

\!\(\*FractionBox[\(\[DifferentialD]Vc[z]\), \(\[DifferentialD]z\)]\) \[Congruent] 4 \[Pi] (\!\(\*SubscriptBox[\(D\), \(c\)]\)[z]\!\(\*SuperscriptBox[\()\), \(2\)]\) \!\(\*FractionBox[\(\[DifferentialD]\(\*SubscriptBox[\(D\), \(c\)]\)[z]\), \(\[DifferentialD]z\)]\) [\!\(\*SuperscriptBox[\(Mpc\), \(3\)]\)],

in the form of an ```InterpolatingFunction```. The default values are (\!\(\*TemplateBox[{\"\\\"Planck18\\\"\", \"https://arxiv.org/pdf/1807.06209\"},\n\"HyperlinkURL\"]\) Table II):
				{\[CapitalOmega]m, \[CapitalOmega]\[Kappa]} = {0.3111, 0},
				H0 = 67.66 \!\(\*FractionBox[\(km/s\), \(Mpc\)]\).
These values can be altered with the options {\"H0\", \"\[CapitalOmega]m\", \"\[CapitalOmega]\[Kappa]\"}.
The maximum value for redshift (the default is zmax=20) can be altered with the option \"zmax\". 

See ```ComovingDistance``` for the definition of \!\(\*SubscriptBox[\(D\), \(c\)]\)[z].";


PowerLawPlusPeak::usage="```PowerLawPlusPeak[m1, q]``` returns the normalized POWER LAW+PEAK distribution (see \!\(\*TemplateBox[{\"\\\"2207.02771\\\"\", \"https://arxiv.org/pdf/2207.02771\"},\n\"HyperlinkURL\"]\) Appendix A.2.) 
at (m1, q). The distribution hyper parameters are set to the values found in GWTC-3 and can be altered with the options: 

{\"\[Lambda]peak\", \"\[Alpha]\", \"m_min\",\"m_max\",\"\[Delta]m\",\"\[Mu]m\", \"\[Sigma]m\", \"\[Beta]q\" , \"qmin\"}.

\"qmin\" sets the minimum value of the mass ratio, the default is 0.1.";


MadauDickinsonProfile::usage="```MadauDickinsonProfile[z]``` returns the Madau-Dickinson profile

	\[Psi](z | \[Alpha]z, \[Beta]z, zp) \[Proportional] \!\(\*FractionBox[SuperscriptBox[\((1 + z)\), \(\[Alpha]z\)], \(1\(\\\ \)\(+\)\(\\\ \)\*SuperscriptBox[\((\*FractionBox[\(1 + z\), \(1 + zp\)])\), \(\[Alpha]z + \[Beta]z\)]\(\\\ \)\)]\) [\!\(\*SuperscriptBox[\(Gpc\), \(-3\)]\) \!\(\*SuperscriptBox[\(yr\), \(-1\)]\)],
	\[Psi](z | \[Alpha]z, \[Beta]z, zp) -> \!\(\*SubscriptBox[\(R\), \(0\)]\) (1+z\!\(\*SuperscriptBox[\()\), \(\[Alpha]z\)]\),  z<<1,

at redshift ```z```. The parameters (\[Alpha]z, \[Beta]z, zp) can be set with the options {\"\[Alpha]z\", \"\[Beta]z\", \"zp\"}, and \!\(\*SubscriptBox[\(R\), \(0\)]\) can be set 
with the option \"R0\" [\!\(\*SuperscriptBox[\(Gpc\), \(-3\)]\) \!\(\*SuperscriptBox[\(yr\), \(-1\)]\)]. Default values, 

	{\"\[Alpha]z\",\"\[Beta]z\",\"zp\"} = {2.7, 3, 2}   \[And]   \"R0\" = 17 [\!\(\*SuperscriptBox[\(Gpc\), \(-3\)]\) \!\(\*SuperscriptBox[\(yr\), \(-1\)]\)],

insure agreement with GWTC-3  (see \!\(\*TemplateBox[{\"\\\"2207.02771\\\"\", \"https://arxiv.org/pdf/2207.02771\"},\n\"HyperlinkURL\"]\) Appendix A.1).";


SpinDistribution::usage="```SpinDistribution[\[Chi]1, \[Chi]2, z1, z2]``` returns the spin distribution \[ScriptP](\[Chi]1, \[Chi]2, z1, z2), where

	\[ScriptP](\[Chi]1, \[Chi]2, z1, z2) = p(\[Chi]1 | \[Alpha]\[Chi], \[Beta]\[Chi]) p(\[Chi]2 | \[Alpha]\[Chi], \[Beta]\[Chi]) p(\!\(\*OverscriptBox[\(z\), \(->\)]\) | \[Zeta], \[Sigma]t),
	p(\[Chi]i | \[Alpha]\[Chi], \[Beta]\[Chi]) = Beta(\[Alpha]\[Chi], \[Beta]\[Chi]),
	p(\!\(\*OverscriptBox[\(z\), \(->\)]\) | \[Zeta], \[Sigma]t) = \[Zeta] \[ScriptCapitalN](\!\(\*OverscriptBox[\(z\), \(->\)]\) | 0, \[Sigma]t) + (1- \[Zeta]) I(\!\(\*OverscriptBox[\(z\), \(->\)]\)),

with zi \[Congruent] Cos[\[Theta]i] where \[Theta] is the angle between spin and orbital angular momenta.

The parameters (\[Alpha]\[Chi], \[Beta]\[Chi]) of the Beta distribution and (\[Sigma]t, \[Zeta]) for the orientation distribution can be 
altered with the options {\"\[Alpha]\[Chi]\",  \"\[Beta]\[Chi]\", \"\[Sigma]t\", \"\[Zeta]\"}. Default values are set to agree with GWTC-3 
(see \!\(\*TemplateBox[{\"\\\"2207.02771\\\"\", \"https://arxiv.org/pdf/2207.02771\"},\n\"HyperlinkURL\"]\) Appendix A.3).";


Begin["Private`"];


(* ::Section::Closed:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*Matrix Inverse*)


CholeskyInverse[m_]/; MatrixQ[m, NumberQ] && SymmetricMatrixQ[m] := Module[
	{u, iu, inverse, d = Length[m], \[Epsilon]},
	
	u =  Quiet[Check[CholeskyDecomposition[m], m], CholeskyDecomposition::posdef];
	
	If[
		u==m,
		inverse = Quiet[m//Inverse, Inverse::luc],
		iu = Quiet[u//Inverse, Inverse::luc]; 
		inverse = iu . iu\[ConjugateTranspose]
	];
	
	\[Epsilon] = Join[
		inverse . m - IdentityMatrix[d],
		m . inverse-IdentityMatrix[d]
	]//Abs//Max;
	
	{inverse, \[Epsilon]}
]

CholeskyInverse[x___] := $Failed


Protect[CholeskyInverse];


(* ::Subsection::Closed:: *)
(*Cosmology standard functions*)


\[Mu][zp_, \[CapitalOmega]m_, \[CapitalOmega]\[Kappa]_] = (\[CapitalOmega]m*(1 + zp)^3 + \[CapitalOmega]\[Kappa] (1+zp)^2 + (1 -\[CapitalOmega]m - \[CapitalOmega]\[Kappa]))^(-1/2);


Options[ComovingDistance] = {
	"H0" -> 67.66,
	"\[CapitalOmega]m" -> 0.3111,
	"\[CapitalOmega]\[Kappa]" -> 0,
	"zmax" -> 20
};


ComovingDistance[OptionsPattern[]] := Module[
	{int, c = 299792.458`(*km/s*), zmax, H0, \[CapitalOmega]m, \[CapitalOmega]\[Kappa]},
	
	{zmax, H0, \[CapitalOmega]m, \[CapitalOmega]\[Kappa]} = OptionValue[ComovingDistance, #]&/@{"zmax", "H0", "\[CapitalOmega]m", "\[CapitalOmega]\[Kappa]"};
	
	NDSolve[
		{y'[zp] == c/H0 \[Mu][zp, \[CapitalOmega]m, \[CapitalOmega]\[Kappa]], y[0] == 0 },
		y,
		{zp, 0, zmax}
	][[-1,-1,-1]]
	
	(*ComovingDistance[z_] = int[z] c/67.74*)
]


Options[DifferentialComovingVolume] = {
	"H0" -> 67.66,
	"\[CapitalOmega]m" -> 0.3111,
	"\[CapitalOmega]\[Kappa]" -> 0,
	"zmax" -> 20
};


DifferentialComovingVolume[OptionsPattern[]] := Module[
	{dDcdz, zmax, H0, \[CapitalOmega]m, \[CapitalOmega]\[Kappa], c = 299792.458` (*km/s*), Dc, d\[Mu]dz},
	
	{zmax, H0, \[CapitalOmega]m, \[CapitalOmega]\[Kappa]} = OptionValue[DifferentialComovingVolume, #]&/@{"zmax", "H0", "\[CapitalOmega]m", "\[CapitalOmega]\[Kappa]"};
	
	Dc[z_?NumberQ] := c/H0 NIntegrate[\[Mu][zp, \[CapitalOmega]m, \[CapitalOmega]\[Kappa]], {zp, 0, z}];
	
	d\[Mu]dz[z_] = D[\[Mu][z, \[CapitalOmega]m, \[CapitalOmega]\[Kappa]], z];
	
	NDSolve[
		{
			dVcdz'[zp] == 4 \[Pi] (2 Dc[zp] (c/H0 \[Mu][zp, \[CapitalOmega]m, \[CapitalOmega]\[Kappa]])^2 + Dc[zp]^2 c/H0  d\[Mu]dz[zp]),
			dVcdz[0] == 0
		},
		dVcdz,
		{zp, 0, zmax}
	][[-1,-1,-1]]	
	
]


Protect[ComovingDistance, DifferentialComovingVolume];


(* ::Subsection::Closed:: *)
(*Probability distributions*)


Options[MadauDickinsonProfile] = {
	"\[Alpha]z" -> 2.7,
	"\[Beta]z" -> 3,
	"R0" -> 17,
	"zp" -> 2
};


MadauDickinsonProfile[z_, OptionsPattern[]] := Module[
		{expr, \[Alpha]z, \[Beta]z, iR0, zp, R0},
		
		{\[Alpha]z, \[Beta]z, R0, zp} = OptionValue[MadauDickinsonProfile, #]&/@{"\[Alpha]z","\[Beta]z","R0","zp"};
		
		
		expr[red_] := (1+red)^\[Alpha]z/(1+((1+red)/(1+zp))^(\[Alpha]z+\[Beta]z)); (*Eq.(A1) 2207.02771v3*)
		
		
		iR0 = R0/expr[0]; (*pag 12 2207.02771v3*)
		
		expr[z] iR0
];


Block[
	{norm},
	(*
	As confusing as it can be: (A4) 2207.02771v3 states power law with index -\[Alpha]. Since power laws
	are tipically x^-\[Alpha] this would imply x^\[Alpha], but that does not seem in accord with Fig 10 and Fig. 11
	in 2111.03634v5, so I will use x^-\[Alpha] anyway. 
	
	Hard to understand what is the difficult in explicitly write 1 line equations and formulas...
	*)
	norm = Integrate[m1^(-\[Alpha]), {m1, mmin, mmax}, Assumptions->\[Alpha]>0 && mmin>0 && mmax>mmin]; 
	
	Clear@\[ScriptCapitalP];
	
	\[ScriptCapitalP][m1_, \[Alpha]_, mmin_, mmax_]  = m1^(-\[Alpha])/norm;
	
]

\[ScriptCapitalN][m1_, \[Mu]m_, \[Sigma]m_] = PDF[NormalDistribution[\[Mu]m, \[Sigma]m], m1]; (*Eq.(A4) 2207.02771v3*)

f[m_, \[Delta]m_] = Exp[\[Delta]m/m+\[Delta]m/(m-\[Delta]m)];

S[m1_, mmin_, \[Delta]m_]/;mmin<=m1<mmin+\[Delta]m := 1/(f[m1-mmin, \[Delta]m]+1);
S[m1_, mmin_, \[Delta]m_]/;mmin+\[Delta]m<=m1 := 1;
S[m1_, mmin_, \[Delta]m_]/;m1 < mmin := 0;


Clear@p
p[m1_, \[Lambda]_, \[Alpha]_, mmin_, \[Delta]m_, mmax_, \[Mu]m_, \[Sigma]m_] := (
		(1-\[Lambda]) \[ScriptCapitalP][m1, \[Alpha],mmin, mmax] + \[Lambda] \[ScriptCapitalN][m1, \[Mu]m, \[Sigma]m]
	)*S[m1, mmin, \[Delta]m];


Clear@\[ScriptP]
\[ScriptP][q_, \[Beta]_, m1_, mmin_, \[Delta]m_] := (q^\[Beta]) S[q m1, mmin, \[Delta]m] (*Eq. A7 2207.02771v3*)


ClearAll@norm


norm[\[Lambda]_, \[Alpha]_, mmin_, \[Delta]m_, mmax_, \[Mu]m_, \[Sigma]m_, \[Beta]q_, qmin_] := ( 
	norm[\[Lambda], \[Alpha], mmin, \[Delta]m, mmax, \[Mu]m, \[Sigma]m, \[Beta]q, qmin] = Block[
		{prob},
		
		prob[m1_?NumberQ, q_] :=  p[m1, \[Lambda], \[Alpha], mmin, \[Delta]m, mmax, \[Mu]m, \[Sigma]m] \[ScriptP][q, \[Beta]q, m1, mmin, \[Delta]m];
		
		NIntegrate[
			prob[m1, q], 
			{m1, mmin, mmax}, {q, qmin, 1}
		]
	]
	
) 


Options[PowerLawPlusPeak] = {
	"\[Lambda]peak" -> 0.039,
	"\[Alpha]" -> 3.4,
	"m_min"-> 5.1,
	"m_max" -> 87,
	"\[Delta]m"->4.8,
	"\[Mu]m" -> 34,
	"\[Sigma]m"-> 3.6,
	"\[Beta]q" -> 1.1,
	"qmin" -> 0.1
};


PowerLawPlusPeak[m1_, q_, OptionsPattern[]] := Module[
	{inorm, \[Lambda],\[Alpha],mmin,\[Delta]m,mmax,\[Mu]m,\[Sigma]m,\[Beta]q,qmin},
	
	{\[Lambda],\[Alpha],mmin,\[Delta]m,mmax,\[Mu]m,\[Sigma]m,\[Beta]q,qmin} = OptionValue[PowerLawPlusPeak, #]&/@{"\[Lambda]peak","\[Alpha]","m_min","\[Delta]m","m_max","\[Mu]m","\[Sigma]m","\[Beta]q","qmin"};
	
	inorm = norm[\[Lambda],\[Alpha],mmin,\[Delta]m,mmax,\[Mu]m,\[Sigma]m,\[Beta]q,qmin];
	
	p[m1, \[Lambda], \[Alpha], mmin, \[Delta]m, mmax, \[Mu]m, \[Sigma]m]*\[ScriptP][q, \[Beta]q, m1, mmin, \[Delta]m]/inorm
]	


Block[
	{beta1, beta2, direction},(*params from Eq. (A9) 2207.02771v3*)
	
	beta1 = PDF[BetaDistribution[\[Alpha]x, \[Beta]x], \[Chi]1][[1,1,1]];  (*Eq. (A8) 2207.02771v3*) 
	beta2 = PDF[BetaDistribution[\[Alpha]x, \[Beta]x], \[Chi]2][[1,1,1]];
	
	direction = (1- \[Zeta])/4 + \[Zeta] 2/(\[Pi] \[Sigma]t^2) Exp[-(((z1-1)^2 + (z2-1)^2)/(2 \[Sigma]t^2))]/Erf[Sqrt[2]/\[Sigma]t]^2; (*Eq. (4) 1704.08370v2*)
	Clear@iSpinDistribution;
	
	
	(*\[Phi]1 \[And]\[Phi]2 are random. {s1x, s1y} = \[Chi]1 z1 {Cos[\[Phi]1], Sin[\[Phi]1]}, ...*)
	iSpinDistribution[\[Chi]1_, \[Chi]2_, z1_, z2_, \[Alpha]x_, \[Beta]x_, \[Sigma]t_, \[Zeta]_] = (beta1 beta2 direction);
]


Options[SpinDistribution] = {
	"\[Alpha]\[Chi]" -> 1.6,
	"\[Beta]\[Chi]"->4.12,
	"\[Sigma]t"-> 1.5,
	"\[Zeta]"->0.66
};


SpinDistribution[\[Chi]1_, \[Chi]2_, z1_, z2_, OptionsPattern[]] := Module[
	{\[Alpha]x, \[Beta]x, \[Sigma]t, \[Zeta]},
	
	{\[Alpha]x, \[Beta]x, \[Sigma]t, \[Zeta]} = OptionValue[SpinDistribution, #]&/@{"\[Alpha]\[Chi]", "\[Beta]\[Chi]", "\[Sigma]t", "\[Zeta]"};
	
	iSpinDistribution[\[Chi]1, \[Chi]2, z1, z2, \[Alpha]x, \[Beta]x, \[Sigma]t, \[Zeta]]
]


Protect[MadauDickinsonProfile, PowerLawPlusPeak, SpinDistribution];


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
