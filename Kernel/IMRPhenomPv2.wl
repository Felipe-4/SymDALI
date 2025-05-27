(* ::Package:: *)

(* ::Section:: *)
(*PackageHeader*)


BeginPackage["FelipeBarbosa`SymDALI`IMRPhenomPv2`"]

Begin["`Private`"];


(* ::Section:: *)
(*Defintions*)


Clear@EradRational0815s
EradRational0815s[eta_, s_] = ((0.055974469826360077` eta+0.5809510763115132` eta^2-0.9606726679372312` eta^3+3.352411249771192` eta^4) (1.` +(-0.0030302335878845507` -2.0066110851351073` eta+7.7050567802399215` eta^2) s))/(1.` +(-0.6714403054720589` -1.4756929437702908` eta+7.304676214885011` eta^2) s)


FinalSpin0815s[\[Eta]_, S_] =\[Eta] (3.4641016151377544` -4.399247300629289` \[Eta]+9.397292189321194` \[Eta]^2-13.180949901606242` \[Eta]^3+S (-0.0850917821418767`+S (0.1014665242971878` -2.0967746996832157` \[Eta])+1.`/\[Eta]-5.837029316602263` \[Eta]+S^3 (-0.8676969352555539` +2.064046835273906` \[Eta])+S^2 (-1.3546806617824356` +4.108962025369336` \[Eta])))


FinalSpinInPlane[m1_, m2_, s1z_, s2z_,  \[Chi]p_] := Module[
	
	{M = m1+m2, \[Eta] = (m1 m2)/(m1+m2)^2, q,afParallel, Sperp},
	
	q = m1/M;
	
	afParallel = FinalSpin0815s[\[Eta], (m1/M)^2 s1z + (m2/M)^2 s2z];
	
	Sperp = \[Chi]p q^2;
	
	Sign[afParallel] Sqrt[Sperp^2+afParallel^2]
]


(*re\[Omega]->rd, im\[Omega]->damp*)
re\[Omega][\[Chi]_] :=(0.05947169566573468` -0.14989771215394762` \[Chi]+0.09535606290986028` \[Chi]^2+0.02260924869042963` \[Chi]^3-0.02501704155363241` \[Chi]^4-0.005852438240997211` \[Chi]^5+0.0027489038393367993` \[Chi]^6+0.0005821983163192694` \[Chi]^7)/(1-2.8570126619966296` \[Chi]+2.373335413978394` \[Chi]^2-0.6036964688511505` \[Chi]^4+0.0873798215084077` \[Chi]^6);
im\[Omega][\[Chi]_] :=(0.014158792290965177` -0.036989395871554566` \[Chi]+0.026822526296575368` \[Chi]^2+0.0008490933750566702` \[Chi]^3-0.004843996907020524` \[Chi]^4-0.00014745235759327472` \[Chi]^5+0.0001504546201236794` \[Chi]^6)/(1-2.5900842798681376` \[Chi]+1.8952576220623967` \[Chi]^2-0.31416610693042507` \[Chi]^4+0.009002719412204133` \[Chi]^6);


(*
Asp1 and Asp2 in convert_spins in IMRPhenomPv2.py
*)
Asp1[m1_, m2_, s1x_, s1y_] := m1^2 Sqrt[s1x^2 + s1y^2] (2+ (3 m2)/(2 m1));
Asp2[m1_, m2_, s2x_, s2y_] := m2^2 Sqrt[s2x^2 + s2y^2] (2+ (3 m1)/(2 m2));


(*
\[Chi]p in convert_spins in IMRPhenomPv2.py
*)
\[Chi]p[Asp1_, Asp2_, A1m1sq_] := If[Asp1 >= Asp2, Asp1/A1m1sq, Asp2/A1m1sq]


\[Omega]rd[m1_, m2_, s1x_, s1y_, s1z_, s2x_, s2y_, s2z_] := Module[
	{chip,aeff, Erad,s , \[Eta]},
	
	chip = \[Chi]p[
		Asp1[m1, m2, s1x, s1y],
		Asp2[m1,m2, s2x, s2y],
		1/2 m1 (4 m1+3 m2) (*m1^2 (2 + (3 * m2) / (2 * m1))//Simplify*)
	];
	
	aeff = FinalSpinInPlane[m1, m2, s1z, s2z, chip];
	
	
	\[Eta] = m1 m2/(m1+m2)^2;
	s =(m1^2 s1z + m2^2 s2z)/(m1^2+m2^2);
	
	Erad = EradRational0815s[\[Eta],s];
	
	re\[Omega][aeff]/(1-Erad)
]

\[Omega]damp[m1_, m2_, s1x_, s1y_, s1z_, s2x_, s2y_, s2z_] := Module[
	{aeff, Erad,s , \[Eta], chip},
	
	chip = \[Chi]p[
		Asp1[m1, m2, s1x, s1y],
		Asp2[m1,m2, s2x, s2y],
		1/2 m1 (4 m1+3 m2) (*m1^2 (2 + (3 * m2) / (2 * m1))//Simplify*)
	];
	
	aeff = FinalSpinInPlane[m1, m2, s1z, s2z, chip];
	
	\[Eta] = m1 m2/(m1+m2)^2;
	
	s =(m1^2 s1z + m2^2 s2z)/(m1^2+m2^2);
	
	Erad = EradRational0815s[\[Eta],s];
	
	im\[Omega][aeff]
]


\[Omega]Peak[\[Gamma]2_, \[Gamma]3_, fRD_, fDamp_] := If[
       \[Gamma]2<=1,
       fRD + (fDamp \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2,
       fRD - fDamp \[Gamma]3/\[Gamma]2
   
];


(* ::Section:: *)
(*Package Footer*)


End[];
EndPackage[];
