#!/usr/bin/env wolframscript
(* ::Package:: *)

(*Basic stuff:*)
$HistoryLength = 1;
PacletDirectoryLoad["/home/davi/felipe_barbosa/"]
<<FelipeBarbosa`SymDALI`

vectorDefs = Module[
	{DALIDir, docDir},
	DALIDir = "/home/davi/felipe_barbosa/SymDALI";
	docDir = FileNameJoin[{DALIDir, "Documentation/English/Tutorials"}];
	Import@FileNameJoin[{docDir, "IMRPhenomDComponentsDefinition.wdx"}]
];

SymRules = <||>;
NRules = <||>;


(* ::Chapter:: *)
(*Amplitude*)


(* ::Section:: *)
(*PPN*)


(* Coefficients *)
A0 = 1;
A1 = 0;
A2 = -323/224 + (451 \[Eta])/168;

A3 = (27 \[Delta] \[Chi]a)/8 + (27/8 - (11 \[Eta])/6) \[Chi]s//.{\[Delta]->Sqrt[1-4 \[Eta]], \[Chi]s -> (\[Chi]1+\[Chi]2)/2, \[Chi]a -> (\[Chi]1-\[Chi]2)/2};

A4 = (-27312085/8128512 - (1975055 \[Eta])/338688 + (105271 \[Eta]^2)/24192 +
     (-81/32 + 8 \[Eta]) \[Chi]a^2 - (81/16) \[Delta] \[Chi]a \[Chi]s +
     (-81/32 + (17 \[Eta])/8) \[Chi]s^2)//.{\[Delta]->Sqrt[1-4 \[Eta]], \[Chi]s -> (\[Chi]1+\[Chi]2)/2, \[Chi]a -> (\[Chi]1-\[Chi]2)/2};
A5 = (-85 \[Pi]/64 + (85 \[Pi] \[Eta])/16 +
     \[Delta] (285197/16128 - (1579 \[Eta])/4032) \[Chi]a +
     (285197/16128 - (15317 \[Eta])/672 - (2227 \[Eta]^2)/1008) \[Chi]s)//.{\[Delta]->Sqrt[1-4 \[Eta]], \[Chi]s -> (\[Chi]1+\[Chi]2)/2, \[Chi]a -> (\[Chi]1-\[Chi]2)/2};
A6 = (-177520268561/8583708672 + 
     ((545384828789/5007163392) - (205 \[Pi]^2)/48) \[Eta] - 
     (3248849057 \[Eta]^2)/178827264 + 
     (34473079 \[Eta]^3)/6386688 +
     (1614569/64512 - (1873643 \[Eta])/16128 + (2167 \[Eta]^2)/42) \[Chi]a^2 +
     (31 \[Pi]/12 - (7 \[Pi] \[Eta])/3) \[Chi]s +
     (1614569/64512 - (61391 \[Eta])/1344 + (57451 \[Eta]^2)/4032) \[Chi]s^2 +
     \[Delta] \[Chi]a (31 \[Pi]/12 + ((1614569/32256) - (165961 \[Eta])/2688) \[Chi]s))//.{\[Delta]->Sqrt[1-4 \[Eta]], \[Chi]s -> (\[Chi]1+\[Chi]2)/2, \[Chi]a -> (\[Chi]1-\[Chi]2)/2};



(* ::Section:: *)
(*Inspiral*)


insVecAmplitude[\[Eta]_, \[Chi]1_,\[Chi]2_] = Sqrt[\[Eta]] {A0, A1, A2, A3, A4, A5, A6};
\[Omega]insVecAmplitude[\[Omega]_] = (\[Pi] \[Omega])^(#/3)&/@Range[0,6];
(*Include \[Omega]^(-7/6)*)
\[Omega]insVecAmplitude[\[Omega]_] = \[Omega]^(-7/6) \[Omega]insVecAmplitude[\[Omega]];


(* ::Text:: *)
(*Light is something*)


{\[ScriptCapitalA]0[\[Eta]_,  M_], \[ScriptCapitalA]1[\[Eta]_, M_], \[ScriptCapitalA]2[\[Eta]_, M_], \[ScriptCapitalA]3[\[Eta]_, M_, \[Chi]1_, \[Chi]2_], \[ScriptCapitalA]4[\[Eta]_, M_, \[Chi]1_, \[Chi]2_], \[ScriptCapitalA]5[\[Eta]_, M_, \[Chi]1_, \[Chi]2_], \[ScriptCapitalA]6[\[Eta]_, M_, \[Chi]1_, \[Chi]2_]} = Block[
	{G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]List },
	\[Omega]List = \[Omega]insVecAmplitude[\[Omega]]//. \[Omega] -> M G;
	
	insVecAmplitude[\[Eta], \[Chi]1, \[Chi]2]*M^2*\[Omega]List
];


insvec  =Sqrt[\[Eta]] (PhenomCoeff[\[Eta], \[Chi]PN, #]&/@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[1;;3]])//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };
\[Omega]insvec = \[Omega]^(-7/6)*(\[Omega]^((#+6)/3)&/@Range[3]);


{\[Rho]1[\[Eta]_, M_, \[Chi]1_, \[Chi]2_],\[Rho]2[\[Eta]_, M_, \[Chi]1_, \[Chi]2_], \[Rho]3[\[Eta]_, M_, \[Chi]1_, \[Chi]2_]}= Block[
	{G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]List},
	\[Omega]List = \[Omega]insvec//.\[Omega]-> G M;
	
	insvec*M^2*\[Omega]List


];


InsAmpExpr  = {
	U\[ScriptCapitalA]0[\[Eta],M],U\[ScriptCapitalA]1[\[Eta],M],U\[ScriptCapitalA]2[\[Eta],M],
	U\[ScriptCapitalA]3[\[Eta],M,\[Chi]1,\[Chi]2],U\[ScriptCapitalA]4[\[Eta],M,\[Chi]1,\[Chi]2],
	U\[ScriptCapitalA]5[\[Eta],M,\[Chi]1,\[Chi]2],U\[ScriptCapitalA]6[\[Eta],M,\[Chi]1,\[Chi]2]} . (\[Omega]insVecAmplitude[\[Omega]]//. {\[Omega] -> f, \[Pi]->1}) + {U\[Rho]1[\[Eta],M,\[Chi]1,\[Chi]2],U\[Rho]2[\[Eta],M,\[Chi]1,\[Chi]2],U\[Rho]3[\[Eta],M,\[Chi]1,\[Chi]2]} . (\[Omega]insvec//.\[Omega]->f);


(* ::Section::Closed:: *)
(*Ringdown and Damping -Phase*)


aeff[\[Eta]_, \[Chi]1_, \[Chi]2_] = Block[
	{S, \[Delta] = Sqrt[1-4 \[Eta]]},
	S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;
	S + 2 Sqrt[3.] \[Eta] + (-0.085` S+0.102` S^2-1.355` S^3-0.868` S^4) \[Eta]-4.399` \[Eta]^2+(-5.837` S-2.097` S^2+4.109` S^3+2.064` S^4) \[Eta]^2+9.397` \[Eta]^3-13.181` \[Eta]^4
]//Simplify;

Erad[\[Eta]_, \[Chi]1_, \[Chi]2_] = Block[
	{EradNS=0.0559745 \[Eta]+0.580951 \[Eta]^2-0.960673 \[Eta]^3+3.35241 \[Eta]^4, S, \[Delta] = Sqrt[1- 4 \[Eta]]},
	S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;
	S = S/(1-2 \[Eta]);
	
	(EradNS (1+(-0.00303023` -2.00661` \[Eta]+7.70506` \[Eta]^2) S))/(1+ (-0.67144` -1.47569` \[Eta] + 7.30468` \[Eta]^2) S)
]//Simplify;

(*re\[Omega] -> interpolation for ringdown and the other is for damping:*)
re\[Omega][\[Chi]_] :=(0.05947169566573468` -0.14989771215394762` \[Chi]+0.09535606290986028` \[Chi]^2+0.02260924869042963` \[Chi]^3-0.02501704155363241` \[Chi]^4-0.005852438240997211` \[Chi]^5+0.0027489038393367993` \[Chi]^6+0.0005821983163192694` \[Chi]^7)/(1-2.8570126619966296` \[Chi]+2.373335413978394` \[Chi]^2-0.6036964688511505` \[Chi]^4+0.0873798215084077` \[Chi]^6);
im\[Omega][\[Chi]_] :=(0.014158792290965177` -0.036989395871554566` \[Chi]+0.026822526296575368` \[Chi]^2+0.0008490933750566702` \[Chi]^3-0.004843996907020524` \[Chi]^4-0.00014745235759327472` \[Chi]^5+0.0001504546201236794` \[Chi]^6)/(1-2.5900842798681376` \[Chi]+1.8952576220623967` \[Chi]^2-0.31416610693042507` \[Chi]^4+0.009002719412204133` \[Chi]^6);

\[Omega]RdDamping[int_, Erad_] = int/(1 - Erad);


\[Gamma]2[\[Eta]_, \[Chi]1_, \[Chi]2_] = (PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[6]])//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };
\[Gamma]3[\[Eta]_, \[Chi]1_, \[Chi]2_] =(PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[7]])//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


(* ::Section::Closed:: *)
(*MR Amplitude*)


\[Gamma]1[\[Eta]_, \[Chi]1_, \[Chi]2_] = Block[
	{dummy = FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[5]], res},
	res = PhenomCoeff[\[Eta], \[Chi]PN, dummy]//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };
	res Sqrt[\[Eta]]
]; 


Block[
	{
		ringdownFrequency, dampingFrequency, \[Gamma]1 = U\[Gamma]1[\[Eta], \[Chi]1, \[Chi]2], \[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]1, \[Chi]2], \[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]1, \[Chi]2],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]
	},
	\[Omega] = f M G;
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta], \[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];
	
	MRAmpExpr = \[Omega]^(-7/6) (M^2 \[Gamma]1 \[Gamma]3 dampingFrequency)/((\[Omega] - ringdownFrequency)^2 + (\[Gamma]3 dampingFrequency)^2) Exp[-((\[Gamma]2 (\[Omega] - ringdownFrequency))/(\[Gamma]3 dampingFrequency))];


]


(* ::Section::Closed:: *)
(*Intermediate*)


Block[
	{\[Omega]},
	\[Omega] = G M f;
	IntAmpfrequency = Sqrt[\[Eta]] M^2 \[Omega]^(-7/6) {1, \[Omega], \[Omega]^2, \[Omega]^3, \[Omega]^4}
]


With[{
	point1 = IntAmpfrequency//.f->f1,
	point2 = IntAmpfrequency//.f->f2,
	point3 = IntAmpfrequency//.f->f3,
	point1D = D[IntAmpfrequency, f]//.f->f1,
	point3D = D[IntAmpfrequency, f]//.f->f3},
	
	sol = LinearSolve[{point1, point2, point3, point1D, point3D}, {v1, f2^(-7/6) v2,v3,d1,d3}]//Simplify;
	sol = (sol//.{f1 -> 0.014/(M G)})//Simplify;
]


{
	\[Delta]0[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]1[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]2[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]3[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]4[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_]
} = sol;


v2[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = With[
	{
		dummy =(PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[4]])//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 +(1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 },
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]
		
	},
	dummy Sqrt[\[Eta]] M^2 (G M)^(-7/6)
];


{
		\[CapitalDelta]0[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]1[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]2[\[Eta]_, M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]3[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]4[\[Eta]_, M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_]
	} = {
		\[Delta]0[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],\[Delta]1[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[Delta]2[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],\[Delta]3[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[Delta]4[\[Eta],M,f2, f3,v1,v2,v3,d1,d3]
	}*((Sqrt[\[Eta]] M^2 \[Omega]^(-7/6))*{1, \[Omega], \[Omega]^2, \[Omega]^3, \[Omega]^4}//. \[Omega]-> G M)//Simplify;
	
{
		\[CapitalDelta]0[\[Eta]_, M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]1[\[Eta]_, M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]2[\[Eta]_, M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]3[\[Eta]_, M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]4[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_]
	} = {
		\[CapitalDelta]0[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]1[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]2[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]3[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]4[\[Eta],M,f2, f3,v1,v2,v3,d1,d3]
	}//.G-> UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]//Simplify;


Block[
	{ G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
		 
		ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]],
		dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]],  peakFrequency, \[Gamma]2, \[Gamma]3,
		f1, f2,f3, v1,v2,v3, d1, d3, \[Omega]
	},
	\[Omega] = M G f;

	\[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]1, \[Chi]2];
	\[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]1, \[Chi]2];
	
	peakFrequency = U\[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2, \[Gamma]3];
	f1 = 0.014/(G M); f3 = peakFrequency/(M G);  f2 =(0.014+peakFrequency)/(2 G M);
	
	v1 = InsAmpExpr//.f-> f1; v3 = MRAmpExpr//.f->f3; v2 = Uv2[\[Eta], M, \[Chi]1, \[Chi]2];
	d1 = D[InsAmpExpr, f]//.f-> f1; d3 = D[MRAmpExpr, f]//.f->f3;
	
	
	IntAmpExpr = {
		U\[CapitalDelta]0[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],U\[CapitalDelta]1[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],
		U\[CapitalDelta]2[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],U\[CapitalDelta]3[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],
		U\[CapitalDelta]4[\[Eta],M,f2,f3,v1,v2,v3,d1,d3]
	} . (f^(-7/6)*{1, f, f^2, f^3, f^4});
]


(* ::Section::Closed:: *)
(*\[ScriptCapitalA]IMR*)


\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins_, \[ScriptCapitalA]Int_, \[ScriptCapitalA]MR_, f_, M_, \[Omega]Peak_] = Block[
	{c = UnitConvert["SpeedOfLight", ("Gigaparsecs")/("Seconds")][[1]], 
	G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]},
	\[Omega] = f M G;
	Sqrt[5/6] c/(2 \[Pi]^(2/3)) G^2 (
		\[ScriptCapitalA]Ins us\[Theta][0.014 - \[Omega]] + \[ScriptCapitalA]Int us\[Theta][(\[Omega] - 0.014) (\[Omega]Peak - \[Omega])] + \[ScriptCapitalA]MR us\[Theta][(0.2 - \[Omega]) (\[Omega]- \[Omega]Peak)]
	)
];


\[ScriptCapitalA]IMRExpr = Block[
	{ 
		G= UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], 
		ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]],
		dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]],  peakFrequency, \[Gamma]2, \[Gamma]3,
		f1, f2,f3, v1,v2,v3, d1, d3, \[Omega]
	},
	\[Omega] = M G f;

	\[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]1, \[Chi]2];
	\[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]1, \[Chi]2];
	
	peakFrequency = U\[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2, \[Gamma]3];
	
	U\[ScriptCapitalA]IMR[InsAmpExpr, IntAmpExpr, MRAmpExpr, f, M, peakFrequency]
];


(* ::Section::Closed:: *)
(*Making Block function*)


AmpDefs = <||>;


varsDefs = HoldForm[{
	\[Delta] = Sqrt[1 - 4 \[Eta]],
	\[Chi]s = (\[Chi]1+\[Chi]2),
	\[Chi]a = (\[Chi]1-\[Chi]2),
	\[Chi]PN = Sqrt[1 - 4 \[Eta]] (\[Chi]1-\[Chi]2)/2  + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2,
	S = 1/4 (1+Sqrt[1 - 4 \[Eta]])^2 \[Chi]1 + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 \[Chi]2,
	Shat = (1/4 (1+Sqrt[1 - 4 \[Eta]])^2 \[Chi]1 + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 \[Chi]2)/(1- 2 \[Eta])
}];


AmpDefs["Ins"] = Block[
	{rule, sym, defs},
	sym = {a0,a1,a2,a3,a4,a5,a6,r1,r2,r3};
	defs = {
		\[ScriptCapitalA]0[\[Eta],M],\[ScriptCapitalA]1[\[Eta],M],\[ScriptCapitalA]2[\[Eta],M],\[ScriptCapitalA]3[\[Eta],M,\[Chi]1,\[Chi]2],\[ScriptCapitalA]4[\[Eta],M,\[Chi]1,\[Chi]2],\[ScriptCapitalA]5[\[Eta],M,\[Chi]1,\[Chi]2],\[ScriptCapitalA]6[\[Eta],M,\[Chi]1,\[Chi]2],
		\[Rho]1[\[Eta],M,\[Chi]1,\[Chi]2],\[Rho]2[\[Eta],M,\[Chi]1,\[Chi]2],\[Rho]3[\[Eta],M,\[Chi]1,\[Chi]2]
	};
	
	rule = MapThread[
		Rule,
		{sym, defs}
	];
	
	
	
	HoldForm[{
		\[ScriptCapitalA]0[\[Eta]_,  M_] = a0, 
		\[ScriptCapitalA]1[\[Eta]_, M_] = a1, 
		\[ScriptCapitalA]2[\[Eta]_, M_] = a2, 
		\[ScriptCapitalA]3[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = a3, 
		\[ScriptCapitalA]4[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = a4, 
		\[ScriptCapitalA]5[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = a5, 
		\[ScriptCapitalA]6[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = a6,
		\[Rho]1[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = r1,
		\[Rho]2[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = r2, 
		\[Rho]3[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = r3
	}]/.rule
];


Block[
	{
		ringdownFrequency, dampingFrequency, \[Gamma]1 = U\[Gamma]1[\[Eta], \[Chi]1, \[Chi]2], \[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]1, \[Chi]2], \[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]1, \[Chi]2],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]
	},
	\[Omega] = f M G;
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta], \[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];
	
	MRAmpExpr = \[Omega]^(-7/6) (\[Gamma]1 \[Gamma]3 dampingFrequency)/((\[Omega] - ringdownFrequency)^2 + (\[Gamma]3 dampingFrequency)^2) Exp[-((\[Gamma]2 (\[Omega] - ringdownFrequency))/(\[Gamma]3 dampingFrequency))];


]


AmpDefs["MR"] = HoldForm[{
	\[Gamma]1[\[Eta]_, \[Chi]1_, \[Chi]2_] = g1,
	\[Gamma]2[\[Eta]_, \[Chi]1_, \[Chi]2_] = g2,
	\[Gamma]3[\[Eta]_, \[Chi]1_, \[Chi]2_] = g3,
	aeff[\[Eta]_, \[Chi]1_, \[Chi]2_] = a ,
	re\[Omega][\[Chi]_] = r,
	im\[Omega][\[Chi]_] = i,
	Erad[\[Eta]_, \[Chi]1_, \[Chi]2_] = erad,
	\[Omega]RdDamping[int_, Erad_] = s,
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ]
}]/.{
	g1 -> \[Gamma]1[\[Eta], \[Chi]1, \[Chi]2], g2 -> \[Gamma]2[\[Eta], \[Chi]1, \[Chi]2], g3 -> \[Gamma]3[\[Eta], \[Chi]1, \[Chi]2],
	a-> aeff[\[Eta], \[Chi]1, \[Chi]2], r -> re\[Omega][\[Chi]], i-> im\[Omega][\[Chi]], erad-> Erad[\[Eta], \[Chi]1, \[Chi]2], s-> \[Omega]RdDamping[int, Erad]
};


AmpDefs["Int"] = Block[
	{rule, defs},
	defs = {
		\[CapitalDelta]0[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],\[CapitalDelta]1[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]2[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]3[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]4[\[Eta],M,f2,f3,v1,v2,v3,d1,d3],
		v2[\[Eta], M, \[Chi]1, \[Chi]2]
	};
	
	rule = MapThread[
		Rule,
		{{D0,D1,D2,D3,D4,V2},defs}
	];

	HoldForm[{
		\[CapitalDelta]0[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] = D0,
		\[CapitalDelta]1[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] =D1,
		\[CapitalDelta]2[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] = D2,
		\[CapitalDelta]3[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] = D3,
		\[CapitalDelta]4[\[Eta]_, M_,f2_,f3_, v1_, v2_, v3_, d1_, d3_] = D4,
		v2[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] = V2
	}]/.rule

];


AmpDefs["\[ScriptCapitalA]IMR"] = HoldForm[{
	\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins_, \[ScriptCapitalA]Int_, \[ScriptCapitalA]MR_, f_, M_, \[Omega]Peak_] = x
}]/.x->\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins, \[ScriptCapitalA]Int, \[ScriptCapitalA]MR, f, M, \[Omega]Peak];


AmpDefs["Commute D"] = HoldForm[{
	(*Commute D with \[Phi]IMR*)
	\[ScriptCapitalA]IMR/: D[\[ScriptCapitalA]IMR[x__], n___] := With[
		{args = D[#, n]&/@({x}[[1;;3]])},
		\[ScriptCapitalA]IMR@@(Join[args, {x}[[4;;6]]])
	]
}];


(*Make rules to eliminate U*)
AllFunctionHeads = (AmpDefs//Values)[[All, 1, All, 1,0]]//Flatten//DeleteDuplicates;
AllFunctionHeads = DeleteElements[AllFunctionHeads, {Symbol}];
LocalVars = varsDefs[[1,All,1]];


Module[{undefined =  ("U" <> #)&/@(ToString/@AllFunctionHeads)},
	
	undefined = ToExpression[undefined,InputForm];
	Headrules = Thread@Rule[undefined, AllFunctionHeads]
];


FinalExpr = Module[
	{dummy},
	dummy = HoldForm[
		{aIMR}
	]/.aIMR->\[ScriptCapitalA]IMRExpr;

	dummy//.Headrules
];


defs = Module[
	{functionDefs = AmpDefs//Values, orderedDefs},
	functionDefs = HoldForm@@@#&/@functionDefs; (*{HoldForm[HoldForm[]], HoldForm[HoldForm[]],...}*)
	orderedDefs = {HoldForm@@@varsDefs, functionDefs, HoldForm@@@FinalExpr}//Flatten; (*{HoldForm[HoldForm[]], HoldForm[HoldForm[]],...}*)
	orderedDefs = Flatten[HoldForm@@orderedDefs]; (*HoldForm[allDefs]*)
	orderedDefs = CompoundExpression@@@(HoldForm[Evaluate@orderedDefs]) 
];


declarations = HoldForm[Evaluate@{AllFunctionHeads, Sequence@@LocalVars}];


$blockexpr = $Block@@@HoldForm[Join[declarations, defs]//Evaluate];


With[
	{list1 =varsDefs[[1, All,1]], list2 = varsDefs[[1, All,2]]},
	KeepDefs = Thread@Rule[list1,list2];
]


expr =  Hold[
	{test, {M, \[Eta], \[Chi]1, \[Chi]2}, {{}}},
	Evaluate@$blockexpr,
	"KeepDefs" -> KeepDefs
]//.HoldForm[X_] :> X;


(* ::Section:: *)
(*Calculating derivatives*)


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {M, \[Eta], \[Chi]1,\[Chi]2};


(*All derivatives up to order 3*)
derivatives = Combinations[vars, 5];


(*We need only 1 \[Delta]p_i at a time, and no more than one derivative because they are linear in the phase*)
derivatives//Length


PrependTo[derivatives, {}]


Combinations[vars, 3]//Length
Combinations[vars, 4]//Length


derivatives[[36;;70]]


expr =  Hold[
	{\[ScriptA]IMR, {f, M, \[Eta], \[Chi]1, \[Chi]2}, {i}},
	Evaluate@$blockexpr,
	"KeepDefs" -> KeepDefs,
	"IncludeZeroDerivative"->False
]//.HoldForm[X_] :> X;


Clear@Ds;
expr2 = HoldForm[EchoTiming[MemoryConstrained[QuietEcho[x], 6 10^9], ToString[i]], Evaluate[{i,derivatives[[36;;70]] }]]/.x-> expr;
expr2 = expr2/.Hold-> DerivativeRules;
Ds = Table@@expr2;
Export["Amplitude_Ds_order_4.wdx", Ds]
