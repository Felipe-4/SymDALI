(* ::Package:: *)

(*SetOptions[EvaluationNotebook[], DefaultNewCellStyle->"Code"];
SetOptions[EvaluationNotebook[], WindowElements->{"MemoryMonitor","VerticalScrollBar","MenuBar", "HorizontalScrollBar"}]
SetDirectory[NotebookDirectory[]];

Get["maTHEMEatica.wl"];

colors=<|
	"background"->RGBColor["#000000"],
	"fontcolor"->RGBColor["#eeeeee"],
	"primary"->RGBColor["#B87333"],
	"variable"->RGBColor["#55f7df"],
	"module"->RGBColor["#e638e9"],
	"block"->RGBColor["#FFFF00"],
	"error"->RGBColor["#FF0000"],
	"headhighlight"->RGBColor["#02584c"]
|>;
SetColors[colors];
CreateStyleSheet[];
ApplyStyleSheet[];*)


SetOptions[EvaluationNotebook[], NotebookAutoSave->True];


(*Basic stuff:*)
$HistoryLength = 1;
PacletDirectoryLoad["/home/cosmo-ufes/Documentos/GitHub"];
<<FelipeBarbosa`SymDALI`

vectorDefs = Module[
	{DALIDir, docDir},
	DALIDir = FindFile["FelipeBarbosa`SymDALI`"]//FileNameDrop//ParentDirectory;
	docDir = FileNameJoin[{DALIDir, "Documentation/English/Tutorials"}];
	Import@FileNameJoin[{docDir, "IMRPhenomDComponentsDefinition.wdx"}]
];

SymRules = <||>;
NRules = <||>;


(* ::Chapter:: *)
(*Phase*)


(* ::Section::Closed:: *)
(*PPN-parameters:*)


(* Define the coefficients *)
\[Phi]0 = 1;
\[Phi]1 = 0;
\[Phi]2[\[Eta]_] = 3715/756 + 55 \[Eta]/9;
\[Phi]3[\[Eta]_, \[Chi]1_, \[Chi]2_] = (-16 \[Pi] + (113 \[Delta] \[Chi]a)/3 + (113/3 - 76 \[Eta]/3) \[Chi]s)//.{\[Chi]s->(\[Chi]1+\[Chi]2)/2, \[Chi]a->(\[Chi]1-\[Chi]2)/2, \[Delta] -> Sqrt[1 - 4 \[Eta]]};

\[Phi]4[\[Eta]_, \[Chi]1_, \[Chi]2_] = (15293365/508032 + 27145 \[Eta]/504 + 3085 \[Eta]^2/72 + (-405/8 + 200 \[Eta]) \[Chi]a^2 - (405/4) \[Delta] \[Chi]a \[Chi]s + (-405/8 + 5 \[Eta]/2) \[Chi]s^2)//.{\[Chi]s->(\[Chi]1+\[Chi]2)/2, \[Chi]a->(\[Chi]1-\[Chi]2)/2, \[Delta] -> Sqrt[1 - 4 \[Eta]]};

\[Phi]5[\[Eta]_, \[Chi]1_, \[Chi]2_] = (1 + Log[\[Pi] M \[Omega]]) * (38645 \[Pi]/756 - 65 \[Pi] \[Eta]/9 + 
    \[Delta] (-732985/2268 - 140 \[Eta]/9) \[Chi]a + 
    (-732985/2268 + 24260 \[Eta]/81 + 340 \[Eta]^2/9) \[Chi]s)//.{\[Chi]s->(\[Chi]1+\[Chi]2)/2, \[Chi]a->(\[Chi]1-\[Chi]2)/2, \[Delta] -> Sqrt[1 - 4 \[Eta]], Log[x_] ->0};

\[Phi]6[\[Eta]_, \[Chi]1_, \[Chi]2_] = (11583231236531/4694215680 - (6848 EulerGamma)/21 - 
   (640 \[Pi]^2)/3 + (-15737765635/3048192 + (2255 \[Pi]^2)/12) \[Eta] + 
   76055 \[Eta]^2/1728 - 127825 \[Eta]^3/1296 - 
   (6848/63) Log[64 \[Pi] M \[Omega]] + (2270/3) \[Pi] \[Delta] \[Chi]a + 
   ((2270 \[Pi])/3 - 520 \[Pi] \[Eta]) \[Chi]s)//.{\[Chi]s->(\[Chi]1+\[Chi]2)/2, \[Chi]a->(\[Chi]1-\[Chi]2)/2, \[Delta] -> Sqrt[1 - 4 \[Eta]], Log[x_]-> Log[64]};

\[Phi]7[\[Eta]_, \[Chi]1_, \[Chi]2_] = (77096675 \[Pi]/254016 + (378515 \[Pi] \[Eta])/1512 - (74045 \[Pi] \[Eta]^2)/756 + 
   \[Delta] (-25150083775/3048192 + (26804935 \[Eta])/6048 - (1985 \[Eta]^2)/48) \[Chi]a + 
   (-25150083775/3048192 + (10566655595 \[Eta])/762048 - 
      (1042165 \[Eta]^2)/3024 + (5345 \[Eta]^3)/36) \[Chi]s)//.{\[Chi]s->(\[Chi]1+\[Chi]2)/2, \[Chi]a->(\[Chi]1-\[Chi]2)/2, \[Delta] -> Sqrt[1 - 4 \[Eta]]};


(* ::Section::Closed:: *)
(*Ins-Phase*)


(* ::Text:: *)
(*In  the  notation  of  the  arXiv : 1903.04467 v3  we  want  a  vector*)
(*{\[CurlyPhi]minus2, \[CurlyPhi]0, \[CurlyPhi]1, \[CurlyPhi]2, \[CurlyPhi]3, \[CurlyPhi]4, \[CurlyPhi]5, \[CurlyPhi]5l, \[CurlyPhi]6, \[CurlyPhi]6l, \[CurlyPhi]7} . In  GR  there  is  no  \[CurlyPhi]minus2  or  \[CurlyPhi]1, we  shall  set  their  values  to  1*3/(128 \[Eta])*)
(*and  multiply  by  a  vector  1 + \[Delta]\[CurlyPhi]  where  the  default  value  of  \[Delta]\[CurlyPhi]minus2  and  \[Delta]\[CurlyPhi]1  is - 1 (to recover GR) and -1+\[Delta]\[CurlyPhi] when we want to sample on them:*)


Block[
	{},
	insVecPhase[\[Eta]_, \[Chi]1_, \[Chi]2_] = 3*{
		1, (*\[CurlyPhi]minus2*)
		1, (*\[CurlyPhi]0*)
		1, (*\[CurlyPhi]1*)
		\[Phi]2[\[Eta]], (*\[CurlyPhi]2*)
		\[Phi]3[\[Eta], \[Chi]1, \[Chi]2], (*\[CurlyPhi]3*)
		\[Phi]4[\[Eta], \[Chi]1, \[Chi]2], (*\[CurlyPhi]4*)
		\[Phi]5[\[Eta], \[Chi]1, \[Chi]2], (*\[CurlyPhi]5*)
		\[Phi]5[\[Eta], \[Chi]1, \[Chi]2], (*\[CurlyPhi]5l*)
		\[Phi]6[\[Eta], \[Chi]1, \[Chi]2], (*\[CurlyPhi]6*)
		-(6848/63), (*\[CurlyPhi]6l*)
		\[Phi]7[\[Eta], \[Chi]1, \[Chi]2] (*\[CurlyPhi]7*)
	}/(128 \[Eta]);
	
	\[Omega]InsVecPhase[\[Omega]_] = Join[
		{(\[Pi] \[Omega])^(-7/3)}, (*-2*)
		(\[Pi] \[Omega])^((#-5)/3)&/@Range[0,5], (*0 to 5*)
		{Log[\[Pi] \[Omega]], (\[Pi] \[Omega])^(1/3), Log[\[Pi] \[Omega]] (\[Pi] \[Omega])^(1/3), (\[Pi] \[Omega])^(2/3)} (*5l, 6, 6l, 7*)
	]
];


(* ::Text:: *)
(*Now, introduce the \[Delta]s:*)


Block[
	{v = 1 + {-1 + \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, -1 + \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}},
	InsVecPhase[\[Eta]_, \[Chi]1_,\[Chi]2_, \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_] = (
		v*insVecPhase[\[Eta], \[Chi]1, \[Chi]2]
	)
];


{
	\[CapitalPhi]minus2[\[Eta]_, \[Delta]\[CurlyPhi]minus2_], \[CapitalPhi]0[\[Eta]_, \[Delta]\[CurlyPhi]0_], \[CapitalPhi]1[\[Eta]_, \[Delta]\[CurlyPhi]1_], \[CapitalPhi]2[\[Eta]_, \[Delta]\[CurlyPhi]2_], \[CapitalPhi]3[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[CurlyPhi]3_], 
	\[CapitalPhi]4[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[CurlyPhi]4_], \[CapitalPhi]5[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[CurlyPhi]5_], \[CapitalPhi]5l[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[CurlyPhi]5l_], \[CapitalPhi]6[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[CurlyPhi]6_], 
	\[CapitalPhi]6l[\[Eta]_, \[Delta]\[CurlyPhi]6l_], \[CapitalPhi]7[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[CurlyPhi]7_]
} = InsVecPhase[\[Eta], \[Chi]1, \[Chi]2, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];


{\[CapitalSigma]1[\[Eta]_, \[Chi]1_, \[Chi]2_], \[CapitalSigma]2[\[Eta]_, \[Chi]1_, \[Chi]2_], \[CapitalSigma]3[\[Eta]_, \[Chi]1_, \[Chi]2_], \[CapitalSigma]4[\[Eta]_, \[Chi]1_, \[Chi]2_]} = \[Eta]^-1*(PhenomCoeff[\[Eta], \[Chi]PN, #]&/@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[8;;11]])//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };


{\[CapitalSigma]1[\[Eta],\[Chi]1,\[Chi]2],\[CapitalSigma]2[\[Eta],\[Chi]1,\[Chi]2],\[CapitalSigma]3[\[Eta],\[Chi]1,\[Chi]2],\[CapitalSigma]4[\[Eta],\[Chi]1,\[Chi]2]} . {\[Omega], 3/4 \[Omega]^(4/3), 3/5 \[Omega]^(5/3), 1/2 \[Omega]^2}


InsExpr = -\[Pi]/4 + {
	U\[CapitalPhi]minus2[\[Eta],\[Delta]\[CurlyPhi]minus2],U\[CapitalPhi]0[\[Eta],\[Delta]\[CurlyPhi]0],U\[CapitalPhi]1[\[Eta],\[Delta]\[CurlyPhi]1],U\[CapitalPhi]2[\[Eta],\[Delta]\[CurlyPhi]2],U\[CapitalPhi]3[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]3],U\[CapitalPhi]4[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]4],U\[CapitalPhi]5[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]5],
	U\[CapitalPhi]5l[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]5l],U\[CapitalPhi]6[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]6],U\[CapitalPhi]6l[\[Eta],\[Delta]\[CurlyPhi]6l],U\[CapitalPhi]7[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]7]
} . \[Omega]InsVecPhase[\[Omega]]  + {U\[CapitalSigma]1[\[Eta],\[Chi]1,\[Chi]2],U\[CapitalSigma]2[\[Eta],\[Chi]1,\[Chi]2],U\[CapitalSigma]3[\[Eta],\[Chi]1,\[Chi]2],U\[CapitalSigma]4[\[Eta],\[Chi]1,\[Chi]2]} . {\[Omega], 3/4 \[Omega]^(4/3), 3/5 \[Omega]^(5/3), 1/2 \[Omega]^2};


(* ::Section::Closed:: *)
(*Int-Phase*)


(*Find vectors for Intermediate Phase and their variables:*)
Position[vectorDefs[[All, 1, All, 0]] /. HoldPattern -> Identity, #] & /@ {IntVecPhase, \[Omega]IntVecPhase}

vectorDefs[[3 ;; 4, 1]]

DownValues[IntVecPhase] = {vectorDefs[[3]]};
DownValues[\[Omega]IntVecPhase] = {vectorDefs[[4]]};


(* ::Text:: *)
(*Now  Introduce \[Delta]\[Beta]:*)


{\[Beta]1[\[Eta]_, \[Chi]1_, \[Chi]2_], \[Beta]2[\[Eta]_, \[Chi]1_, \[Chi]2_], \[Beta]3[\[Eta]_, \[Chi]1_, \[Chi]2_]} = \[Eta] IntVecPhase[\[Eta], \[Chi]PN]//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };


IntVecPhase[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Beta]2_, \[Delta]\[Beta]3_] = (1 + {0, \[Delta]\[Beta]2, \[Delta]\[Beta]3})*{\[Beta]1[\[Eta], \[Chi]1, \[Chi]2], \[Beta]2[\[Eta], \[Chi]1, \[Chi]2], \[Beta]3[\[Eta], \[Chi]1, \[Chi]2]};


{\[CapitalBeta]1[\[Eta]_, \[Chi]1_, \[Chi]2_], \[CapitalBeta]2[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Beta]2_], \[CapitalBeta]3[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Beta]3_]} = (IntVecPhase[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Beta]2,\[Delta]\[Beta]3]/\[Eta])//Simplify;


IntExpr = {U\[CapitalBeta]1[\[Eta], \[Chi]1, \[Chi]2], U\[CapitalBeta]2[\[Eta], \[Chi]1, \[Chi]2, \[Delta]\[Beta]2],U\[CapitalBeta]3[\[Eta], \[Chi]1, \[Chi]2, \[Delta]\[Beta]3]} . \[Omega]IntVecPhase[\[Omega]];


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


\[Gamma]2[\[Eta]_, \[Chi]PN_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[6]];
\[Gamma]3[\[Eta]_, \[Chi]PN_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[7]];


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


(* ::Section::Closed:: *)
(*MR-Phase*)


Position[vectorDefs[[All, 1, All, 0]]/.HoldPattern->Identity, #]&/@{MRVecPhase, \[Omega]MRVecPhase}


vectorDefs[[5;;6,1]]
DownValues[MRVecPhase] = {vectorDefs[[5]]};
DownValues[\[Omega]MRVecPhase] = {vectorDefs[[6]]};


{\[Alpha]1[\[Eta]_, \[Chi]1_, \[Chi]2_], \[Alpha]2[\[Eta]_, \[Chi]1_, \[Chi]2_], \[Alpha]3[\[Eta]_, \[Chi]1_, \[Chi]2_], \[Alpha]4[\[Eta]_, \[Chi]1_, \[Chi]2_]} = \[Eta] MRVecPhase[\[Eta], \[Chi]PN]//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };


\[Alpha]5[\[Eta]_, \[Chi]1_, \[Chi]2_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[19]]//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };


(* ::Text:: *)
(*Introduce  \[Delta]\[Alpha] :*)


MRVecPhase[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Alpha]2_, \[Delta]\[Alpha]3_, \[Delta]\[Alpha]4_] = (1 + {0, \[Delta]\[Alpha]2, \[Delta]\[Alpha]3, \[Delta]\[Alpha]4})*{\[Alpha]1[\[Eta], \[Chi]1, \[Chi]2], \[Alpha]2[\[Eta], \[Chi]1, \[Chi]2], \[Alpha]3[\[Eta], \[Chi]1, \[Chi]2], \[Alpha]4[\[Eta], \[Chi]1, \[Chi]2]};


{\[CapitalAlpha]1[\[Eta]_, \[Chi]1_, \[Chi]2_], \[CapitalAlpha]2[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Alpha]2_], \[CapitalAlpha]3[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Alpha]3_], \[CapitalAlpha]4[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Alpha]4_]} =  MRVecPhase[\[Eta], \[Chi]1, \[Chi]2, \[Delta]\[Alpha]2, \[Delta]\[Alpha]3, \[Delta]\[Alpha]4]/\[Eta]//Simplify;


MRExpr = Block[
	{ringdownFrequency, dampingFrequency, \[Alpha]5},
	\[Alpha]5 = U\[Alpha]5[\[Eta], \[Chi]1, \[Chi]2];
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];
	 {U\[CapitalAlpha]1[\[Eta],\[Chi]1,\[Chi]2],U\[CapitalAlpha]2[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]2],U\[CapitalAlpha]3[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]3],U\[CapitalAlpha]4[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]4]} . \[Omega]MRVecPhase[\[Omega], ringdownFrequency, dampingFrequency, \[Alpha]5]
]//Simplify;


(* ::Section::Closed:: *)
(*C(1)*)


C1[beta0_, beta1\[Omega]_, alpha0_, alpha1\[Omega]_, \[Omega]_, ringdownFrequency_] = (
			(beta0 + beta1\[Omega])*us\[Theta][(ringdownFrequency/2 - \[Omega]) (\[Omega]-0.018)] + 
			(alpha0 +alpha1\[Omega])us\[Theta][(\[Omega] - ringdownFrequency/2) (0.2-\[Omega])]
);


Module[{
	ringdownFrequency, dampingFrequency, \[Alpha]5
	},
	\[Alpha]5 = U\[Alpha]5[\[Eta], \[Chi]1, \[Chi]2];
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]];

	\[Delta]\[Beta]1[d\[Phi]Ins_, d\[Phi]Int_] = d\[Phi]Ins - d\[Phi]Int;
	\[Delta]\[Beta]0[\[CapitalPhi]Ins_, \[CapitalPhi]Int_, \[Delta]\[Beta]1_] = \[CapitalPhi]Ins - \[CapitalPhi]Int - \[Delta]\[Beta]1 0.018;

	\[Delta]\[Alpha]1[d\[CapitalPhi]Int_, d\[CapitalPhi]MR_, \[Delta]\[Beta]1_] = d\[CapitalPhi]Int + \[Delta]\[Beta]1 - d\[CapitalPhi]MR;
	
	\[Delta]\[Alpha]0[\[CapitalPhi]Int_, \[CapitalPhi]MR_, \[Delta]\[Beta]0_, \[Delta]\[Beta]1\[Omega]rd_, \[Delta]\[Alpha]1\[Omega]rd_] := \[CapitalPhi]Int + \[Delta]\[Beta]1\[Omega]rd/2 + \[Delta]\[Beta]0 - \[CapitalPhi]MR - \[Delta]\[Alpha]1\[Omega]rd/2;
	
	Block[
		{dphiIns, dphiInt1,dphiInt2,dphiMR, phiIns, phiInt1, phiInt2,phiMR, beta0, beta1, alpha0, alpha1},
		
		(*derivative functions ###################################################### *)
		dphiIns = D[InsExpr, \[Omega]]//.\[Omega]->0.018;

		dphiInt1 = D[IntExpr, \[Omega]]//.\[Omega]->0.018;
		
		dphiInt2 = D[IntExpr, \[Omega]]//.\[Omega]->ringdownFrequency/2;
		
		dphiMR = D[MRExpr, \[Omega]]//.\[Omega]->ringdownFrequency/2;
		
		(*normal functions ########################################################################*)
		phiIns = InsExpr//.\[Omega]->0.018;

		phiInt1 = IntExpr//.\[Omega]->0.018;
		
		phiInt2 =IntExpr//.\[Omega]-> ringdownFrequency/2;
		
		phiMR = MRExpr//.\[Omega]-> ringdownFrequency/2;
	
		
		beta1 = U\[Delta]\[Beta]1[dphiIns, dphiInt1];
		beta0 = U\[Delta]\[Beta]0[phiIns, phiInt1, beta1];
		
		alpha1 = U\[Delta]\[Alpha]1[dphiInt2, dphiMR, beta1];
		alpha0 = U\[Delta]\[Alpha]0[phiInt2, phiMR, beta0, beta1*ringdownFrequency, alpha1*ringdownFrequency];
		
		
		
		exprC1 = UC1[beta0, beta1*\[Omega], alpha0, alpha1*\[Omega], \[Omega], ringdownFrequency];
	]
];


(* ::Section::Closed:: *)
(*\[Phi]IMR*)


\[Phi]IMR[\[Phi]Ins_, \[Phi]Int_, \[Phi]MR_, c1_, ringdownFrequency_, \[Omega]_] := (
	\[Phi]Ins*us\[Theta][-\[Omega]+0.018] +
	\[Phi]Int*us\[Theta][(ringdownFrequency/2 - \[Omega]) (\[Omega]-0.018)]  +
	\[Phi]MR*us\[Theta][(\[Omega] - ringdownFrequency/2) (0.2-\[Omega])] + 
	c1
)


Module[
	{ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]],
	\[Phi]Ins, \[Phi]Int, \[Phi]MR},
	
	\[Phi]Ins = InsExpr;
	
	\[Phi]Int = IntExpr;
	
	\[Phi]MR = MRExpr;
	
	\[Phi]IMRexpr = U\[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, exprC1, ringdownFrequency, \[Omega]];
	\[Phi]IMRexprRef = U\[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, exprC1, ringdownFrequency, \[Omega]]//.\[Omega]->\[Omega]ref;
]//Simplify;


(* ::Section::Closed:: *)
(*t0*)


\[CapitalDelta]t0[\[Omega]_, d\[Phi]MR_, \[Omega]ref_] = d\[Phi]MR (\[Omega] - \[Omega]ref);

Module[
	{d\[Phi]MR, ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]],
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Eta],\[Chi]1, \[Chi]2]], UErad[\[Eta], \[Chi]1, \[Chi]2]], \[Alpha]5, peakFrequency, \[Gamma]2, \[Gamma]3},
	\[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]PN];
	\[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]PN];
	
	peakFrequency = U\[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2, \[Gamma]3];
	
	\[Alpha]5 = U\[Alpha]5[\[Eta], \[Chi]1, \[Chi]2];
	
	
	d\[Phi]MR =D[MRExpr, \[Omega]]//.\[Omega]->peakFrequency;
		
	\[CapitalDelta]t0expr = U\[CapitalDelta]t0[\[Omega], d\[Phi]MR, \[Omega]ref]
];


(* ::Text:: *)
(*Pos[expri, hj] ={POS1, POS2, ...},  j = 1, ..., N*)


(* ::Section::Closed:: *)
(*Making Block function*)


(* ::Text:: *)
(*First all variable declarations:*)


varsDefs = HoldForm[{
	\[Delta] = Sqrt[1 - 4 \[Eta]],
	\[Chi]s = (\[Chi]1+\[Chi]2)/2,
	\[Chi]a = (\[Chi]1-\[Chi]2)/2,
	\[Chi]PN = Sqrt[1 - 4 \[Eta]] (\[Chi]1-\[Chi]2)/2  + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2,
	S = 1/4 (1+Sqrt[1 - 4 \[Eta]])^2 \[Chi]1 + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 \[Chi]2,
	Shat = (1/4 (1+Sqrt[1 - 4 \[Eta]])^2 \[Chi]1 + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 \[Chi]2)/(1- 2 \[Eta])
}];


(* ::Text:: *)
(*Now, all function definitions, starting at Primitives and ending on \[Phi]IMR:*)


FunctionDefs = <||>;


FunctionDefs["Ins+Int"]=With[{
	rules =Thread@Rule[
		{-2,0,1,2,3,4,5,5l,6,6l,7},
		{
			\[CapitalPhi]minus2[\[Eta],\[Delta]\[CurlyPhi]minus2],\[CapitalPhi]0[\[Eta],\[Delta]\[CurlyPhi]0],\[CapitalPhi]1[\[Eta],\[Delta]\[CurlyPhi]1],\[CapitalPhi]2[\[Eta],\[Delta]\[CurlyPhi]2],\[CapitalPhi]3[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]3],\[CapitalPhi]4[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]4],
			\[CapitalPhi]5[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]5],\[CapitalPhi]5l[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]5l],\[CapitalPhi]6[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]6],\[CapitalPhi]6l[\[Eta],\[Delta]\[CurlyPhi]6l],\[CapitalPhi]7[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]7]
		}
	]},

 
   HoldForm[{
	\[CapitalPhi]minus2[\[Eta]_, \[Delta]\[CurlyPhi]minus2_]=-2,
	\[CapitalPhi]0[\[Eta]_,\[Delta]\[CurlyPhi]0_]=0,
	\[CapitalPhi]1[\[Eta]_,\[Delta]\[CurlyPhi]1_]=1,
	\[CapitalPhi]2[\[Eta]_,\[Delta]\[CurlyPhi]2_]=2,
	\[CapitalPhi]3[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]3_]=3,
	\[CapitalPhi]4[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]4_]=4,
	\[CapitalPhi]5[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]5_]=5,
	\[CapitalPhi]5l[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]5l_]=5l,
	\[CapitalPhi]6[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]6_]=6,
	\[CapitalPhi]6l[\[Eta]_,\[Delta]\[CurlyPhi]6l_]=6l,
	\[CapitalPhi]7[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]7_]=7,
	\[CapitalSigma]1[\[Eta]_, \[Chi]1_, \[Chi]2_] = A, 
	\[CapitalSigma]2[\[Eta]_, \[Chi]1_, \[Chi]2_] = B, 
	\[CapitalSigma]3[\[Eta]_, \[Chi]1_, \[Chi]2_] = c, 
	\[CapitalSigma]4[\[Eta]_, \[Chi]1_, \[Chi]2_] = d,
	\[CapitalBeta]1[\[Eta]_, \[Chi]1_, \[Chi]2_] = x,
	\[CapitalBeta]2[\[Eta]_,\[Chi]1_, \[Chi]2_, \[Delta]\[Beta]2_] = y,
	\[CapitalBeta]3[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Beta]3_] =z
	
}]/.Join[
	rules,
	{x-> \[CapitalBeta]1[\[Eta],\[Chi]1, \[Chi]2], y-> \[CapitalBeta]2[\[Eta],\[Chi]1, \[Chi]2, \[Delta]\[Beta]2],z-> \[CapitalBeta]3[\[Eta],\[Chi]1, \[Chi]2, \[Delta]\[Beta]3]},
	{A ->\[CapitalSigma]1[\[Eta], \[Chi]1, \[Chi]2], B -> \[CapitalSigma]2[\[Eta], \[Chi]1, \[Chi]2], c -> \[CapitalSigma]3[\[Eta], \[Chi]1, \[Chi]2], d -> \[CapitalSigma]4[\[Eta], \[Chi]1, \[Chi]2] }
]

];


FunctionDefs["MR"] = HoldForm[{
	aeff[\[Eta]_, \[Chi]1_, \[Chi]2_] := x,
	Erad[\[Eta]_, \[Chi]1_, \[Chi]2_] := y,
	
	re\[Omega][\[Chi]_] := z,
	im\[Omega][\[Chi]_] := w,
	\[Omega]RdDamping[int_, Erad_] := \[ScriptX],
	\[Gamma]2[\[Eta]_, \[Chi]PN_] := \[ScriptY],
	\[Gamma]3[\[Eta]_, \[Chi]PN_] := \[ScriptZ],
	
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] := \[ScriptW],
	\[Alpha]5[\[Eta]_, \[Chi]1_, \[Chi]2_] := \[ScriptA],
	
	\[CapitalAlpha]1[\[Eta]_,\[Chi]1_,\[Chi]2_] = 1,
	\[CapitalAlpha]2[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[Alpha]2_] =2 ,
	\[CapitalAlpha]3[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[Alpha]3_] = 3,
	\[CapitalAlpha]4[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[Alpha]4_] = 4
	
}]/.{
	x -> aeff[\[Eta],\[Chi]1, \[Chi]2], y-> Erad[\[Eta],\[Chi]1, \[Chi]2], z-> re\[Omega][\[Chi]], w-> im\[Omega][\[Chi]],
	\[ScriptX]-> \[Omega]RdDamping[int, Erad], \[ScriptY]-> \[Gamma]2[\[Eta], \[Chi]PN], \[ScriptZ]-> \[Gamma]3[\[Eta], \[Chi]PN], \[ScriptW]-> \[Omega]Peak[\[Omega]RD, \[Omega]DAMP, \[Gamma]2, \[Gamma]3],\[ScriptA]-> \[Alpha]5[\[Eta], \[Chi]1, \[Chi]2],
	mr-> MRVecPhase[\[Eta], \[Chi]PN, \[Delta]\[Alpha]2, \[Delta]\[Alpha]3, \[Delta]\[Alpha]4], \[Omega]mr-> \[Omega]MRVecPhase[\[Omega], \[Omega]RD, \[Omega]DAMP, \[Alpha]5],
	1 -> \[CapitalAlpha]1[\[Eta],\[Chi]1,\[Chi]2], 2-> \[CapitalAlpha]2[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]2], 3 ->\[CapitalAlpha]3[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]3], 4 -> \[CapitalAlpha]4[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]4]

};


FunctionDefs["C1"] = HoldForm[{
	
	C1[beta0_, beta1\[Omega]_, alpha0_, alpha1\[Omega]_, \[Omega]_,ringdownFrequency_]  := w,
	
	\[Delta]\[Beta]1[d\[Phi]Ins_, d\[Phi]Int_] := d\[Phi]Ins - d\[Phi]Int,
	
	\[Delta]\[Beta]0[\[CapitalPhi]Ins_, \[CapitalPhi]Int_, \[Delta]\[Beta]1_] := \[CapitalPhi]Ins - \[CapitalPhi]Int - \[Delta]\[Beta]1 0.018,

	\[Delta]\[Alpha]1[d\[CapitalPhi]Int_, d\[CapitalPhi]MR_, \[Delta]\[Beta]1_] := d\[CapitalPhi]Int + \[Delta]\[Beta]1 - d\[CapitalPhi]MR,
	
	\[Delta]\[Alpha]0[\[CapitalPhi]Int_, \[CapitalPhi]MR_, \[Delta]\[Beta]0_, \[Delta]\[Beta]1\[Omega]rd_, \[Delta]\[Alpha]1\[Omega]rd_] := \[CapitalPhi]Int + \[Delta]\[Beta]1\[Omega]rd/2 + \[Delta]\[Beta]0 - \[CapitalPhi]MR - \[Delta]\[Alpha]1\[Omega]rd/2
}]//.{
	w -> C1[beta0, beta1\[Omega], alpha0, alpha1\[Omega], \[Omega],ringdownFrequency]
};


FunctionDefs["\[Phi]IMR+t0"] = HoldForm[{
	\[Phi]IMR[\[Phi]Ins_, \[Phi]Int_, \[Phi]MR_, c1_, ringdownFrequency_, \[Omega]_] := x,
	\[CapitalDelta]t0[\[Omega]_, d\[Phi]MR_, \[Omega]ref_] := d\[Phi]MR (\[Omega] - \[Omega]ref)
}]//.{
	x-> \[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, c1, ringdownFrequency, \[Omega]]
};


FunctionDefs["Commute D"] = HoldForm[{
	(*distribute derivatives*)
	Unprotect[D],
	DownValues[D] = {},
	D[\[Phi]IMR[x__] - \[Phi]IMR[y__] - \[CapitalDelta]t0[z__], n___] := D[\[Phi]IMR[x],n] - D[\[Phi]IMR[y],n] - D[\[CapitalDelta]t0[z], n],
	Protect[D],

	(*Commute D with \[Phi]IMR*)
	\[Phi]IMR/: D[\[Phi]IMR[x__], n___] := With[
	{args = D[#, n]&/@({x}[[1;;4]])},
	\[Phi]IMR@@(Join[args, {x}[[5;;6]]])
	],
	
	(*Commute D with C1*)
	C1/: D[C1[x__], n___] := With[
		{args = D[#, n]&/@({x}[[1;;4]])},
		C1@@(Join[args, {x}[[5;;6]]])
	],
	(*Commute with \[Delta]\[Beta]0, \[Delta]\[Beta]1, \[Delta]\[Alpha]0*)
	\[Delta]\[Beta]0/: D[\[Delta]\[Beta]0[x__], n___] := \[Delta]\[Beta]0@@(D[#, n]&/@{x}),
	\[Delta]\[Beta]1/: D[\[Delta]\[Beta]1[x__], n___] := \[Delta]\[Beta]1@@(D[#,n]&/@{x}),
	\[Delta]\[Alpha]0/: D[\[Delta]\[Alpha]0[x__], n___] := \[Delta]\[Alpha]0@@(D[#,n]&/@{x})
	
}];


(* ::Text:: *)
(*Make rules to eliminate the U from expressions:*)


(*Make rules to eliminate U*)
AllFunctionHeads = (FunctionDefs//Values)[[All, 1, All, 1,0]]//Flatten//DeleteDuplicates;
AllFunctionHeads = DeleteElements[AllFunctionHeads, {D, DownValues, Symbol}];
LocalVars = varsDefs[[1,All,1]];



Module[{undefined =  ("U" <> #)&/@(ToString/@AllFunctionHeads)},
	
	undefined = ToExpression[undefined,InputForm];
	Headrules = Thread@Rule[undefined, AllFunctionHeads]
];


(* ::Text:: *)
(*Make  the  final  expression :*)


FinalExpr = Module[
	{dummy},
	dummy = HoldForm[
		{phiIMR - phiIMRref - deltat}
	]/.{phiIMR -> \[Phi]IMRexpr, phiIMRref -> \[Phi]IMRexprRef, deltat -> \[CapitalDelta]t0expr};

	dummy//.Headrules
];


(* ::Text:: *)
(*The $Block  structure : *)
(*$Block[*)
(*	{{heads__}, localvars__},*)
(*	localvardefs__;*)
(*	functiondefs__;*)
(*	expr*)
(*] ,*)
(*following that:*)


defs = Module[
	{functionDefs = FunctionDefs//Values, orderedDefs},
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
	{test, {\[Omega], \[Eta], \[Chi]1, \[Chi]2}, 1},
	Evaluate@$blockexpr,
	"KeepDefs" -> KeepDefs
]//.HoldForm[X_] :> X;


<<FelipeBarbosa`SymDALI`


res = EchoTiming[DerivativeRules@@expr];


res[[All,1]]


res[[1,2]]


DGrad[f_, vars_] := Module[
	{h = 1. 10^-6, dummy, Point1, Point2},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = dummy[[i]]+ h;
		dummy,
		{i, Length@vars}
	];

	
	Point2 = ConstantArray[vars, Length[vars]];
	
	
	(f@@@Point1 - f@@@Point2)/h
]


Block[
	{testF, \[Phi]1, \[Phi]2, ND\[Omega],\[Omega], \[Eta], \[Chi]1, \[Chi]2, D\[Omega],
	\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4,
	D\[Eta], D\[Chi]1, D\[Chi]2, \[Omega]ref = 0.0001, M, G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
	nd1, nd2, nd3, nd4, grad, ngrad},
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = ConstantArray[0, 16];
	
	testF[\[Omega]_, \[Eta]_, \[Chi]1_, \[Chi]2_] = res[[1,2]]; DownValues[testF] = DownValues[testF]//. HoldForm[x_]:> x;
	
	D\[Omega] = res[[2,2]]; OwnValues[D\[Omega]] = OwnValues[D\[Omega]]//. HoldForm[x_]:> x;
	D\[Eta] = res[[3,2]]; OwnValues[D\[Eta]] = OwnValues[D\[Eta]]//. HoldForm[x_]:> x;
	D\[Chi]1 = res[[4,2]]; OwnValues[D\[Chi]1] = OwnValues[D\[Chi]1]//. HoldForm[x_]:> x;
	D\[Chi]2 = res[[5,2]]; OwnValues[D\[Chi]2] = OwnValues[D\[Chi]2]//. HoldForm[x_]:> x;
	
	\[Eta] = RandomReal[{0.1, 0.24}];
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2];
	M = RandomReal[{20,100}];
	\[Omega] = 20 M G;
	(*numerical gradient:*)
	ngrad = DGrad[testF, {\[Omega], \[Eta], \[Chi]1, \[Chi]2}];
	grad = {D\[Omega], D\[Eta], D\[Chi]1, D\[Chi]2};
	
	{grad, ngrad}
]//Quiet


(*Compiler`$CCompilerOptions = {
    "ShellCommandFunction" -> Print, 
    "SystemCompileOptions" -> "-march=native -finline-functions -funroll-loops -flto -fPIC -O3 -fno-fast-math -ftree-vectorize -fvect-cost-model=dynamic"
};
*)

(*Compiler`$CCompilerOptions={
	"ShellCommandFunction"->Print, 
	"SystemCompileOptions"->"-march=native -finline-functions -funroll-loops -flto -fPIC -O3"};*)

(*Compiler`$CCompilerOptions =Compiler`$CCompilerOptions={
	"ShellCommandFunction"->Print, 
	"SystemCompileOptions"->" -fPIC -O2"};*)


c = Hold[
	{{\[Omega], _Real,1},\[Eta],\[Chi]1,\[Chi]2,\[Omega]ref, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4},
	Evaluate[res[[1,2]]//N], 
	CompilationTarget->"C",
	RuntimeOptions->"Speed"(*,
	RuntimeAttributes->{Listable},
	Parallelization->True*)
	
]//.{HoldForm[x_]:> x, us\[Theta]-> UnitStep};


c= Compile@@c;


c2 = c;
c2[[7]] = Function[{\[Omega],\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4},
	"Total_Phase"
];


(* ::Section::Closed:: *)
(*Checking against Ripple*)


python = StartExternalSession[
	{

		"System" -> "Python"
    }
];


ExternalEvaluate[python, {
	   "import numpy as np",
            "from ripplegw.waveforms import IMRPhenomD as IMRD",
            "from ripplegw.waveforms import IMRPhenomD_utils as IMRD_utils"
       }]


(* ::Subsection::Closed:: *)
(*Ripple Functions*)


helperCoeffs = ExternalFunction[python, "def Coeffs(theta):
	a = IMRD_utils.get_coeffs(theta)
	return np.array(a)"
];

RippleCoeffs[\[Theta]_] := helperCoeffs[\[Theta]]//Normal


helperTransitionFrequencies = ExternalFunction[python, "def transitionfrequencies(theta, gamma2, gamma3):
	a = IMRD_utils.get_transition_frequencies(theta, gamma2, gamma3)
	return np.array(a)
"];

RippleTransitionFrequencies[\[Theta]_, \[Gamma]2_, \[Gamma]3_] := helperTransitionFrequencies[\[Theta], \[Gamma]2, \[Gamma]3]//Normal


helperInspiralPhase = ExternalFunction[python, "def InspiralPhase(f, theta, coeffs):
	a = np.array(f)
	return np.array(IMRD.get_inspiral_phase(a, theta, coeffs) )"];

RippleInspiralPhase[f_, \[Theta]_, coeffs_] := helperInspiralPhase[f, \[Theta], coeffs]//Normal


helperIntPhase = ExternalFunction[python, "def IntPhase(f, theta, coeffs):
	a = np.array(f)
	return np.array(IMRD.get_IIa_raw_phase(a, theta, coeffs) )"];
	
RippleIntPhase[f_, \[Theta]_, coeffs_] := helperIntPhase[f, \[Theta], coeffs]//Normal


helperMRPhase = ExternalFunction[python, "def MRPhase(f, theta, coeffs, fRD, fDAMP):
	a = np.array(f)
	return np.array(IMRD.get_IIb_raw_phase(a, theta, coeffs, fRD, fDAMP) )"];
	
RippleMRPhase[f_, \[Theta]_, coeffs_, fRD_,fDAMP_] := helperMRPhase[f, \[Theta], coeffs, fRD, fDAMP]//Normal


helperTotalPhase = ExternalFunction[python, "def totalPhase(f, theta, coeffs, transition_frequencies):
	a = np.array(f)
	return np.array(IMRD.Phase(a, theta, coeffs, transition_frequencies) )"];
	
RippleTotalPhase[f_, \[Theta]_, coeffs_, transition_] := helperTotalPhase[f, \[Theta], coeffs, transition]//Normal


helperCoeffs = ExternalFunction[python, "def Coeffs(theta):
	a = IMRD_utils.get_coeffs(theta)
	return np.array(a)"
];

RippleCoeffs[\[Theta]_] := helperCoeffs[\[Theta]]//Normal

helperArg = ExternalFunction[python, "def h0(f, \[Theta]in, \[Theta]extr, coeffs,  fref):
	b = np.array(f)
	a = IMRD._gen_IMRPhenomD(b, \[Theta]in, \[Theta]extr, coeffs,  fref)
	return np.array(a)"
];
(*the function takes the arguments: Mc, eta, chi1, chi2, dist_mpc, tc, phic, inclination*)
Argument[f_, \[Theta]in_,\[Theta]ex_, coeffs_, fref_] := Arg[helperArg[f, \[Theta]in, \[Theta]ex, coeffs, fref]//Normal]


(* ::Subsection::Closed:: *)
(*Comparison*)


Attributes[RelativeDiff] = {Listable};
RelativeDiff[x_, y_] := With[
	{diff = x-y},
	Max[Abs@{diff/x, diff/y}]
]


test := Module[
	{f, M, \[Chi]1, \[Chi]2,pos, m1, m2,\[Eta],tc, \[Phi]c, Ripple, MMA, \[Omega], \[Theta]in, \[Theta]ex, coeffs,  
	G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],diff, transition,
	ringdown,\[Omega]ref,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4},

	f = Range[20, 2048,0.25];
	{m1, m2} = ReverseSort@RandomReal[{10,120},2];
	
	M = (m1+m2);
	
	\[Eta] = (m1 m2)/M^2;
	tc =0;  
	\[Phi]c = 0; 
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2];
	\[Omega] = G M f;
	\[Omega]ref = 20 G M;
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = ConstantArray[0,16];

	\[Theta]ex = {1, tc, \[Phi]c};
	\[Theta]in = {m1, m2, \[Chi]1, \[Chi]2};
	coeffs = RippleCoeffs[\[Theta]in];

	Ripple = Argument[f, \[Theta]in, \[Theta]ex, coeffs, 20];
	Ripple = ResourceFunction["PhaseUnwrap"][Ripple];
	
	MMA = -c2[\[Omega],\[Eta],\[Chi]1,\[Chi]2,\[Omega]ref,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4];
	
	diff = RelativeDiff@@{MMA, Ripple};
	transition = RippleTransitionFrequencies[\[Theta]in,  Sequence@@coeffs[[6;;7]] ];
	ringdown = transition[[-2]] M G;
	
	
	Ripple = Riffle[\[Omega], Ripple]//Partition[#,2]&;
	MMA = Riffle[\[Omega],MMA]//Partition[#,2]&;

	pos = FirstPosition[\[Omega], x_/; x>=0.19]//Last; (*0.2 is the upper cutoff for IMRPhenomD. *)
	diff = Riffle[\[Omega], diff]//Partition[#,2]&;
	{
		ListLinePlot[
			Take[diff, pos], 
			GridLines->{{{0.018,Red}, {ringdown/2, Red}}, None}, 
			PlotRange->All, ImageSize->Medium
		],
		ListLinePlot[
			{Take[Ripple, pos], Take[MMA, pos]}, GridLines->{{{0.018,Red}, {ringdown/2, Red}}, None},PlotRange->All,
			PlotLegends->{"Python", "MMA"}, ImageSize->Medium
		]
	}
	

]


test//Quiet


(* ::Section::Closed:: *)
(*Derivatives and compilation*)


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {\[Omega],\[Eta], \[Chi]1, \[Chi]2(*, \[Omega]ref*)(*, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4*)};


(*All derivatives up to order 3*)
derivatives = Combinations[vars, 5];


Table[Length@Combinations[vars, i], {i, 1 5}]


derivatives//Length


PrependTo[derivatives, {}];


derivatives//Length


Clear@expr;

With[{
	vars = {\[Omega],\[Eta], \[Chi]1, \[Chi]2, \[Omega]ref, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}
	},
	expr =  Hold[
	{\[CapitalPhi]IMR,  vars, {i}},
	Evaluate@$blockexpr,
	"KeepDefs" -> KeepDefs,
	"IncludeZeroDerivative"->False
]//.HoldForm[X_] :> X;
	expr = Hold[Evaluate@expr];
	
	expr = DerivativeRules@@@expr;

]


ParallelTable@@Join[expr, Hold[{i, derivatives[[1;;5]]}, "Method" -> "CoarsestGrained"]]//QuietEcho;


DistributeDefinitions@@(KeepDefs[[All,1]])
DistributeDefinitions[DerivativeRules]
DistributeDefinitions["FelipeBarbosa`SymDALI`*"]


Clear@Ds1to4
Ds4to5 = MemoryConstrained[EchoTiming[DerivativeRules@@expr], 6 10^9];


Export["Phase_Ds_order_5_part1.wdx", Ds4to5]


Clear@expr;

With[{
	vars = {\[Omega],\[Eta], \[Chi]1, \[Chi]2, \[Omega]ref, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}
	},
	expr =  HoldForm[
	{\[CapitalPhi]IMR,  vars, {i}},
	Evaluate@$blockexpr,
	"KeepDefs" -> KeepDefs,
	"IncludeZeroDerivative"->False
]//.HoldForm[X_] :> X;
]


Clear@Ds4to5

expr2 = Hold[Evaluate@expr, Evaluate@{i, derivatives[[101;;126]]}]//.HoldForm-> DerivativeRules;

D5part2 = MemoryConstrained[ Table@@expr2, 5 10^9];

Export["Phase_Ds_order_5_part2.wdx", D5part2]


$D[{x__}, \[CapitalPhi]IMR][y__]/;Total[{x}[[6;;-1]]] === 1 := Module[
	{pos = Position[{x}[[6;;-1]], 1], ds, args},
	ds = Join[
		{x}[[1;;5]],
		ConstantArray[0, 16]
	];
	args = Join[
		{y}[[1;;5]],
		ReplacePart[{y}[[6;;-1]], pos -> 1]
	];
	
	$D[ds, \[CapitalPhi]IMR]@@args
];

$D[{x__}, \[CapitalPhi]IMR][y__]/; Total[{x}[[6;;-1]]] > 1 := 0


(*
	Consider that \[CapitalPsi] = \[CapitalPhi]IMR[\[Omega]] - \[CapitalPhi][\[Omega]ref] - t0 (\[Omega]-\[Omega]ref) and that \[Omega] variables has to be a vector.
	The derivative with respect to \[Omega]ref is minus the derivative of \[Omega] evaluated at the value of \[Omega]ref passed as a vector
*)

$D[{x__}, \[CapitalPhi]IMR][y__]/; {x}[[1]] === 0 && {x}[[5]] > 0 := Module[
	{ds, args},
	ds = Join[
		{x}[[{5}]], 
		{x}[[2;;4]], 
		{0}, 
		{x}[[6;;-1]]
	];
	
	args = Join[
		{y}[[{5}]], 
		{y}[[2;;-1]]
	];
	
	-Last[
		$D[ds, \[CapitalPhi]IMR]@@args
	]
]


$D[{x__}, \[CapitalPhi]IMR][y__]/; ({x}[[1]] > 0 && {x}[[5]] > 0)  := 0


(* ::Subsection:: *)
(*Compiling Phase terms:*)


Module[{ds = Import["Phase_Ds_order_1_to_4.wdx"]}, Ds = ds[[2]]];


(*Set hold pattern on the lfhs of the rule*)
Ds = MapAt[HoldPattern, Ds, {All, 1}];


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis[x_HoldForm] := Module[
	{dummy, vars},
	vars = {{\[Omega],  _Real,  1}, \[Eta], \[Chi]1, \[Chi]2, \[Omega]ref, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4};
	dummy = Hold[
		Evaluate[vars], 
		Evaluate[x],
		RuntimeOptions->{
			"CatchMachineOverflow"->False,
			"CatchMachineIntegerOverflow"->False,
			"EvaluateSymbolically"->False,
			"RuntimeErrorHandler"->None,
			"WarningMessages"->True
		}
	]//.{HoldForm[y_] :> y, us\[Theta]->UnitStep, Rational->Divide, Piecewise[{{val_, cond_}}, else_]-> If[cond, val, else]};
	
	Compile@@dummy
]


<<CompiledFunctionTools`


compiledDs = MapAt[
	compileThis,
	Ds, 
	{All, 2}
];


<<CCompilerDriver`
<<CCodeGenerator`


$CCompilerDefaultDirectory = "/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/";


list = MapIndexed[
	LibraryGenerate[#1[[2]], "phi" <> ToString[#2//First]]&,
	compiledDs
];


(* ::Subsection:: *)
(*Defining phase terms for RosettaStone*)


With[
	{d =FileNames["phi*", {"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/"}] },
	list = SortBy[(StringReplace[FileBaseName[#], "phi"->""]//ToExpression)&]@d
];


phis = Block[
	{Cvariables},
	Cvariables= ConstantArray[{Real, 0, "Constant"},20];
	Cvariables = Join[{{Real,1,"Constant"}}, Cvariables];

MapThread[
	(#1 -> libraryFunction[
		#2, 
		FileBaseName[#2], 
		Cvariables,
		{Real, 1}
	])&,
	{Ds[[2, All,1]],list }
]
];


phis = phis//.libraryFunction[
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/phi6.so",
	"phi6",
	{{Real,1,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"}},
	{Real,1}] -> libraryFunction[
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/phi6.so",
	"phi6",
	{{Real,1,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"},{Real,0,"Constant"}},
	{Real,0}];


phis = phis//.libraryFunction->LibraryFunction;


(* ::Chapter:: *)
(*Amplitude*)


(* ::Section::Closed:: *)
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



(* ::Section::Closed:: *)
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


(* ::Section:: *)
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
	{test, {M, \[Eta], \[Chi]1, \[Chi]2}, 1},
	Evaluate@$blockexpr,
	"KeepDefs" -> KeepDefs
]//.HoldForm[X_] :> X;


<<FelipeBarbosa`SymDALI`


res = EchoTiming[DerivativeRules@@expr]//QuietEcho;


DGrad[f_, vars_] := Module[
	{h = 1. 10^-6, dummy, Point1, Point2},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = dummy[[i]]+ h;
		dummy,
		{i, Length@vars}
	];

	
	Point2 = ConstantArray[vars, Length[vars]];
	
	
	(f@@@Point1 - f@@@Point2)/h
]


Block[
	{testF, \[Phi]1, \[Phi]2, ND\[Omega],\[Omega], \[Eta], \[Chi]1, \[Chi]2, DM,
	D\[Eta], D\[Chi]1, D\[Chi]2, M, G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
    grad, ngrad, f},
    f = RandomReal[{20, 1024}];
	
	testF[M_, \[Eta]_, \[Chi]1_, \[Chi]2_] = res[[1,2]]; DownValues[testF] = DownValues[testF]//. HoldForm[x_]:> x;
	
	DM = res[[2,2]]; OwnValues[DM] = OwnValues[DM]//. HoldForm[x_]:> x;
	D\[Eta] = res[[3,2]]; OwnValues[D\[Eta]] = OwnValues[D\[Eta]]//. HoldForm[x_]:> x;
	D\[Chi]1 = res[[4,2]]; OwnValues[D\[Chi]1] = OwnValues[D\[Chi]1]//. HoldForm[x_]:> x;
	D\[Chi]2 = res[[5,2]]; OwnValues[D\[Chi]2] = OwnValues[D\[Chi]2]//. HoldForm[x_]:> x;
	
	\[Eta] = RandomReal[{0.1, 0.24}];
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2];
	M = RandomReal[{20,100}];
	(*numerical gradient:*)
	ngrad = DGrad[testF, {M, \[Eta], \[Chi]1, \[Chi]2}];
	grad = {DM, D\[Eta], D\[Chi]1, D\[Chi]2};
	
	{grad, ngrad}
]//Quiet


c = Hold[
	{{f, _Real, 1}, M, \[Eta], \[Chi]1,\[Chi]2},
	Evaluate[res[[1,2]]//N], 
	CompilationTarget->"C",
	RuntimeOptions->"Speed"(*,
	RuntimeAttributes->{Listable},
	Parallelization->True*)
	
]//.{HoldForm[x_]:> x, us\[Theta]-> UnitStep};


Compiler`$CCompilerOptions =Compiler`$CCompilerOptions={
	"ShellCommandFunction"->Print, 
	"SystemCompileOptions"->" -fPIC -O2", 
	"CleanIntermediate"->True};


c= Compile@@c;


c2 = c;
c2[[7]] = Function[{f,M,\[Eta],\[Chi]1,\[Chi]2},
	"Total_Amp"
];


Block[
	{expr1,expr2, f=Range[20, 2048, 0.25], \[Eta]=0.2, \[Chi]1, \[Chi]2, M=60},
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1},2];
	
	expr2 = res[[1,2]]; 
	OwnValues[expr2] = OwnValues[expr2]//.HoldForm[x_]:> x;
	
	expr1 = c2[f,M, \[Eta], \[Chi]1, \[Chi]2];
	
	(*{expr1;//AbsoluteTiming, expr2;//AbsoluteTiming}*)
	
	Short[expr1]//AbsoluteTiming
	

]


(* ::Section::Closed:: *)
(*Testing Against Ripple*)


python = StartExternalSession[
	{

		"System" -> "Python"
    }
];


ExternalEvaluate[python, {
	   "import numpy as np",
            "from ripplegw.waveforms import IMRPhenomD as IMRD",
            "from ripplegw.waveforms import IMRPhenomD_utils as IMRD_utils"
       }]


helperCoeffs = ExternalFunction[python, "def Coeffs(theta):
	a = IMRD_utils.get_coeffs(theta)
	return np.array(a)"
];

RippleCoeffs[\[Theta]_] := helperCoeffs[\[Theta]]//Normal

helperArg = ExternalFunction[python, "def h0(f, \[Theta]in, \[Theta]extr, coeffs,  fref):
	b = np.array(f)
	a = IMRD._gen_IMRPhenomD(b, \[Theta]in, \[Theta]extr, coeffs,  fref)
	return np.array(a)"
];


Amplitude[f_, \[Theta]in_,\[Theta]ex_, coeffs_, fref_] := Abs[helperArg[f, \[Theta]in, \[Theta]ex, coeffs, fref]//Normal]


helperTransitionFrequencies = ExternalFunction[python, "def transitionfrequencies(theta, gamma2, gamma3):
	a = IMRD_utils.get_transition_frequencies(theta, gamma2, gamma3)
	return np.array(a)
"];

RippleTransitionFrequencies[\[Theta]_, \[Gamma]2_, \[Gamma]3_] := helperTransitionFrequencies[\[Theta], \[Gamma]2, \[Gamma]3]//Normal


Attributes[RelativeDiff] = {Listable};
RelativeDiff[x_, y_] := With[
	{diff = x-y},
	Max[Abs@{diff/x, diff/y}]
]


test\[ScriptCapitalA] := Module[
	{f, M, \[Chi]1, \[Chi]2,pos, m1, m2,\[Eta], \[Theta], Ripple, MMA, \[Omega], \[Theta]in, \[Theta]ex, coeffs,  ripplefreqs,\[Omega]p,posPeak,diff,
	G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]
	},

	f = Range[20, 2048,0.25];
	{m1, m2} = ReverseSort@RandomReal[{10,120},2];
	M = (m1+m2);
	\[Eta] = (m1 m2)/M^2;
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2];
	\[Omega] = G M f;


	\[Theta]ex = {1, 0, 0};
	\[Theta]in = {m1,m2, \[Chi]1, \[Chi]2};
	coeffs = RippleCoeffs[\[Theta]in];
	ripplefreqs = RippleTransitionFrequencies[\[Theta]in, coeffs[[6]], coeffs[[7]]];
	\[Omega]p = ripplefreqs[[4]] G M;
	

	Ripple = Amplitude[f, \[Theta]in, \[Theta]ex, coeffs,  20];
	MMA = 10^3 c2[f, M, \[Eta], \[Chi]1, \[Chi]2];
	diff = RelativeDiff[Ripple, MMA]//Quiet;

	Ripple = Riffle[\[Omega], Ripple]//Partition[#,2]&;
		(*\[ScriptCapitalA] ~ 1/dL, MMA unities are GPC and dL is set to 1 Mpc*)
	MMA = Riffle[\[Omega],MMA]//Partition[#,2]&;
	diff = Riffle[\[Omega],diff]//Partition[#,2]&;

	pos = FirstPosition[\[Omega], x_/; x>=0.2]//Last; (*0.2 is the upper cutoff for IMRPhenomD. *)
	posPeak = FirstPosition[\[Omega], x_/; x>=\[Omega]p]//Last;

	{
		ListLinePlot[
			Take[diff, pos], 
			GridLines->{{{0.014,Red}, {\[Omega]p, Red}}, None}, 
			PlotRange->All, ImageSize->Medium, Background->White, ScalingFunctions->"Log10"
		],
		
		ListLinePlot[
			{MMA, Ripple}[[All, 1;;pos]],  PlotLegends->{"MMA", "Python"}, ImageSize->Medium,PlotRange->All,
			GridLines -> {{{0.014, Red}, {\[Omega]p, Red}}, None},
			Background->White
		]
	}
]



test\[ScriptCapitalA]


(* ::Section::Closed:: *)
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


PrependTo[derivatives, {}];


Combinations[vars, 3]//Length
Combinations[vars, 4]//Length


derivatives[[36;;70]]


expr =  Hold[
	{\[ScriptA]IMR, {f, M, \[Eta], \[Chi]1, \[Chi]2}, {i}},
	Evaluate@$blockexpr,
	"KeepDefs" -> KeepDefs,
	"IncludeZeroDerivative"->False
]//.HoldForm[X_] :> X;


<<FelipeBarbosa`SymDALI`


Clear@Ds;
expr2 = HoldForm[MemoryConstrained[x, 6 10^9], Evaluate[{i,derivatives[[1;;5]] }]]/.x-> expr;
expr2 = expr2/.Hold-> DerivativeRules;
Ds = Table@@expr2//QuietEcho;
(*Export["Amplitude_Ds_order_4.wdx", Ds]*)


(* ::Subsection:: *)
(*Importing amplitude terms*)


Module[{ds = Import["Amplitude_Ds_order_1_to_3.wdx"]}, Ds = ds[[2]]];


compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{f,_Real,1}, M, \[Eta], \[Chi]1, \[Chi]2}, 
		Evaluate[x],
		RuntimeOptions->{
			"CatchMachineOverflow"->False,
			"CatchMachineIntegerOverflow"->False,
			"EvaluateSymbolically"->False,
			"RuntimeErrorHandler"->None,
			"WarningMessages"->True
		}
	]//.{HoldForm[y_] :> y, us\[Theta]->UnitStep};
	
	Compile@@dummy
]


compileThis[Ds[[-1, 2]]]


<<CompiledFunctionTools`
CompilePrint[]


compiledDs = MapAt[
	compileThis,
	Ds[[2]], 
	{All, 2}
];


<<CCompilerDriver`


$CCompilerDefaultDirectory = "/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/";


Needs["CCodeGenerator`"]


MapIndexed[
	LibraryGenerate[#1[[2]], "a" <> ToString[#2//First]]&,
	compiledDs
]


list\[ScriptCapitalA] = {
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/a1.so",
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/a2.so",
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/a3.so",
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/a4.so",
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/a5.so"
};


\[ScriptCapitalA]s = MapThread[
	(#1 -> LibraryFunction[
		#2, 
		FileBaseName[#2], 
		{{Real, 1, "Constant"}, {Real, 0, "Constant"}, {Real, 0, "Constant"}, {Real, 0, "Constant"}, {Real, 0, "Constant"}},
		{Real, 1}
	])&,
	{Ds\[ScriptCapitalA][[2, All,1]], list\[ScriptCapitalA]}
]


NRules = <||>;


NRules["\[CapitalPhi]IMR"] = phis;
NRules["\[ScriptA]IMR"] = \[ScriptCapitalA]s;


Export["/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/NRules/RosettaStone.wdx", NRules];
