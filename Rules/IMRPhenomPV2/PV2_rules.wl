(* ::Package:: *)

(*

SetOptions[EvaluationNotebook[], WindowElements->{"MemoryMonitor","VerticalScrollBar","MenuBar", "HorizontalScrollBar"}]
SetDirectory[NotebookDirectory[]];
SetOptions[EvaluationNotebook[], DefaultNewCellStyle->"Code"]

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


Quit


(*Basic stuff:*)
$HistoryLength = 1;
PacletDirectoryLoad[ParentDirectory[NotebookDirectory[], 3]];
<<FelipeBarbosa`SymDALI`

vectorDefs = Module[
	{DALIDir, docDir},
	DALIDir = FindFile["FelipeBarbosa`SymDALI`"]//FileNameDrop//ParentDirectory;
	docDir = FileNameJoin[{DALIDir, "Documentation/English/Tutorials"}];
	Import@FileNameJoin[{docDir, "IMRPhenomDComponentsDefinition.wdx"}]
];

SymRules = <||>;
NRules = <||>;


(* ::Chapter::Closed:: *)
(*Phase*)


(* ::Section::Closed:: *)
(*IMRPhenomD C(1) Phase*)


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
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

InsExpr = InsExpr//.{\[Chi]1-> s1z, \[Chi]2->s2z};


InsExpr


(* ::Subsubsection::Closed:: *)
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

IntExpr = IntExpr//.{\[Chi]1->s1z, \[Chi]2->s2z};


(* ::Subsubsection::Closed:: *)
(*Ringdown and Damping -Phase*)


Clear@\[Chi]p

Block[
	{Asp1, Asp2,A1m1sq, q = m1/(m1+m2)},
	Asp1[m1_, m2_, s1x_, s1y_] = m1^2 Sqrt[s1x^2 + s1y^2] (2+ (3 m2)/(2 m1));
	Asp2[m1_, m2_, s2x_, s2y_] = m2^2 Sqrt[s2x^2 + s2y^2] (2+ (3 m1)/(2 m2));
	A1m1sq = 1/2 m1 (4 m1+3 m2);
	
	\[Chi]p[m1_, m2_, s1x_, s1y_, s2x_, s2y_] = If[
		Asp1[m1, m2, s1x, s1y]>=Asp2[m1, m2, s2x, s2y], 
		Evaluate[Asp1[m1, m2, s1x, s1y]/A1m1sq//Simplify],
		Evaluate[Asp2[m1, m2, s2x, s2y]/A1m1sq//Simplify]
	];
	
	Sperp[m1_, m2_, s1x_, s1y_, s2x_, s2y_] = If[
		Asp1[m1, m2, s1x, s1y]>=Asp2[m1, m2, s2x, s2y], 
		Evaluate[q^2 Asp1[m1, m2, s1x, s1y]/A1m1sq//Simplify],
		Evaluate[q^2 Asp2[m1, m2, s2x, s2y]/A1m1sq//Simplify]
	]

]



(*In Ripple this is ```FinalSpin0815s```, S \[Congruent] (m1/M)^2 s1z + (m2/M)^2 s2z*)
Clear@aParallel
aParallel[m1_, m2_, s1z_, s2z_] = Block[
	{S, M = m1+m2, \[Eta]},
	\[Eta] = (m1 m2)/M^2; 
	S=(m1/M)^2 s1z + (m2/M)^2 s2z;
	\[Eta] (
	3.4641016151377544` -
	4.399247300629289` \[Eta]+
	9.397292189321194` \[Eta]^2-
	13.180949901606242` \[Eta]^3+
	S (
		-0.0850917821418767`+
		S (0.1014665242971878` -2.0967746996832157` \[Eta])+
		1.`/\[Eta]-5.837029316602263` \[Eta]+
		S^3 (-0.8676969352555539` +2.064046835273906` \[Eta])+
		S^2 (-1.3546806617824356` +4.108962025369336` \[Eta])
	)
)//Simplify
];


aeff[aParallel_, Sperp_] = If[aParallel>=0,Sqrt[Sperp^2 + aParallel^2],-Sqrt[Sperp^2 + aParallel^2] ];


(*EradRational0815s in ripple, s\[Congruent] (m1^2 s1z + m2^2 s2z)/(m1^2+m2^2)*)
Clear@Erad
Erad[m1_, m2_, s1z_, s2z_] = Block[
	{s = (m1^2 s1z + m2^2 s2z)/(m1^2+m2^2), eta = m1 m2/(m1+m2)^2},
	((0.055974469826360077` eta + 0.5809510763115132` eta^2-0.9606726679372312` eta^3+3.352411249771192` eta^4) (1.` +(-0.0030302335878845507` -2.0066110851351073` eta+7.7050567802399215` eta^2) s))/(1.` +(-0.6714403054720589` -1.4756929437702908` eta+7.304676214885011` eta^2) s)//Simplify
];


(*re\[Omega] -> interpolation for ringdown and the other is for damping: I think this is frm the PhenomHM paper*)
re\[Omega][\[Chi]_] :=(0.05947169566573468` -0.14989771215394762` \[Chi]+0.09535606290986028` \[Chi]^2+0.02260924869042963` \[Chi]^3-0.02501704155363241` \[Chi]^4-0.005852438240997211` \[Chi]^5+0.0027489038393367993` \[Chi]^6+0.0005821983163192694` \[Chi]^7)/(1-2.8570126619966296` \[Chi]+2.373335413978394` \[Chi]^2-0.6036964688511505` \[Chi]^4+0.0873798215084077` \[Chi]^6);
im\[Omega][\[Chi]_] :=(0.014158792290965177` -0.036989395871554566` \[Chi]+0.026822526296575368` \[Chi]^2+0.0008490933750566702` \[Chi]^3-0.004843996907020524` \[Chi]^4-0.00014745235759327472` \[Chi]^5+0.0001504546201236794` \[Chi]^6)/(1-2.5900842798681376` \[Chi]+1.8952576220623967` \[Chi]^2-0.31416610693042507` \[Chi]^4+0.009002719412204133` \[Chi]^6);

\[Omega]RdDamping[int_, Erad_] = int/(1 - Erad);


\[Gamma]2[\[Eta]_, \[Chi]PN_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[6]];
\[Gamma]3[\[Eta]_, \[Chi]PN_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[7]];


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


(* ::Subsubsection::Closed:: *)
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
	{ringdownFrequency, dampingFrequency, \[Alpha]5, aeff},
	
	
	\[Alpha]5 = U\[Alpha]5[\[Eta], \[Chi]1, \[Chi]2];
	
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	 
	  
	{U\[CapitalAlpha]1[\[Eta],\[Chi]1,\[Chi]2],U\[CapitalAlpha]2[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]2],U\[CapitalAlpha]3[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]3],U\[CapitalAlpha]4[\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[Alpha]4]} . \[Omega]MRVecPhase[
		\[Omega],
		ringdownFrequency, 
		dampingFrequency, 
		\[Alpha]5
	]
]//Simplify;


MRExpr = MRExpr//.{\[Chi]1->s1z, \[Chi]2->s2z};


(* ::Subsubsection::Closed:: *)
(*C(1)*)


C1[beta0_, beta1\[Omega]_, alpha0_, alpha1\[Omega]_, \[Omega]_, ringdownFrequency_] = (
			(beta0 + beta1\[Omega])*us\[Theta][(ringdownFrequency/2 - \[Omega]) (\[Omega]-0.018)] + 
			(alpha0 +alpha1\[Omega])us\[Theta][(\[Omega] - ringdownFrequency/2) (0.2-\[Omega])]
);


Module[{
	ringdownFrequency, dampingFrequency, \[Alpha]5, aeff
	},
	\[Alpha]5 = U\[Alpha]5[\[Eta], \[Chi]1, \[Chi]2];
	
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][aeff], UErad[m1, m2, s1z, s2z]];

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


exprC1 = exprC1//.{\[Chi]1->s1z, \[Chi]2->s2z};


(* ::Subsubsection::Closed:: *)
(*\[Phi]IMR*)


\[Phi]IMR[\[Phi]Ins_, \[Phi]Int_, \[Phi]MR_, c1_, ringdownFrequency_, \[Omega]_] := (
	\[Phi]Ins*us\[Theta][-\[Omega]+0.018] +
	\[Phi]Int*us\[Theta][(ringdownFrequency/2 - \[Omega]) (\[Omega]-0.018)]  +
	\[Phi]MR*us\[Theta][(\[Omega] - ringdownFrequency/2) (0.2-\[Omega])] + 
	c1
)


Module[
	{ringdownFrequency, aeff,
	\[Phi]Ins, \[Phi]Int, \[Phi]MR},
	
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	\[Phi]Ins = InsExpr;
	
	\[Phi]Int = IntExpr;
	
	\[Phi]MR = MRExpr;
	
	\[Phi]IMRexpr = U\[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, exprC1, ringdownFrequency, \[Omega]];
	\[Phi]IMRexprRef = U\[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, exprC1, ringdownFrequency, \[Omega]]//.\[Omega]->\[Omega]ref;
]//Simplify;


(* ::Subsubsection::Closed:: *)
(*t0*)


\[CapitalDelta]t0[\[Omega]timesd\[Phi]MR_] = \[Omega]timesd\[Phi]MR;

(*Module[
	{d\[Phi]MR,
	ringdownFrequency, aeff},
	
	
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	d\[Phi]MR =D[MRExpr, \[Omega]]//.\[Omega]->ringdownFrequency;
		
	\[CapitalDelta]t0expr = U\[CapitalDelta]t0[\[Omega], d\[Phi]MR, \[Omega]ref]
];*)


Module[{
	ringdownFrequency, dampingFrequency, \[Alpha]5, aeff
	},
	\[Alpha]5 = U\[Alpha]5[\[Eta], \[Chi]1, \[Chi]2];
	
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][aeff], UErad[m1, m2, s1z, s2z]];

	\[Delta]\[Beta]1[d\[Phi]Ins_, d\[Phi]Int_] = d\[Phi]Ins - d\[Phi]Int;
	\[Delta]\[Beta]0[\[CapitalPhi]Ins_, \[CapitalPhi]Int_, \[Delta]\[Beta]1_] = \[CapitalPhi]Ins - \[CapitalPhi]Int - \[Delta]\[Beta]1 0.018;

	\[Delta]\[Alpha]1[d\[CapitalPhi]Int_, d\[CapitalPhi]MR_, \[Delta]\[Beta]1_] = d\[CapitalPhi]Int + \[Delta]\[Beta]1 - d\[CapitalPhi]MR;
	
	\[Delta]\[Alpha]0[\[CapitalPhi]Int_, \[CapitalPhi]MR_, \[Delta]\[Beta]0_, \[Delta]\[Beta]1\[Omega]rd_, \[Delta]\[Alpha]1\[Omega]rd_] := \[CapitalPhi]Int + \[Delta]\[Beta]1\[Omega]rd/2 + \[Delta]\[Beta]0 - \[CapitalPhi]MR - \[Delta]\[Alpha]1\[Omega]rd/2;
	
	Block[
		{dphiIns, dphiInt1,dphiInt2,dphiMR, phiIns, phiInt1, phiInt2,phiMR, beta0, beta1,  alpha1, d\[Phi]MR},
		
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
		
		
		(*########################### This is Specifically for t0:*)
		
		d\[Phi]MR =D[MRExpr, \[Omega]]//.\[Omega]->ringdownFrequency;
		
		\[CapitalDelta]t0expr = U\[CapitalDelta]t0[\[Omega] (d\[Phi]MR + alpha1)]
		
	]
];


(* ::Text:: *)
(*Pos[expri, hj] ={POS1, POS2, ...},  j = 1, ..., N*)


(* ::Section::Closed:: *)
(*\[Phi]_JSF \[And] \[Epsilon]*)


(* ::Subsubsection::Closed:: *)
(*\[Phi]_JSF*)


\[Phi]Jsf[m1_,m2_, s1x_, s1y_, s2x_, s2y_] = ArcTan[
	m1^2 s1x + m2^2 s2x,
	m1^2 s1y + m2^2 s2y
];


(* ::Subsubsection::Closed:: *)
(*\[Epsilon] PPN coefficients*)


angcoeffs = <|
    "epsiloncoeff1" -> -0.18229166666666666 - (5 dm)/(64.0 m2),
    
    "epsiloncoeff2" -> (-15 dm m2 chil)/(128.0 mtot2 eta) - (35 m2^2 chil)/(128.0 mtot2 eta),
    
    "epsiloncoeff3" -> -1.7952473958333333 - (4555 dm)/(7168.0 m2) - (515 eta)/384.0 - (15 dm^2 eta)/(256.0 m2^2) - (175 dm eta)/(256.0 m2),
    
    "epsiloncoeff4" -> -(35 Pi)/48.0 - (5 dm Pi)/(16.0 m2) + (5 dm^2 chil)/(16.0 mtot2) + (5 dm m2 chil)/(3.0 mtot2) 
        + (2545 m2^2 chil)/(1152.0 mtot2) + (2035 dm m2 chil)/(21504.0 mtot2 eta) + (2995 m2^2 chil)/(9216.0 mtot2 eta),
    
    "epsiloncoeff5" -> 4.318908476114694 + (27895885 dm)/(2.1676032*10^7 m2) + (39695 eta)/86016.0 
        + (1615 dm^2 eta)/(28672.0 m2^2) - (265 dm eta)/(14336.0 m2) + (955 eta2)/576.0 + (15 dm^3 eta2)/(1024.0 m2^3) 
        + (35 dm^2 eta2)/(256.0 m2^2) + (2725 dm eta2)/(3072.0 m2) - (15 dm m2 Pi chil)/(16.0 mtot2 eta) 
        - (35 m2^2 Pi chil)/(16.0 mtot2 eta) + (375 dm^2 m2^2 chil2)/(256.0 mtot4 eta) + (1815 dm m2^3 chil2)/(256.0 mtot4 eta) 
        + (1645 m2^4 chil2)/(192.0 mtot4 eta)
|>;


Block[
	{m2,m1,dm,mtot,eta,eta2,eta3,eta4,mtot2,mtot4,mtot6,mtot8,chil2,chip2,chip4,dm2,dm3, 
	rules, association},
	
	rules  = { (*All defs Ripple does in IMRPhenomPv2_utils.py*)
		m2->q/(1.` +q),m1->1.`/(1.` +q),dm->m1-m2,mtot->1.`,eta->m1 m2,eta2->eta eta,
		eta3->eta2*eta,eta4->eta3*eta,mtot2->mtot*mtot,mtot4->mtot2*mtot2,
		mtot6->mtot4*mtot2,mtot8->mtot6*mtot2,chil2-> chil*chil,chip2->chip*chip,
		chip4->chip2*chip2,dm2->dm*dm,dm3->dm2*dm
	};
	
	association = Simplify/@(angcoeffs//.rules);
	
	(*When they call this function they don't pass q<=1 but rather q>=1, So lets change q->1/q to have
	the usual variable:*)
	
	association = association//.q->qp^-1;
	association = Simplify/@(association//.qp->q);
	NewDictionary = association;
]


Clear[{"Global`\[CapitalEpsilon]*"}]


(*Define the functions now.
All functions of chil, q, chip, with q<=1
*)

\[CapitalEpsilon]1[q_] = NewDictionary["epsiloncoeff1"];
\[CapitalEpsilon]2[q_, chil_] = NewDictionary["epsiloncoeff2"];
\[CapitalEpsilon]3[q_] = NewDictionary["epsiloncoeff3"];
\[CapitalEpsilon]4[q_, chil_] = NewDictionary["epsiloncoeff4"];
\[CapitalEpsilon]5[q_, chil_] = NewDictionary["epsiloncoeff5"];


(* ::Subsubsection::Closed:: *)
(*\[Chi]l*)


(*
chil in gen_IMRPhenomPv2 in IMRPhenomPv2.py
*)
Module[
	{qf = m2/m1},
	((1+q)/q (m1 s1z + m2 s2z)/(m1+m2)//.q->qp^-1)//.qp->qf (*just so we work with q <= 1*)
]//Simplify


Clear@\[Chi]l
\[Chi]l[m1_, m2_, s1z_, s2z_]= s1z + m2/m1 s2z; (*q<=1*)


(* ::Subsubsection::Closed:: *)
(*\[Epsilon] function*)


\[Epsilon][term1_, term2_, term3_, term4_, term5_] := (
	term1+term2+term3+term4+term5
)


Module[{\[Chi]l = U\[Chi]l[m1, m2, s1z, s2z]}, 
	
	\[Epsilon]Expression = U\[Epsilon][
		U\[CapitalEpsilon]1[q] (\[Omega]^-1-\[Omega]ref^-1),
		U\[CapitalEpsilon]2[q, \[Chi]l] (\[Omega]^(-2/3)-\[Omega]ref^(-2/3)),
		U\[CapitalEpsilon]3[q] (\[Omega]^(-1/3) - \[Omega]ref^(-1/3)),
		U\[CapitalEpsilon]4[q, \[Chi]l] Log[f/fref],
		U\[CapitalEpsilon]5[q, \[Chi]l] (\[Omega]^(1/3) - \[Omega]ref^(1/3))
	]//.{\[Omega]->\[Omega]p, \[Omega]ref->\[Omega]pref};
]


Block[
	{
	chil,q = m2/m1, \[Epsilon]Expression, \[Omega]p = (m1+m2) UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]] \[Pi] f, 
	\[Omega]pref =  (m1+m2) UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]] \[Pi] fref
	}, 
	chil = \[Chi]l[m1, m2, s1z, s2z];
	\[Epsilon]Expression = \[Epsilon][
		\[CapitalEpsilon]1[q] (\[Omega]^-1-\[Omega]ref^-1),
		\[CapitalEpsilon]2[q, chil] (\[Omega]^(-2/3)-\[Omega]ref^(-2/3)),
		\[CapitalEpsilon]3[q] (\[Omega]^(-1/3) - \[Omega]ref^(-1/3)),
		\[CapitalEpsilon]4[q, chil] Log[\[Omega]/\[Omega]ref],
		\[CapitalEpsilon]5[q, chil] (\[Omega]^(1/3) - \[Omega]ref^(1/3))
	]//.{\[Omega]->\[Omega]p, \[Omega]ref->\[Omega]pref};
	
	ManualGrad\[Epsilon][f_, fref_, m1_, m2_, s1z_, s2z_] = 2 D[\[Epsilon]Expression, {{m1, m2}}];
]


Module[
	{f, fref, m1, m2, s1z, s2z, G=UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], derivativeTools, manual},
	
	{m1, m2} = ReverseSort@RandomReal[{20, 100}, 2];
	
	f = RandomReal[{1, 0.2/(G (m1+m2))}];
	
	fref = RandomReal[{1, f}];
	
	{s1z, s2z} = RandomReal[{-1,1},2];
	
	derivativeTools = TestGrad\[CapitalPsi][f, fref, m1, m2, 0, 0,s1z, 0, 0, s2z][[1;;2]];
	manual = ManualGrad\[Epsilon][f, fref, m1, m2, s1z, s2z];
	
	(*Echo[{"DT:", derivativeTools}];
	Echo[{"M:", manual}];*)
	RelativeDiff@@{derivativeTools, manual}
	
	
]


(* ::Section::Closed:: *)
(*Making Block function*)


varsDefs = HoldForm[{
	\[Omega] = G (m1+m2) f,
	\[Omega]p =G (m1+m2) \[Pi] f,
	\[Omega]ref = G (m1+m2) fref,
	\[Omega]pref = G (m1+m2) \[Pi] fref,
	\[Eta] = m2 m1/(m2+m1)^2,
	q= m2/m1,
	\[Delta] = Sqrt[1 - 4 \[Eta]],
	\[Chi]s = (s1z+s2z)/2,
	\[Chi]a = (s1z-s2z)/2,
	\[Chi]PN = Sqrt[1 - 4 \[Eta]] (s1z-s2z)/2  + (1-76 \[Eta]/113) (s1z+s2z)/2,
	S = 1/4 (1+Sqrt[1 - 4 \[Eta]])^2 s1z + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 s2z,
	Shat = (1/4 (1+Sqrt[1 - 4 \[Eta]])^2 s1z + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 s2z)/(1- 2 \[Eta])
}];


(* ::Text:: *)
(*Now, all function definitions, starting at Primitives and ending on \[Phi]IMR :*)


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
	\[CapitalPhi]minus2[\[Eta]_, \[Delta]\[CurlyPhi]minus2_]:=-2,
	\[CapitalPhi]0[\[Eta]_,\[Delta]\[CurlyPhi]0_] :=0,
	\[CapitalPhi]1[\[Eta]_,\[Delta]\[CurlyPhi]1_]:=1,
	\[CapitalPhi]2[\[Eta]_,\[Delta]\[CurlyPhi]2_]:=2,
	\[CapitalPhi]3[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]3_]:=3,
	\[CapitalPhi]4[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]4_]:=4,
	\[CapitalPhi]5[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]5_]:=5,
	\[CapitalPhi]5l[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]5l_]:=5l,
	\[CapitalPhi]6[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]6_]:=6,
	\[CapitalPhi]6l[\[Eta]_,\[Delta]\[CurlyPhi]6l_]:=6l,
	\[CapitalPhi]7[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[CurlyPhi]7_]:=7,
	\[CapitalSigma]1[\[Eta]_, \[Chi]1_, \[Chi]2_] := A, 
	\[CapitalSigma]2[\[Eta]_, \[Chi]1_, \[Chi]2_] := B, 
	\[CapitalSigma]3[\[Eta]_, \[Chi]1_, \[Chi]2_] := c, 
	\[CapitalSigma]4[\[Eta]_, \[Chi]1_, \[Chi]2_] := d,
	\[CapitalBeta]1[\[Eta]_, \[Chi]1_, \[Chi]2_] := x,
	\[CapitalBeta]2[\[Eta]_,\[Chi]1_, \[Chi]2_, \[Delta]\[Beta]2_] := y,
	\[CapitalBeta]3[\[Eta]_, \[Chi]1_, \[Chi]2_, \[Delta]\[Beta]3_] :=z
	
}]/.Join[
	rules,
	{x-> \[CapitalBeta]1[\[Eta],\[Chi]1, \[Chi]2], y-> \[CapitalBeta]2[\[Eta],\[Chi]1, \[Chi]2, \[Delta]\[Beta]2],z-> \[CapitalBeta]3[\[Eta],\[Chi]1, \[Chi]2, \[Delta]\[Beta]3]},
	{A ->\[CapitalSigma]1[\[Eta], \[Chi]1, \[Chi]2], B -> \[CapitalSigma]2[\[Eta], \[Chi]1, \[Chi]2], c -> \[CapitalSigma]3[\[Eta], \[Chi]1, \[Chi]2], d -> \[CapitalSigma]4[\[Eta], \[Chi]1, \[Chi]2] }
]

];


FunctionDefs["MR"] = HoldForm[{
	aParallel[m1_, m2_, s1z_, s2z_] := x,
	
	aeff[aParallel_, Sperp_] := x0,
	Erad[m1_, m2_, s1z_, s2z_] := y,
	
	Sperp[m1_, m2_, s1x_, s1y_, s2x_, s2y_] := y0,
	re\[Omega][\[Chi]_] := z,
	im\[Omega][\[Chi]_] := w,
	\[Omega]RdDamping[int_, Erad_] := \[ScriptX],
	\[Gamma]2[\[Eta]_, \[Chi]PN_] := \[ScriptY],
	\[Gamma]3[\[Eta]_, \[Chi]PN_] := \[ScriptZ],
	
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] := \[ScriptW],
	\[Alpha]5[\[Eta]_, \[Chi]1_, \[Chi]2_] := \[ScriptA],
	
	\[CapitalAlpha]1[\[Eta]_,\[Chi]1_,\[Chi]2_] := 1,
	\[CapitalAlpha]2[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[Alpha]2_] :=2 ,
	\[CapitalAlpha]3[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[Alpha]3_] := 3,
	\[CapitalAlpha]4[\[Eta]_,\[Chi]1_,\[Chi]2_,\[Delta]\[Alpha]4_] := 4
	
}]/.{
	x -> aParallel[m1, m2, s1z, s2z], x0-> aeff[aParallel, Sperp],
	y-> Erad[m1, m2, s1z, s2z], 
	y0->Sperp[m1, m2, s1x, s1y, s2x, s2y],
	z-> re\[Omega][\[Chi]], w-> im\[Omega][\[Chi]],
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
	\[CapitalDelta]t0[\[Omega]timesd\[Phi]MR_] := \[Omega]timesd\[Phi]MR
}]//.{
	x-> \[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, c1, ringdownFrequency, \[Omega]]
};


FunctionDefs["PV2_exclusive"] = HoldForm[{
	\[CapitalEpsilon]1[q_] := x1,
	\[CapitalEpsilon]2[q_, chil_] := x2,
	\[CapitalEpsilon]3[q_] := x3,
	\[CapitalEpsilon]4[q_, chil_] := x4,
	\[CapitalEpsilon]5[q_, chil_] := x5,
	
	\[Phi]Jsf[m1_,m2_, s1x_, s1y_, s2x_, s2y_] := ArcTan[m1^2 s1x + m2^2 s2x,m1^2 s1y + m2^2 s2y],
	\[Chi]l[m1_, m2_, s1z_, s2z_] := s1z + m2/m1 s2z,
	\[Epsilon][term1_, term2_, term3_, term4_, term5_] := term1+term2+term3+term4+term5
}]/.{x1 ->\[CapitalEpsilon]1[q], x2->\[CapitalEpsilon]2[q, chil], x3->\[CapitalEpsilon]3[q], x4->\[CapitalEpsilon]4[q, chil], x5->\[CapitalEpsilon]5[q, chil]};


FunctionDefs["Commute D"] = HoldForm[{
	(*distribute derivatives*)
	Unprotect[D],
	DownValues[D] = {},
	D[\[Phi]IMR[x__] - \[CapitalDelta]t0[z__] + 2 \[Phi]Jsf[y__] + 2\[Epsilon][w__], n___] := D[\[Phi]IMR[x],n] - D[\[CapitalDelta]t0[z], n] + 2 D[\[Phi]Jsf[y], n] + 2 D[\[Epsilon][w], n],
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
	\[Delta]\[Alpha]0/: D[\[Delta]\[Alpha]0[x__], n___] := \[Delta]\[Alpha]0@@(D[#,n]&/@{x}),
	
	(*commute with \[Epsilon]:*)
	\[Epsilon]/: D[\[Epsilon][x__], n__] := \[Epsilon]@@(D[#, n]&/@{x})
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


FinalExpr = Module[
	{dummy},
	dummy = HoldForm[
		{phiIMR  - deltat + 2 phiJSF + 2 epsilon}
	]/.{phiIMR -> \[Phi]IMRexpr, deltat -> \[CapitalDelta]t0expr, phiJSF-> U\[Phi]Jsf[m1,m2, s1x, s1y, s2x, s2y],epsilon-> \[Epsilon]Expression };

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


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {m1, m2, s1x, s1y, s1z, s2x, s2y, s2z};


(*All derivatives up to order 3*)
derivatives = Combinations[vars, 1];
PrependTo[derivatives, {}];


With[{
	vars = {
		f, fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
	}
	},
	expr =  Hold[
		{test, vars, derivatives},
		Evaluate@$blockexpr,
		"KeepDefs" -> KeepDefs,
		"IncludeZeroDerivative"->False
	]//.HoldForm[X_] :> X;
]


<<FelipeBarbosa`SymDALI`


res = MemoryConstrained[
	DerivativeRules@@expr,
	8 10^9];


Block[
	{testF,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,f, fref,
	\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4,
	 G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]},
	
	testF[f_, fref_, m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_] = res[[1,2]]; 
	DownValues[testF] = DownValues[testF]//. HoldForm[x_]:> x;
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = ConstantArray[0, 16];
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100}, 2]];
	{s1x,s1y,s1z}  = RandomReal[{-1,1},3];
	{s2x,s2y,s2z}  = RandomReal[{-1,1},3];
	f = 50;
	fref=10;

	
	testF[f,fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z]
]


SetDirectory[NotebookDirectory[]]


Export["Phase_Ds_order_0_to_1.wdx", res]


(* ::Section::Closed:: *)
(*Testing the function*)


Clear@RelativeDiff

Attributes[RelativeDiff] = {Listable};
RelativeDiff[x_,y_]/;x==0 &&y==0 := 0
RelativeDiff[0, y_]/; y!=0 := 1
RelativeDiff[x_, 0]/; x!=0 := 1

RelativeDiff[x_, y_]/; x!=0 &&y!=0 := With[
	{diff = x-y},
	Max[Abs@{diff/x, diff/y}]
]


SetDirectory[NotebookDirectory[]];
res = Import["Phase_Ds_order_0_to_1.wdx"];


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	rule = MapThread[
		Rule,
		{{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}, ConstantArray[0,16]}
	];
	ClearAll[Test\[CapitalPsi]];
	Test\[CapitalPsi][f_,fref_, m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_] = res[[1,2]]//.Join[rule, {G -> g}];

]
DownValues[Test\[CapitalPsi]] = DownValues[Test\[CapitalPsi]]//.HoldForm[x_]:> x;


(* ::Subsection::Closed:: *)
(*Testing the Phase against Ripple*)


DeleteObject/@ExternalSessions[]

Clear@python

python = StartExternalSession["Python"];
ExternalEvaluate[python,"
import numpy as np

import jax
import jax.numpy as jnp

from jax import grad, vmap
from functools import partial

from ripplegw.waveforms import IMRPhenomPv2
from ripplegw import get_match_arr, get_eff_pads
from ripplegw import ms_to_Mc_eta
from ripplegw.constants import MSUN, gt
"]


R\[CapitalPsi] = ExternalFunction[python,"
def Ripple_hp(f, f_ref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dL, tc,iota, phi_ref):
    f_array = jnp.array(f)

    Mc = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    eta = m1 * m2 / ((m1 + m2)**2)

    theta = jnp.array([Mc, eta, s1x, s1y, s1z, s2x, s2y, s2z, dL, tc, phi_ref, iota])

    hp, hc, t2m, tm2m, zeta, epsilon, phi_Jsf, t0, a, psi, alpha = IMRPhenomPv2.gen_IMRPhenomPv2_hphc(f_array, theta, f_ref)
    
    result = psi + 2*phi_Jsf - (-2*np.pi*t0)*f_array + 2*epsilon
    
    return result.tolist()


"]


Clear@Test

Test := Module[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA,Ripplehp, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	
	{s1x, s1y, s1z} = RandomReal[{-1,1}, 3];
	{s2x, s2y, s2z} = RandomReal[{-1,1}, 3];
	\[Chi]s = (s1z+s2z)/2; \[Chi]a = (s1z-s2z)/2;
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	
	MMA = Test\[CapitalPsi][f, 10, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z ];
	
	Ripplehp = R\[CapitalPsi][f,10, m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,1,0,\[Iota],\[Phi]Ref];
	
	diff = RelativeDiff@@{MMA, Ripplehp};
	
	{ListLinePlot[
		{MMA, Ripplehp},
		Frame->True,
		PlotLegends->{"MMA", "Ripple"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.018/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	],
	ListLinePlot[
		diff,
		Frame->True,
		PlotLegends->{"Relative Difference"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.018/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	]}
	
]


Test


(* ::Subsection::Closed:: *)
(*Comparing numerical x Symbolic derivatives*)


NGrad//Clear
NGrad[f_, vars_, n_] := Module[
	{h = 1. 10^-6, dummy, Point1, Point2},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = vars[[i]] + h vars[[i]];
		dummy,
		{i, n+1, Length@vars}
	];
	
	
	Point2 = ConstantArray[vars, (Length[vars] - n)];
	
	
	(f@@@Point1 - f@@@Point2)/(h vars[[n+1;;-1]])
]


{SymRules, NRules} = DerivativeRulesLoad["IMRPhenomPv2"];



(*GRAD WITH COMPILED FUNCTIONS:*)
Block[{\[Delta]s, args}, 
	
	\[Delta]s = ConstantArray[0, 16];
	args = Join[{f,fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z}, \[Delta]s];
	
	ClearAll[TestGrad\[CapitalPsi]];
	
	TestGrad\[CapitalPsi][f_,fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_] = Table[
	(NRules["\[CapitalPhi]IMR"])[[i, 3]]@@args,
	{i, 2, 9}
	]

]
DownValues[TestGrad\[CapitalPsi]] = DownValues[TestGrad\[CapitalPsi]]//.HoldForm[x_]:> x;


(*BELOW IS USEFUL FOR THE BLOCK FUNCTIONS*)



(*Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	rule = MapThread[
		Rule,
		{{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}, ConstantArray[0,16]}
	];
	
	ClearAll[TestGrad\[CapitalPsi]];
	
	TestGrad\[CapitalPsi][f_,fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_] = res[[2;;-1, 2]]//.Join[rule, {G -> g}];

]
DownValues[TestGrad\[CapitalPsi]] = DownValues[TestGrad\[CapitalPsi]]//.HoldForm[x_]:> x;*)


Clear@Test
Test := Module[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, f, Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	
	{s1x, s1y, s1z} = RandomReal[{-1,1}, 3];
	{s2x, s2y, s2z} = RandomReal[{-1,1}, 3];
	
	f = RandomReal[{10., 0.2/(G (m1+m2))}]//List;
	fref = RandomReal[10];
	
	Symbolic = Last/@TestGrad\[CapitalPsi][f, fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z ];
	vars = {f[[1]], fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z};
	Numeric = NGrad[Test\[CapitalPsi], vars, 2];
	
	RelativeDiff@@{Symbolic, Numeric}

	
	

	
]


(*Table[Test//Round, {30}]//MatrixForm*)
Table[Round@Test, {20}]//MatrixForm


(* ::Section::Closed:: *)
(*Compiling*)


Module[{ds = Import["Phase_Ds_order_0_to_1.wdx"]}, Ds = ds//.G->UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]];


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis[x_HoldForm] := Module[
	{dummy, vars},
	vars = {
		{f,  _Real,  1}, fref, 
		m1,m2,s1x,s1y,s1z,s2x,s2y,s2z, 
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
	};
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


CompilePrint[compiledDs[[-1,2]]]


<<CCompilerDriver`
<<CCodeGenerator`


ParentDirectory[NotebookDirectory[], 2]


$CCompilerDefaultDirectory = FileNameJoin[{
	ParentDirectory[NotebookDirectory[], 2],
	"/LibraryResources/",
	$SystemID,
	"/DerivativeRules/IMRPhenomPv2/NRules"
}]


list = MapIndexed[
	LibraryGenerate[#1[[2]], "phi" <> ToString[#2//First]]&,
	compiledDs
];


(* ::Chapter::Closed:: *)
(*Amplitude*)


(* ::Section::Closed:: *)
(*IMRPhenomD C(1)*)


(* ::Subsection::Closed:: *)
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



(* ::Subsection::Closed:: *)
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


InsAmpExpr = InsAmpExpr//.{\[Chi]1->s1z, \[Chi]2->s2z};


(* ::Subsection::Closed:: *)
(*Ringdown and Damping -Phase*)


Clear@\[Chi]p

Block[
	{Asp1, Asp2,A1m1sq, q = m1/(m1+m2)},
	Asp1[m1_, m2_, s1x_, s1y_] = m1^2 Sqrt[s1x^2 + s1y^2] (2+ (3 m2)/(2 m1));
	Asp2[m1_, m2_, s2x_, s2y_] = m2^2 Sqrt[s2x^2 + s2y^2] (2+ (3 m1)/(2 m2));
	A1m1sq = 1/2 m1 (4 m1+3 m2);
	
	\[Chi]p[m1_, m2_, s1x_, s1y_, s2x_, s2y_] = If[
		Asp1[m1, m2, s1x, s1y]>=Asp2[m1, m2, s2x, s2y], 
		Evaluate[Asp1[m1, m2, s1x, s1y]/A1m1sq//Simplify],
		Evaluate[Asp2[m1, m2, s2x, s2y]/A1m1sq//Simplify]
	];
	
	Sperp[m1_, m2_, s1x_, s1y_, s2x_, s2y_] = If[
		Asp1[m1, m2, s1x, s1y]>=Asp2[m1, m2, s2x, s2y], 
		Evaluate[q^2 Asp1[m1, m2, s1x, s1y]/A1m1sq//Simplify],
		Evaluate[q^2 Asp2[m1, m2, s2x, s2y]/A1m1sq//Simplify]
	]

]



(*In Ripple this is ```FinalSpin0815s```, S \[Congruent] (m1/M)^2 s1z + (m2/M)^2 s2z*)
Clear@aParallel
aParallel[m1_, m2_, s1z_, s2z_] = Block[
	{S, M = m1+m2, \[Eta]},
	\[Eta] = (m1 m2)/M^2; 
	S=(m1/M)^2 s1z + (m2/M)^2 s2z;
	\[Eta] (
	3.4641016151377544` -
	4.399247300629289` \[Eta]+
	9.397292189321194` \[Eta]^2-
	13.180949901606242` \[Eta]^3+
	S (
		-0.0850917821418767`+
		S (0.1014665242971878` -2.0967746996832157` \[Eta])+
		1.`/\[Eta]-5.837029316602263` \[Eta]+
		S^3 (-0.8676969352555539` +2.064046835273906` \[Eta])+
		S^2 (-1.3546806617824356` +4.108962025369336` \[Eta])
	)
)//Simplify
];


aeff[aParallel_, Sperp_] = If[aParallel>=0,Sqrt[Sperp^2 + aParallel^2],-Sqrt[Sperp^2 + aParallel^2] ];


(*EradRational0815s in ripple, s\[Congruent] (m1^2 s1z + m2^2 s2z)/(m1^2+m2^2)*)
Clear@Erad
Erad[m1_, m2_, s1z_, s2z_] = Block[
	{s = (m1^2 s1z + m2^2 s2z)/(m1^2+m2^2), eta = m1 m2/(m1+m2)^2},
	((0.055974469826360077` eta + 0.5809510763115132` eta^2-0.9606726679372312` eta^3+3.352411249771192` eta^4) (1.` +(-0.0030302335878845507` -2.0066110851351073` eta+7.7050567802399215` eta^2) s))/(1.` +(-0.6714403054720589` -1.4756929437702908` eta+7.304676214885011` eta^2) s)//Simplify
];


(*re\[Omega] -> interpolation for ringdown and the other is for damping: I think this is frm the PhenomHM paper*)
re\[Omega][\[Chi]_] :=(0.05947169566573468` -0.14989771215394762` \[Chi]+0.09535606290986028` \[Chi]^2+0.02260924869042963` \[Chi]^3-0.02501704155363241` \[Chi]^4-0.005852438240997211` \[Chi]^5+0.0027489038393367993` \[Chi]^6+0.0005821983163192694` \[Chi]^7)/(1-2.8570126619966296` \[Chi]+2.373335413978394` \[Chi]^2-0.6036964688511505` \[Chi]^4+0.0873798215084077` \[Chi]^6);
im\[Omega][\[Chi]_] :=(0.014158792290965177` -0.036989395871554566` \[Chi]+0.026822526296575368` \[Chi]^2+0.0008490933750566702` \[Chi]^3-0.004843996907020524` \[Chi]^4-0.00014745235759327472` \[Chi]^5+0.0001504546201236794` \[Chi]^6)/(1-2.5900842798681376` \[Chi]+1.8952576220623967` \[Chi]^2-0.31416610693042507` \[Chi]^4+0.009002719412204133` \[Chi]^6);

\[Omega]RdDamping[int_, Erad_] = int/(1 - Erad);


Clear[\[Gamma]2, \[Gamma]3]
\[Gamma]2[\[Eta]_, \[Chi]1_, \[Chi]2_] = (PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[6]])//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };
\[Gamma]3[\[Eta]_, \[Chi]1_, \[Chi]2_] =(PhenomCoeff[\[Eta], \[Chi]PN, #]&@FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[7]])//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


(* ::Subsection::Closed:: *)
(*MR Amplitude*)


\[Gamma]1[\[Eta]_, \[Chi]1_, \[Chi]2_] = Block[
	{dummy = FelipeBarbosa`SymDALI`IMRPhenomD`Private`PhenomDTableV[[5]], res},
	res = PhenomCoeff[\[Eta], \[Chi]PN, dummy]//.{\[Chi]PN->Sqrt[1-4 \[Eta]] (\[Chi]1-\[Chi]2)/2 + (1-76 \[Eta]/113) (\[Chi]1+\[Chi]2)/2 };
	res Sqrt[\[Eta]]
]; 


Block[
	{
		ringdownFrequency, dampingFrequency, \[Gamma]1 = U\[Gamma]1[\[Eta], \[Chi]1, \[Chi]2], \[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]1, \[Chi]2], \[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]1, \[Chi]2],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega], aeff
	},
	\[Omega] = f M G;
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	MRAmpExpr = \[Omega]^(-7/6) (M^2 \[Gamma]1 \[Gamma]3 dampingFrequency)/((\[Omega] - ringdownFrequency)^2 + (\[Gamma]3 dampingFrequency)^2) Exp[-((\[Gamma]2 (\[Omega] - ringdownFrequency))/(\[Gamma]3 dampingFrequency))];


]


MRAmpExpr = MRAmpExpr//.{\[Chi]1->s1z, \[Chi]2->s2z};


(* ::Subsection::Closed:: *)
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
		 
		ringdownFrequency,dampingFrequency,aeff,  peakFrequency, \[Gamma]2, \[Gamma]3,
		f1, f2,f3, v1,v2,v3, d1, d3, \[Omega]
	},
	\[Omega] = M G f;

	\[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]1, \[Chi]2];
	\[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]1, \[Chi]2];
	
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
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


IntAmpExpr = IntAmpExpr//.{\[Chi]1->s1z, \[Chi]2->s2z};


(* ::Subsection::Closed:: *)
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
		ringdownFrequency, dampingFrequency,aeff,  peakFrequency, \[Gamma]2, \[Gamma]3,
		f1, f2,f3, v1,v2,v3, d1, d3, \[Omega]
	},
	\[Omega] = M G f;
	
	\[Gamma]2 = U\[Gamma]2[\[Eta], \[Chi]1, \[Chi]2];
	\[Gamma]3 = U\[Gamma]3[\[Eta], \[Chi]1, \[Chi]2];
	
	aeff = Uaeff[UaParallel[m1, m2, s1z, s2z], USperp[m1, m2, s1x, s1y, s2x, s2y]];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][aeff], UErad[m1, m2, s1z, s2z]];
	
	peakFrequency = U\[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2, \[Gamma]3];
	
	U\[ScriptCapitalA]IMR[InsAmpExpr, IntAmpExpr, MRAmpExpr, f, M, peakFrequency]
];


\[ScriptCapitalA]IMRExpr = \[ScriptCapitalA]IMRExpr//.{\[Chi]1->s1z, \[Chi]2->s2z};


(* ::Section::Closed:: *)
(*PhenomPV2 exclusive*)


(* ::Subsection::Closed:: *)
(*c\[Beta] \[And] s\[Beta]*)


(*
v==\[Omega]^(1/3) (convert_spins in IMRPhenomPv2_utils.py)
\[Omega] == (G M \[Pi])/c^3*f, [f] = Hz (convert_spins in IMRPhenomPv2_utils.py)
*)

Block[
	{x = \[Omega]p^(2/3)}, (*Here \[Omega] \[Congruent] (G M \[Pi])/c^3*f, [f] = Hz*)
	
	Clear@L2PNR;
	
	L2PNR[\[Omega]p_, \[Eta]_] = \[Eta]/Sqrt[x] (1+ (3/2+\[Eta]/6) x + (3.375 - (19 \[Eta])/8 - \[Eta]^2/24) x^2)//Simplify; 
]


(*
	v==\[Omega]^(1/3) (convert_spins in IMRPhenomPv2_utils.py)
	L0-> Stands for L2PNR*M^2
	Sp -> \[Chi]p m1^2; m1 >= m2 (PhenomPCoreTwistUp in IMRPhenomPV2_utils.py)
	SL== s1z m1^2 + s2z m2^2 (PhenomPCoreTwistUp in IMRPhenomPV2_utils.py and 
	convert_spins in IMRPhenomPv2_utils.py )
	
*)



Cos\[Beta][L0_, Sp_, SL_] := Module[
	{s, cosBeta, cosBetaHalf, sinBetaHalf}, 
    
    s = Sp / (L0 + SL); 
    
    cosBeta = (1 + s^2)^(-1/2); 
    cosBetaHalf = Sqrt[(1 + cosBeta)/2]; 
    (*sinBetaHalf = Sqrt[(1 - cosBeta)/2]; *)
    
    cosBetaHalf
]

Sin\[Beta][L0_, Sp_, SL_] := Module[
	{s, cosBeta, cosBetaHalf, sinBetaHalf}, 
    
    s = Sp / (L0 + SL); 
    
    cosBeta = (1 + s^2)^(-1/2); 
    (*cosBetaHalf = Sqrt[(1 + cosBeta)/2];*) 
    sinBetaHalf = Sqrt[(1 - cosBeta)/2]; 
    
    sinBetaHalf
]



C\[Beta]Expression = UCos\[Beta][
	(m1+m2)^2 UL2PNR[\[Omega]p, \[Eta]],
	m1^2 U\[Chi]p[m1,m2, s1x, s1y, s2x, s2y],
	m1^2 s1z + m2^2 s2z
]

S\[Beta]Expression = USin\[Beta][
	(m1+m2)^2 UL2PNR[\[Omega]p, \[Eta]],
	m1^2 U\[Chi]p[m1,m2, s1x, s1y, s2x, s2y],
	m1^2 s1z + m2^2 s2z
]


(* ::Subsection::Closed:: *)
(*\[Alpha]*)


Clear@\[Chi]l
\[Chi]l[q_, s1z_, s2z_]= s1z + q s2z; (*q<=1*)


(* ::Subsubsection::Closed:: *)
(*cos\[Theta]Jsf \[And] \[Phi]Jsf*)


Clear[Cos\[Theta]Jsf, \[Phi]Jsf]
(*  (in convert_spins in IMRPhenomPv2_utils)
	J0x_sf== m1^2 s1x+ m2^2 s2x
	J0y_sf== m1^2 s1y + m2^2 s2y
	J0z_sf== M^2 L2PNR[\[Omega]ref, \[Eta]] + m1^2 s1z + m2^2 s2z
	thetaJ_sf
	phiJ_sf
*)


Cos\[Theta]Jsf[J0x_, J0y_, J0z_] := 1/Sqrt[1 + ((J0x^2 + J0y^2)/(J0z^2)) ];

\[Phi]Jsf[J0x_, J0y_] := ArcTan[J0x, J0y];


cos\[Theta]JExpression = UCos\[Theta]Jsf[
	m1^2 s1x +m2^2 s2x,
	m1^2 s1y + m2^2 s2y,
	m1^2 s1z + m2^2 s2z +(m1+m2)^2 UL2PNR[\[Omega]pref, \[Eta]]
];


\[Phi]JExpression = U\[Phi]Jsf[
	m1^2 s1x + m2^2 s2x,
	m1^2 s1y + m2^2 s2y
];


(* ::Subsubsection::Closed:: *)
(*\[Alpha]0 \[And] Cos\[Theta]JN*)


Block[
	{NxSf,NySf,NzSf,tmpX,tmpY,tmpZ,kappa, ROTATEY, ROTATEZ, NxJf,NzJf},
	(* Helper functions for LALtoPhenomP *)
	
	(*in IMRPhenomPv2_utils*)
	ROTATEZ[angle_, x_, y_, z_] := {
	    x * Cos[angle] - y * Sin[angle],
		x * Sin[angle] + y * Cos[angle],
		z
	};

	(*in IMRPhenomPv2_utils*)
	ROTATEY[angle_, x_, y_, z_] := {
		x * Cos[angle] + z * Sin[angle],
		y,
		-x * Sin[angle] + z * Cos[angle]
	};


	(*Nx_sf,Ny_sf, Nz_sf in convert_spins in IMRPhenomPv2_utils*)
	NxSf = Sin[incl] * Cos[Pi/2 - phiRef];
	NySf = Sin[incl] * Sin[Pi/2 - phiRef];
	NzSf = Cos[incl];

	(*tmp_x, tmp_y, tmp_z in convert_spins in IMRPhenomPv2_utils*)
	{tmpX, tmpY, tmpZ} = {NxSf, NySf, NzSf};
	{tmpX, tmpY, tmpZ} = ROTATEZ[-phiJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEY[-thetaJSf, tmpX, tmpY, tmpZ];
	
	(*kappa in convert_spins in IMRPhenomPv2_utils*)
	kappa = -ArcTan[tmpX, tmpY];

	
	(*alpha0 in convert_spins in IMRPhenomPv2_utils*)
	{tmpX, tmpY, tmpZ} = {0, 0, 1};
	{tmpX, tmpY, tmpZ} = ROTATEZ[-phiJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEY[-thetaJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEZ[kappa, tmpX, tmpY, tmpZ];
	(*change y<-> x with respect to Ripple because ArcTan syntax is inverse of jnp.arctan2*)
	alpha0 = ArcTan[tmpX, tmpY];
	
	
	(*thetaJN in convert_spins in IMRPhenomPv2_utils*)
	{tmpX, tmpY, tmpZ} = {NxSf, NySf, NzSf};
	{tmpX, tmpY, tmpZ} = ROTATEZ[-phiJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEY[-thetaJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEZ[kappa, tmpX, tmpY, tmpZ];
	{NxJf, NzJf}= {tmpX, tmpZ};
	
	Clear[Cos\[Theta]JN];
	Cos\[Theta]JN[incl_,Cos\[Theta]JSf_, phiJSf_, phiRef_] = (NzJf//Simplify)//.{
	Cos[thetaJSf] -> Cos\[Theta]JSf,
	Sin[thetaJSf] -> Sqrt[1-Cos\[Theta]JSf^2]
	};
]


Block[{res},
	res = (alpha0//FullSimplify)//.{
		Cos[thetaJSf] -> Cos\[Theta]JSf,
		Sin[thetaJSf] -> Sqrt[1-Cos\[Theta]JSf^2]
	};
	
	
	Cases[res, x_Symbol, Infinity]//DeleteDuplicates//Echo;
	
	\[Alpha]0[Cos\[Theta]JSf_, phiJSf_, phiRef_, incl_] = res;
]//QuietEcho


\[Alpha]0Expr = U\[Alpha]0[
	cos\[Theta]JExpression,
	\[Phi]JExpression,
	\[Phi]ref,
	\[Iota]
]


Cos\[Theta]JNExpression = UCos\[Theta]JN[
	\[Iota],
	cos\[Theta]JExpression,
	\[Phi]JExpression,
	\[Phi]ref
]


(* ::Subsubsection::Closed:: *)
(*PPN terms*)


angcoeffs = <|
    "alphacoeff1" -> -0.18229166666666666 - (5 dm)/(64.0 m2),
    
    "alphacoeff2" -> (-15 dm m2 chil)/(128.0 mtot2 eta) - (35 m2^2 chil)/(128.0 mtot2 eta),
    
    "alphacoeff3" -> -1.7952473958333333 - (4555 dm)/(7168.0 m2) - (15 chip2 dm m2^3)/(128.0 mtot4 eta2) 
        - (35 chip2 m2^4)/(128.0 mtot4 eta2) - (515 eta)/384.0 - (15 dm^2 eta)/(256.0 m2^2) - (175 dm eta)/(256.0 m2),
    
    "alphacoeff4" -> -(35 \[Pi])/48.0 - (5 dm \[Pi])/(16.0 m2) + (5 dm^2 chil)/(16.0 mtot2) 
        + (5 dm m2 chil)/(3.0 mtot2) + (2545 m2^2 chil)/(1152.0 mtot2) - (5 chip2 dm m2^5 chil)/(128.0 mtot6 eta3) 
        - (35 chip2 m2^6 chil)/(384.0 mtot6 eta3) + (2035 dm m2 chil)/(21504.0 mtot2 eta) + (2995 m2^2 chil)/(9216.0 mtot2 eta),
    
    "alphacoeff5" -> 4.318908476114694 + (27895885 dm)/(2.1676032*10^7 m2) - (15 chip4 dm m2^7)/(512.0 mtot8 eta4) 
        - (35 chip4 m2^8)/(512.0 mtot8 eta4) - (485 chip2 dm m2^3)/(14336.0 mtot4 eta2) + (475 chip2 m2^4)/(6144.0 mtot4 eta2) 
        + (15 chip2 dm^2 m2^2)/(256.0 mtot4 eta) + (145 chip2 dm m2^3)/(512.0 mtot4 eta) + (575 chip2 m2^4)/(1536.0 mtot4 eta) 
        + (39695 eta)/86016.0 + (1615 dm^2 eta)/(28672.0 m2^2) - (265 dm eta)/(14336.0 m2) + (955 eta2)/576.0 
        + (15 dm^3 eta2)/(1024.0 m2^3) + (35 dm^2 eta2)/(256.0 m2^2) + (2725 dm eta2)/(3072.0 m2) 
        - (15 dm m2 Pi chil)/(16.0 mtot2 eta) - (35 m2^2 Pi chil)/(16.0 mtot2 eta) 
        + (15 chip2 dm m2^7 chil2)/(128.0 mtot8 eta4) + (35 chip2 m2^8 chil2)/(128.0 mtot8 eta4) 
        + (375 dm^2 m2^2 chil2)/(256.0 mtot4 eta) + (1815 dm m2^3 chil2)/(256.0 mtot4 eta) + (1645 m2^4 chil2)/(192.0 mtot4 eta)
|>;


Block[
	{m2,m1,dm,mtot,eta,eta2,eta3,eta4,mtot2,mtot4,mtot6,mtot8,chil2,chip2,chip4,dm2,dm3, 
	rules, association},
	
	rules  = { (*All defs Ripple does in IMRPhenomPv2_utils.py*)
		m2->q/(1.` +q),m1->1.`/(1.` +q),dm->m1-m2,mtot->1.`,eta->m1 m2,eta2->eta eta,
		eta3->eta2*eta,eta4->eta3*eta,mtot2->mtot*mtot,mtot4->mtot2*mtot2,
		mtot6->mtot4*mtot2,mtot8->mtot6*mtot2,chil2-> chil*chil,chip2->chip*chip,
		chip4->chip2*chip2,dm2->dm*dm,dm3->dm2*dm
	};
	
	association = Simplify/@(angcoeffs//.rules);
	
	(*When they call this function they don't pass q<=1 but rather q>=1, So lets change q->1/q to have
	the usual variable:*)
	
	association = association//.q->qp^-1;
	association = Simplify/@(association//.qp->q);
	NewDictionary = association;
]


Clear[{"Global`\[CapitalAlpha]*"}]


(*Define the functions now.
All functions of chil, q, chip, with q<=1
*)
\[CapitalAlpha]1[q_] = NewDictionary["alphacoeff1"];
\[CapitalAlpha]2[q_, chil_] = NewDictionary["alphacoeff2"];
\[CapitalAlpha]3[q_, chip_] = NewDictionary["alphacoeff3"];
\[CapitalAlpha]4[q_, chil_, chip_] = NewDictionary["alphacoeff4"];
\[CapitalAlpha]5[q_,chil_, chip_] = NewDictionary["alphacoeff5"];


(* ::Subsubsection::Closed:: *)
(*\[Alpha] expr*)


\[Alpha][term1_, term2_, term3_, term4_, term5_, term0_] := (
	term1 + term2 + term3 + term4 + term5 + term0 
)


Module[
	{\[Alpha]0 = \[Alpha]0Expr, \[Chi]l = \[Chi]l[q, s1z, s2z], \[Chi]p = \[Chi]p[m1, m2, s1x, s1y, s2x, s2y]},
	\[Alpha]Expression = U\[Alpha][
		U\[CapitalAlpha]1[q] (\[Omega]^-1 - \[Omega]ref^-1),
		U\[CapitalAlpha]2[q, \[Chi]l] (\[Omega]^(-2/3)-\[Omega]ref^(-2/3)),
		U\[CapitalAlpha]3[q, \[Chi]p] (\[Omega]^(-1/3) - \[Omega]ref^(-1/3)),
		U\[CapitalAlpha]4[q, \[Chi]l, \[Chi]p] Log[\[Omega]/\[Omega]ref],
		U\[CapitalAlpha]5[q, \[Chi]l, \[Chi]p] (\[Omega]^(1/3)-\[Omega]ref^(1/3)), 
		\[Alpha]0
	]
];

\[Alpha]Expression = \[Alpha]Expression//.{\[Omega]->\[Omega]p, \[Omega]ref->\[Omega]pref};


(* ::Subsection::Closed:: *)
(*T2m \[And] Tm2m*)


T2m[alpha_, cos\[Beta]_, sin\[Beta]_, cos\[Theta]_] = 1/8 E^(-2 I alpha) Sqrt[5/\[Pi]] (cos\[Beta] Sqrt[1-cos\[Theta]] E^(I alpha)-  Sqrt[1+cos\[Theta]] sin\[Beta])^4


Tm2m[alpha_, cos\[Beta]_, sin\[Beta]_, cos\[Theta]_] = 1/8 E^(-2 I alpha) Sqrt[5/\[Pi]] (cos\[Beta] Sqrt[1 + cos\[Theta]] E^(I alpha) +Sqrt[1-cos\[Theta]] sin\[Beta])^4


(* ::Text:: *)
(*Notice that Tm2m[\[Alpha], cos\[Beta], sin\[Beta], cos\[Theta]] == T2m[\[Alpha], cos\[Beta], -sin\[Beta], -cos\[Theta]], So we need only 1 function:*)


T2mExpression = UT2m[
	\[Alpha]Expression,
	C\[Beta]Expression,
	S\[Beta]Expression,
	Cos\[Theta]JNExpression
];

Tm2mExpression = UT2m[
	\[Alpha]Expression,
	C\[Beta]Expression,
	-S\[Beta]Expression,
	-Cos\[Theta]JNExpression
];


(* ::Subsection::Closed:: *)
(*\[Zeta]*)


Block[
	{
		NxSf,NySf,NzSf,tmpX,tmpY,tmpZ,kappa, ROTATEY, ROTATEZ, NxJf,NzJf,
		Xxsf,Xysf,Xzsf,PArunxJf,PArunyJf,PArunzJf,
		QArunxJf,QArunyJf,QArunzJf,XdotPArun,XdotQArun
	},
	(* Helper functions for LALtoPhenomP *)
	ROTATEZ[angle_, x_, y_, z_] := {
	    x * Cos[angle] - y * Sin[angle],
		x * Sin[angle] + y * Cos[angle],
		z
	};

	ROTATEY[angle_, x_, y_, z_] := {
		x * Cos[angle] + z * Sin[angle],
		y,
		-x * Sin[angle] + z * Cos[angle]
	};

(* Compute components of N in the source frame *)
	NxSf = Sin[incl] * Cos[Pi/2 - phiRef];
	NySf = Sin[incl] * Sin[Pi/2 - phiRef];
	NzSf = Cos[incl];

(* Apply rotations to N *)
	{tmpX, tmpY, tmpZ} = {NxSf, NySf, NzSf};
	{tmpX, tmpY, tmpZ} = ROTATEZ[-phiJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEY[-thetaJSf, tmpX, tmpY, tmpZ];

	kappa = -ArcTan[tmpX, tmpY];

(* Compute alpha0 by rotating LN *)
	{tmpX, tmpY, tmpZ} = {0, 0, 1};
	{tmpX, tmpY, tmpZ} = ROTATEZ[-phiJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEY[-thetaJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEZ[kappa, tmpX, tmpY, tmpZ];
	
	
	(*Finally we determine thetaJ, by rotating N*)
	{tmpX, tmpY, tmpZ} = {NxSf, NySf, NzSf};
	{tmpX, tmpY, tmpZ} = ROTATEZ[-phiJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEY[-thetaJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEZ[kappa, tmpX, tmpY, tmpZ];
	{NxJf, NzJf}= {tmpX, tmpZ};
	
	
	
	(*\[Zeta] stuff: *)
	
	Xxsf = -Cos[incl] Sin[phiRef];
	Xysf = -Cos[incl] Cos[phiRef];
	Xzsf = Sin[incl];

	{tmpX, tmpY, tmpZ} = {Xxsf, Xysf, Xzsf};
	{tmpX, tmpY, tmpZ} = ROTATEZ[-phiJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEY[-thetaJSf, tmpX, tmpY, tmpZ];
	{tmpX, tmpY, tmpZ} = ROTATEZ[kappa, tmpX, tmpY, tmpZ];

	(* Components of X in the J frame *)

	PArunxJf = 0.0;
	PArunyJf = -1.0;
	PArunzJf = 0.0;

	(* Q = NxP *)
	QArunxJf = NzJf;
	QArunyJf = 0.0;
	QArunzJf = -NxJf;

	(* Compute dot products *)
	XdotPArun = tmpX * PArunxJf + tmpY * PArunyJf + tmpZ * PArunzJf;
	XdotQArun = tmpX * QArunxJf + tmpY * QArunyJf + tmpZ * QArunzJf;
	
	Clear@\[Zeta];
	\[Zeta][Cos\[Theta]JSf_,phiJSf_,phiRef_,incl_] = (ArcTan[XdotPArun, XdotQArun]//FullSimplify)//.{
		Cos[thetaJSf] -> Cos\[Theta]JSf,
		Sin[thetaJSf] -> Sqrt[1-Cos\[Theta]JSf^2]
	};		
]


\[Zeta]Expression = U\[Zeta][
	cos\[Theta]JExpression,
	\[Phi]JExpression,
	\[Phi]ref,
	\[Iota]
]


(* ::Section::Closed:: *)
(*Detector Pattern Functions*)


f1[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_, D33_] = Module[
	{expr},
	expr = (
		(D22 - D11 Cos[\[Theta]]^2) Cos[\[Phi]]^2
		- D33 Sin[\[Theta]]^2
		+ Sin[2 \[Theta]] (D13 Cos[\[Phi]] + D23 Sin[\[Phi]])
		+ (D11 - D22 Cos[\[Theta]]^2) Sin[\[Phi]]^2
		- D12 Sin[2 \[Phi]] (1+Cos[\[Theta]]^2)
	);
	
	expr = (expr)//.{Cos[d_]^2 :> (1+ Cos[2 d])/2, Sin[d_]^2 :> (1- Cos[2 d])/2 };
	
	
	expr//Simplify
] 


f2[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_] = (
	2 Sin[\[Theta]] (D13 Sin[\[Phi]] - D23 Cos[\[Phi]]) +
	Cos[\[Theta]] (2 D12 Cos[2 \[Phi]] + (D22-D11) Sin[2 \[Phi]])
)//Simplify


Fplus[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	Cos[2 \[Psi]] f1[\[Theta], \[Phi], D11, D12, D13, D22, D23, D33] + 
	Sin[2 \[Psi]] f2[\[Theta], \[Phi], D11, D12, D13, D22, D23]
);

Fx[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	- Sin[2 \[Psi]] f1[\[Theta], \[Phi], D11, D12, D13, D22, D23, D33] +
	Cos[2 \[Psi]] f2[\[Theta], \[Phi], D11, D12, D13, D22, D23]
);


FpPlusFc[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	Fplus[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33] + 
	I Fx[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33]
)//FullSimplify;

FpMinusFc[\[Theta]_, \[Phi]_, \[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	Fplus[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33] -
	I Fx[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33]
)//FullSimplify;


\[CapitalDelta]t2[p1_, p2_, p3_, \[Theta]_, \[Phi]_] = Module[
	{\[Delta], r, c = UnitConvert["SpeedOfLight"][[1]], detectorposition},
	
	r = FromSphericalCoordinates[{1, \[Theta], \[Phi]}];
	
	detectorposition = {p1,p2,p3};
	
	(detectorposition . r/c)
]


(* ::Section::Closed:: *)
(*Making Block Function*)


varsDefs = HoldForm[{
	M = m1+m2,
	\[Omega] = G (m1+m2) f,
	\[Omega]p =G (m1+m2) \[Pi] f,
	\[Omega]ref = G (m1+m2) fref,
	\[Omega]pref = G (m1+m2) \[Pi] fref,
	\[Eta] = m2 m1/(m2+m1)^2,
	q= m2/m1,
	\[Delta] = Sqrt[1 - 4 \[Eta]],
	\[Chi]s = (s1z+s2z)/2,
	\[Chi]a = (s1z-s2z)/2,
	\[Chi]PN = Sqrt[1 - 4 \[Eta]] (s1z-s2z)/2  + (1-76 \[Eta]/113) (s1z+s2z)/2,
	S = 1/4 (1+Sqrt[1 - 4 \[Eta]])^2 s1z + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 s2z,
	Shat = (1/4 (1+Sqrt[1 - 4 \[Eta]])^2 s1z + 1/4 (1-Sqrt[1 - 4 \[Eta]])^2 s2z)/(1- 2 \[Eta])
}];


AmpDefs = <||>;


AmpDefs["Ins"] = Block[
	
	{rule, sym, defs},
	sym = {a0,a1,a2,a3,a4,a5,a6,r1,r2,r3};
	
	defs = {
		\[ScriptCapitalA]0[\[Eta],M],\[ScriptCapitalA]1[\[Eta],M],\[ScriptCapitalA]2[\[Eta],M],\[ScriptCapitalA]3[\[Eta],M,s1z,s2z],\[ScriptCapitalA]4[\[Eta],M,s1z,s2z],\[ScriptCapitalA]5[\[Eta],M,s1z,s2z],\[ScriptCapitalA]6[\[Eta],M,s1z,s2z],
		\[Rho]1[\[Eta],M,s1z,s2z],\[Rho]2[\[Eta],M,s1z,s2z],\[Rho]3[\[Eta],M,s1z,s2z]
	};
	
	rule = MapThread[
		Rule,
		{sym, defs}
	];
	
	
	
	HoldForm[{
		\[ScriptCapitalA]0[\[Eta]_,  M_] := a0, 
		\[ScriptCapitalA]1[\[Eta]_, M_] := a1, 
		\[ScriptCapitalA]2[\[Eta]_, M_] := a2, 
		\[ScriptCapitalA]3[\[Eta]_, M_, s1z_, s2z_] := a3, 
		\[ScriptCapitalA]4[\[Eta]_, M_, s1z_, s2z_] := a4, 
		\[ScriptCapitalA]5[\[Eta]_, M_, s1z_, s2z_] := a5, 
		\[ScriptCapitalA]6[\[Eta]_, M_, s1z_, s2z_] := a6,
		\[Rho]1[\[Eta]_, M_, s1z_, s2z_] := r1,
		\[Rho]2[\[Eta]_, M_, s1z_, s2z_] := r2, 
		\[Rho]3[\[Eta]_, M_, s1z_, s2z_] := r3
	}]/.rule
];


AmpDefs["MR"] = HoldForm[{
	\[Gamma]1[\[Eta]_, s1z_, s2z_] := g1,
	\[Gamma]2[\[Eta]_, s1z_, s2z_] := g2,
	\[Gamma]3[\[Eta]_, s1z_, s2z_] := g3,
	
	aParallel[m1_, m2_, s1z_, s2z_] := x,
	
	aeff[aParallel_, Sperp_] := a,
	Erad[m1_, m2_, s1z_, s2z_] := y,
	
	Sperp[m1_, m2_, s1x_, s1y_, s2x_, s2y_] := y0,
	
	
	re\[Omega][\[Chi]_] := r,
	im\[Omega][\[Chi]_] := i,
	(*Erad[\[Eta]_, \[Chi]1_, \[Chi]2_] = erad,*)
	\[Omega]RdDamping[int_, Erad_] := s,
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] := If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ]
}]/.{
	g1 -> \[Gamma]1[\[Eta], s1z, s2z], g2 -> \[Gamma]2[\[Eta], s1z, s2z], g3 -> \[Gamma]3[\[Eta], s1z, s2z], x-> aParallel[m1, m2, s1z, s2z],
	a-> aeff[aParallel, Sperp], r -> re\[Omega][\[Chi]], i-> im\[Omega][\[Chi]], y-> Erad[m1, m2, s1z, s2z], s-> \[Omega]RdDamping[int, Erad],
	y0-> Sperp[m1, m2, s1x, s1y, s2x, s2y]
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
		\[CapitalDelta]0[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] := D0,
		\[CapitalDelta]1[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] :=D1,
		\[CapitalDelta]2[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] := D2,
		\[CapitalDelta]3[\[Eta]_, M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] := D3,
		\[CapitalDelta]4[\[Eta]_, M_,f2_,f3_, v1_, v2_, v3_, d1_, d3_] := D4,
		v2[\[Eta]_, M_, \[Chi]1_, \[Chi]2_] := V2
	}]/.rule

];


AmpDefs["\[ScriptCapitalA]IMR"] = HoldForm[{
	\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins_, \[ScriptCapitalA]Int_, \[ScriptCapitalA]MR_, f_, M_, \[Omega]Peak_] := x
}]/.x->\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins, \[ScriptCapitalA]Int, \[ScriptCapitalA]MR, f, M, \[Omega]Peak];


AmpDefs["PV2Exclusive"] = HoldForm[{
	L2PNR[\[Omega]p_, \[Eta]_] := l2,
	\[Chi]p[m1_, m2_, s1x_, s1y_, s2x_, s2y_] := Chip,
	Cos\[Beta][L0_, Sp_, SL_] := cbeta,
	Sin\[Beta][L0_, Sp_, SL_] := sbeta,
	Cos\[Theta]Jsf[J0x_, J0y_, J0z_] := 1/Sqrt[1 + ((J0x^2 + J0y^2)/(J0z^2)) ],
	\[Phi]Jsf[J0x_, J0y_] := ArcTan[J0x, J0y],
	Cos\[Theta]JN[incl_,Cos\[Theta]JSf_, phiJSf_, phiRef_] := c\[Theta],
	\[Alpha]0[Cos\[Theta]JSf_, phiJSf_, phiRef_, incl_] := res,
	\[CapitalAlpha]1[q_] := a1,
	\[CapitalAlpha]2[q_, chil_] := a2,
	\[CapitalAlpha]3[q_, chip_] := a3,
	\[CapitalAlpha]4[q_, chil_, chip_] := a4,
	\[CapitalAlpha]5[q_,chil_, chip_] := a5,
	\[Alpha][term1_, term2_, term3_, term4_, term5_, term0_] := (term1 + term2 + term3 + term4 + term5 + term0 ),
	T2m[alpha_, cos\[Beta]_, sin\[Beta]_, cos\[Theta]_] := 1/8 E^(-2 I alpha) Sqrt[5/\[Pi]] (cos\[Beta] Sqrt[1-cos\[Theta]] E^(I alpha)-  Sqrt[1+cos\[Theta]] sin\[Beta])^4,
	\[Zeta][Cos\[Theta]JSf_,phiJSf_,phiRef_,incl_] :=zeta
}]/.{
 l2-> L2PNR[\[Omega]p, \[Eta]], cbeta -> Cos\[Beta][L0, Sp, SL], sbeta-> Sin\[Beta][L0, Sp, SL], c\[Theta]-> Cos\[Theta]JN[incl,Cos\[Theta]JSf, phiJSf, phiRef],
 res -> \[Alpha]0[Cos\[Theta]JSf, phiJSf, phiRef, incl], a1-> \[CapitalAlpha]1[q], a2->\[CapitalAlpha]2[q, chil], a3-> \[CapitalAlpha]3[q, chip], a4-> \[CapitalAlpha]4[q, chil, chip],
 a5-> \[CapitalAlpha]5[q, chil, chip], zeta-> \[Zeta][Cos\[Theta]JSf,phiJSf,phiRef,incl],
 Chip-> \[Chi]p[m1,m2,s1x,s1y,s2x,s2y]

};


AmpDefs["FpFc"] = HoldForm[{
	FpPlusFc[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] := FPPLUSFC,
	FpMinusFc[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] := FPMINUSFC,
	\[CapitalDelta]t2[p1_, p2_, p3_, \[Theta]_, \[Phi]_] := timeDelay
}]/.{
	FPPLUSFC -> FpPlusFc[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33],
	FPMINUSFC->FpMinusFc[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33],
	timeDelay-> \[CapitalDelta]t2[p1,p2, p3, \[Theta], \[Phi]]
};


(*AmpDefs["Commute D"] = HoldForm[{
	(*Commute D with \[Phi]IMR*)
	\[ScriptCapitalA]IMR/: D[\[ScriptCapitalA]IMR[x__], n___] := With[
		{args = D[#, n]&/@({x}[[1;;3]])},
		\[ScriptCapitalA]IMR@@(Join[args, {x}[[4;;6]]])
	]
}];*)


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
	dummy = HoldForm[{
		aIMR/2 (Exp[2 I zeta] FplusPlusFcross t2m  + Exp[-2 I zeta] FplusMinusFcross tm2m) Exp[-I 2 \[Pi] f \[Delta]t]
	}]/.{
		aIMR->\[ScriptCapitalA]IMRExpr, zeta-> \[Zeta]Expression, FplusPlusFcross-> UFpPlusFc[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33],t2m ->T2mExpression,
		FplusMinusFcross-> UFpMinusFc[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33], tm2m-> Tm2mExpression, \[Delta]t-> U\[CapitalDelta]t2[p1,p2,p3,\[Theta],\[Phi]],
		alpha->\[Alpha]Expression, aMR-> MRAmpExpr
	};

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


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]ref, \[Iota], \[Theta], \[Phi], \[Psi]};


(*All derivatives up to order 3*)
derivatives = Combinations[vars, 1];
PrependTo[derivatives, {}];


derivatives


With[{
	vars = {f, fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]ref, \[Iota], \[Theta], \[Phi], \[Psi], p1, p2, p3, D11, D12, D13, D22, D23, D33}
	},
	expr =  Hold[
		{test, vars, derivatives},
		Evaluate@$blockexpr,
		"KeepDefs" -> KeepDefs,
		"IncludeZeroDerivative"->False
	]//.HoldForm[X_] :> X;
]


<<FelipeBarbosa`SymDALI`


Ds = MemoryConstrained[
	DerivativeRules@@expr,
	8 10^9]//EchoTiming;


Block[
	{testF,
	m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]ref, \[Iota], \[Theta], \[Phi], \[Psi],
	f, fref,p1, p2, p3, D11, D12, D13, D22, D23, D33,
	 M, G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]
	},
	
	{p1, p2, p3, D11, D12, D13, D22, D23, D33} = ConstantArray[0, 9];
	D22 = 1;
	
	testF[f_, fref_, m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_,\[Phi]ref_,\[Iota]_,\[Theta]_,\[Phi]_,\[Psi]_] = Ds[[1,2]]; 
	DownValues[testF] = DownValues[testF]//. HoldForm[x_] :> x;
	
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100}, 2]];
	{s1x,s1y,s1z}  = RandomReal[{-1,1},3];
	{s2x,s2y,s2z}  = RandomReal[{-1,1},3];
	
	{\[Phi]ref, \[Phi]}= RandomReal[{0, 2 \[Pi]}, 2];
	{\[Theta], \[Iota], \[Psi]} = RandomReal[{0, \[Pi]}, 3];
	
	f = RandomReal[{10, 0.2/(G (m1+m2))}];
	fref=10;

	
	testF[f,fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Phi]ref,\[Iota],\[Theta],\[Phi],\[Psi]]
](*//Cases[#, x_Symbol, Infinity]&//DeleteDuplicates*)


SetDirectory[NotebookDirectory[]]


Export["Amp_Ds_order_0_to_1.wdx", Ds]


(* ::Section::Closed:: *)
(*Testing function*)


Clear@RelativeDiff

Attributes[RelativeDiff] = {Listable};
RelativeDiff[x_,y_]/;x==0 && y==0 := 0
RelativeDiff[x_, y_]/; x==0&&y!=0 := 1
RelativeDiff[x_, y_]/; y==0&&x!=0 := 1

RelativeDiff[x_, y_]/; x!=0 &&y!=0 := With[
	{diff = x-y},
	Max[Abs@{diff/x, diff/y}]
]


SetDirectory[NotebookDirectory[]];
Ds = Import["Amp_Ds_order_0_to_1.wdx"];


Block[{g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	ClearAll[Test\[ScriptCapitalA]];
	Test\[ScriptCapitalA][f_, fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_,\[Phi]ref_,\[Iota]_,\[Theta]_,\[Phi]_,\[Psi]_] = Ds[[1,2]]//.G -> g;

]
DownValues[Test\[ScriptCapitalA]] = DownValues[Test\[ScriptCapitalA]]//.HoldForm[x_]:> x;


(* ::Subsection::Closed:: *)
(*Testing the Amplitude against Ripple*)


DeleteObject/@ExternalSessions[]

Clear@python

python = StartExternalSession["Python"];
ExternalEvaluate[python,"
import numpy as np

import jax
import jax.numpy as jnp

from jax import grad, vmap
from functools import partial

from ripplegw.waveforms import IMRPhenomPv2
from ripplegw import get_match_arr, get_eff_pads
from ripplegw import ms_to_Mc_eta
from ripplegw.constants import MSUN, gt
"]


R\[ScriptCapitalA] = ExternalFunction[python,"
def Ripple_hp(f, f_ref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dL, tc,iota, phi_ref):
    f_array = jnp.array(f)

    Mc = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    eta = m1 * m2 / ((m1 + m2)**2)

    theta = jnp.array([Mc, eta, s1x, s1y, s1z, s2x, s2y, s2z, dL, tc, phi_ref, iota])

    hp, hc, t2m, tm2m, zeta, epsilon, phi_Jsf, t0, a, psi, alpha = IMRPhenomPv2.gen_IMRPhenomPv2_hphc(f_array, theta, f_ref)
    
    result = (jnp.exp(1j*2*zeta)*t2m + jnp.exp(-1j*2*zeta)*tm2m )*a/2
    #result = tm2m
    return result.tolist()


"]


(* ::Text:: *)
(*Note that *)


Fplus[0,0,0,0,0,0,1,0,0]
Fx[0,0,0,0,0,0,1,0,0]


Clear@Test

Test := Block[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA,Ripplehp, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diffRe,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
		p1, p2, p3, D11, D12, D13, D22, D23, D33, \[Theta], \[Phi], \[Psi], ReMMA, ReRipple, ImMMA, ImRipple, diffIm
	},
	
	
	{p1, p2, p3, D11, D12, D13, D22, D23, D33, \[Theta], \[Phi], \[Psi]} = ConstantArray[0, 12];
	D22 = 1;
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	
	{s1x, s1y, s1z} = RandomReal[{-1,1}, 3];
	{s2x, s2y, s2z} = RandomReal[{-1,1}, 3];
	
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	
	MMA = 10^3/(2 Sqrt[5/(64 \[Pi])] ) Test\[ScriptCapitalA][f, 10, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref,\[Iota], \[Theta], \[Phi],\[Psi]];
	
	Ripplehp = R\[ScriptCapitalA][f,10, m1,m2,s1x,s1y,s1z,s2x,s2y,s2z, 1, 0, \[Iota],\[Phi]Ref];
	
	ReMMA = Re[MMA];
	ReRipple = Re[Ripplehp];
	
	ImMMA = Im[MMA];
	ImRipple = Im[Ripplehp];
	
	diffRe = RelativeDiff@@{ReMMA, ReRipple};
	
	diffIm = RelativeDiff@@{ImMMA, ImRipple};
	
	
	{ListLinePlot[
		{ReMMA, ReRipple},
		Frame->True,
		PlotLegends->{"MMA", "Ripple"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	],
	ListLinePlot[
		diffRe,
		Frame->True,
		PlotLegends->{"Relative Difference"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	],
	
	ListLinePlot[
		{ImMMA, ImRipple},
		Frame->True,
		PlotLegends->{"MMA", "Ripple"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	],
	ListLinePlot[
		diffIm,
		Frame->True,
		PlotLegends->{"Relative Difference"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	]}
	
]


Test


(* ::Subsection::Closed:: *)
(*Testing Symbolic vs Numerical derivatives*)


NGrad//Clear
NGrad[f_, vars_, n_] := Module[
	{h = 1. 10^-6, dummy, Point1, Point2},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = vars[[i]] + h vars[[i]];
		dummy,
		{i, n+1, Length@vars}
	];
	
	
	Point2 = ConstantArray[vars, (Length[vars] - n)];
	
	
	(f@@@Point1 - f@@@Point2)/(h vars[[n+1;;-1]])
]


{SymRules, NRules} =  DerivativeRulesLoad["IMRPhenomPv2"];


NRules["\[ScriptA]IMR"][[1]]


(*Symbolic grad from compiled functions*)

Block[
	{p1, p2, p3, D11, D12, D13, D22, D23, D33, args},
	
	{D12, D13, D23} = {0,0,0};
	{D11, D22, D33} = {1,0,1};
	{p1, p2, p3} = {1,1,1};
	
	args = Join[{f,fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Phi]ref,\[Iota],\[Theta],\[Phi],\[Psi]}, {p1, p2, p3, D11, D12, D13, D22, D23, D33}];

	ClearAll[TestGrad\[ScriptCapitalA]];
	
	TestGrad\[ScriptCapitalA][f_, fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_,\[Phi]ref_,\[Iota]_,\[Theta]_,\[Phi]_,\[Psi]_] = Table[
		NRules["\[ScriptA]IMR"][[i, 3]]@@args, 
		{i, 2, 14}
	
	]
];
DownValues[TestGrad\[ScriptCapitalA]] = DownValues[TestGrad\[ScriptCapitalA]]//.HoldForm[x_]:> x;



(*BELOW IS THE DEFINITION FOR BLOCK FUNCTIONS*)


(*Block[{g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	
	ClearAll[TestGrad\[ScriptCapitalA]];
	
	TestGrad\[ScriptCapitalA][f_, fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_,\[Phi]ref_,\[Iota]_,\[Theta]_,\[Phi]_,\[Psi]_] = Ds[[2;;-1, 2]]//.G->g;

]
DownValues[TestGrad\[ScriptCapitalA]] = DownValues[TestGrad\[ScriptCapitalA]]//.HoldForm[x_]:> x;*)


RandomSpin[] := Module[
	{norm=2, v},
	While[norm>1, v = RandomReal[{-1,1}, 3]; norm = Norm[v]];
	v
]


Clear@Test
Test := Block[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, f, Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,\[Phi]ref, \[Iota],
		\[Theta], \[Phi], \[Psi], G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		D12, D13, D23, D11, D22, D33, p1, p2, p3, Rediff, Imdiff
	},
	
	
	
	{D12, D13, D23} = {0,0,0};
	{D11, D22, D33} = {1,0,1};
	{p1, p2, p3} = {1,1,1};
	
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	
	{s1x, s1y, s1z} = RandomSpin[];
	{s2x, s2y, s2z} = RandomSpin[];
	
	f = RandomReal[{10., 0.2/(G (m1+m2))}]//List;
	fref = RandomReal[10];
	{\[Phi]ref, \[Phi]} = RandomReal[{0, 2 \[Pi]},2];
	{\[Theta], \[Iota], \[Psi]}  = RandomReal[{0, \[Pi]}, 3];
	
	
	Symbolic = Last/@TestGrad\[ScriptCapitalA][f, fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]ref, \[Iota], \[Theta], \[Phi], \[Psi]];
	vars = {f[[1]], fref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]ref, \[Iota], \[Theta], \[Phi], \[Psi]};
	Numeric = NGrad[Test\[ScriptCapitalA], vars, 2];
	
	Rediff = RelativeDiff@@{Abs@Symbolic, Abs@Numeric};
	Imdiff = RelativeDiff@@{Arg@Symbolic, Arg@Numeric};
	
	{
		Rediff,
		Imdiff
	}//Transpose//Round//MatrixForm

	
	

	
]


Table[Test, {50}]


(* ::Section:: *)
(*Compiling*)


Module[{ds = Import["Amp_Ds_order_0_to_1.wdx"]}, Ds = ds//.G->UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]];


compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{
			{f, _Real, 1}, fref,
			m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Phi]ref,\[Iota],
			\[Theta],\[Phi],\[Psi],
			p1,p2,p3,D11,D12,D13,D22,D23,D33
		
		}, 
		Evaluate[x],
		RuntimeOptions->{
			"CatchMachineOverflow"->False,
			"CatchMachineIntegerOverflow"->False,
			"EvaluateSymbolically"->False,
			"RuntimeErrorHandler"->None,
			"WarningMessages"->True
		}
	]//.{HoldForm[y_] :> y, us\[Theta]->UnitStep, Rational[any_, any2_] :> Divide[any, any2], Complex[a_, b_] :> a +I b};
	
	Compile@@dummy
]


<<CompiledFunctionTools`
compileThis[Ds[[-1, 2]]]//CompilePrint


compiledDs = MapAt[
	compileThis,
	Ds, 
	{All, 2}
];


<<CCompilerDriver`


$CCompilerDefaultDirectory = FileNameJoin[{
	ParentDirectory[NotebookDirectory[], 2],
	"/LibraryResources/",
	$SystemID,
	"/DerivativeRules/IMRPhenomPv2/NRules"
}]


Needs["CCodeGenerator`"]


Clear@Ds


MapIndexed[
	LibraryGenerate[#1[[2]], "a" <> ToString[#2//First]]&,
	compiledDs
];


(* ::Chapter::Closed:: *)
(*Testing Compiled functions for the Waveform*)


(* ::Input:: *)
(**)


Clear@RelativeDiff

Attributes[RelativeDiff] = {Listable};
RelativeDiff[x_,y_]/;x==0 && y==0 := 0
RelativeDiff[x_, y_]/; x==0&&y!=0 := 1
RelativeDiff[x_, y_]/; y==0&&x!=0 := 1

RelativeDiff[x_, y_]/; x!=0 &&y!=0 := With[
	{diff = x-y},
	Max[Abs@{diff/x, diff/y}]
]


NotebookDirectory[]//ParentDirectory[#, 2]&


SetDirectory[FileNameJoin[{
	NotebookDirectory[]//ParentDirectory[#, 2]&,
	"/LibraryResources", 
	$SystemID,
	"DerivativeRules/IMRPhenomPv2/NRules"
}]]


{\[CapitalPsi], \[ScriptCapitalA]} = Module[
	{\[Psi], a, direc, \[Psi]name, \[Psi]args, aname, aArgs},
	direc =  FileNameJoin[{
		NotebookDirectory[]//ParentDirectory[#, 2]&,
		"/LibraryResources", 
		$SystemID,
		"DerivativeRules/IMRPhenomPv2/NRules"
	}];
	
	(*CHANGE .dylib to something else if you are in linux or Windows*)
	\[Psi]name = FileNameJoin[{direc, "phi1.so"}];
	aname = FileNameJoin[{direc, "a1.so"}];
	
	aArgs = Join[{{Real,1}}, ConstantArray[Real, 23]];
	\[Psi]args = Join[{{Real,1}}, ConstantArray[Real, 25]];
	
	\[Psi] = LibraryFunctionLoad[\[Psi]name, "phi1", \[Psi]args, {Real,1}];
	a = LibraryFunctionLoad[aname, "a1", aArgs, {Complex,1}];
	
	{\[Psi], a}

]


Clear[hp, hx]
hp[f_,fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_,\[Phi]ref_,\[Iota]_] := Module[
	{\[Theta], \[Phi], \[Psi], p1, p2, p3, D11, D12, D13, D22, D23, D33, aArgs, \[Psi]Args},
	
	
	
	{\[Theta], \[Phi], \[Psi], p1, p2, p3} = ConstantArray[0, 6];
	
	{D11, D12, D13, D22, D23, D33} = ConstantArray[0,6];
	
	(*Fplus=1, Fx=0*)
	D22=1; 
	
	aArgs = {f, fref, m1,m2, s1x,s1y,s1z,s2x,s2y,s2z,\[Phi]ref,\[Iota],\[Theta], \[Phi],\[Psi], p1,p2,p3, D11, D12, D13, D22, D23, D33};
	
	\[Psi]Args = {f, fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z, Sequence@@ConstantArray[0, 16]};
	
	(\[ScriptCapitalA]@@aArgs)*Exp[-I*(\[CapitalPsi]@@\[Psi]Args)]
]

hx[f_,fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_,\[Phi]ref_,\[Iota]_] := Module[
	{\[Theta], \[Phi], \[Psi], p1, p2, p3, D11, D12, D13, D22, D23, D33, aArgs, \[Psi]Args},
	
	
	
	{\[Theta], \[Phi], \[Psi], p1, p2, p3} = ConstantArray[0, 6];
	
	{D11, D12, D13, D22, D23, D33} = ConstantArray[0,6];
	
	(*Fplus=0, Fx=1*)
	D12=1/2; 
	
	aArgs = {f, fref, m1,m2, s1x,s1y,s1z,s2x,s2y,s2z,\[Phi]ref,\[Iota],\[Theta], \[Phi],\[Psi], p1,p2,p3, D11, D12, D13, D22, D23, D33};
	
	\[Psi]Args = {f, fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z, Sequence@@ConstantArray[0, 16]};
	
	(\[ScriptCapitalA]@@aArgs)*Exp[-I*(\[CapitalPsi]@@\[Psi]Args)]
]


DeleteObject/@ExternalSessions[]

Clear@python

python = StartExternalSession["Python"];
ExternalEvaluate[python,"
import numpy as np

import jax
import jax.numpy as jnp

from jax import grad, vmap
from functools import partial

from ripplegw.waveforms import IMRPhenomPv2
from ripplegw import get_match_arr, get_eff_pads
from ripplegw import ms_to_Mc_eta
from ripplegw.constants import MSUN, gt
"]


Rh = ExternalFunction[python,"
def Ripple_hp(f, f_ref, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dL, tc,iota, phi_ref):
    f_array = jnp.array(f)

    Mc = (m1 * m2)**(3/5) / (m1 + m2)**(1/5)
    eta = m1 * m2 / ((m1 + m2)**2)

    theta = jnp.array([Mc, eta, s1x, s1y, s1z, s2x, s2y, s2z, dL, tc, phi_ref, iota])

    hp, hc = IMRPhenomPv2.gen_IMRPhenomPv2_hphc(f_array, theta, f_ref)
    
    result = hp + hc
    
    return result.tolist()
"]


Clear@Test

Test := Block[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA,Rippleh, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diffRe,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
		ReMMA, ReRipple, ImMMA, ImRipple, diffIm
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	
	{s1x, s1y, s1z} = RandomReal[{-1,1}, 3];
	{s2x, s2y, s2z} = RandomReal[{-1,1}, 3];
	
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	(*f_,fref_,m1_,m2_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_,\[Phi]ref_,\[Iota]_*)
	MMA = 10^3/(2 Sqrt[5/(64 \[Pi])] ) (
		hp[f, 10, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota]] + 
		hx[f, 10, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota]]
	);
	
	Rippleh = Rh[f,10, m1,m2,s1x,s1y,s1z,s2x,s2y,s2z, 1, 0, \[Iota],\[Phi]Ref];
	
	ReMMA = Re[MMA];
	ReRipple = Re[Rippleh];
	
	ImMMA = Im[MMA];
	ImRipple = Im[Rippleh];
	
	diffRe = RelativeDiff@@{ReMMA, ReRipple};
	
	diffIm = RelativeDiff@@{ImMMA, ImRipple};
	
	
	{ListLinePlot[
		{ReMMA, ReRipple},
		Frame->True,
		PlotLegends->{"MMA", "Ripple"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	],
	ListLinePlot[
		diffRe,
		Frame->True,
		PlotLegends->{"Relative Difference"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	],
	
	ListLinePlot[
		{ImMMA, ImRipple},
		Frame->True,
		PlotLegends->{"MMA", "Ripple"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	],
	ListLinePlot[
		diffIm,
		Frame->True,
		PlotLegends->{"Relative Difference"},
		PlotRange->All,
		ImageSize->Medium,
		GridLines->{{0.014/(G (m2+m1))}, None},
		GridLinesStyle->Directive[Red, 13, Dashed],
		Background->White,
		DataRange->{10, 0.2/(G (m1+m2))}
	]}
	
]


Test


(* ::Chapter:: *)
(*Rosetta Stone*)


(* ::Section:: *)
(*Phase*)


ParentDirectory[NotebookDirectory[], 2]


Module[
	{d, direc = ParentDirectory[NotebookDirectory[], 2], dummy},
	
	dummy = FileNameJoin[{
		direc,
		"/LibraryResources",
		$SystemID,
		"DerivativeRules/IMRPhenomPv2/NRules"
	}];
	
	d =FileNames["phi*", {
		dummy
	}]; 
	
	list = SortBy[(StringReplace[FileBaseName[#], "phi"->""]//ToExpression)&]@d;
	
	(*THIS IS IMPORTANT IT TAKES FROM THE FILE NAME EVERETHING BEFORE "SymDALI/...":*)
	list = FileNameDrop[#, 5]&/@list
];


SetDirectory[NotebookDirectory[]]


Dterms =Block[
	{ Ds},
	
	Ds = Import["Phase_Ds_order_0_to_1.wdx"];
	
	Ds[[All,1]]
];

Dterms = Dterms//.test->\[CapitalPhi]IMR;


phis = Block[
	{Cvariables},
	
	Cvariables= ConstantArray[{Real, 0},25];
	
	Cvariables = Join[{{Real,1}}, Cvariables];

MapThread[
	(#1 -> LF[
		#2, 
		FileBaseName[#2], 
		Cvariables,
		{Real, 1}
	])&,
	{Dterms, list}
]
];


(*put the rules in the correct format*)
phis2 = phis//.{($D[{n__}, \[CapitalPhi]IMR][x__] -> LF[y__]) :>  TagRule[\[CapitalPhi]IMR, $D[{n}, \[CapitalPhi]IMR][x],  LF[y]]};


(*The derivative of the function with respect to any \[Delta]p (at first order) is the function with the \[Delta]p in question replaced by 1*)

$D[{x__}, \[CapitalPhi]IMR][y__]/;Total[{x}[[11;;-1]]] === 1 -> Aux\[CapitalPhi]IMR1;

AuxRule1 = Aux\[CapitalPhi]IMR1[x__][y__] :>  Module[
	{pos = Position[{x}[[11;;-1]], 1], ds, args},
	
	ds = Join[
		{x}[[1;;10]],
		ConstantArray[0, 16]
	];
	
	args = Join[
		{y}[[1;;10]],
		ReplacePart[{y}[[11;;-1]], pos -> 1]
	];
	
	If[
		DeleteDuplicates[ds] === {0}, 
		\[CapitalPhi]IMR@@args,
		$D[ds, \[CapitalPhi]IMR]@@args	
	
	]
]


(*more than 1 derivative in \[Delta]pi is zero because each appears linearly*)
sym1 = $D[{x__}, \[CapitalPhi]IMR]/; Total[{x}[[11;;-1]]] > 1 -> 0


phis3 = Join[
	phis2,
	{
		TagRule[\[CapitalPhi]IMR, $D[{x__}, \[CapitalPhi]IMR]/;Total[{x}[[11;;-1]]] === 1, Aux\[CapitalPhi]IMR1[x]]
	}
];


(* ::Section:: *)
(*Amplitude*)


Module[
	{d, direc = ParentDirectory[NotebookDirectory[], 2],  dummy},
	dummy = FileNameJoin[{direc,"/LibraryResources", $SystemID,"/DerivativeRules/IMRPhenomPv2/NRules"}];
	d =FileNames["a*", {
		dummy
	}];
	
	list\[ScriptCapitalA] = SortBy[(StringReplace[FileBaseName[#], "a"->""]//ToExpression)&]@d;
	
	list\[ScriptCapitalA] = FileNameDrop[#,5]&/@list\[ScriptCapitalA];
];


Dterms\[ScriptCapitalA] =Module[ {Ds = Import["Amp_Ds_order_0_to_1.wdx"]}, Ds[[All,1]]];


Dterms\[ScriptCapitalA] = Dterms\[ScriptCapitalA]//.test->\[ScriptA]IMR;


Dterms\[ScriptCapitalA][[1]]


\[ScriptCapitalA]s = MapThread[
	(#1 -> LF[
		#2, 
		FileBaseName[#2], 
		Join[{{Real, 1}}, ConstantArray[{Real,0}, 23]],
		{Complex, 1}
	])&,
	
	{Dterms\[ScriptCapitalA], list\[ScriptCapitalA]}
];


\[ScriptCapitalA]s2 = \[ScriptCapitalA]s//.{($D[{n__}, \[ScriptA]IMR][x__] -> LF[y__]) :>  TagRule[\[ScriptA]IMR, $D[{n}, \[ScriptA]IMR][x],  LF[y]]};


(* ::Section:: *)
(*Exporting files*)


NRules = <||>;


NRules["\[CapitalPhi]IMR"] = phis3;
NRules["\[ScriptA]IMR"] = \[ScriptCapitalA]s2;
NRules["Aux\[CapitalPhi]IMR1"] = {AuxRule1};


ParentDirectory[NotebookDirectory[],2]


Module[
	{name = ParentDirectory[NotebookDirectory[],2]},
	
	name = FileNameJoin[{name,"LibraryResources/", $SystemID,  "/DerivativeRules/IMRPhenomPv2/NRules/RosettaStone.wdx"}];
	
	Export[name, NRules]
]


SymRules = <||>;

SymRules["\[CapitalPhi]IMR"] = {sym1};


Module[
	{name = ParentDirectory[NotebookDirectory[],2]},
	
	name = FileNameJoin[{name,"LibraryResources/",$SystemID,  "/DerivativeRules/IMRPhenomPv2/SymRules/file.wdx"}];
	
	Export[name, SymRules]
]
