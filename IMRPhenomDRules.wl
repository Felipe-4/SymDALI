(* ::Package:: *)

(*Basic stuff:*)
$HistoryLength=0;

PacletDirectoryLoad["/Users/felipe/Documents/GitHub"];
<<FelipeBarbosa`SymDALI`

vectorDefs = Module[
	{DALIDir, docDir},
	DALIDir = FindFile["FelipeBarbosa`SymDALI`"]//FileNameDrop//ParentDirectory;
	docDir = FileNameJoin[DALIDir, "Documentation/English/Tutorials"];
	Import@FileNameJoin[docDir, "IMRPhenomDComponentsDefinition.wdx"]
];

SymRules = <||>;
NRules = <||>;


(* ::Section:: *)
(*Ins-Phase*)


(*Find vectors for Inspiral Phase and their variables:*)

Position[vectorDefs[[All, 1,All, 0]]/.HoldPattern->Identity, #]&/@{InsVecPhase, \[Omega]InsVecPhase}

vectorDefs[[1;;2, 1]]
DownValues[InsVecPhase] = {vectorDefs[[1]]};
DownValues[\[Omega]InsVecPhase] = {vectorDefs[[2]]};


(*We have to include the \[Delta]\[CurlyPhi]'s in the vector, following arXiv:1903.04467v3*)
extraInsVecPhase[\[Eta]_, \[Chi]1_, \[Chi]2_] = Module[
	{\[Chi]s, \[Chi]a, \[CurlyPhi]4, \[CurlyPhi]5l, \[CurlyPhi]6, \[CurlyPhi]6l, \[CurlyPhi]7},
	
	\[Chi]a = (\[Chi]1-\[Chi]2)/2; \[Chi]s =( \[Chi]1+\[Chi]2)/2;
	
	(*following the idea of \[CurlyPhi] -> (1 + \[Delta]\[CurlyPhi])\[CurlyPhi] we check that -
	InsVecPhase[\[Eta], \[Chi]1, \[Chi]2][[#]]&/@{1,3,4} === 3/(128 \[Eta])  {\[CurlyPhi]0, \[CurlyPhi]2, \[CurlyPhi]3};
	*)
	\[CurlyPhi]4 = 15293365/508032+(27145 \[Eta])/504+(3085 \[Eta]^2)/72+(-(405/8)+200 \[Eta]) \[Chi]a^2-405/4 Sqrt[1-4 \[Eta]] \[Chi]a \[Chi]s+(-(405/8)+(5 \[Eta])/2) \[Chi]s^2;
	\[CurlyPhi]5l = ((38645 \[Pi])/756-(65 \[Pi] \[Eta])/9+(-(732985/2268)-(140 \[Eta])/9) Sqrt[1-4 \[Eta]] \[Chi]a+(-(732985/2268)+(24260 \[Eta])/81+(340 \[Eta]^2)/9) \[Chi]s);
	\[CurlyPhi]6 = 11583231236531/4694215680-(6848 EulerGamma)/21-(640 \[Pi]^2)/3+(-(15737765635/3048192)+(2255 \[Pi]^2)/12) \[Eta]+(76055 \[Eta]^2)/1728-(127825 \[Eta]^3)/1296+2270/3 \[Pi] Sqrt[1-4 \[Eta]] \[Chi]a+((2270 \[Pi])/3-520 \[Pi] \[Eta]) \[Chi]s-6848/63 Log[64];
	\[CurlyPhi]6l = -(6848/63);
	\[CurlyPhi]7 = (77096675 \[Pi])/254016+(378515 \[Pi] \[Eta])/1512-(74045 \[Pi] \[Eta]^2)/756+Sqrt[1-4 \[Eta]] (-(25150083775/3048192)+(26804935 \[Eta])/6048-(1985 \[Eta]^2)/48) \[Chi]a+(-(25150083775/3048192)+(10566655595 \[Eta])/762048-(1042165 \[Eta]^2)/3024+(5345 \[Eta]^3)/36) \[Chi]s;
	
	Join[ (*We use this vector to make 
		Dot[extraphase[\[Eta],\[Chi]1, \[Chi]2]*deltas[\[Delta]\[CurlyPhi]__],  extra\[Omega]vec[\[Omega]]]
		*)
		{3/(128 \[Eta])},
		{InsVecPhase[\[Eta], \[Chi]1, \[Chi]2][[1]]},
		{3/(128 \[Eta])},
		InsVecPhase[\[Eta], \[Chi]1, \[Chi]2][[#]]&/@{3,4},
		3*{\[CurlyPhi]4, \[CurlyPhi]5l, \[CurlyPhi]6, \[CurlyPhi]6l, \[CurlyPhi]7}/(128 \[Eta])
	]
];

InsVecDeltas[\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] = {\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7};

extra\[Omega]InsVecPhase[\[Omega]_] =Join[
	(\[Pi] \[Omega])^((#-5)/3)&/@{-2, 0,1,2,3,4},
	{Log[\[Pi] \[Omega]], (\[Pi] \[Omega])^(1/3), (\[Pi] \[Omega])^(1/3) Log[\[Pi] \[Omega]], (\[Pi] \[Omega])^(2/3)} (*{5l,6, 6l, 7}*)
];

Extra\[Phi]Ins[\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_, \[Omega]_, \[Eta]_, \[Chi]1_, \[Chi]2_] = Dot[
	UextraInsVecPhase[\[Eta], \[Chi]1, \[Chi]2]*UInsVecDeltas[\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7],
	Uextra\[Omega]InsVecPhase[\[Omega]]	
];


(*Extend Defs with -\[Pi]/4:*)
InsVecPhase[\[Eta]_, \[Chi]1_, \[Chi]2_] = Prepend[InsVecPhase[\[Eta],\[Chi]1, \[Chi]2], -\[Pi]/4];
\[Omega]InsVecPhase[\[Omega]_] = Prepend[\[Omega]InsVecPhase[\[Omega]], 1];


(*U \[Congruent] undefined. You will define UInsVecPhase to be InsVecPhase latter and so on*)
\[Phi]Ins[\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_, \[Omega]_, \[Eta]_, \[Chi]1_, \[Chi]2_] = (
	UInsVecPhase[\[Eta], \[Chi]1, \[Chi]2] . U\[Omega]InsVecPhase[\[Omega]] + 
	Extra\[Phi]Ins[\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7, \[Omega], \[Eta], \[Chi]1, \[Chi]2]
);


(* ::Section:: *)
(*Int-Phase*)


(*Find vectors for Intermediate Phase and their variables:*)
Position[vectorDefs[[All, 1, All, 0]] /. HoldPattern -> Identity, #] & /@ {IntVecPhase, \[Omega]IntVecPhase}

vectorDefs[[3 ;; 4, 1]]

DownValues[IntVecPhase] = {vectorDefs[[3]]};
DownValues[\[Omega]IntVecPhase] = {vectorDefs[[4]]};


extraIntVecPhase[\[Eta]_, \[Chi]Pn_] = IntVecPhase[[3;;4]];
IntVecDeltas[\[Delta]\[Beta]2_, \[Delta]\[Beta]3_] = {\[Delta]\[Beta]2, \[Delta]\[Beta]3};

extra\[Phi]Int[\[Omega]_, \[Eta]_, \[Chi]PN_, \[Delta]\[Beta]2_, \[Delta]\[Beta]3_] = Dot[
	UextraIntVecPhase[\[Eta], \[Chi]PN]*UIntVecDeltas[\[Delta]\[Beta]2, \[Delta]\[Beta]3],
	{Log[\[Omega]], -1/3 \[Omega]^-3}
];


\[Phi]Int[\[Omega]_, \[Eta]_, \[Chi]PN_, \[Delta]\[Beta]2_, \[Delta]\[Beta]3_] = UIntVecPhase[\[Eta], \[Chi]PN] . U\[Omega]IntVecPhase[\[Omega]] + extra\[Phi]Int[\[Omega],\[Eta],\[Chi]PN,\[Delta]\[Beta]2,\[Delta]\[Beta]3];


(* ::Section:: *)
(*Ringdown and Damping -Phase*)


(* ::Section:: *)
(*MR-Phase*)
