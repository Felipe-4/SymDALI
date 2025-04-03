(* ::Package:: *)

h[f_, fref_, M_,\[Eta]_,\[Chi]1_,\[Chi]2_, \[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_] = Block[
	{\[Omega], G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],\[Omega]ref},
	\[Omega] = f M G; 
	\[Omega]ref = fref M G;
	
	\[ScriptA]IMR[f,M, \[Eta], \[Chi]1, \[Chi]2] Exp[-I \[CapitalPhi]IMR[\[Omega], \[Eta],\[Chi]1,\[Chi]2, \[Omega]ref, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4]]
];


DownValues[h]


ParentDirectory[NotebookDirectory[], 2]


Module[
	{direc = ParentDirectory[NotebookDirectory[], 2]},
	
	Export[
		FileNameJoin[{direc, "LibraryResources", $SystemID, "DerivativeRules/IMRPhenomD/Defs.wdx"}],
		DownValues[h][[1]]//.RuleDelayed->Rule
	]
]



(*Export[
	"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/IMRPhenomD/Defs.wdx",
	DownValues[h][[1]]//.RuleDelayed->Rule
]*)
