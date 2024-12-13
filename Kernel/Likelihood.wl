(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`Likelihood`"];
iGW\[ScriptCapitalL]::usage = "iGW\[ScriptCapitalL][{h__}, detecs_Integer, {vars__List}, {f0_, f1_, \[CapitalDelta]f_}, {PSD__List}, {strainData__}]
calculates the GW Likelihood.";


Begin["Private`"];


(* ::Section:: *)
(*Definitions*)


iGW\[ScriptCapitalL][{h__}, detecs_Integer, {vars__List}, {f0_, f1_, \[CapitalDelta]f_}, {PSD__List}, {strainData__}] := Module[
	{
		detecH, frequency = Range[f0, f1, \[CapitalDelta]f], ObspointsDetecs, DetectorArrays, RemainingHeads, RemainingObsPoints, 
		RemainingArray, \[CapitalDelta]s
	},
	
	(*Firs element in {h__} is always the detector head*)
	detecH = {h}[[1]];
	(*Create the variables for all detectors:*)
	ObspointsDetecs = Join[#, {frequency}]&/@({vars}[[1;;detecs]]);
	
	(*Just map the detector head to all obspoints:*)
	DetectorArrays = detecH@@@ObspointsDetecs; 
	
	Clear[ObspointsDetecs];
	
	(*Make a list of other heads and other ObsPoints:*)
	RemainingObsPoints = Join[#, {frequency}]&/@({vars}[[detecs+1;;-1]]);
	RemainingHeads = {h}[[2;;-1]];
	
	(*Make the array that accounts for the wavefor apart of the detector info:*)
	RemainingArray = Times@@MapThread[
		#1@@#2&,
		{RemainingHeads, RemainingObsPoints}
	];
	Clear[RemainingObsPoints];
	
	\[CapitalDelta]s = MapThread[
		(#1 - RemainingArray*#2)&,
		{{strainData}, DetectorArrays}	
	]//Abs;
	
	- 2 \[CapitalDelta]f Total/@(\[CapitalDelta]s^2/{PSD})
]


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
