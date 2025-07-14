(* ::Package:: *)

Quit


SetDirectory[NotebookDirectory[]]


PacletDirectoryLoad[ParentDirectory[NotebookDirectory[], 3]];
<<FelipeBarbosa`SymDALI`


(* ::Section::Closed:: *)
(*Defs*)


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


Remove[Fx]


Fplus[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	Cos[2 \[Psi]] f1[\[Theta], \[Phi], D11, D12, D13, D22, D23, D33] + 
	Sin[2 \[Psi]] f2[\[Theta], \[Phi], D11, D12, D13, D22, D23]
);

Fx[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	- Sin[2 \[Psi]] f1[\[Theta], \[Phi], D11, D12, D13, D22, D23, D33] +
	Cos[2 \[Psi]] f2[\[Theta], \[Phi], D11, D12, D13, D22, D23]
);


\[CapitalDelta]t2[p1_, p2_, p3_, \[Theta]_, \[Phi]_] := Module[
	{\[Delta], r, c = UnitConvert["SpeedOfLight"][[1]], detectorposition},
	
	r = FromSphericalCoordinates[{1, \[Theta], \[Phi]}];
	
	detectorposition = {p1,p2,p3};
	
	(detectorposition . r/c)
]


FpFc[\[Theta]_, \[Phi]_, \[Psi]_,  cos\[Iota]_, f_ , Dij_, pi_] := Module[
	{\[Delta]t, D11, D12,D13,D22,D23,D33, comps, p1, p2, p3},
	
	{p1, p2, p3} = pi;
	
	{D11,D12,D13,D22,D23,D33} = Dij;
	
	\[Delta]t = \[CapitalDelta]t2[p1,p2,p3, \[Theta], \[Phi]];
	
	{
		Fplus[\[Theta], \[Phi],\[Psi],D11,D12,D13,D22,D23,D33],
		Fx[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33]
	} Exp[- I 2 \[Pi] f \[Delta]t]
]


DetectorTensor["H1"] = Module[
    {nx= {-0.2239, 0.7998, 0.5569}, ny = {-0.9140, 0.0261, -0.4049}, d},
    d= (nx\[TensorProduct]nx - ny\[TensorProduct]ny)/2;
    Extract[d, SymmetrizedIndependentComponents[{3,3}, Symmetric[All]]]
];

Vertex["H1"] = {-2.16141492636 10^6,  -3.83469517889 10^6   , 4.60035022664 10^6};


(* ::Section:: *)
(*Calculating Derivatives:*)


Remove[Fplus]


expr = HoldForm[
$Block[
	{{f1,f2, Fplus, Fcross}},
	
	f1[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_, D33_] := \[ScriptF]1;
	f2[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_] := \[ScriptF]2;
	
	Fplus[f1_, f2_, \[Psi]_] := Cos[2 \[Psi]] f1 + Sin[2 \[Psi]] f2;
	Fcross[f1_, f2_, \[Psi]_] := -Sin[2 \[Psi]] f1 + Cos[2 \[Psi]] f2;
	
	{
		Fplus[f1[\[Theta],\[Phi],D11,D12,D13,D22,D23,D33], f2[\[Theta],\[Phi],D11,D12,D13,D22,D23], \[Psi]], 
		Fcross[f1[\[Theta],\[Phi],D11,D12,D13,D22,D23,D33], f2[\[Theta],\[Phi],D11,D12,D13,D22,D23], \[Psi]]
	}*Exp[-I 2 \[Pi] f \[Delta]t]
]]//.{
	\[ScriptF]1 -> f1[\[Theta],\[Phi],D11,D12,D13,D22,D23,D33], 
	\[ScriptF]2 -> f2[\[Theta],\[Phi],D11,D12,D13,D22,D23], 
	\[Delta]t-> \[CapitalDelta]t2[p1,p2,p3, \[Theta], \[Phi]]
};


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


derivatives = Combinations[{\[Theta], \[Phi], \[Psi]}, 3];


PrependTo[derivatives, {}];


vars = {f, \[Theta], \[Phi], \[Psi], pi, Dij};


Clear[FpFc]


expr2 = HoldForm[Evaluate[{FpFc, vars, derivatives}], Evaluate@expr, "IncludeZeroDerivative"->False]//.HoldForm[x_] :> x;


expr2


Remove["$x*"]
res = DerivativeRules@@expr2;
(*Export["Detector_Ds_order_0_to_3.wdx", Ds];*)


SetDirectory[NotebookDirectory[]]


Export["Detector_Ds_order_0_to_3.wdx", res];


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


Block[
	{
		rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
		(*p1,p2,p3,D11, D12, D13, D22, D23, D33*)
	}, 
	
	ClearAll[TestFpFc];
	
	TestFpFc[f_, \[Theta]_, \[Phi]_, \[Psi]_, p1_, p2_, p3_, D11_,D12_,D13_,D22_,D23_,D33_] = res[[1,2]]//.Join[{G -> g}];

]
DownValues[TestFpFc] = DownValues[TestFpFc]//.HoldForm[x_]:> x;


(* ::Subsection::Closed:: *)
(*Testing the Phase against lal*)


DeleteObject/@ExternalSessions[];
Clear@python


DeleteObject/@ExternalSessions[];

python = StartExternalSession["Python"];

ExternalEvaluate[python, "
from lal import antenna
import numpy as np
"]

lalFpFc = ExternalFunction[python, "

def fpfc(ra, dec, psi):
	res = antenna.AntennaResponse('H1', ra, dec, psi=psi, times=630696086.1999)
	return np.array([res.plus, res.cross])
"]


ExternalEvaluate["Python","
from astropy.time import Time

gps_time = 630696086.1999

# Convert GPS time to UTC and then to GMST in radians
time = Time(gps_time, format='gps', scale='utc')
gmst_radians = time.sidereal_time('mean', 'greenwich').radian

print(\"GMST in radians:\", gmst_radians)


import lal

# GPS time as float (e.g., time of GW150914)
gps_time_float = 630696086.1999

# Step 1: Convert to LIGOTimeGPS object
gps_time = lal.LIGOTimeGPS(gps_time_float)

# Step 2: Compute GMST in radians
gmst_rad = lal.GreenwichMeanSiderealTime(gps_time)

print(\"GMST (radians):\", gmst_rad)

"]


Clear@test
test := Module[
	{r1, r2, \[Theta], \[Phi],\[Alpha], \[Delta], \[Psi], MMA},
	\[Delta] = RandomReal[{-\[Pi],\[Pi]}];
	\[Psi] = RandomReal[{0, \[Pi]}];
	\[Alpha] = RandomReal[{0, 2 \[Pi]}];
	
	\[Theta] = \[Pi]/2-\[Delta]; 
	\[Phi] = \[Alpha] - (-2.821265599576199 10^-6);
	
	MMA = TestFpFc[
		10,
		\[Theta], \[Phi],\[Psi],
		0,0,0, (*H1 DATA AND 0 TIME DELAY*)
		-0.392632395`,-0.07760990999999999`,-0.247384255`,0.31949941499999995`,0.22798825499999997`,0.07309679999999999`
	];
	
	r1 = RelativeDiff[
		MMA[[1]],
		(lalFpFc[\[Alpha], \[Delta], \[Psi]]//Normal)[[1]]
	];
	
	r2 = RelativeDiff[
		MMA[[2]],
		(lalFpFc[\[Alpha], \[Delta], \[Psi]]//Normal)[[2]]
	];
	
	{r1,r2}
	
	
]


test//ScientificForm


(* ::Text:: *)
(*If you get some error from the Python function just run again. I don't know why, but sometimes python external functions generate errors inside ```Do``` or ```Table```*)


list = Table[test, 10^4];


(* ::Text:: *)
(*The points where the transition from negative to positive happens have the highest relative differences, bcs one implementation may go to zero slightly faster than the other*)


list[[All,1]]//Sort//ListPlot[#, PlotRange->All, ScalingFunctions->"Log10"]&
list[[All,2]]//Sort//ListPlot[#, PlotRange->All, ScalingFunctions->"Log10"]&


(* ::Subsection::Closed:: *)
(*Comparing numerical x Symbolic derivatives*)


NGrad//Clear

NGrad[f_, vars_, ni_, nf_] := Module[
	{h = 1. 10^-6, dummy, Point1, Point2, denominator},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = If[vars[[i]]==0, h, vars[[i]] + h vars[[i]]];
		dummy,
		{i, ni, nf}
	];
	
	
	Point2 = ConstantArray[vars, nf-ni+1];
	
	denominator = Table[If[vars[[i]]==0, h, h vars[[i]]], {i, ni, nf}];

	(f@@@Point1 - f@@@Point2)/denominator
]


(*{SymRules, NRules} = DerivativeRulesLoad["IMRPhenomD"];*)


(*(*GRAD WITH COMPILED FUNCTIONS:*)
Block[{\[Delta]s, args}, 
	
	\[Delta]s = ConstantArray[0, 16];
	args = Join[{f,fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z}, \[Delta]s];
	
	ClearAll[TestGrad\[CapitalPsi]];
	
	TestGrad\[CapitalPsi][f_,fref_,m1_,m2_,s1x_,s1y_,s1z_,s2x_,s2y_,s2z_] = Table[
	(NRules["\[CapitalPhi]IMR"])[[i, 3]]@@args,
	{i, 2, 9}
	]

]
DownValues[TestGrad\[CapitalPsi]] = DownValues[TestGrad\[CapitalPsi]]//.HoldForm[x_]:> x;*)


(*BELOW IS USEFUL FOR THE BLOCK FUNCTIONS*)



Block[{}, 
	
	ClearAll[TestGradFpFc];
	
	TestGradFpFc[
		f_, \[Theta]_, \[Phi]_, \[Psi]_,
		p1_, p2_, p3_, D11_,D12_,D13_,D22_,D23_,D33_
	] = res[[2;;4, 2]];

]
DownValues[TestGradFpFc] = DownValues[TestGradFpFc]//.HoldForm[x_]:> x;


Clear@Test
Test := Module[
	{
		\[Theta], \[Phi], \[Psi], f, vars, Symbolic, Numeric, r1, r2
	},
	
	f = RandomReal[{10.,1024}];
	{\[Theta], \[Psi]}=RandomReal[{0, \[Pi]},2];
	\[Phi] = RandomReal[{0, 2 \[Pi]}];
	
	Symbolic = TestGradFpFc[
		f, \[Theta], \[Phi], \[Psi],
		Sequence@@Vertex["H1"],
		Sequence@@DetectorTensor["H1"]
	];
	vars = {f, \[Theta], \[Phi], \[Psi],Sequence@@Vertex["H1"],Sequence@@DetectorTensor["H1"]};
	
	Numeric = NGrad[TestFpFc, vars, 2, 4];
	
	r1 = RelativeDiff@@{Symbolic[[All, 1]], Numeric[[All,1]]}; (*diff between the Fp derivatives*)
	r2 = RelativeDiff@@{Symbolic[[All, 2]], Numeric[[All,2]]}; (*diff between the Fc derivatives*)
	
	Max/@{r1,r2} (*get only the max diffs*)
	

	
]


Test//ScientificForm


(*there are errors at \[Eta] = 0.25 bcs numerical derivative will try evaluation at 0.25 + 10^-6 0.25*)


l = Table[Test, {5 10^4}];


ListPlot[l[[All,1]]//Sort, PlotRange->All, ScalingFunctions->"Log10"]
ListPlot[l[[All,2]]//Sort, PlotRange->All, ScalingFunctions->"Log10"]


(* ::Section:: *)
(*Compiling*)


Ds = Import["Detector_Ds_order_0_to_3.wdx"];


Ds[[1]]


compileThis[x_HoldForm] := Module[
	{dummy, vars},
	vars = {{f,  _Real,  1}, \[Theta], \[Phi], \[Psi], {pi, _Real, 1}, {Dij, _Real, 1}};
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
	]//.{
		HoldForm[y_] :> Block[
			{p1,p2,p3, D11, D12,D13,D22,D23,D33}, 
			{p1,p2,p3} = pi; 
			{D11, D12, D13, D22, D23, D33} = Dij;
			y
		], 
		us\[Theta]->UnitStep, 
		Rational->Divide, 
		Complex[any_,any2_] :> any + I any2 };
	
	Compile@@dummy
]


<<CompiledFunctionTools`


compiledDs = MapAt[
	compileThis,
	Ds, 
	{All, 2}
];


<<CCompilerDriver`


FileNameJoin[{
	ParentDirectory[NotebookDirectory[], 2],
	"/LibraryResources/",
	$SystemID, 
	"/DerivativeRules/Detectors/NRules/"

}]


$CCompilerDefaultDirectory = FileNameJoin[{
	ParentDirectory[NotebookDirectory[], 2],
	"/LibraryResources/",
	$SystemID, 
	"/DerivativeRules/Detectors/NRules/"

}]


Needs["CCodeGenerator`"]


MapIndexed[
	LibraryGenerate[#1[[2]], "D" <> ToString[#2//First], {
	"SystemCompileOptions" -> "-march=native -O3 -finline-functions -funroll-loops -flto -ftree-vectorize -fno-fast-math -fPIC"}]&,
	compiledDs
];


(* ::Section:: *)
(*SymRules and NRules*)


Module[
	{d, direc = ParentDirectory[NotebookDirectory[], 2]},
	
	direc = FileNameJoin[{direc, "/LibraryResources",$SystemID, "DerivativeRules/Detectors/NRules/"}];
	d = FileNames["D*", {direc}];
	
	list = SortBy[(StringReplace[FileBaseName[#], "D"->""]//ToExpression)&]@d;
	
	list = FileNameDrop[#,5]&/@list;
]


Dterms =Module[ {Ds = Import["Detector_Ds_order_0_to_3.wdx"]}, Ds[[All,1]] ];


{{f,  _Real,  1}, \[Theta], \[Phi], \[Psi], {pi, _Real, 1}, {Dij, _Real, 1}};


Detectors = Block[
	{Cvariables},(*f, \[Theta], \[Phi], \[Psi], pi, Dij*)
	Cvariables= {{Real,1}, {Real,0}, {Real,0}, {Real,0}, {Real,1}, {Real,1}};
	

MapThread[
	(#1 -> LF[
		#2, 
		FileBaseName[#2], 
		Cvariables,
		{Complex, 2}
	])&,
	{Dterms, list}
]
];


Detectors2 = Detectors//.{($D[{n__}, S][x__] -> LF[z__]):> TagRule[S, $D[{n}, S][x], LF[z]]};


NRules = <||>;
NRules["FpFc"] = Detectors2;

Module[
	{name = ParentDirectory[NotebookDirectory[], 2]},
	
	name = FileNameJoin[{name, "/LibraryResources",$SystemID, "DerivativeRules/Detectors/NRules/RosettaStone.wdx"}];
	Export[name, NRules]
]
