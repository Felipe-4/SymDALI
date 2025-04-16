(* ::Package:: *)

Quit


SetDirectory[NotebookDirectory[]]


PacletDirectoryLoad[ParentDirectory[NotebookDirectory[], 3]];
<<FelipeBarbosa`SymDALI`


(* ::Section::Closed:: *)
(*New*)


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


\[CapitalDelta]t2[p1_, p2_, p3_, \[Theta]_, \[Phi]_] := Module[
	{\[Delta], r, c = UnitConvert["SpeedOfLight"][[1]], detectorposition},
	
	r = FromSphericalCoordinates[{1, \[Theta], \[Phi]}];
	
	detectorposition = {p1,p2,p3};
	
	(detectorposition . r/c)
]


S[\[Theta]_,\[Phi]_, \[Psi]_, cos\[Iota]_, f_ , Dij_, p1_, p2_,p3_] := Module[
	{\[Delta]t, D11, D12,D13,D22,D23,D33, comps},
	
	\[Delta]t = \[CapitalDelta]t2[p1,p2,p3, \[Theta], \[Phi]];
	comps = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]];
	
	{D11,D12,D13,D22,D23,D33} = Extract[Dij, comps];
	
	(
		1/2 Fplus[\[Theta], \[Phi],\[Psi],D11,D12,D13,D22,D23,D33] (1 + cos\[Iota]^2) -
		I cos\[Iota] Fx[\[Theta],\[Phi],\[Psi],D11,D12,D13,D22,D23,D33]
	) Exp[- I 2 \[Pi] f \[Delta]t]
]


NRules = Module[
	{ass = Import["/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/RosettaStone.wdx"]},
	
	
	LibraryFunctionLoad@@@((ass//Values)[[1,All, 2]])
]


functions = Module[
	{ass = Import["/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/RosettaStone.wdx"]},
	
	
	((ass//Values)[[1,All, 1]])
]


NRules[[1]]


test := Module[
	{sin\[Delta], \[Theta] = RandomReal[{0,\[Pi]}], \[Phi] = RandomReal[{-\[Pi], \[Pi]}], cos\[Iota] = RandomReal[{-1,1}], 
	\[Psi]  = RandomReal[{0, 2 \[Pi]}], old, new},
	
	sin\[Delta] = Sin[\[Pi]/2-\[Theta]];
	
	old = NRules[[1]]@@{\[Phi], sin\[Delta], cos\[Iota],\[Psi], 0, Vertex["H1"], DetectorTensor["H1"], 20};
	new = S[\[Theta], \[Phi], \[Psi], cos\[Iota], 20, DetectorTensor["H1"], Vertex["H1"]];
	SetPrecision[new,8] ==  SetPrecision[old,8]
]


DetectorTensor["H1"] = With[
    {nx= {-0.2239, 0.7998, 0.5569}, ny = {-0.9140, 0.0261, -0.4049}},
    (nx\[TensorProduct]nx - ny\[TensorProduct]ny)/2
];

Vertex["H1"] = {-2.16141492636 10^6,  -3.83469517889 10^6   , 4.60035022664 10^6};


(* ::Section::Closed:: *)
(*Calculating Derivatives:*)


expr = HoldForm[$Block[
	{{f1,f2, Fplus, Fcross}},
	
	f1[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_, D33_] = \[ScriptF]1;
	f2[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_] = \[ScriptF]2;
	
	Fplus[f1_, f2_, \[Psi]_] := Cos[2 \[Psi]] f1 + Sin[2 \[Psi]] f2;
	Fcross[f1_, f2_, \[Psi]_] := -Sin[2 \[Psi]] f1 + Cos[2 \[Psi]] f2;
	
	(
		Fplus[f1[\[Theta],\[Phi],D11,D12,D13,D22,D23,D33], f2[\[Theta],\[Phi],D11,D12,D13,D22,D23], \[Psi]] (1+cos\[Iota]^2)/2 - 
		Fcross[f1[\[Theta],\[Phi],D11,D12,D13,D22,D23,D33], f2[\[Theta],\[Phi],D11,D12,D13,D22,D23], \[Psi]] I cos\[Iota]
	)*Exp[-I 2 \[Pi] f \[Delta]t]
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


derivatives = Combinations[{\[Theta], \[Phi], \[Psi], cos\[Iota]}, 3];


derivatives = DeleteElements[derivatives, {{cos\[Iota], cos\[Iota], cos\[Iota]}}];


PrependTo[derivatives, {}];


vars = {f, \[Theta], \[Phi], \[Psi], cos\[Iota], p1, p2, p3, D11, D12, D13, D22, D23, D33};


Clear[S]


expr2 = Hold[Evaluate[{S, vars, derivatives}], Evaluate@expr, "IncludeZeroDerivative"->False]//.HoldForm[x_] :> x;


Ds = DerivativeRules@@expr2;
Export["Detector_Ds_order_0_to_3.wdx", Ds];


(* ::Section:: *)
(*Compiling*)


Ds = Import["Detector_Ds_order_0_to_3.wdx"];


compileThis[x_HoldForm] := Module[
	{dummy, vars},
	vars = {{f,  _Real,  1}, \[Theta], \[Phi], \[Psi], cos\[Iota], p1,p2,p3, D11, D12,D13,D22,D23,D33};
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
	]//.{HoldForm[y_] :> y, us\[Theta]->UnitStep, Rational->Divide, Complex[any_,any2_] :> any + I any2 };
	
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
	LibraryGenerate[#1[[2]], "D" <> ToString[#2//First]]&,
	compiledDs
];


(* ::Section:: *)
(*SymRules and NRules*)


Quit


d = FileNames["D*", {"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/"}] 


Module[
	{d, direc = ParentDirectory[NotebookDirectory[], 2]},
	
	direc = FileNameJoin[{direc, "/LibraryResources",$SystemID, "DerivativeRules/Detectors/NRules/"}];
	d = FileNames["D*", {direc}];
	
	list = SortBy[(StringReplace[FileBaseName[#], "D"->""]//ToExpression)&]@d;
	
	list = FileNameDrop[#,5]&/@list;
]


Dterms =Module[ {Ds = Import["Detector_Ds_order_0_to_3.wdx"]}, Ds[[All,1]] ];


Detectors = Block[
	{Cvariables},
	Cvariables= ConstantArray[{Real, 0}, 13];
	Cvariables = Join[{{Real,1}}, Cvariables];

MapThread[
	(#1 -> LF[
		#2, 
		FileBaseName[#2], 
		Cvariables,
		{Complex, 1}
	])&,
	{Dterms, list}
]
];


Detectors2 = Detectors//.{($D[{n__}, S][x__] -> LF[z__]):> TagRule[S, $D[{n}, S][x], LF[z]]};


(*More than 2 derivatives in cos\[Iota] -> 0  automatically.*)
$D[{n__}, S][y__]/; {n}[[5]] > 2 -> 0;


NRules = <||>;
NRules["S"] = Detectors2;

Module[
	{name = ParentDirectory[NotebookDirectory[], 2]},
	
	name = FileNameJoin[{name, "/LibraryResources",$SystemID, "DerivativeRules/Detectors/NRules/RosettaStone.wdx"}];
	Export[name, NRules]
]


SymRules = <||>;
SymRules["S"] = {$D[{n__}, S]/; {n}[[5]] > 2 -> 0};
Module[
	{name = ParentDirectory[NotebookDirectory[], 2]},
	
	name = FileNameJoin[{name, "/LibraryResources",$SystemID, "DerivativeRules/Detectors/SymRules/file.wdx"}];
	Export[name, SymRules]
]
