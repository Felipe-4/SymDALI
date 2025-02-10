(* ::Package:: *)

PacletDirectoryLoad["/home/cosmo-ufes/Documentos/GitHub/"];
<<FelipeBarbosa`SymDALI`


(* ::Section::Closed:: *)
(*old*)


polarizationTensors[\[Alpha]_, sin\[Delta]_, GMST_] := Module[
	{\[Theta], \[Delta], \[Phi] = \[Alpha] - GMST, eplus, ecross, R1, R2, T, D, eplusGeoFrame, ecrossGeoFrame},
	(*Polarization frame tensors*)
	\[Theta]= \[Pi]/2 - \[Delta];

    eplus = {{1,0,0}, {0,-1,0}, {0,0,0}};
    ecross = {{0,1,0}, {1,0,0}, {0,0,0}};
    
    R1 = List[
        {-Sin[\[Phi]], Cos[\[Phi]], 0}, (*\hat{\[Phi]} decomposed in \hat{e}_i*)
        {Cos[\[Theta]] Cos[\[Phi]], Cos[\[Theta]] Sin[\[Phi]], -Sin[\[Theta]]}, (*\hat{\[Theta]} decomposed in \hat{e}_i*)
        {Sin[\[Theta]] Cos[\[Phi]], Sin[\[Theta]] Sin[\[Phi]], Cos[\[Theta]]} (*\hat{r}decomposed in \hat{e}_i*)
    ];
	
	R1 = (R1//.{Sin[\[Delta]] -> sin\[Delta], Cos[\[Delta]]-> Sqrt[1- sin\[Delta]^2]})//FullSimplify;


    R2 = RotationMatrix[-\[Psi], {0,0,1}];
    T = R2 . R1//Simplify;
    
	eplusGeoFrame = (T\[Transpose] . eplus . T);
    ecrossGeoFrame = (T\[Transpose] . ecross . T);
    
    {eplusGeoFrame, ecrossGeoFrame}

]


(* ::Text:: *)
(*Meaning that our Function can be written as *)


\[CapitalDelta]t[detectorposition_, \[Alpha]_, sin\[Delta]_, gmst_] := Module[
	{\[Delta], \[Phi] = \[Alpha]-gmst,\[Theta], r, c = UnitConvert["SpeedOfLight"][[1]]},
	
	\[Theta] = \[Pi]/2 - \[Delta];
	
	r = FromSphericalCoordinates[{1, \[Theta], \[Phi]}]//.{Sin[\[Delta]] -> sin\[Delta], Cos[\[Delta]] -> Sqrt[1- sin\[Delta]^2]};
	
	-(detectorposition . r/c)
]


(*make some trace operator that commutes with D:*)
Unprotect@Tr
Tr/: D[Tr[x_], y___] := Tr[D[x,y]];
Tr[0] := 0 (*This shoudn't happen but the derivatives of order 3 or higher in cos\[Iota] give you 0 straightfoward, 
instead of Dot[something, 0, somethingelse] in any case you can do this with trace so that it is accounted for.
*)

{SymRules["detectors"], NRules["detectors"]} = Module[
	{\[Delta]t, \[Delta], \[Theta], \[Phi] = \[Alpha]-gmst, matrix1, matrix2, delta},

	\[Theta] = \[Pi]/2 - \[Delta];

	delta = \[CapitalDelta]t[detectorPosition, \[Alpha], sin\[Delta], gmst];
	
	matrix1 = List[
        {-Sin[\[Phi]], Cos[\[Phi]], 0}, 
        {Cos[\[Theta]] Cos[\[Phi]], Cos[\[Theta]] Sin[\[Phi]], -Sin[\[Theta]]},
        {Sin[\[Theta]] Cos[\[Phi]], Sin[\[Theta]] Sin[\[Phi]], Cos[\[Theta]]} 
    ];
	
	matrix1 = (matrix1//.{Sin[\[Delta]] -> sin\[Delta], Cos[\[Delta]]-> Sqrt[1- sin\[Delta]^2]})//FullSimplify;

	matrix2 = RotationMatrix[-\[Psi], {0,0, 1 }];

	Unevaluated@DerivativeRules[
		{S1, {\[Alpha], sin\[Delta], cos\[Iota], \[Psi]}, 1},

		$Block[
		{{m1, m2,m3,m4, \[CapitalDelta]}},

		
		(*NOTE THAT WE ONLY DECLARE VARIABLES THAT APPEAR IN DERIVATIVES*)
		\[CapitalDelta][\[Alpha]_, sin\[Delta]_] := DELTA;
		m1[\[Alpha]_, sin\[Delta]_] := MATRIX1;
		m2[\[Psi]_] := MATRIX2;
		m3[\[Alpha]_, sin\[Delta]_] := TMATRIX1;
		m4[\[Psi]_] := TMATRIX2;

		(*T = m2.m1:*)
		Dot[
			(m3[\[Alpha], sin\[Delta]] . m4[\[Psi]] . (1/2 (1+cos\[Iota]^2) eplus - I cos\[Iota] ecross) . m2[\[Psi]] . m1[\[Alpha], sin\[Delta]]) Exp[-I 2 \[Pi] f \[CapitalDelta][\[Alpha], sin\[Delta]]],
			Transpose[detectorTensor]
		]//Tr
	],
	"KeepDefs" -> {eplus -> {{1,0,0}, {0,-1,0}, {0,0,0}}, ecross -> {{0,1,0}, {1,0,0}, {0,0,0}}}
]//.{DELTA ->delta, MATRIX1 -> matrix1, MATRIX2 -> matrix2, TMATRIX1-> Transpose[matrix1], TMATRIX2 -> Transpose[matrix2]}

]; (*Insert the correct trace operator*)


newSymRules = <||>; newNRules = <||>;

{newSymRules["detectors"], newNRules["detectors"]} = Module[
	{addVars,headsNR, headsSR},

	(*Make a function to change the heads:*)
	addVars[x_]/; (AtomQ[x[[0]]]) := Join[
		x, 
		x[[0]][gmst_, detectorPosition_, detectorTensor_, f_]
	];
	
	addVars[x_] := Module[
		{newHead = x[[0]]}, 
		newHead[[1]] = Join[newHead[[1]], {0,0,0,0}];
		Join[
			newHead@@x, 
			newHead[gmst_, detectorPosition_, detectorTensor_, f_]
		]
	];
	
	headsNR = addVars/@(NRules["detectors"][[All, 1]]);
	headsSR = addVars/@(SymRules["detectors"][[All, 1]]);

	(*Return the new lists of rules:*)
	{
		Thread@Rule[headsSR, SymRules["detectors"][[All,2]]],
		Thread@Rule[headsNR, NRules["detectors"][[All,2]]]
	}
];


newNRules["detectors"] = newNRules["detectors"]//.Tr[x_] :> Sum[x[[i,i]], {i,1,3}];


newNRules["detectors"][[1,2]]


compileTHis[x_] := Module[
	{ivars, dummy},
	ivars = {\[Alpha], sin\[Delta], cos\[Iota], \[Psi], gmst,{detectorPosition, _Real, 1},{detectorTensor, _Real, 2}, f};
	
	dummy = Hold[
		Evaluate[ivars],
		Evaluate[x],
		(*CompilationTarget->"C",*)
		RuntimeOptions->{
			"CatchMachineOverflow"->False,
			"CatchMachineIntegerOverflow"->False,
			"EvaluateSymbolically"->False,
			"RuntimeErrorHandler"->None,
			"WarningMessages"->True
		},
		RuntimeAttributes->{Listable},
		Parallelization->True
	]/.HoldForm[y_]:>y;
	
	Compile@@dummy
]


{H1Tensor, H1pos} = With[{
	nx = {-0.22389266154, 0.79983062746,0.55690487831}, ny = {-0.91397818574, 0.02609403989, -0.40492342125},
	pos = {-2.16141492636 10^6, -3.83469517889 10^6,4.60035022664 10^6}
},
	{
		0.5 (nx\[TensorProduct]nx - ny\[TensorProduct]ny),
		pos
	}

];


compiledDs = MapAt[
	compileTHis,
	newNRules["detectors"],
	{All, 2}
];


<<CCompilerDriver`
$CCompilerDefaultDirectory = "/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/";
Needs["CCodeGenerator`"]


list = MapIndexed[
	LibraryGenerate[#1[[2]], "D" <> ToString[#2//First]]&,
	compiledDs
];


newNRules["detectors"][[1,1]]


Dets = Block[{ivars},
	ivars = ConstantArray[{Real, 0}, 5];
	ivars = Join[
		ivars,
		{{Real, 1}, {Real, 2}, {Real, 0}}
	];
	Echo[ivars];
	
	MapThread[
		(#1 -> lf[
			#2, 
			FileBaseName[#2], 
			ivars,
			Complex
		])&,
		{compiledDs[[All,1]], list}
	]



];


ass = <||>;
ass["detectors"] = Dets;
Export["/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/RosettaStone.wdx", ass]


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
	
	-(detectorposition . r/c)
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


(* ::Section::Closed:: *)
(*Compiling*)


Ds = Import["Detector_Ds_order_0_to_3.wdx"];


Ds[[1,1]]


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


$CCompilerDefaultDirectory = "/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/";


Needs["CCodeGenerator`"]


MapIndexed[
	LibraryGenerate[#1[[2]], "D" <> ToString[#2//First]]&,
	compiledDs
];


(* ::Section:: *)
(*SymRules and NRules*)


With[
	{d =FileNames["D*", {"/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/"}] },
	list = SortBy[(StringReplace[FileBaseName[#], "D"->""]//ToExpression)&]@d
];


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
Export["/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/RosettaStone.wdx", NRules]


SymRules = <||>;
SymRules["S"] = {$D[{n__}, S]/; {n}[[5]] > 2 -> 0};
Export["/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/SymRules/file.wdx", SymRules];
