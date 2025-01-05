(* ::Package:: *)

Quit


PacletDirectoryLoad["/home/cosmo-ufes/Documentos/GitHub/"];
<<FelipeBarbosa`SymDALI`


(*\[CurlyPhi]: latitude*)
(*\[Lambda]: longitude*)

(*DetectorPosition[h_, \[CurlyPhi]_, \[Lambda]_] := *)Module[
	{R, a =6378137 , b = 6356752.314, X,Y,Z},
	
	R = a^2/Sqrt[a^2 Cos[\[CurlyPhi]]^2 + b^2 (Sin[\[CurlyPhi]]^2) ];
	
	X = (R+h) Cos[\[CurlyPhi]] Cos[\[Lambda]];
	Y = (R+h) Cos[\[CurlyPhi]] Sin[\[Lambda]];
	Z =  (0.993306 R +h) Sin[\[CurlyPhi]];
	
	{X,Y,Z}
]


\[CapitalDelta]t[DetectorPosition[h, \[CurlyPhi], \[Lambda]], \[Alpha], Sin[\[Delta]], gmst]


Module[
	{ \[Phi]  = \[Alpha]-gmst,\[Theta], r, c = UnitConvert["SpeedOfLight"][[1]]},
	\[Theta] = \[Pi]/2 - \[Delta];
	
	r = FromSphericalCoordinates[{1,\[Theta], \[Phi]}];
	
	\[CapitalDelta]texpr = -(DetectorPosition[h, \[CurlyPhi], \[Lambda]] . r/c)
]


Block[
	{Fp, Fc},
	Fp = Sin[\[Zeta]]  Sin[\[Zeta]] (a[\[Delta],\[Alpha], \[Gamma],\[Lambda],\[CurlyPhi], gmst] Cos[2 \[Psi]] + b[\[Delta],\[Alpha], \[Gamma],\[Lambda],\[CurlyPhi], gmst] Sin[2 \[Psi]]);
	Fc = Sin[\[Zeta]] (b[\[Delta],\[Alpha], \[Gamma],\[Lambda],\[CurlyPhi], gmst] Cos[2 \[Psi]] - a[\[Delta],\[Alpha], \[Gamma],\[Lambda],\[CurlyPhi], gmst] Sin[2 \[Psi]]);
	S1expr = (Fp (1+cos\[Iota]^2)/2 - Fc I cos\[Iota]) E^(-I 2 \[Pi] f \[CapitalDelta]texpr)//Simplify
]


Block[
	{b ,\[Delta], \[Phi] = \[Alpha] - gmst,\[Theta]},
	\[Theta] = \[Pi]/2-\[Delta];
	
	b := (Cos[2 \[Gamma]] Sin[\[Lambda]] Cos[\[Theta]] Cos[2 (\[Phi] - \[CurlyPhi])] 
  + (1/4) Sin[2 \[Gamma]] (3 - Cos[2 \[Lambda]]) Cos[\[Theta]] Sin[2 (\[Phi] - \[CurlyPhi])] 
  + Cos[2 \[Gamma]] Cos[\[Lambda]] Sin[\[Theta]] Cos[\[Phi] - \[CurlyPhi]] 
  + (1/2) Sin[2 \[Gamma]] Sin[2 \[Lambda]] Sin[\[Theta]] Sin[\[Phi] - \[CurlyPhi] ])//.{\[Lambda] ->"\[CurlyPhi]", \[CurlyPhi]->"\[Lambda]"}; (*GWFAST uses the opposite def for \[CurlyPhi] and \[Lambda]*)


	
	b//.{"\[CurlyPhi]"->\[CurlyPhi], "\[Lambda]"->\[Lambda]}
]


Block[
	{a ,\[Delta], \[Phi] = \[Alpha] - gmst,\[Theta]},
	\[Theta] = \[Pi]/2-\[Delta];
	
	a := ((1/16) Sin[2 \[Gamma]] (3 - Cos[2 \[Lambda]]) (3 + Cos[2 \[Theta]]) Cos[2 (\[Phi] - \[CurlyPhi])] 
  - (1/4) Cos[2 \[Gamma]] Sin[\[Lambda]] (3 + Cos[2 \[Theta]]) Sin[2 (\[Phi] - \[CurlyPhi])] 
  + (1/4) Sin[2 \[Gamma]] Sin[2 \[Lambda]] Sin[2 \[Theta]] Cos[\[Phi] - \[CurlyPhi]] 
  - (1/2) Cos[2 \[Gamma]] Cos[\[Lambda]] Sin[2 \[Theta]] Sin[\[Phi] - \[CurlyPhi]] 
  + (3/4) Sin[2 \[Gamma]] Cos[\[Lambda]]^2 Sin[\[Theta]]^2)//.{\[Lambda] ->"\[CurlyPhi]", \[CurlyPhi]->"\[Lambda]"}; (*GWFAST uses the opposite def for \[CurlyPhi] and \[Lambda]*)


	
	a//.{"\[CurlyPhi]"->\[CurlyPhi], "\[Lambda]"->\[Lambda]}
]


$blockExpr = Hold[
	{{DetectorPosition, a, b}},
	
	DetectorPosition[h_, \[CurlyPhi]_, \[Lambda]_] := {
		Cos[\[Lambda]] Cos[\[CurlyPhi]] (h+40680631590769/Sqrt[40680631590769 Cos[\[CurlyPhi]]^2+4.040829998154436`*^13 Sin[\[CurlyPhi]]^2]),
		Cos[\[CurlyPhi]] Sin[\[Lambda]] (h+40680631590769/Sqrt[40680631590769 Cos[\[CurlyPhi]]^2+4.040829998154436`*^13 Sin[\[CurlyPhi]]^2]),
		Sin[\[CurlyPhi]] (h+4.040831544290039`*^13/Sqrt[40680631590769 Cos[\[CurlyPhi]]^2+4.040829998154436`*^13 Sin[\[CurlyPhi]]^2])
	};
	
	a[\[Delta]_,\[Alpha]_, \[Gamma]_,\[Lambda]_,\[CurlyPhi]_, gmst_] := (
		3/4 Cos[\[Delta]]^2 Cos[\[CurlyPhi]]^2 Sin[2 \[Gamma]]+
		1/16 (3+Cos[2 (\[Pi]/2-\[Delta])]) Cos[2 (-gmst+\[Alpha]-\[Lambda])] (3-Cos[2 \[CurlyPhi]]) Sin[2 \[Gamma]]+
		1/2 Cos[2 \[Gamma]] Cos[\[CurlyPhi]] Sin[2 (\[Pi]/2-\[Delta])] Sin[gmst-\[Alpha]+\[Lambda]]-
		1/4 Cos[2 \[Gamma]] (3+Cos[2 (\[Pi]/2-\[Delta])]) Sin[2 (-gmst+\[Alpha]-\[Lambda])] Sin[\[CurlyPhi]]+
		1/4 Cos[gmst-\[Alpha]+\[Lambda]] Sin[2 \[Gamma]] Sin[2 (\[Pi]/2-\[Delta])] Sin[2 \[CurlyPhi]]
	);

	b[\[Delta]_,\[Alpha]_, \[Gamma]_,\[Lambda]_,\[CurlyPhi]_, gmst_] := (
		Cos[2 \[Gamma]] Cos[\[Delta]] Cos[gmst-\[Alpha]+\[Lambda]] Cos[\[CurlyPhi]]+
		1/4 (3-Cos[2 \[CurlyPhi]]) Sin[2 \[Gamma]] Sin[\[Delta]] Sin[2 (-gmst+\[Alpha]-\[Lambda])]+
		Cos[2 \[Gamma]] Cos[2 (-gmst+\[Alpha]-\[Lambda])] Sin[\[Delta]] Sin[\[CurlyPhi]]-
		1/2 Cos[\[Delta]] Sin[2 \[Gamma]] Sin[gmst-\[Alpha]+\[Lambda]] Sin[2 \[CurlyPhi]]
	);
	
	x
]/.x -> S1expr;

$blockExpr  = HoldForm[Evaluate@$blockExpr];
$blockExpr = $Block@@@$blockExpr;


DerivativeCombinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


DS = DerivativeCombinations[{\[Alpha], \[Delta], cos\[Iota], \[Psi]}, 3];
PrependTo[DS, {}]


(* ::Text:: *)
(*where the detector position is given in meters following the LIGO convention. All of these considerations lead to the definition: *)


ds  = Hold[
		{S1, {\[Alpha], \[Delta], cos\[Iota], \[Psi], \[CurlyPhi], \[Lambda], h, \[Gamma], \[Zeta], gmst, f}, DS},
		Evaluate[$blockExpr],
		"IncludeZeroDerivative"->False
]/.HoldForm[x_]:>x;


res = DerivativeRules@@ds;


Block[
	{\[Alpha],\[Delta],cos\[Iota],\[Psi],\[CurlyPhi],\[Lambda],h,\[Gamma],\[Zeta],gmst, test, f},
	Thread@Set[
		{\[Alpha],\[Delta],cos\[Iota],\[Psi],\[CurlyPhi],\[Lambda],h,\[Gamma],\[Zeta],gmst}, 
		{3.4461599999999994`,0.4080839999999999`,-0.8272916044289224`,2.3`,0.53342313506`,-1.58430937078`,-6.574`,1.2566370614359172`,\[Pi]/2,1.1870088824`*^9}
	];
	f = Range[20, 1024, 0.125];
	
	test = res[[2,1,2]]; OwnValues[test] = OwnValues[test]//.HoldForm[y_]:>y;
	test;
	(*\[Alpha], \[Delta], cos\[Iota], \[Psi], \[CurlyPhi], \[Lambda], h, \[Gamma], \[Zeta], gmst, {f, _Real, 1}*)
	
	
		(*c[\[Alpha],\[Delta],cos\[Iota],\[Psi],\[CurlyPhi],\[Lambda],h,\[Gamma],\[Zeta],gmst, f]*)	
	

]


compileTHis[x_] := Module[
	{ivars, dummy},
	ivars = {\[Alpha], \[Delta], cos\[Iota], \[Psi], \[CurlyPhi], \[Lambda], h, \[Gamma], \[Zeta], gmst, {f, _Real, 1}};
	
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
		}
	]/.HoldForm[y_]:>y;
	
	Compile@@dummy
]


compiledDs = MapAt[
	compileTHis,
	res[[2]],
	{All, 2}
];


<<CCompilerDriver`
$CCompilerDefaultDirectory = "/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/";
Needs["CCodeGenerator`"]


list = MapIndexed[
	LibraryGenerate[#1[[2]], "D" <> ToString[#2//First]]&,
	compiledDs
];


Dets = Block[{ivars},
	ivars = ConstantArray[{Real, 0}, 10];
	ivars = Join[
		ivars,
		{{Real, 1}}
	];
	
	MapThread[
		(#1 -> LibraryFunction[
			#2, 
			FileBaseName[#2], 
			ivars,
			{Complex, 1}
		])&,
		{compiledDs[[All,1]], list}
	]



];


ass = <||>;
ass["detectors"] = Dets;
Export["/home/cosmo-ufes/Documentos/GitHub/SymDALI/LibraryResources/Linux-x86-64/DerivativeRules/Detectors/NRules/RosettaStone.wdx", ass]
