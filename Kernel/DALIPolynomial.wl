(* ::Package:: *)

(*Needs["maTHEMEatica`"]
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
SetColors[colors]
CreateStyleSheet[]
ApplyStyleSheet[]*)


(* ::Section:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`DALIPolynomial`"];


CompiledPolynomial::usage="CompiledPolynomial[daliList, fiducialPoint ] gives the DALI polynomial associated to the coefficients in daliList expanded around fiducialPoint";
PermutationsNumber::usage = "PermutationsNumber[list] gives the number of all possible permutations of the elements in list.";
StanPolynomial::usage="StanPolynomial[DALITensors, fiducialPoint_]
DALITensors_List: list of DALI tensors as outputed by GWDALICoefficients or DALICoefficients;
fiducialPoint_List: list of vars and their expansion point, for instance {{x1,0}, {x2,1}, ...};
Output: Stan code to be compilled for this model.
Example:

>>Block[{v1,v2,l, list},
	l = Length@SymmetrizedIndependentComponents[{3,3}, Symmetric[All]];
	v1 = RandomReal[{1,2}, 3];
	v2 = RandomReal[{1,2}, l];
	list = {{Flatten[v1\[TensorProduct]v1]}, {Flatten[v2\[TensorProduct]v1], Flatten[v2\[TensorProduct]v2]}};
	
	StanPolynomial[list, {{x,2}, {y,3}, {z,0}}]
]

>>> \"parameters{real x; real y; real z;} model{ target +=-97.59653871634103 + y*(106.6330469515925 + y*(-41.0198425626578 + y*(6.497242207656599 - 0.36233383566103594*y - 1.4622384260186643*z) + (18.25722660084694 - 2.0303585485253195*z)*z) + z*(-70.14591515184834 + (15.362704701457346 - 1.1200882138725219*z)*z)) + x*(87.0145343531812 + y*(-74.63414678790865 + y*(20.105662735773176 - 1.6511359723966994*y - 4.491156243292682*z) + (32.82130742191058 - 3.6044006532284394*z)*z) + x*(-26.671093619701924 + y*(16.055078079792068 - 2.324374825711876*y - 3.536429842538419*z) + x*(3.24359961426922 - 0.13561424012605597*x - 1.0101387721380386*y - 0.7093544931687271*z) + (11.632000203844619 - 1.2672040196842107*z)*z) + z*(-56.686806403151195 + (12.296688555050064 - 0.888176252628384*z)*z)) + z*(84.47371744496604 + z*(-27.389210417613445 + (3.9426922526206813 - 0.21260670717869937*z)*z));}\"
";
SymbolicVector::usage="SymbolicVector[listofLIComponents, head] applies head to the list of Linear Independent components of a tensor listOfLIComponents"
PreprocessDALItensors
TaylorForm


Begin["`Private`"]


(* ::Section:: *)
(*Definitions*)


PermutationsNumber[list_List] := With[

{factorialList = Factorial/@(Tally[list][[All,2]])},
		
		Divide[
		(Length@list)!,
		Times@@factorialList
	]
]

PermutationsNumber[x___] := Throw[$Failed, failTag[PermutationsNumber]]


(*Compile[{x, y}, Module[{h}, h[a_] := 3]]*)


c[i_,j_] := -1/(i! j!)
c[i_,j_]/;i==j := -1/(2 (i!)^2)


(*vector = Range[5]*)


(*Clear@v*)


(*i[dim_, order_] := Module[
	{IndComp, internalFunction},
	
	(IndComp[#] = SymmetrizedIndependentComponents[ConstantArray[dim, #], Symmetric[All]])&/@Range[order];
	(IndComp[#] = ArrayReshape[IndComp[#], {Length[IndComp[#]]*order, 1}])&/@Range[order];
	internalFunction = HoldComplete[
		{matrixVector, iList, k},
		
		matrixVector = Table[
		temp[vector, LIComponents, order,  Length[LIComponents]],
			{LIComponents,  x}
		];
		matrixVector
		
		
		(*Do[
			KroneckerProduct[matrixVector[[j]], matrixVector[[i]]]//Flatten,
			{i, 1, order},
			{j, i, order}
		]*)
		
	]/.{x-> IndComp/@Range[order]}
]*)


(*Clear@vector*)


(*c = SymmetrizedIndependentComponents[ConstantArray[11, 3], Symmetric[All]];*)


(*temp = Hold[{{vector, _Real, 1}, {LIComponents, _Real, 2}, {l}, {L}}, 
	Block[
		{res}, 
		res = Extract[vector, LIComponents];
		res = ArrayReshape[res, {L, l}];
		Times@@@res
	],
	CompilationTarget->"C",
	RuntimeOptions->"Speed"
];

temp = Compile@@temp*)


(*ct = Compile[{{v1, _Real, 1}, {v2, _Real,1}}, 
	Table[
		v1[[i]]*v2[[j]],
		{i, 1 , Length[v1]},
		{j, 1, Length[v2]}
	]//Flatten,
	CompilationTarget->"C",
	RuntimeOptions->"Speed"
]*)


(*f[vector_] = i[11,2]; DownValues[f] = DownValues[f]/.HoldComplete->Module;*)


(*Hold[{{vector, _Real, 1}}, Evaluate[i[11,4]], CompilationTarget->"C",RuntimeOptions->"Speed"]/.HoldComplete->Module;
cf = Compile@@%*)


(*v = RandomReal[{0,1},11]*)


(*cf[v];//AbsoluteTiming*)


(*F=Function[
{Typed[vector,"PackedArray"::["Real64", 1]], Typed[LIComponents,"PackedArray"::["Real64", 2]]}, 
	Block[{res},
	res = Map[
		Part[vector, #]&,
		LIComponents,
		{2}
	];
	
	Map[(Times@@#)&, res]
]
	
	
	]*)


(*FunctionCompile[F]*)


(*?CompilerOptions*)


SymbolicVector[LIComponents_, head_Symbol]/;MatrixQ[LIComponents, NumericQ] := Times@@@Map[
	head,
	LIComponents,
	{2}
]

SymbolicVector[x___] := Throw[$Failed, failTag[SymbolicVector]]


(*Multiply the LI components by their multiplicity and the c[i,j] from Taylor expansion*)

PreprocessDALItensors[DALIlist_List, dimension_Integer] := Module[
	{n = Length@DALIlist, LIComponents, factorialmultiplicity, factorialmultiplicityTensor},
	
	(LIComponents[#] = SymmetrizedIndependentComponents[ConstantArray[dimension, #], Symmetric[All]])&/@Range[n];
	
	(factorialmultiplicity[#] = PermutationsNumber/@LIComponents[#])&/@Range[n];

	factorialmultiplicityTensor = Table[
		
		(Flatten@KroneckerProduct[factorialmultiplicity[i], factorialmultiplicity[j]])*c[i,j],
		{i, 1, n}, {j,1,i}
	
	];
	
	factorialmultiplicityTensor*DALIlist

]

PreprocessDALItensors[x___] := Throw[$Failed, failTag[PreprocessDALItensors]]


(*(LIComponents[#] = SymmetrizedIndependentComponents[ConstantArray[12, #], Symmetric[All]])&/@Range[10];*)


\[CapitalDelta]p[i_] := ToExpression["p"<>ToString[i]]


(*fiducial point = {{x1,value}, {x2,value},...}
output is the polynomial*)

TaylorForm[DALIlist_List, fiducialPoint_?MatrixQ] := Module[

	{n = Length@DALIlist, dimension = Length@fiducialPoint, p, LIComponents, \[CapitalDelta]p, dummyvar, symtensor, dalilist, rules},
	
	(LIComponents[#] = SymmetrizedIndependentComponents[ConstantArray[dimension, #], Symmetric[All]])&/@Range[n];
	
	\[CapitalDelta]p/@Range[dimension]//Evaluate = Unique[ConstantArray[p, dimension]];
	
	dalilist = PreprocessDALItensors[DALIlist, dimension];
	
	dummyvar = 0;
	(*Need to prioritize memory here, so nested loop and old style*)
	
	Do[
		
		Do[
			
			symtensor = (Flatten@KroneckerProduct[#1,#2])&@@{
				SymbolicVector[LIComponents[i], \[CapitalDelta]p],
				SymbolicVector[LIComponents[j], \[CapitalDelta]p]
			}; (*This is a redundancy, but these are fast to calculate and I think memory will be a problem before 
				performance
			*)
			
			dummyvar = (dummyvar + Plus@@(symtensor*dalilist[[i,j]]))//Simplify;,
			{j,1,i}
		],
		
		{i,1,n}
	];
	
	rules = MapThread[
			Rule,
			{\[CapitalDelta]p/@Range[dimension], Subtract@@@fiducialPoint}
		];
	(*Put in the HornerForm:*)
		
		
	HornerForm[dummyvar]//.rules
	
]

TaylorForm[x___] := Throw[$Failed, failTag[TaylorForm]]


StanParser[a_] := Module[
	{res},
	res = HornerForm[a];
	res = Block[{Power=pow}, res];
	CForm[res]//ToString
]

StanParser[x___] := Throw[$Failed, failTag[StanParser]]


MakeStanCode[expression_, vars_List] := Module[
	{iVars = ToString[("real "<> ToString[#])&/@vars], iExpression},
	
	iExpression = StanParser[expression];	
	
	iVars = StringReplace[
		iVars, 
		{
			"," -> ";",
			"}"  ->  ";}"
		}
	];
	
	StringJoin["parameters", iVars, " model{ target +=" , iExpression , ";}"]
]

MakeStanCode[x___] := Throw[$Failed, failTag[MakeStanCode]]


StanPolynomial[DALIOutput_, fiducialPoint_] := Module[
	{TaylorPolynomial, vars},
	
	Catch[
        TaylorPolynomial = TaylorForm[DALIOutput, fiducialPoint];
        vars = fiducialPoint[[All,1]];
        MakeStanCode[TaylorPolynomial, vars]
      ]
]

StanPolynomial[x___] := Throw[$Failed, failTag[StanPolynomial]] 


CompiledPolynomial::fail = "The function failed. The failure occured in function `1`"


CompiledPolynomial[GetDALIOutput_, fiducialPoint_?MatrixQ] := Module[
    {TaylorPolynomial, vars},
    
    
    Catch[
        TaylorPolynomial = TaylorForm[GetDALIOutput, fiducialPoint];
        vars = fiducialPoint[[All,1]];
        Compile[
            Evaluate@vars,
            Evaluate@(TaylorPolynomial),
            CompilationTarget->"C",
            RuntimeOptions->{"CatchMachineOverflow"->False, "CatchMachineIntegerOverflow"->False, "EvaluateSymbolically" ->False},
            RuntimeAttributes->{Listable}, Parallelization->True
        ],
        _failTag,
		(Message[CompiledPolynomial::fail, Style[First@#2, Red]];
      #1)&
  ]
    
]

CompiledPolynomial[x___] := Throw[$Failed, failTag[CompiledPolynomial]]


(* ::Section:: *)
(*Package Footer*)


End[];


EndPackage[];
