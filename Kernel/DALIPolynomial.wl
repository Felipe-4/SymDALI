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
(*SymbolicVector::usage="SymbolicVector[listofLIComponents, head] applies head to the list of Linear Independent components of a tensor listOfLIComponents"
PreprocessDALItensors
TaylorForm*)


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


c[i_,j_] := -1/(i! j!)
c[i_,j_]/;i==j := -1/(2 (i!)^2)


SymbolicVector[LIComponents_, head_Symbol]/;MatrixQ[LIComponents, NumericQ] := Times@@@Map[
	head,
	LIComponents,
	{2}
]

SymbolicVector[x___] := Throw[$Failed, failTag[SymbolicVector]]


(*Multiply the LI components by their multiplicity and the c[i,j] from Taylor expansion*)

PreprocessDALItensors[DALIlist_List, dimension_Integer]/; MatrixQ[DALIlist, NumericQ] := Module[
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


(*fiducial point = {{x1,value}, {x2,value},...}
output is the polynomial*)

TaylorForm[DALIlist_List, fiducialPoint_?MatrixQ]/;MatrixQ[DALIlist, NumericQ] := Module[

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
		
	dummyvar//.rules
	
]

TaylorForm[x___] := Throw[$Failed, failTag[TaylorForm]]


CompiledPolynomial::fail = "The function failed. The failure occured in function `1`"


CompiledPolynomial[GetDALIOutput_?MatrixQ, fiducialPoint_?MatrixQ] := Module[
    {TaylorPolynomial, vars},
    
    
    Catch[
        TaylorPolynomial = TaylorForm[GetDALIOutput, fiducialPoint];
        vars = fiducialPoint[[All,1]];
        Compile[
            Evaluate@vars,
            Evaluate@(HornerForm[TaylorPolynomial]//N),
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
