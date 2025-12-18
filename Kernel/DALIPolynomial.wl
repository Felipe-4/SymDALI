(* ::Package:: *)

(* ::Section::Closed:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`DALIPolynomial`"];


Unprotect@TaylorForm;

TaylorForm//ClearAll

TaylorForm::usage="TaylorForm[DALITensors_List]
DALITensors: the exact List that comes out from ```DALITensors```

returns the TaylorForm of the derivative expansion, that is, a list of tensors 
to be contracted with {\!\(\*SuperscriptBox[\(\[CapitalDelta]p\), \(2\)]\), \!\(\*SuperscriptBox[\(\[CapitalDelta]p\), \(3\)]\), ..., \!\(\*SuperscriptBox[\(\[CapitalDelta]p\), \(2  n\)]\)}, where ```n``` is the derivative order of 
the expansion.

obs.:The tensors are represented by their LI components, supllemented of their
multiplicities, ready to be contracted with the LI components of \!\(\*SuperscriptBox[\(\[CapitalDelta]p\), \(i\)]\).";


Unprotect@ProcessDALITensors;

ProcessDALITensors//ClearAll

ProcessDALITensors::usage="ProcessDALITensors[DALI_List]
DALI_List: the exact List that comes out from ```DALITensors```
returns the list with only the LI components of the tensors and the multiplicity of each component, to be
contracted with \!\(\*SuperscriptBox[\(\[CapitalDelta]p\), \(i\)]\)s directly.";


Begin["`Private`"]


(* ::Section:: *)
(*Definitions*)


(* ::Subsection::Closed:: *)
(*Old usage messages*)


(*CompiledPolynomial::usage="CompiledPolynomial[daliList, fiducialPoint ] gives the DALI polynomial associated to the coefficients in daliList expanded around fiducialPoint";
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
TaylorForm*)


(* ::Subsection::Closed:: *)
(*PermutationsNumber*)


ClearAll@PermutationsNumber


PermutationsNumber::usage=" PermutationNumber[{x,y,...,z}]
returns the number of distinct permutations of (x,y,..z)

>>PermutationsNumber[{1,2,3}] ==3!
>>True

>>PermutationsNumber[{1,2,2}]== 3
>>True
";


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


(* ::Subsection::Closed:: *)
(*DALIComponents*)


DALIComponents//ClearAll


DALIComponents::usage="DALIComponents[dim, order]
generates a list of the LI components of the gradients used to build the DALI structure.
>>DALIComponents[3, 2]
>>{
	{{1},{2},{3}},
	{{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
}";


DALIComponents[dim_, order_] := SymmetrizedIndependentComponents[
	ConstantArray[dim, #],
	Symmetric[All]
]&/@Range[order]


(* ::Subsection::Closed:: *)
(*ProcessDALITensors*)


ProcessDALITensors[DALITensors_List] := Module[
	{
		order = Length[DALITensors], dim = Sqrt[DALITensors[[1,1]]//Length], 
		iTensorList, vectors, SymbolicDALI, highOrderLI, HOdim, HOLI, HOM, HOMultiplicity, HOLITensors 
	},
	
	(*First you add the c[i,j] contribution in front of the tensors*)
	iTensorList = Table[
		c[i,j]DALITensors[[i,j]],
		{i,order},
		{j,i}
	];
	
	(*Make a new DALI_List whose tensor components will match the multiplicity of each component in the original List*)
	vectors = DALIComponents[dim, order];(*gradients LI components*)
	vectors = Map[PermutationsNumber, vectors, {2}]; (*adding multiplicities of each component*)
	
	(*remember the order: {
		{(1,1)}, {(1,2), (2,2)}, {(1,3), (2,3), (3,3)}, ...
	}*)
	SymbolicDALI = Do[
		Sow[#, j]&@(Flatten[vectors[[j]]\[TensorProduct]vectors[[i]]]),
		{i, 1, order},
		{j, i, order}
	]//Reap//Last;     
	
	
	
	
	(*add these multiplicities to iTensorList:*)
	iTensorList = SymbolicDALI*iTensorList;
	
	(*
		Now we just extract the LI components of the high order terms: (1,1), (2,2), (3,3)
		and include their multiplicities:
	*)
	
	
	HOLITensors = Table[
		HOdim = Sqrt[(iTensorList[[i,i]]//Length)]; (*higher order dimensions*)
		HOLI = SymmetrizedIndependentComponents[{HOdim, HOdim}, Symmetric[All]]; (*HO LIs*)
		HOMultiplicity = PermutationsNumber/@HOLI; (*Multiplicities of these LI components*)
		HOM = ArrayReshape[iTensorList[[i,i]], {HOdim, HOdim}]; (*Matrix Format*)
		Extract[HOM, HOLI]*HOMultiplicity (*extract the LI components and add the multiplicityes*), 
		{i, 1, order}
	];
	
	Do[
		iTensorList[[i,i]] = HOLITensors[[i]], (*modify iTensorList with the LI components and multiplicities*)
		{i,1,order}
	];
	
	iTensorList
	
]


ProcessDALITensors//Protect;


(* ::Subsection::Closed:: *)
(*SymbolicVector*)


(*SymbolicVector[LIComponents_, head_Symbol]/;MatrixQ[LIComponents, NumericQ] := Times@@@Map[
	head,
	LIComponents,
	{2}
]

SymbolicVector[x___] := Throw[$Failed, failTag[SymbolicVector]]*)


(* ::Subsection::Closed:: *)
(*PreprocessDALITensors*)


(*(*Multiply the LI components by their multiplicity and the c[i,j] from Taylor expansion*)

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

PreprocessDALItensors[x___] := Throw[$Failed, failTag[PreprocessDALItensors]]*)


(*\[CapitalDelta]p[i_] := ToExpression["p"<>ToString[i]]*)


(* ::Subsection::Closed:: *)
(*TaylorForm*)


ClearAll@DALISymmetry

DALISymmetry[1, 1] = Symmetric[{}];

DALISymmetry[2, 1] = Symmetric[{1,2}];
DALISymmetry[2, 2] = {Symmetric[{1,2}], Symmetric[{3,4}]};

DALISymmetry[3,1] = Symmetric[{1,2,3}];
DALISymmetry[3,2] = {Symmetric[{1,2,3}], Symmetric[{4,5}]};
DALISymmetry[3,3] = {Symmetric[{1,2,3}], Symmetric[{4,5,6}]};

DALISymmetry[4,1] = Symmetric[{1,2,3,4}];
DALISymmetry[4,2] = {Symmetric[{1,2,3,4}], Symmetric[{5,6}]};
DALISymmetry[4,3] = {Symmetric[{1,2,3,4}], Symmetric[{5,6,7}]};
DALISymmetry[4,4] = {Symmetric[{1,2,3,4}], Symmetric[{5,6,7,8}]};


TaylorForm[DALIlist_List] := Module[
	{
		LIC, rules, DALITensor, dim  = Sqrt[Length[DALIlist[[1,1]]]], order= Length@DALIlist,
		list, dummy, multiplicities
	},
	(*LI components of the DALI tensors according to their position in the DALI list*)
	LIC[i_,j_] := SymmetrizedIndependentComponents[ConstantArray[dim, i+j], DALISymmetry[i,j]];
	
	Do[
		(*Make the rules: LI component -> value*)
		rules[i,j] = MapThread[Rule, {LIC[i,j], Flatten[DALIlist[[i,j]]]}],
		{i, order},
		{j, 1, i}
	];
	
	Do[
		(*Create the Symmetrized Array Object to represent the Tensors:*)
		DALITensor[i,j] = SymmetrizedArray[rules[i,j], ConstantArray[dim,i+j], DALISymmetry[i,j]];
		(*
			Introduce c[i,j] to account for the coefficients of the expansion in the tensors
			and Make these tensors totally symmetric:
		*)
		DALITensor[i,j] = Symmetrize[
			c[i,j]*DALITensor[i,j],
			Symmetric[All]
		],
		{i, order},
		{j, 1, i}
	];
	
	
	(*Organize Tensors by their rank and  summ different tensors with the same rank*)
	list = Table[
		dummy[i,j],
		{i,order},
		{j, 1, i}
	]//Flatten;
	
	list = SortBy[list, (List@@#//Total)&]; (*This will organize ranks into increasing order*)
	
	(*This puts different ranks into different sublists and sums dummys in the same sublist with each other*)
	list = Plus@@@GatherBy[list, (List@@#//Total)&];
	
	(*Use the actual tensors:*)
	list = list//.dummy->DALITensor;
	
	(*Extract the LI components of every rank:*)
	Clear@LIC;
	(LIC[#] =  SymmetrizedIndependentComponents[ConstantArray[dim, #], Symmetric[All]])&/@Range[2, 2*order];
	list = Table[
		Extract[list[[i]], LIC[i+1]],
		{i, 1, 2*order-1}
	];
	
	(*now just introduce the multiplicities for every rank:*)
	(multiplicities[#] = PermutationsNumber/@LIC[#])&/@Range[2, 2*order];
	Table[
		list[[i]]*multiplicities[i+1],
		{i, 1, 2*order-1}
	]
]


Protect@TaylorForm;


(* ::Subsection:: *)
(*Legacy Stan and Compiled Polynomial*)


(*StanParser[a_] := Module[
	{res},
	res = HornerForm[a];
	res = Block[{Power=pow}, res];
	CForm[res]//ToString
]

StanParser[x___] := Throw[$Failed, failTag[StanParser]]*)


(*MakeStanCode[expression_, vars_List] := Module[
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

MakeStanCode[x___] := Throw[$Failed, failTag[MakeStanCode]]*)


(*StanPolynomial[DALIOutput_, fiducialPoint_] := Module[
	{TaylorPolynomial, vars},
	
	Catch[
        TaylorPolynomial = TaylorForm[DALIOutput, fiducialPoint];
        vars = fiducialPoint[[All,1]];
        MakeStanCode[TaylorPolynomial, vars]
      ]
]

StanPolynomial[x___] := Throw[$Failed, failTag[StanPolynomial]] *)


(*CompiledPolynomial::fail = "The function failed. The failure occured in function `1`"*)


(*CompiledPolynomial[GetDALIOutput_, fiducialPoint_?MatrixQ] := Module[
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

CompiledPolynomial[x___] := Throw[$Failed, failTag[CompiledPolynomial]]*)


(* ::Section::Closed:: *)
(*Package Footer*)


End[];


EndPackage[];
