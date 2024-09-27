(* ::Package:: *)

(* ::Section:: *)
(*Dark Mode and some other settings*)


(*
SetOptions[EvaluationNotebook[], DefaultNewCellStyle->"Code"]
Needs["maTHEMEatica`"]
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


(* ::Section::Closed:: *)
(*BeginPackage and Summary*)


BeginPackage["FelipeBarbosa`SymDALI`DALICoefficients`"];


DALICoefficients::usage="
Assuming
              ln(\[ScriptCapitalL]) = - \!\(\*FractionBox[\(1\), \(2\)]\) (\!\(\*SuperscriptBox[\(m\), \(i\)]\) - \!\(\*SuperscriptBox[\(\[Mu]\), \(i\)]\)) (\!\(\*SuperscriptBox[\(C\), \(\:207b1\)]\)\!\(\*SubscriptBox[\()\), \(ij\)]\) (\!\(\*SuperscriptBox[\(m\), \(j\)]\) - \!\(\*SuperscriptBox[\(\[Mu]\), \(j\)]\))    and    \!\(\*SuperscriptBox[\(\[Mu]\), \(i\)]\) \[Congruent] \[Mu](\!\(\*SubscriptBox[\(p\), \(\[Alpha]\)]\), \!\(\*SubscriptBox[OverscriptBox[\(x\), \(\[RightVector]\)], \(i\)]\)),
the call                                                                  
            DALICoefficients[\[Mu], {{\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, {\!\(\*SubscriptBox[\(p0\), \(1\)]\), \!\(\*SubscriptBox[\(p0\), \(2\)]\),...}, n}, {\!\(\*SuperscriptBox[\(C\), \(-1\)]\), False}, ObservationPoints]
generates the list of coefficients for the DALI expansion with respect to {\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, around {\!\(\*SubscriptBox[\(p0\), \(1\)]\),\!\(\*SubscriptBox[\(p0\), \(2\)]\),..}
to order n in derivatives, where \!\(\*SuperscriptBox[\(C\), \(-1\)]\) is a non-diagonal matrix and ObservationPoints = {\!\(\*SubscriptBox[OverscriptBox[\(x\), \(\[RightVector]\)], \(1\)]\), \!\(\*SubscriptBox[OverscriptBox[\(x\), \(\[RightVector]\)], \(2\)]\), ...}.

For a diagonal matrix, such that:
                            ln(\[ScriptCapitalL]) = - \!\(\*FractionBox[\(1\), \(2\)]\) \!\(\*FractionBox[SuperscriptBox[\((\*SuperscriptBox[\(m\), \(i\)] - \*SuperscriptBox[\(\[Mu]\), \(i\)])\), \(2\)], SuperscriptBox[\((\*SubscriptBox[\(\[Sigma]\), \(i\)])\), \(2\)]]\)    and     \!\(\*SuperscriptBox[\(\[Mu]\), \(i\)]\) \[Congruent] \[Mu](\!\(\*SubscriptBox[\(p\), \(\[Alpha]\)]\), \!\(\*SubscriptBox[OverscriptBox[\(x\), \(\[RightVector]\)], \(i\)]\)),
the call
             DALICoefficients[\[Mu], {{\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, {\!\(\*SubscriptBox[\(p0\), \(1\)]\), \!\(\*SubscriptBox[\(p0\), \(2\)]\),...}, n}, {\[Sigma], True}, ObservationPoints],
generates the list of coefficients for the DALI expansion with respect to {\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, around {\!\(\*SubscriptBox[\(p0\), \(1\)]\),\!\(\*SubscriptBox[\(p0\), \(2\)]\),..}
to order n in derivatives, where \[Sigma] = {\!\(\*SubscriptBox[\(\[Sigma]\), \(1\)]\), \!\(\*SubscriptBox[\(\[Sigma]\), \(2\)]\), ...} and ObservationPoints = {\!\(\*SubscriptBox[OverscriptBox[\(x\), \(\[RightVector]\)], \(1\)]\), \!\(\*SubscriptBox[OverscriptBox[\(x\), \(\[RightVector]\)], \(2\)]\), ...}.
";


GWDALICoefficients::usage="
Assuming
                      ln(\[ScriptCapitalL]) = - \!\(\*FractionBox[\(1\), \(2\)]\)(g-h, g-h)   and  (u,\[Nu]) \[Congruent] 4 \[ScriptCapitalR] \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)  \!\(\*FractionBox[\(\*OverscriptBox[\(u\), \(~\)]\[Conjugate]  \((\[ScriptF])\)\\\ \*OverscriptBox[\(\[Nu]\), \(~\)] \((\[ScriptF])\)\), \(\*SubscriptBox[\(S\), \(n\)] \((\[ScriptF])\)\)]\) d\[ScriptF],
the \!\(\*
StyleBox[\"call\", \"Code\"]\)
             GWDALICoefficients[\!\(\*OverscriptBox[\(h\), \(~\)]\), {{\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, {\!\(\*SubscriptBox[\(p0\), \(1\)]\), \!\(\*SubscriptBox[\(p0\), \(2\)]\),...}, n}, Sn, {fmin, fmax, \[CapitalDelta]f}],
generates the list of coefficients for the DALI expansion with respect to {\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, around {\!\(\*SubscriptBox[\(p0\), \(1\)]\),\!\(\*SubscriptBox[\(p0\), \(2\)]\),..}
to order n in derivatives, where Sn = {\!\(\*SubscriptBox[\(S\), \(n\)]\)(\!\(\*SubscriptBox[\(f\), \(1\)]\)), \!\(\*SubscriptBox[\(S\), \(n\)]\)(\!\(\*SubscriptBox[\(f\), \(2\)]\)), ...} and \!\(\*OverscriptBox[\(h\), \(~\)]\) = \!\(\*OverscriptBox[\(h\), \(~\)]\)(\!\(\*SubscriptBox[\(p\), \(\[Alpha]\)]\),f).
";


NDALICoefficients::usage = "Generate the list of DALI tensors to a given derivative order, for likelihoods defined numerically"


Begin["`Private`"];


(* ::Section:: *)
(*Shared Functions*)


Needs["Developer`"]
ToSymbol[a_String,i_Integer] := ToExpression[ a <> ToString[i]]
ToSymbol[a_String] := ToExpression[a]
SymbolConcatenate[sym_Symbol, i_Integer] := ToExpression[ToString[sym] <>  ToString[i]]


(* ::Subsection::Closed:: *)
(*Utilities to implement SymRules and NRules*)


(*
	SymRulesValues: list of values from SymRules, without the lfhs of the rules, i.e. Flatten[Values@SymRules][[All,2]]
	Output: list with distinct numbers and a head for each number, i.e. {{0,1.2,3., ...}, {h1,h2, h3,...}}
*)

CreateHeads[SymRulesValues_]/; (VectorQ[SymRulesValues, NumberQ]) := Module[
	{list, heads},
	list = DeleteDuplicates[SymRulesValues];
	
	heads = Unique@(("h" <> ToString[#])&/@Range[Length@list]);
	
	(*Return heads and their corresponding values*)
	{list, heads}
]

CreateHeads[x___] := Throw[$Failed, failTag[CreateHeads]]


(*
	heads: flat list of all heads in SymRules, i.e. Flatten[Values@SymRules][[All, 1, 0, -1]]
	lefthandsides: all lefthandsides in SymRules without pattern vars, i.e. Flatten[Values@SymRules][[All, 1, 0]]
	righthandsides: list of heads, replacing corresponding numbers, i.e AuxFunc/@(Flatten[Values@SymRules][[All, 2]]),
	where AuxFunc[xi] := headi, where headi is defined in "CreateHeads" and headi[x__] := xi.
	
	Output: None. It defines the following kind of UpValues  f/: $D[{1,0,0}, f] = hi, for all expressions in SymRules.
*)

SymUpValues[heads_List, lefthandsides_List, righthandsides_List]/;(
	DeleteDuplicates[heads[[All,0]]] === {Symbol} && DeleteDuplicates[righthandsides[[All,0]]] === {Symbol} &&
	DeleteDuplicates[lefthandsides[[All, 0]]] === {$D}
) := MapThread[
	TagSet,
	{heads, lefthandsides, righthandsides}
]

SymUpValues[x___] := Throw[$Failed, failTag[SymUpValues]]


(*
	lhside: normal expression of the form head[x1_, x2_,...]
	rhside: CompiledFunction object CompiledFunction[...]
	Output: None. It defines the assignement head[x1_,x2_,...] := CompiledFunction[...][x1,x2,...]
	* The CompiledFunction must have the RuntimeOption "EvaluateSymbolically" ->False
*)

NDownValue[lhside_, rhside_CompiledFunction] := With[
	{lhs = lhside, rhs = rhside},
	Quiet[Inactive[SetDelayed][lhs, rhs@@(lhs[[All,1]])]//Activate]
]

NDownValue[x___] := Throw[$Failed, failTag[NDownValue]]


(*
	lhsides: all left hand sides of NRules with unique heads, i.e. Flatten[Values@NRules][[All, 1]]/.heads -> uniqueheads
	rhsides: all rhsides of NRules, i.e. Flatten[Values@NRules][[All,2]] 
	tags: all unique heads in the lhside of NRules, i. e., Last[#,#]&/@(Flatten[Values@NRules][[All, 1, 0]])/. heads -> uniqueheads
	* Last[#,#]& is to return tag in "$D[{1,1,1}, tag]" and tag also in  "tag" 
	
	Output: None. It performs the following kind of assignements 
	uniquehead[x__] := CompiledFunction[...][x] 
	uniquehead/: $D[{1,0,0}, uniquehead] = CompiledFunction[...]
*)

NUpValues[tags_List, lhsides_List, rhsides_List]/;(
	DeleteDuplicates[tags[[All,0]]] === {Symbol} && DeleteDuplicates[rhsides[[All,0]]] === {CompiledFunction}
) :=Module[
	{atomsPos, itags, ilhsides, irhsides},
	
	atomsPos = Position[AtomQ/@(lhsides[[All,0]]), True];
	
	MapThread[
		NDownValue,
		{Extract[lhsides, atomsPos], Extract[rhsides, atomsPos]}
	];
	
	itags = Delete[tags, atomsPos];
	ilhsides = Delete[lhsides, atomsPos][[All,0]];
	irhsides = Delete[rhsides, atomsPos];
	
	MapThread[
		TagSet,
		{itags, ilhsides, irhsides}
	];
]

NUpValues[x___] := Throw[$Failed, failTag[NUpValues]]


(* ::Subsection:: *)
(*Step2*)


(* ::Subsubsection::Closed:: *)
(*GenDaliTerm*)


(*These are now lists of vectors carrying LI components only*)

GenDaliTerm[grads1_List, grads2_List,matrix_List, True]/;(
	MatrixQ[grads1, NumberQ] &&MatrixQ[grads2, NumberQ]
) := Flatten[Transpose[grads1] . grads2]


GenDaliTerm[grads1_, grads2_, matrix_, False]/;(
	MatrixQ[grads1, NumberQ] &&MatrixQ[grads2, NumberQ]
) := Flatten[Transpose[grads1] . matrix . grads2]

GenDaliTerm[x___] := Throw[$Failed, failTag[GenDaliTerm]];


(*TEST:*)

(*Create lists of symmetrized arrays*)


(*list1 = Table[
SymmetrizedArray[
    MapThread[Rule,{SymmetrizedIndependentComponents[{4,4,4}, Symmetric[All]], RandomReal[{10,100}, {20}]}],
    {4,4,4}, Symmetric[All]
],
{10}
]; 

list2 = Table[
SymmetrizedArray[
    MapThread[Rule,{SymmetrizedIndependentComponents[{4,4,4,4}, Symmetric[All]], RandomReal[{10,100}, {35}] }],
    {4,4,4,4}, Symmetric[All]
],
{10}
];*)

(*

(*Take diagonal tensor products*)
a = Total@MapThread[TensorProduct, {list1, list2} ];

(*and LI components of the result*)

aValues = Values[a["ArrayRules"]];
(*Now take the tensor product among LI components only*)
values1 = (Values[#["ArrayRules"]]&)/@list1;
values2 = (Values[#["ArrayRules"]]&)/@list2;
values3 = GenDaliTerm[values1, values2,  True, {1,2}];
(*Comparison:*)
values3 == aValues
*)

(*Now for a non diagonal case:*)
(*matrix = RandomReal[{10,100}, {10,10}];

TensorProduct[list1[[1]], list2[[1]]]*2.

term  = Sum[matrix[[i,j]]*TensorProduct[list1[[i]], list2[[j]]], {i, 10}, {j, 10}];

(*termValues = Values[term["ArrayRules"]];*)

(*termValues == GenDaliTerm[values1, values2, False, matrix]*)

Clear[list1, list2, a, aValues, values1, values2, values3, matrix,  termValues]*)


(* ::Subsubsection::Closed:: *)
(*RecoverSymmetric*)


(*Note that for the Fisher Matrix and any product of the same tensors you include the components for perms between
pairs of indices, so it is best to leave it that way:*)

IndComponents[dim_Integer, i_,j_] := Module[
    {list1, list2}, 
    
    list1 = SymmetrizedIndependentComponents[ConstantArray[dim, i], Symmetric[All]];
    list2 = SymmetrizedIndependentComponents[ConstantArray[dim, j],Symmetric[All]];
    
    Outer[{#1,#2}&, list1, list2, 1]//Flatten[#,1]&
]


(*TEST: JUST format*)

(*IndComponents[3, 2, 1]
IndComponents[3,1,1]*)


RecoverSymmetric[vectorOfLIcomponents_?VectorQ, dim_?IntegerQ,  {i_,  j_} ] := Module[
    {listOfRules,  listIndcomp, \[Alpha], \[Beta], \[Alpha]list, \[Beta]list, joinedList},
    listIndcomp = IndComponents[dim,i,j];
    
    listOfRules = Association[
        Rule@@@(Riffle[listIndcomp, vectorOfLIcomponents]//Partition[#,2]&)
    ];
    
    \[Alpha]list = SymbolConcatenate[\[Alpha],#]&/@Range[i]; 
    \[Beta]list = SymbolConcatenate[\[Beta], #]&/@Range[j];
    joinedList = Join[\[Alpha]list, \[Beta]list];
    
    
    Table[
        listOfRules[{Sort[\[Alpha]list], Sort[\[Beta]list]}], 
        
        Evaluate[
            Sequence@@(Riffle[joinedList, ConstantArray[dim, i+j]]//Partition[#,2]&)
        ]
    ]
]


(*TEST: Just the Format for now*)

(*DEFINE YOUR TENSORS FOR THE DaliList*)


(*fm = SymmetrizedArray[
    Rule@@@(Riffle[
         SymmetrizedIndependentComponents[{4,4}, Symmetric[All]],
         RandomReal[{0,10}, {10}]
    ]//Partition[#,2]&), 
    {4,4},
    Symmetric[All]
];

tensor1 = SymmetrizedArray[
    Rule@@@(Riffle[
         SymmetrizedIndependentComponents[{4,4,4}, {{Symmetric[{1,2}]},  {Symmetric[{3}] } } ],
         RandomReal[{0,10}, {40}]
    ]//Partition[#,2]&), 
    {4,4, 4},
    {{Symmetric[{1,2}]},  {Symmetric[{3}]}}
];

tensor2 = SymmetrizedArray[
    Rule@@@(Riffle[
         SymmetrizedIndependentComponents[{4,4,4,4}, {{Symmetric[{1,2}]},  {Symmetric[{3,4}] } } ],
         RandomReal[{0,10}, {100}]
    ]//Partition[#,2]&), 
    {4,4, 4,4},
    {{Symmetric[{1,2}]},  {Symmetric[{3,4}]}}
];

(*Create the DALI list:*)
daliList = Map[Normal, {{fm}, {tensor1, tensor2}}, {2}];

(*Now create the values and the corresponding DALIlist of LI components*)
fmValues = Values[fm["ArrayRules"]];
tensor1Values = Values[tensor1["ArrayRules"]];
tensor2Values = Values[tensor2["ArrayRules"]];

daliListValues = {{(Normal@fm)//Flatten},{tensor1Values, tensor2Values}};

(*Compare RecoverSymmetric against Normal:*)
MapIndexed[RecoverSymmetric[#1, 4, #2]&, daliListValues, {2}] == daliList

Clear[fm, tensor1, tensor2, daliList, fmValues, tensor1Values, tensor2Values, daliListValues]*)



(* ::Subsubsection:: *)
(*GenDaliList*)


GenDaliList[list_Symbol, order_Integer, diag_?TrueQ, matrix_List, dim_] := Module[
    {dummy, clearList, DaliList},
    
    clearList[i_, order] := ReleaseHold@ToExpression[ToString[list] <> ToString[i] <> "=.", InputForm, Hold];
    
    DaliList = Do[
        dummy = Divide[SymbolConcatenate[list,i], matrix^2]; 
        
        Do[
            Sow[#,k]&@GenDaliTerm[SymbolConcatenate[list,k], dummy, matrix, diag]; 
            clearList[i, k],
            {k, i, order}  (*Not a fan of nested Do, but I want to generate dummy just once for each i *)
         ],
       
       {i, 1, order}
    ]//Reap//Last;
    
    DaliList
]


GenDaliList[list_Symbol, order_Integer, diag_/;(diag==False), matrix_List, dim_] := Module[
    {clearList, DaliList},
    
    clearList[i_, order] := ReleaseHold@ToExpression[ToString[list] <> ToString[i] <> "=.", InputForm, Hold];
    
    DaliList = Do[
        Sow[#,k]&@GenDaliTerm[SymbolConcatenate[list,k], SymbolConcatenate[list,i], matrix, diag]; 
        clearList[i,k],
        {i, 1, order},{k, i, order}
    ]//Reap//Last;
    
    DaliList
    
]


(*TEST:*)

(*Let's make fake lists of gradients and symmetric tensors.*)
(*list1 = Table[RandomReal[{1,10}, {4}], {15}];

list2 = Table[
SymmetrizedArray[
    MapThread[Rule,{SymmetrizedIndependentComponents[{4,4}, Symmetric[All]], RandomReal[{10,100}, {10}]}],
    {4,4}, Symmetric[All]
],
{15}
];

list3 = Table[
SymmetrizedArray[
    MapThread[Rule,{ SymmetrizedIndependentComponents[{4,4,4}, Symmetric[All]], RandomReal[{10,100}, {20}] }],
    {4,4,4}, Symmetric[All]
],
{15}
];

(*Take the lists of LI components*)
Values1 = list1;
Values2 = (Values[#["ArrayRules"]]&)/@list2;
Values3 = (Values[#["ArrayRules"]]&)/@list3;*)

(*From here on comment either the Digonal or nonDiagonal case according to what you want to test:*)


(*DiagonalCASE*)

(*Create a sigma vector*)
(*testsigma = RandomReal[{1.,2.}, {15}];

(*Make the DALI list manually from the symmetric tensors*)
manualDALIlist = Table[
	Total@(Normal/@MapThread[TensorProduct, {ToSymbol["list", k]/(testsigma^2), ToSymbol["list", j]}]),
     {k,1,3},
     {j,1,k}
];

(*Feed the values to the automatic generator*)
automaticDAliLIst = GenDaliList["Values", 3, True, testsigma, 4];

automaticDAliLIst == manualDALIlist
*)


(*NON DIAGONALCASE*)

(*testmatrix = RandomReal[{1,2}, {15,15}];

(*Make the dali list manually:*)
manualDALIlist = Table[
	Normal@Sum[testmatrix[[p,m]]*TensorProduct[
	    (ToSymbol["list", k])[[p]], 
	    (ToSymbol["list", j])[[m]]
        ], 
      {p,15}, {m,15}],
     {k,1,3},
     {j,1,k}
];

automaticDAliLIst = GenDaliList["Values", 3, False, testmatrix, 4];

manualDALIlist == automaticDAliLIst*)


(* ::Section::Closed:: *)
(*DALICoefficients*)


(* ::Subsection::Closed:: *)
(*step1*)


(*
I could also replace the full symbolic array by an association, with this I would probably not need an full array
	functionhead: head of function that you want to take the gradient
	dummyvariables: Matrix of variables of the form
	 {
		{p1, ..., pn}, {Subscript[s, 1], p2, ..., pn}, {Subscript[s, 1],Subscript[s, 2], p3, ...,pn}, ...,{Subscript[s, 1],...,Subscript[s, n-1],pn}, {p1,...,pn, X1, ..., XM}
	},
	* Assuming your function depends on (p1,...,pn,X1,...,XM) and your gradients are with respect to pi. 
	(Subscript[s, i] just ensure 0's in LD components of the arrays)
	tensorRank: tensor rank of the function ```functionhead```.
	Output: the gradient of functionhead[p1,...,pn, X1,...,XM], i.e. an array of rank (tensorRank+1), 
	with zeros in the LD components of the array. 
*)

TakeGrad[functionhead_Symbol,  dummyvariables_List, tensorRank_Integer]/;(
	functionhead[[0]] === Symbol && Positive[tensorRank] && DeleteDuplicates[Flatten[dummyvariables][[All,0]]] === {Symbol}
) := MapIndexed[
		D[#1, {dummyvariables[[#2//Last]]}]&,
		functionhead@@(dummyvariables//Last),
		{tensorRank}
]


(*Same thing as before, for an scalar*)
TakeGrad[functionhead_Symbol, dummyvariables_List, tensorRank_Integer]/;(
	 tensorRank==0 && DeleteDuplicates[Flatten[dummyvariables][[All, 0]]] === {Symbol}
) := D[
	functionhead@@(dummyvariables//Last), 
	{dummyvariables[[1]]}
]

TakeGrad[x___] := Throw[$Failed, failTag[TakeGrad]]


(*
	head: head such that head[p1,...,pn, X1,...,Xm] gives your array.
	Orighs: heads that appear in the explicit expression head[p1,...,pn, X1,...,Xm].
	uniquehs: uniqueHeads for which UpValues were defined from NRules. 
	OpsPoints: set of values to pass to the function in the form of {pf1,pf2,..., Overscript[X, ->]1, Overscript[X, ->]2, ...},
	pfi are numbers and Overscript[X, ->]i are lists of numbers.
	
	Output: List of gradients (only the LI components of head[p1,...,pn, X1,...,Xm]) evaluated at 
	{
		{pf1,pf2,...,Subscript[( Overscript[X, ->]1), 1], Subscript[(Overscript[X, ->]2), 1], ...},
		{pf1,pf2,...,Subscript[( Overscript[X, ->]1), 2], Subscript[(Overscript[X, ->]2), 2], ...},
		{pf1,pf2,...,Subscript[( Overscript[X, ->]1), 3], Subscript[(Overscript[X, ->]2), 3], ...},
		...
	}
*)

CalculateGrads[head_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List,tensorRank_Integer]/;(
	Length[{Orighs}] === Length[{uniquehs}] && Positive@tensorRank && 
	Length[ObsPoints] === Length@FunctionVariables[head]
) := Module[
    {dummyhead, LIComponents, x, dummyvariables, result, dim = Count[ObsPoints, _?NumberQ]},
    
    (*Construct pattern variables*)
    dummyvariables = Unique@ConstantArray[x, Length[ObsPoints]];
    
    LIComponents = SymmetrizedIndependentComponents[ConstantArray[dim, tensorRank], Symmetric[All]];
    
    Inactive[Set][(*make dummyhead[x1_, x2_,..., head1_, head2_,...]*)
		dummyhead@@(Pattern[#,_]&/@Join[dummyvariables, {Orighs}]),
		Extract[head@@dummyvariables, LIComponents];
	]//Activate;
    
    result = dummyhead@@Join[ObsPoints, {uniquehs}];
    
    (*
    assuming the components of the function are listable, this is a list of lists, where each list is a 
    LI component evaluated at several points.
    Transpose to return a list of lists where each sublist is all LI components evaluated at a single point
    *)
    
    Transpose@result
]


CalculateGrads[x___] := Throw[$Failed, failTag[CalculateGrads]]


(*TEST:*)

(*SetAttributes[sym2, Listable];
sym2[x_, y_, z_, w_] =  { {2 Sin[x y z] \[ExponentialE]^(w+y), Cos[x+y] }, {0, Sin[x z] } };
testlist = Transpose[RandomReal[{0.,1.},{10,4}]];
CalculateGrads[sym2, testlist, test, $MachinePrecision, 2, 2]

test == (   Delete[#, 3]&/@(Flatten/@(sym2@@testlist) ))

ClearAll[sym2, testlist, test]*)


FunctionVariables[head_] := Module[
	{s},
	s = (DownValues[head][[1,1]]);
	(Hold@@s)/. Hold[h_[x__]] :> {x}[[All,1]]
]


(*
	head: head of the function to take gradients
	vars: {p1, ..., pn, X1, ..., Xm}; 
	obspoints: set of values to pass to the function in the form of {pf1,pf2,..., Overscript[X, ->]1, Overscript[X, ->]2, ...}
	order: desired derivative order
	Output: {list1, list2,...}, where 
	Subscript[list, i] = {\[Del]^ihead[pf1, pf2...,( Overscript[X, ->]1) ] }
*)
GenGrads[head_Symbol, vars_List, obspoints_List, {Orighs___Symbol}, {uniquehs___Symbol}, order_Integer]/;(
	Positive[order] && Length[{Orighs}] === Length[{uniquehs}] && 
	Length[obspoints] === Length@FunctionVariables[head] === Length[vars]
 ) := Module[
	{dummyvariables, dim = Count[obspoints, ?NumberQ], h1, h2, Auxh},
    
    (*Create dummyvariables for TakeGrad:*)
     dummyvariables = Module[
        {result}, (*Maybe there is a functional way for this*)
        Reap[
            result = Sow@vars[[1;;dim]];
            Do[
                result = Sow@ReplacePart[result, i -> Unique["s"]], 
                {i, 1, dim -1}
            ];
        Sow[vars];
        
        ]//Last//Last
    ];
        
    (*Create auxiliar heads for TakeGrad*)
    Auxh[x_?OddQ] := h1;
    Auxh[x_?EvenQ] := h2;
    
    (*Define a function of the form f[p1,...,pn,X1, ...,Xm] to act on:*)
    Inactive[Set][
        h1@@(Pattern[#, _]&/@vars),
        head@@FunctionVariables[head]
    ]//Activate;
    
    (*Implementing the loop:*)
    Table[
        Inactive[Set][
			Auxh[i+1]@@(Pattern[#,_]&/@vars),
			TakeGrad[Auxh[i], dummyvariables, i-1]
		]//Activate;
       
       Clear[Evaluate[Auxh[i]]];
    
       CalculateGrads[Auxh[i+1], {Orighs}, {uniquehs}, obspoints, i],
       {i, 1, order} 
    ]
]


GenGrads[x___] := Throw[$Failed, failTag[GenGrads]]


(*TEST, Just output format, not its correctness, for now: *)

(*sym2[x_, y_, z_, w_] =  Sin[x y z] \[ExponentialE]^(w+y);
points = Transpose@RandomReal[{0.,1.}, {10, 4}];

GenGrads[sym2[x,y,z,w], {x,y,z}, 3, points, 3, list, $MachinePrecision, {}, {}]

Table[TensorRank/@ToSymbol["list",i], {i,1,3}]
Table[Length/@ToSymbol["list",i], {i,1,3}]
ClearAll[list1,list2,list3,points, sym2 ]*)


(* ::Subsection::Closed:: *)
(*GetDALITensors*)


SetAttributes[DALICoefficients, HoldFirst];

DALICoefficients[FTheory_, parameters_List, {\[Sigma]_List, diagonal_?BooleanQ}, ObservationPoints_List,wp_:$MachinePrecision, symbolicRules_:{}, numericalRules_:{}]/;And[(*Conditions on the Input:*)
        Head/@(List@@Unevaluated[FTheory]) == ConstantArray[Symbol,  Length@(List@@Unevaluated[FTheory]) ],
        Length@parameters==3,
        Head/@parameters == {List, List, Integer},
        parameters[[-1]] > 0, 
        Length[ parameters[[1]] ] ==  Length[ parameters[[2]] ],
        Length[\[Sigma]] == Length[ObservationPoints]
    ]:= Module[

{\[Sigma]list, obspoints, DaliList, list, dim, DOrder},

(*##############################################  DECLARING BASIC VARIABLES ################################################################*)

dim = parameters[[1]]//Length;
DOrder = parameters//Last;
\[Sigma]list =\[Sigma];

obspoints = Transpose[ (*REDUNDANT BECAUSE LISTABILITY ALLOWS FOR f[x1,x2, list1, list2]*)
         ({parameters[[2]], #}//Flatten)&/@ObservationPoints
     ];

(*############################################EXECUTING FUNCTIONS ##########################################################################*)

GenGrads[FTheory, parameters[[1]], dim, obspoints, DOrder, list, wp, symbolicRules,  numericalRules];

GenDaliList[list, Max@DOrder, diagonal, \[Sigma]list,  dim]
]


(*HEADTEST[x_, y_, z_, w_] =  Sin[x y z]*E^(w+y);

points = RandomReal[{0.,1.}, 10];
\[Sigma] = RandomReal[{0.,1.}, 10];

test = DALICoefficients[HEADTEST[x,y,z,w], {{x,y,z}, {1,2,5}, 4}, {\[Sigma], True}, points];
Map[TensorRank, test, {2}]
Clear[HEADTEST, points, \[Sigma], test]*)


(* ::Section::Closed:: *)
(*GWDALICoefficients*)


(*Likelihood def. eq. 42 of https://arxiv.org/pdf/1809.02293*)
GWGenDaliTerm[gradlist1_, gradlist2_, SensitivityVector_, \[CapitalDelta]f_] := Module[
    {complexSum },
    
    complexSum = Flatten@Total@MapThread[ 
         KroneckerProduct, 
         {Conjugate[gradlist1],  Divide[gradlist2, SensitivityVector]}
     ];
    4 \[CapitalDelta]f Re[complexSum]  (*I think this 4 \[CapitalDelta]f can be just absorbed in the normalization...*)
]


(*list1 =  Flatten/@RandomComplex[{0, 1+I}, {10,3,3}];
list2 = Flatten/@RandomComplex[{0,1+I}, {10, 3,3,3}];
sensitivity = RandomReal[{0.1 ,2.}, 10];

test = GWGenDaliTerm[list1, list2, sensitivity, 0.2];
(*TensorRank@test
ArrayDepth@test*)
test

Clear[list1,list2,sensitivity, frequency, test]*)


GWGenDaliList[varsPrefix_Symbol, DOrder_Integer, Sensitivity_List, dim_, \[CapitalDelta]f_ ] := Module[
    {dummy, clearList},
    
    clearList[i_, DOrder] := ReleaseHold@ToExpression[ToString[varsPrefix] <> ToString[i] <> "=.", InputForm, Hold];
    
    dummy = Do[
    
        Sow[#,k]&@GWGenDaliTerm[SymbolConcatenate[varsPrefix,k], SymbolConcatenate[varsPrefix, i], Sensitivity, \[CapitalDelta]f]; 
        clearList[i, k],
        {i, 1, DOrder}, {k, i, DOrder}
    
    ]//Reap//Last;
    
 dummy
]


SetAttributes[GWDALICoefficients, HoldFirst]

GWDALICoefficients[
    h_, 
    parameters_List,
    SensitivityVector_, 
    ObsPoints_, 
    wp_,
    symbolicRules_, 
    numericalRules_
    
    ]/;(
        Head/@(List@@Unevaluated[h]) == ConstantArray[Symbol, Length[parameters[[1]]] + 1 ]&&
        Length[parameters]==3 &&
        Head/@parameters== {List, List, Integer} &&
        parameters[[-1]] > 0 && 
        Length[parameters[[1]]] == Length[ parameters[[2]] ] &&
        Length[ObsPoints] == 3 &&
        Length[SensitivityVector] == Round[(ObsPoints[[2]] - ObsPoints[[1]])/ObsPoints[[3]] + 1]
    ):= Module[
    {\[Sigma]list, obspoints, DaliList, dim,  list, Dorder},
    
    dim = parameters[[1]]//Length;
    Dorder = parameters//Last;
    
    \[Sigma]list = SensitivityVector;
    
    obspoints = Transpose[
        ({parameters[[2]], #}//Flatten)&/@(Range@@ObsPoints)
    ];
    
    
    EchoTiming[GenGrads[h, parameters[[1]], dim, obspoints, Max@Dorder, list, wp, symbolicRules, numericalRules], "GenGrads"];
    
    
       
    DaliList = EchoTiming[GWGenDaliList[list, Dorder, SensitivityVector, dim, ObsPoints[[3]] ], "GWGenDaliList"]; 
    
    DaliList
    
]


(* ::Section:: *)
(*NDALICoefficients*)


(* ::Subsection::Closed:: *)
(*step1 part a)*)


(*
	 h: function head
	 var: symbol to serve as dependent variable
*)

NTakeDerivative[h_Symbol, var_Symbol] := D[
	h[var], var
]//Simplify[#, ComplexityFunction->ByteCount, TimeConstraint->0.1]&

NTakeDerivative[x___] := Throw[$Failed, failTag[NTakeDerivative]]


(*
	h: function head
	Points: Vector of values
	
*)

NCalculateDerivatives[h_Symbol, Points_List] := Module[
    {dummyHead,x},
    
    SetAttributes[dummyHead, Listable];
    
    dummyHead[x_] = h[x];
    
    dummyHead[Points]
]

NCalculateDerivatives[x___] := Throw[$Failed, failTag[NCalculateDerivatives]]


(*TEST*)

(*f[x_] = x^2;
Map[f,Range[6]]
NCalculateDerivatives[f,Range[6]]
Clear[f]*)


(*
	h: head of the function
	Points: List of points at wich you want to calculate the function and its derivatives
	DOrder: maximum derivative order
	
	Output: Matrix where each row is one derivative of the function calculated at all required points, starts at 1 and 
	ends at d^n/d x^n
*)

NGenDerivatives[h_Symbol, Points_, n_Integer] := Module[
    {Auxh, x, h1, h2},
    
    Auxh[x_?OddQ] := h1;
    Auxh[x_?EvenQ] := h2;
    
    h1[x_]  = h[x];
    
    Table[
       Evaluate[Auxh[i+1][x_]] =  NTakeDerivative[Auxh[i], x];
       
       Clear[Evaluate@Auxh[i]];
       
       NCalculateDerivatives[Auxh[i+1], Points],
        {i, 1, n}
    ]
]


NGenDerivatives[x___] := Throw[$Failed, failTag[NGenDerivatives]]


(*TEST:*)

(*testlist = RandomReal[{0.,1.}, 10];
f[x_] = 3 x^2;
NGenDerivatives[f, testlist, 3]//MatrixForm

Clear[testlist, f]*)


(* ::Subsection::Closed:: *)
(*step1 part b)*)


(*memoization because the same value will be called several times*)
NIndComponents[dim_, tensorRank_] :=  NIndComponents[dim,tensorRank] = SymmetrizedIndependentComponents[ConstantArray[dim,tensorRank], Symmetric[All] ];			


(*Not sure if it is optimized for vectors*)

GenerateSymmetric[LIComponents_, dim_, tensorRank_, LILabels_Symbol] := Module[
	{rules},
	
	rules = MapThread[
		Rule,
		{LILabels[tensorRank], LIComponents}
	];

	Normal@SymmetrizedArray[rules, ConstantArray[dim, tensorRank], Symmetric[All]]

]


(*TEST:*)

(*test = { {1,2,3}, {0,5,6},{0,0,1}};
test2 = GenerateSymmetric[test];
test//MatrixForm
test2//MatrixForm
test3 = {1,2,3,4};
GenerateSymmetric[test3]
Clear[test,test2,test3 ]*)


(*
	h: function head, such that h[p1,..., pn, X] gives the full array (pi are the parameter variables)
	obspoints: {pf1,pf2, ..., Overscript[X, ->]}
	tensorRank: TensorRank[h[p1,..., pn, X]]
	
	Output:  = ( Overscript[Subscript[\[Integral], 0], #] dz \[Del]^tensorRank h[pf1,pf2,..., z])&/@Overscript[X, ->]
*)

NCalculateGrads[h_Symbol, obspoints_List, tensorRank_Integer]/;(
	Positive[tensorRank] && Length[obspoints] === Length@FunctionVariables[h]
) := Module[
    {y, int,LIComponents, vars, dim = Count[obspoints, _?NumberQ], obsVarsnumber = Count[obspoints, _List],
    numberLI
    },
    
    
    (*Make variables for X1, ..., Xm*)
    vars = Unique[
        ("x" <> ToString[#])&/@Range[obsVarsnumber]
    ];
    
    LIComponents = SymmetrizedIndependentComponents[ConstantArray[dim, tensorRank], Symmetric[All]];
    numberLI = Length[LIComponents];
    
    (*Extract the LI components with pi -> pfi \[And] Xi -> dummyVariablei*)
    LIComponents = Extract[
        h@@Join[obspoints[[1;;dim]], vars], 
        LIComponents
    ];
    
    
    int = NDSolve[
        {   (*This only makes sense if you have just one dummy variable*)
            Derivative[1][y]@@vars == LIComponents, 
            y[0] == ConstantArray[0, numberLI] 
        }, 
        y,
        
        {vars[[1]], 0, obspoints[[-1,-1]]},
        PrecisionGoal->8
   ]//Flatten//Last//Last;
   
   int[obspoints[[-1]]]
]


(*
Overload for tensorRank===0
p\[Alpha] stands for grad variables
*)
NCalculateGrads[h_Symbol, obspoints_List, 0, p\[Alpha]_List]/;(
	 Length[obspoints] === Length@FunctionVariables[h]
) := Module[
    {y, int, vars,  obsVarsnumber = Count[obspoints, _List], dummy},
    
    
    (*Make variables for X1, ..., Xm*)
    vars = Unique[
        ("x" <> ToString[#])&/@Range[obsVarsnumber]
    ];
    
    (*make a head with the variables properly ordered*)
    With[
		{ivars = DeleteDuplicates[Join[p\[Alpha], FunctionVariables[h]]]},
		
		Evaluate[dummy@@(Pattern[#,_]&/@ivars)] = h@@FunctionVariables[h]
    ];
   
    int = NDSolve[
        {   (*This only makes sense if you have just one dummy variable*)
            Derivative[1][y]@@vars == dummy@@Join[obspoints[[1;;-2]], vars], 
            y[0] == 0
        }, 
        y,
        
        {vars[[1]], 0, obspoints[[-1,-1]]},
        PrecisionGoal->8
   ]//Flatten//Last//Last;
   
   int[obspoints[[-1]]]
]

NCalculateGrads[x___] := Throw[$Failed, failTag[NCalculateGrads]]


(*
	h: head of the function to take gradients
	vars: {p1, ..., pn, X}; 
	obspoints: set of values to pass to the function in the form of {pf1,pf2,..., Overscript[X, ->]}
	n: desired derivative order
	Output: {list1, list2,...}, where 
	
	Subscript[list, i] =( Overscript[Subscript[\[Integral], 0], #] dz \[Del]^i h[pf1,pf2,..., z])&/@Overscript[X, ->] 
*)

NGenGrads[h_Symbol, {vars__Symbol},  obspoints_, n_Integer]/;(
	Length[obspoints] === Length@FunctionVariables[h] === Length[{vars}] && Positive[n] 
) := Module[
    {dummyvariables, dim = Length[{vars}]-1,Auxh, h1,h2  },
    
    (*Creating dummyvariables for TakeGrad:*)
	dummyvariables = Module[ 
        {result}, 
        Reap[
            result = Sow@({vars}[[1;;dim]]);
            Do[
                result = Sow@ReplacePart[result, i -> Unique["s"] ], 
                {i, 1, dim -1}
            ];
        Sow[{vars}];
        
        ]//Last//Last
    ];
    
        
    (*Create auxiliar heads for TakeGrad*)
    Auxh[x_?OddQ] = h1;
    Auxh[x_?EvenQ] = h2;
    
    
    Inactive[Set][
		h1@@(Pattern[#,_]&/@{vars}),
		h[vars]
	]//Activate;
    
   
    Table[
        Inactive[Set][
            Auxh[i+1]@@(Pattern[#,_]&/@{vars}),
            TakeGrad[Auxh[i], dummyvariables, i-1]
        ]//Activate;
        
        Clear[Evaluate[Auxh[i]]];
        
        
        NCalculateGrads[Auxh[i+1], obspoints, i],
        {i, 1, n} 
    ]
];


NGenGrads[x___] := Throw[$Failed, failTag[NGenGrads]]


(*TEST, Just the form right now.*)

(*sym2T[x_, y_, z_, w_] =  Sin[x y z] \[ExponentialE]^(w+y);
obspoints = RandomReal[{0,1}, {10,4}, WorkingPrecision->30];
test = NGenGrads[sym2T, 3, obspoints,  4, 30];

Map[TensorRank,test,{2}]

ClearAll[sym2T, obspoints, test]*)


(* ::Subsection::Closed:: *)
(*step1 part c)*)


(*
	x: {n1,n2,n3, n1, n4,......}
	Output: (n1+n2+...)!/(n1! n2! ... 2! r2! ...),
	
	* ri \[Congruent] how many times the number ni appears in "x", for instance:
	SymmetryCoeff[{4,4,9,9,9,1}] == (Plus@@{2,2,9,9,9,1})!/((4!)^2 (9!)^3 1! 2! 3!) , where 2! and 3! appear because 4 and 9 appear 2 and 3 times.
	
	*The function calculates how many distinct permutations of the object:  "Subscript[x, \[Alpha]1...\[Alpha]n1] Subscript[x, \[Beta]1...\[Beta]n2] ... "
	you have to sum so in order to obtain a completely symmetric output. Assuming that the labels of each "x" are 
	totally symmetric and exchange \[Alpha]i by \[Beta]i in Subscript[x, \[Alpha]1...\[Alpha]n1] Subscript[x, \[Beta]1...\[Beta]n2] does not change anything if n1==n2
*)

SymmetryCoeff[x_List] := With[
    {multiplicities = Times@@(#!&/@Extract[Tally@x, {All,2}]) },
    
    Divide[Multinomial@@x, multiplicities]
]


(*
	derivatives: {f'[x1], f''[x1], f'''[x1], ..., (f^(n))[x1] }
	gradients: list of all gradients {\[Del] \[Mu], \[Del]^2 \[Mu], \[Del]^3\[Mu], ..., \[Del]^n \[Mu]};      \[Del]^2 \[Mu] \[Congruent] \[Del] \[TensorProduct] \[Del] \[Mu],   \[Del]^3 \[Mu] \[Congruent] \[Del] \[TensorProduct] \[Del] \[TensorProduct] \[Del] \[Mu], and so on
	*gradients here are full arrays, not a vector of LIComponents.
	partitions: IntegerPartitions[n]
	permNumber: number of distinct perms for each combination: SymmetryCoeff/@partitionlist;
	
	Output: \[Del]^n f[\[Mu]]
*)

ConvertToGrad[derivatives_List, gradients_List, partitions_List, permNumber_List]/;(

	VectorQ[derivatives, NumberQ] && Length[gradients] === Length[derivatives]  && VectorQ[permNumber, NumberQ]
) := Module[
	{numberOfderivatives, tensorList, iderivatives, result},
	
	numberOfderivatives = Length/@partitions;
	
	(*make the list {Subscript[x, \[Alpha]\[Beta]\[Gamma]],  Subscript[x, \[Alpha]\[Beta]] Subscript[x, \[Gamma]] , Subscript[x, \[Alpha]],Subscript[x, \[Beta]],Subscript[x, \[Gamma]]}*)
	tensorList = Outer[Times, Sequence@@Part[gradients, #]]&/@partitions;
	
	(*make the list: {f1'[x1] Subscript[x, \[Alpha]\[Beta]\[Gamma]],f''[x1] Subscript[x, \[Alpha]\[Beta]] Subscript[x, \[Gamma]] , f'''[x1] Subscript[x, \[Alpha]],Subscript[x, \[Beta]],Subscript[x, \[Gamma]]} *)
	iderivatives = derivatives[[#]]&/@numberOfderivatives;
	tensorList = MapThread[
        Times[#1, #2]&, 
        {iderivatives, tensorList}
    ];
    
	Symmetrize[
        Plus@@(tensorList*permNumber), 
        Symmetric[All]
    ]//Normal
]

ConvertToGrad[x___] := Throw[$Failed, failTag[ConvertToGrad]]


(*TEST:*)

(*f[x_] = Log[x];
g[x_,y_,z_,w_] = Sin[x y z w];

test = D[f[g[x,y,z,w]], {{x,y,z}, 6}]/.{x->1., y->3.,z->4.,w->5.};
list1 = (D[f[x],{x, #}]/.x->g[1.,3.,4.,5.])&/@Range[6];
list2 = (D[ g[x,y,z,w], {{x,y,z}, #} ]/.{x->1., y->3.,z->4.,w->5.})&/@Range[6];

ConvertToGrad[list1, list2] == test

Clear[g,f,test,list1,list2]*)


(* ::Text:: *)
(*ConvertToGradIterate:*)
(* *)
(*transforms the lists of primary ingredients for grads into, that lists of derivatives of "f" and gradients of the several variables gradients.*)
(**)
(*ConvertToGradIterate starts getting the highest order grads and then goes to the lowest order, cleaning the highest order grads once it doesn't need them anymore. The output is a list of gradients with the prefix "OutVarsPrefix".*)


(*
	derivativeMatrix: Matrix of derivatives:
	{
		{f'[x1], f''[x1], ...},
		{f'[x2], f''[x2], ...},
		...	
	}
	
	gradMatrix: Matrix of gradients(full arrays not LI components):
	{
		{\[Del]\[Mu][x1], \[Del]^2\[Mu][x1], ... },
		{\[Del]\[Mu][x2], \[Del]^2\[Mu][x2], ...},
		
		...
	}
	
	Output: Matrix of gradients:
	{
		{\[Del]f[\[Mu][x1]], \[Del]f[\[Mu][x2]], \[Del]f[\[Mu][x3]], ...},
		{\[Del]^2f[\[Mu][x1]], \[Del]^2f[\[Mu][x2]], \[Del]^2f[\[Mu][x3]], ...},
		...
	}
*)

ConvertToGradIterate[derivativeMatrix_,  gradMatrix_] := Module[
    {partitions, n = Length[derivativeMatrix[[1]] ], permNumbers, result},
    
     (partitions[#] = IntegerPartitions[#])&/@Range[n];
     (permNumbers[#] = SymmetryCoeff/@partitions[#])&/@Range[n];
     
	Table[
		
		MapThread[
			ConvertToGrad[#1, #2, partitions[i], permNumbers[i]]&,
			{derivativeMatrix[[All,1;;i]], gradMatrix[[All, 1;;i]] }
		],
		
		{i,1,n}
	]
]



(*TEST:*)
(*matr1 = RandomReal[{0,1}, {2,2}];
matr2 = {{ {1,2}, RandomReal[{0,1},{2,2}]},  { {1,2}, RandomReal[{0,1},{2,2}]}};
ConvertToGradIterate[matr1, matr2, "test", 2]*)


(* ::Subsection:: *)
(*NGetDALITensors*)


(* ::Subsubsection:: *)
(*NDALICoefficients*)


(*Module[
	{l1,l2,l3,l4,M},
	l1 = RandomReal[{0,1}, {10,3}];
	l2 = RandomReal[{0,1}, {10,3,3}];
	l3 = RandomReal[{0,1}, {10, 3,3,3}];
	l4 = RandomReal[{0,1}, {10,3,3,3,3}];
	M = {l1,l2,l3,l4};
	
	Map[
		Dimensions,
		Transpose[M],
		{2}
	]
]*)


NDALICoefficients::fail = "The function failed. The failure occured in function `1`";


NDALICoefficients[{f_Symbol, \[Mu]_Symbol}, {{vars__Symbol}, {fp__?NumberQ}, n_Integer}, {\[Sigma]_, diag_?BooleanQ}, obsPoints_]/;(
	Length@FunctionVariables[f] === 1 && Length@FunctionVariables[\[Mu]] == Length[{vars}] + 1 && Positive[n] &&
	Length[{vars}] === Length[{fp}] && VectorQ[obsPoints,NumberQ] && Length[obsPoints] === Length[\[Sigma]] &&
	(
		(VectorQ[\[Sigma], NumberQ] && diag===True ) || (MatrixQ[\[Sigma], NumberQ] && diag===False)
	)

) := Module[
	{i\[Sigma], iobsPoints, dpoints, derivativeMatrix, gradientMatrix, ivars, LIComponents, list},	   
	
	i\[Sigma] = ToPackedArray@N@\[Sigma];
	iobsPoints = {fp,  ToPackedArray@N@obsPoints};
	ivars = Union[{vars}, FunctionVariables[\[Mu]]];
	(LIComponents[#] = SymmetrizedIndependentComponents[
		ConstantArray[Length[{vars}], #], Symmetric[All]
	])&/@Range[n];
	
	
	Catch[
		dpoints = NCalculateGrads[\[Mu], iobsPoints, 0, {vars}];
		derivativeMatrix = NGenDerivatives[f, dpoints, n];
		
		EchoTiming[gradientMatrix = NGenGrads[\[Mu], ivars,  iobsPoints, n], "\[CapitalDelta]t for all the gradients:"];
		
		(*We need to recover full Arrays and transpose the matrix in order to pass to ConvertToGradIterate:*)
		EchoTiming[gradientMatrix = MapIndexed[
			GenerateSymmetric[#1, Length[{vars}], First[#2], LIComponents]&,
			gradientMatrix,
			{2}
		];, "Total \[CapitalDelta]t to recover full array structure:"];
		
		gradientMatrix = Transpose[gradientMatrix];
		
		gradientMatrix = Map[
			ToPackedArray,
			gradientMatrix,
			{2}
		];
		
		(*Also this one*)
		derivativeMatrix = ToPackedArray@Transpose[derivativeMatrix];
		
		EchoTiming[
			gradientMatrix = ConvertToGradIterate[derivativeMatrix,  gradientMatrix];,
			"\[CapitalDelta]t to combine gradients and derivatives:"];
			
		(*Extract the LI Components and define the "list" symbol:*)
		 Inactive[Set][
			SymbolConcatenate[list,#]&/@Range[n],
			ToPackedArray/@MapIndexed[
				Extract[#1, LIComponents[First[#2]]]&,
				gradientMatrix,
				{2}
			]
		]//Activate;
		
		
		
		
		Clear[gradientMatrix];
		
		EchoTiming[GenDaliList[list, n, diag, i\[Sigma], Length[{vars}]], "GenDALIList:"],
		_failTag, 
		(Message[NDALICoefficients::fail, Style[First@#2, Red]];#1) &
	] 
]


(*Module[
	{Data,m, \[Sigma] },
	Import["/Users/felipe/Documents/GitHub/DALI/Data/DALI-1000SNe.csv", "Data"][[1]];
	Data = Drop[Import["/Users/felipe/Documents/GitHub/DALI/Data/DALI-1000SNe.csv", "Data"], 1];
	m = Data[[All,2]];
	z = Data[[All,1]];
	\[Sigma] = Data[[All,3]]; (*REALLY REDUNDANT BCS ALL UNCERTAINTIES ARE THE SAME NUMBER*)

	M = (IdentityMatrix[1000] -ConstantArray[1/1000, {1000, 1000}])/0.15^2;
]
\[Mu][\[CapitalOmega]m_, w0_, wa_, zp_] = (\[CapitalOmega]m*(1 +zp )^3 + Exp[-((3*wa*zp)/(1+zp))]*(1-\[CapitalOmega]m)*(1+zp)^(3*(1+w0+wa)))^(-1/2);

g[x_] = 5*Log[10, x];
obspoints = {0.285, -1.0, 0., z};*)


(*NDALICoefficients[
	{g, \[Mu]},
	{ {\[CapitalOmega]m,w0, wa}, {0.285, -1., 0}, 3},
	{M, False},
	z
];//AbsoluteTiming*)


(* ::Subsubsection::Closed:: *)
(*Checks of NGetDALITensors and GetDALITensors*)


(*DIAGONAL CROSS CHECKING*)

(*MACHINE PRECISION AND NDSOLVE AUTOMATIC IS REALLY BAD FOR HIGH ORDERS*)
(*
analytic[x_,y_,z_,w_] = Log[Integrate[Sin[x y z] E^(k+y), {k,0,w}]];
head2[x_,y_,z_,k_] = Sin[x y z] E^(k+y);

obspoints = Table[{1,2,3,i},{i,15}];
\[Sigma] = RandomReal[{1.,2.}, 15];

dummyLog[x_] = Log[x];


ANalytic = EchoTiming[DALICoefficients[analytic[a1,a2,a3,a4], {{a1,a2,a3},{1,2,3}, 4 },  {\[Sigma], True},  Range[15],  4]];


(*numerical = EchoTiming[NGetDALITensors[dummyLog, head2, 3, {\[Sigma], True}, obspoints, 4 ]];*)
Map[TensorRank, ANalytic, {2}]
(*Map[TensorRank, numerical, {2}]*)
(*ANalytic==numerical*) (*NDSolve Introduces Considerable Numerical Errors:*)

(*Table[((ANalytic[[i,j]]- numerical[[i,j]])/(ANalytic[[i,j]]+numerical[[i,j]]))//Max, {i, 4 }, {j,i}]*)

Do[
    gradlist[i] = Table[
        D[analytic[x,y,z,w], {{x,y,z}, i}]//.MapThread[Rule, {{x,y,z,w}, obspoints[[j]]}],
        {j,1,15}
    ]/\[Sigma],
    {i,1,4}
]

resultest = Table[
    MapThread[TensorProduct, {gradlist[i], gradlist[j]}]//Total,
    {i,1,4},
    {j,1,i}
];

Map[TensorRank, resultest, {2}]
resultest == ANalytic
Clear[analytic, head2, obspoints, \[Sigma], dummyLog, ANalytic, numerical, gradlist, resultest]*)


(*NON-DIAGONAL CROSS CHECKING*)

(*analytic[x_,y_,z_,w_] = Log[Integrate[Sin[x y z] \[ExponentialE]^(k+y), {k,0,w}]];
head2[x_,y_,z_,k_] = Sin[x y z] \[ExponentialE]^(k+y);

obspoints = Table[{1,2,3,i},{i,15}];
\[Sigma]list = RandomReal[{1.,2.}, {15,15}];

dummyLog[x_] = Log[x];

ANalytic = GetDALITensors[analytic, 3,  {\[Sigma]list, False},  obspoints,  4];
numerical = NGetDALITensors[dummyLog, head2, 3, {\[Sigma]list, False}, obspoints, 4 ];
Map[TensorRank, ANalytic, {2}]
Map[TensorRank, numerical, {2}]

Table[((ANalytic[[i,j]]- numerical[[i,j]])/(ANalytic[[i,j]]+numerical[[i,j]]))//Max, {i, 4 }, {j,i}]

Do[
    gradlist[i] = Table[
        D[analytic[x,y,z,w], {{x,y,z}, i}]//.MapThread[Rule, {{x,y,z,w}, obspoints[[j]]}],
        {j,1,15}
    ],
    {i,1,4}
]

resultest = Table[
    Sum[\[Sigma]list[[k,l]]*TensorProduct[(gradlist[i])[[k]], (gradlist[j])[[l]]], {k,1,15}, {l,1,15}],
    {i,1,4},
    {j,1,i}
];

Map[TensorRank, resultest, {2}]

resultest==ANalytic

Clear[analytic, head2, obspoints, \[Sigma]list, dummyLog, ANalytic, numerical, gradlist, resultest]*)


(* ::Section::Closed:: *)
(*EndPackage*)


End[];
EndPackage[];
