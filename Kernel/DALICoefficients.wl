(* ::Package:: *)

(* ::Section:: *)
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`DALICoefficients`", {"FelipeBarbosa`SymDALI`DerivativeTools`"}];


Unprotect@DALITensors;
ClearAll@DALITensors


Begin["`Private`"];


(*Get["/Users/felipe/Documents/GitHub/SymDALI/Kernel/DerivativeTools.wl"]
PacletDirectoryLoad["/Users/felipe/Documents/GitHub"];*)


(* ::Section:: *)
(*Definitions*)


Unprotect[
	SymRules, NRules,iFpFc,ihphcIMRPhenomD, ihphcIMRPhenomPv2, \[ScriptA]IMR, \[CapitalPhi]IMR, FpFc
];


ClearAll[
	SymRules, NRules,iFpFc, ihphcIMRPhenomD, ihphcIMRPhenomPv2, \[ScriptA]IMR, \[CapitalPhi]IMR, FpFc
]


<<Developer`
ToSymbol[a_String,i_Integer] := ToExpression[a <> ToString[i]]
ToSymbol[a_String] := ToExpression[a]
SymbolConcatenate[sym_Symbol, i_Integer] := ToExpression[ToString[sym] <>  ToString[i]]
FunctionVariables[head_] := Module[
	{s},
	s = (DownValues[head][[1,1]]);
	(Hold@@s)/. Hold[h_[x__]] :> {x}[[All,1]]
]
Attributes[GenMessage] = {HoldRest};
GenMessage[True, mess_] := True;
GenMessage[False, mess_] := With[{}, Message[mess]; False]


(* ::Subsection::Closed:: *)
(*Make Gradients*)


(* ::Subsubsection::Closed:: *)
(*Function Set-Up*)


CreateHeads::usage="CreateHeads[{SymRules___Rule}] 

{SymRules}: Flat SymRules List, i.e. Flatten[Values@SymRules]
Output: list with distinct numbers and a unique head for each number
Example: 
>>> With[{d = <|
	\"f1 \" -> {$D[{2,0}, f1]-> 2., $D[{2,1}, f1] -> 0.}, 
	\"f2 \" -> {$D[{3}, f2] -> 3., $D[{4}, f2] -> 0.}
|>},
	CreateHeads[Flatten[Values@d]]
]
>>>{{h11,h12,h13},{2.`,0.`,3.`}}";


CreateHeads[{SymRules___Rule}] := Module[
	{list, heads, SymRulesValues},
	
	SymRulesValues = {SymRules}[[All,2]];
	
	list = DeleteDuplicates[SymRulesValues];
	
	heads = Unique@ConstantArray["h", Length[list]];
	
	(*Return heads and their corresponding values*)
	{heads, list}
]

CreateHeads[x___] := Throw[$Failed, failTag[CreateHeads]]


(*
Assuming the lhs is either Condition[$D[{n__}, h], cond] -> something or $D[{n__}, h]-> something
it returns "h"
*)
GetSymRuleHead//ClearAll
GetSymRuleHead[def_Rule]/;def[[1,0]] === Condition := def[[1, 1,-1]]
GetSymRuleHead[def_Rule]/; def[[1,0]] === $D := def[[1,-1]]
GetSymRuleHead[x___] :=Throw[$Failed, failTag[GetSymRuleHead]]


ClearAll@SymUpValues
SymUpValues::usage="SymUpValues[{SymRules___Rule}, AuxHead_]
{SymRules}: Flat SymRules List, where all the rules belong to the same head.
AuxHead: head such that AuxHead[xi] = headi, where xi and headi are corresponding head-number pairs from \"CreateHeads\".
Output: None. It defines the following kind of UpValues  f/: $D[n__, f] = hi, for all expressions in SymRules.

Example:
>>>Block[{d = {$D[{2,0}, f1] -> 2., $D[{n1_,n2_}, f1]/;n2>n1 -> 0.}, heads, aux},
	
	heads = Echo[CreateHeads[d]];
	
	MapThread[
		(aux[#2] = #1)&,
		heads
	];
	Print[{aux[2.], aux[0.]}];
	
	SymUpValues[d, aux];
	Print[UpValues/@{f1,f2}];
	Clear[f1]
]

>>>{{h72,h73},{2.`,0.`}}
{
	{HoldPattern[$D[{2,0},f1]]\[RuleDelayed]h72,HoldPattern[$D[{n1_,n2_},f1]/;n2>n1]\[RuleDelayed]h73},
	{}
}";


SymUpValues[{SymRules__Rule}, AuxHead_Symbol] := Module[
	{iList, dummyFunction,  head = GetSymRuleHead[{SymRules}[[1]]]},

	(*Make a list where the rhs is the head corresponding to the number:*)
	iList = MapAt[AuxHead, {SymRules}, {All, 2}];

	(*place HoldPattern on the lhs and Rule-> RuleDelayed*)
	iList = RuleDelayed@@@(MapAt[HoldPattern, iList, {All, 1}]);
	
	dummyFunction[x_] := (UpValues[x] = iList);
	(*define UpValues:*)
	dummyFunction[head];
]

SymUpValues[x___] := Throw[$Failed, failTag[SymUpValues]]


ProcessSymRules//ClearAll

ProcessSymRules::usage="ProcessSymRules[SymRules_Association]
SymRules: Association with SymRules for all heads
Output: list of all heads such that head[x___] = number_i, for some number_i in the SymRules";


ProcessSymRules[<||>] := 1


ProcessSymRules[SymRules_Association]/;Length[SymRules]>0 := Module[
	{headNumberPair, aux},
	
	(*Make all the head-value pairs*)
	headNumberPair = CreateHeads[Flatten[Values@SymRules]];
	
	(*Make the heads evaluate to corresponding numbers*)
	MapThread[(#1[x___] := #2)&, headNumberPair];
	
	(*Just so you can replace numbers by heads in order to define UpValues*)
	MapThread[
		(aux[#2] = #1)&,
		headNumberPair
	];
	
	(*Make the UpValues for each key in the association*)	
	SymUpValues[#, aux]&/@Values[SymRules];
	
	ClearAll@aux; Remove[aux];
	
	headNumberPair[[1]]
]

ProcessSymRules[x___] := Throw[$Failed, faiTag[ProcessSymRules]]


GCSymRules[heads_List] := Module[
	
	{uniquehs}, 
	
	uniquehs = (UpValues/@heads)//Flatten;
	uniquehs = uniquehs[[All,2]]//DeleteDuplicates;
	
	(UpValues[#] = {})&/@heads; Remove[Evaluate[uniquehs]];
]


newHead[Rule] := SetDelayed
newHead[RuleDelayed] := SetDelayed
newHead[TagRule] := TagSetDelayed
newHead[x___] := Throw[$Failed, failTag[newHead]]


ChangeNRulesHead::usage="ChangingNRulesHead[{rules__}]
{rules}: list of rules of a particular head. The elements are either Rule[Condition[...], rhs] or 
TagRule[tag, Condition[...], rhs]
Output: none, it changes Rule -> SetDelayed and TagRule -> TagSetDelayed";

ChangeNRulesHead[{rules__}] := Module[
	{}, 
	MapAt[
		newHead,
		{rules},
		{All, 0}
	]
]

ChangeNRulesHead[x___] := Throw[$Failed, failTag[ChangeNRulesHead]]


ProcessNRules::usage="Given NRules as an association, 
it will make the relevant definitions in for the lists";

ProcessNRules[NRules_Association] := Module[
	{},
	ChangeNRulesHead/@NRules;
]

ProcessNRules[x___] := Throw[$Failed, failTag[ProcessNRules]]


Parser[expr_, originalHeads_List, newHeads_List] := Module[
	{vars, defs, FullExpr, x},
	vars =  HoldComplete[x]/.x -> originalHeads;
	
	defs = MapThread[set, {originalHeads, newHeads}];
	defs = HoldComplete[x]/. x-> Join[defs, {expr}];
	defs = CompoundExpression@@@defs;
	
	FullExpr = Join[vars, defs];
	
	Block[
		{set=Set}, 
		Block@@FullExpr
	]
]

Parser[x___] := Throw[$Failed, failTag[Parser]]


MakeDefs::usage="MakeDefs[SymRules_Association, NRules_Association]
SymRules: SymRules association
NRules: NRules association
Output: None. It defines DownValues and UpValues for the heads in SymRules and NRules.

Example:
>>>Block[
	{SymRules, NRules, f1, f2},
	
	SymRules = <|
		\"f1\" -> {$D[{1,1}, f1][x_, y_] -> 2., $D[{2, 1}, f1][x_, y_] -> 0.}, 
		\"f2\" -> {$D[{3}, f2][x_] -> 3., $D[{4}, f2][x_] -> 0.}
	|>;
	
	NRules = <|
		\"f1\" -> {
			f1[x_,y_] -> CompiledFunction[\"1\"], 
			$D[{1,0}, f1][x_,y_] -> CompiledFunction[\"2\"], 
			$D[{0,1}, f1][x_,y_] -> CompiledFunction[\"3\"]
		},
		\"f2\" -> {
			f2[x_] ->CompiledFunction[\"4\"], 
			$D[{1}, f2][x_] -> CompiledFunction[\"5\"],
			$D[{2}, f2][x_]-> CompiledFunction[\"6\"] 
		}
	|>;
	
	MakeDefs[SymRules, NRules];
	{
		Comap[{UpValues, DownValues},f1],
		Comap[{UpValues, DownValues},f2]
	}
]//Quiet

{
	{
		{HoldPattern[$D[{1,0},f1]]\[RuleDelayed]CompiledFunction[\"2\"],HoldPattern[$D[{1,1},f1]]\[RuleDelayed]h56,HoldPattern[$D[{2,1},f1]]\[RuleDelayed]h57},
		{HoldPattern[f1[x_,y_]]\[RuleDelayed]CompiledFunction[\"1\"][x,y]}
	},
	{
		{HoldPattern[$D[{1},f2]]\[RuleDelayed]CompiledFunction[\"5\"],HoldPattern[$D[{3},f2]]\[RuleDelayed]h58,HoldPattern[$D[{4},f2]]\[RuleDelayed]h57},
		{HoldPattern[f2[x_]]\[RuleDelayed]CompiledFunction[\"4\"][x]}
	}
}";

MakeDefs[SymRules_Association, NRules_Association] := Module[
	{},
	
	ProcessSymRules[SymRules];
	ProcessNRules[NRules];

]

MakeDefs[x___] := Throw[$Failed, failTag[MakeDefs]]


(* ::Subsubsection::Closed:: *)
(*Calculate symbolic gradient*)


TakeGrad::usage="TakeGrad[functionhead_Symbol,  dummyvariables_List, tensorRank_Integer]
functionhead: Head. Assuming \"functionhead[p1,...,pn,X1,...,XM]\" \[And] gradients with respect to pi
dummyvariables: Matrix of variables of the form
{
	{p1, ..., pn}, 
	{s1, p2, ..., pn}, 
	{s1,s2, p3, ...,pn}, 
	...,
	{s1,...,sn-1,pn}, 
	{p1,...,pn, X1, ..., XM}
},
tensorRank: tensor rank of \"functionhead[p1,...,pn,X1,...,XM]\".
Output: gradient of \"functionhead[p1,...,pn, X1,...,XM]\"  (array of rank \"1 + tensorRank\").
* \"si\" ensures 0's in LD components of arrays
* \"tensorRank===-1\" returns \"functionhead[p1,...,pn, X1,...,XM]\"

Example:
>>>Block[
	{h, dummy = {{x,y}, {s1,x}, {x,y,z}}}, 
	h[x_,y_, z_] = {Exp[x+y] z, Sin[z +x/y]};
	{
		TakeGrad[h, dummy, -1],
		TakeGrad[h, dummy, 1]
	}
]

{
	{\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(x + y\)]\) z,Sin[\!\(\*FractionBox[\(x\), \(y\)]\)+z]},
	{{\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(x + y\)]\) z,\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(x + y\)]\) z},{0,\!\(\*FractionBox[\(Cos[\*FractionBox[\(x\), \(y\)] + z]\), \(y\)]\)}}
}";

TakeGrad//Clear
TakeGrad[functionhead_Symbol,  dummyvariables_List, tensorRank_Integer]/; Positive[tensorRank] := MapIndexed[
		
		D[#1, {dummyvariables[[#2//Last]]}]&,
		
		functionhead@@(dummyvariables//Last),
		
		{tensorRank}
]

(*Same thing as before, for an scalar*)
TakeGrad[functionhead_Symbol, dummyvariables_List, 0]/;(
	DeleteDuplicates[Flatten[dummyvariables][[All, 0]]] === {Symbol}
) := D[
	functionhead@@(dummyvariables//Last), 
	{dummyvariables[[1]]}
]

(*Overload -1 here should take no derivative, necessary for gradients in GW's. *)
TakeGrad[functionhead_Symbol, dummyvariables_List, -1]/;(
	DeleteDuplicates[Flatten[dummyvariables][[All, 0]]] === {Symbol}
) := functionhead@@(dummyvariables//Last)


TakeGrad[x___] := Throw[$Failed, failTag[TakeGrad]]


(* ::Subsubsection::Closed:: *)
(*iCalculateGrads*)


iCalculateGrads//ClearAll


iCalculateGrads::usage="iCalculateGrads[head_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List,tensorRank_Integer, n]

head: function head, such that head[p1,...,pn, X1,...,Xm] gives the array.
Orighs: heads that appear in the explicit expression head[p1,...,pn, X1,...,Xm].
uniquehs: uniqueHeads for which UpValues were defined from NRules. 
OpsPoints: set of values to pass to the function in the form of {pf1,pf2,..., XV1, XV2, ...}
tensorRank: tensorRank of head[p1,...,pn, X1,...,Xm].
n: number of pi's with repect to which you take the gradients

Output: LI components of the gradient. 

* \"tensorRank===0\" evaluates the scalar at ObsPoints {...}.
*\"iCalculateGrads\" assumes that Attributes[h] = {Listable}, i.e., if any of the XVi happen to be a vector the \"head[..., XVi,...]\"automatically
distributes over the values of XVi, implying that the result is a matrix where each row represents 1 LI component
evaluated at several different points.


Example1:
>>>Block[
	{h, obs = {2,3,4,Range[1,5]}},
	h[x_,y_,z_,w_] = {{(x+y+z) w, (x+y+z) Power[w,2]}, {0, (x+y+z) Sqrt[w]}};
	Print[2+3+4];
	iCalculateGrads[h, {},{}, obs, 2, 2]
]
>>>9
>>>{{9,18,27,36,45},{9,36,81,144,225},{9,9 \!\(\*SqrtBox[\(2\)]\),9 \!\(\*SqrtBox[\(3\)]\),18,9 \!\(\*SqrtBox[\(5\)]\)}}

Example2:
>>>Block[
	{h, obs = {2,3,4, Range[0,5]}},
	h[x_,y_,z_,w_] = (x+y+z) w;
	Print[2+3+4];
	iCalculateGrads[h, {},{}, obs,0, 3]
]
>>>9
>>>{0,9,18,27,36,45}";


iCalculateGrads[head_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List, tensorRank_Integer, dim_, rule_, string_]/;(
	Length[{Orighs}] === Length[{uniquehs}] && Positive@tensorRank && 
	Length[ObsPoints] === Length@FunctionVariables[head]
) := Module[
    {dummyhead, LIComponents, x, dummyvariables, result},
    
    (*Construct pattern variables*)
    dummyvariables = Unique@ConstantArray[x, Length[ObsPoints]];
    
    LIComponents = SymmetrizedIndependentComponents[ConstantArray[dim, tensorRank], Symmetric[All]];
    
    
    (*
        make dummyhead[x1_, x2_,..., head1_, head2_,...] = head@@[x1, x2, ...,xf],
        (for higher dimension gradients you need only LI components of "head@@[x1, x2, ...,xf]" )
    *)
    Set@@(set[
		dummyhead@@(Pattern[#,_]&/@Join[dummyvariables, {Orighs}]),
		Block[
			{FelipeBarbosa`GWFORECAST`Private`\[CapitalPhi]IMR},
			FelipeBarbosa`GWFORECAST`Private`\[CapitalPhi]IMR[x__] = 0; (*Optimization that works for PhenomD only at first*)
			Extract[head@@dummyvariables, LIComponents]
		]
	]);
    
    
    result = dummyhead@@Join[ObsPoints, {uniquehs}];
    
   (*
   If head===hphc then each component of the gradient is a matrix with dimension {Length[fvec], 2} (bcs of the auxiliar rule for \[ScriptA]imr)
   then Dimension[result]==={LI, Length[fvec], 2}.
   For the other functions it is convenient that 
   Dimension[result] === {LI,  2, Length[fvec]}
   *)
   If[string==="hphc", result = Transpose/@result];
   
   
    (*Clean the dummy variables and return result:*)
    ClearAll[Evaluate[dummyvariables]]; Remove[Evaluate[dummyvariables]];
    
    Association@MapThread[
		Rule,
		{LIComponents, result}
    ]
    
]//Check[#, Throw["TakeGradFail at order" <> ToString[tensorRank]]]&


iCalculateGrads[h_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List, 0, dim_, rule_, string_]/;(
	Length[{Orighs}] === Length[{uniquehs}]
) := Module[
	{dummyhead, dummyvariables, result,x},
    
    (*Construct pattern variables including {Orighs}*)
    dummyvariables = Unique@ConstantArray[x, Length[ObsPoints]];
    
    (*
        make dummyhead[x1_, x2_,..., head1_, head2_,...] = head@@[x1, x2, ...,xf],
        note the absence of Extract now, bcs the function is not a gradient
    *)
    Set@@(set[
		dummyhead@@(Pattern[#,_]&/@Join[dummyvariables, {Orighs}]),
		Block[
			{FelipeBarbosa`GWFORECAST`Private`\[CapitalPhi]IMR},
			FelipeBarbosa`GWFORECAST`Private`\[CapitalPhi]IMR[x__] = 0; (*Optimization that works for PhenomD only at first*)
			h@@dummyvariables
		]
	]);
    
   result = dummyhead@@Join[ObsPoints, {uniquehs}];
   (*again we want the dimensions {2, Length[fvec]} in the output*)
   If[string==="hphc", result = Transpose[result]];
   
   (*Clean the dummy variables and return result:*)
    ClearAll[Evaluate[dummyvariables]]; Remove[Evaluate[dummyvariables]];
   
   <|{0} -> result|>
   
]//Check[#, Throw["TakeGradFail at order" <> ToString[tensorRank]]]&

iCalculateGrads[x___] := Throw[$Failed, failTag[iCalculateGrads]]


(* ::Subsubsection::Closed:: *)
(*iGenGrads*)


iGenGrads//ClearAll


iGenGrads::usage="iGenGrads[head, dim, obspoints, {Orighs___}, {uniquehs___}, n_, ni_]
head: head of the function to take gradients
** Assuming that the function is head[p1,...pdim, X1, ...,Xm]
dim: gradient dimension
obspoints: values to pass to the function in the form of {pf1,pf2,..., XV1,XV2, ...}
n: DALI order
ni: \"ni===0\" -> gradients from order 0 to n;   \"ni===1\" -> gradients from order 1 to n
Output: List of gradients, with gradients beeing represented by matrices where each row is  LI component calculated at several points.

Example:
>>>Block[
	{h,obs = {2,3,4, SeedRandom[1234];RandomReal[{0,1}, {2,2}], Range[5]}},
	
	h[x_, y_, z_, m_, f_] := Exp[x-y+z] Total[m, 2] f;
	
	iGenGrads[h, 2, obs, {}, {}, 2, 1]
]

>>>{
	{
		{37.413512887965325`,74.82702577593065`,112.24053866389598`},
		{-37.413512887965325`,-74.82702577593065`,-112.24053866389598`}
	},
	{
		{37.413512887965325`,74.82702577593065`,112.24053866389598`},
		{-37.413512887965325`,-74.82702577593065`,-112.24053866389598`},
		{37.413512887965325`,74.82702577593065`,112.24053866389598`}}
}";


iGenGrads::negativeOrder="Requested derivative order is not positive";
iGenGrads::UniqueHeads="Length[{Orighs}] does not match Length[{Uniquehs}]";
iGenGrads::obsPoints = "Length of obsPoints is not bigger than parameter space dimension";


iGenGrads[head_Symbol, dim_Integer, obspoints_List, {Orighs___Symbol}, {uniquehs___Symbol}, order_Integer, ni_, rule_, string_]/;(
	GenMessage[Length[{Orighs}] === Length[{uniquehs}], iGenGrads::UniqueHeads ]&& 
	GenMessage[Length[obspoints] > dim, iGenGrads::obsPoints]
 ) := Module[
	{dummyvariables, h1, h2, Auxh, ivars, iiresult},
    
    (*Create dummyvariables for TakeGrad:*)
    ivars = ConstantArray["x", Length@obspoints]//Unique;
    
     dummyvariables = Module[
        {result}, (*Maybe there is a functional way for this*)
        Reap[
            result = Sow@ivars[[1;;dim]];
            Do[
                result = Sow@ReplacePart[result, i -> Unique["s"]], 
                {i, 1, dim -1}
            ];
        Sow[ivars];
        
        ]//Last//Last
    ];
        
    (*Create auxiliar heads for TakeGrad*)
    Auxh[x_?OddQ] := h1;
    Auxh[x_?EvenQ] := h2;
    
    (*Define a function of the form f[p1,...,pn,X1, ...,Xm], which evaluates to head@@ivars*)
    Inactive[Set][
        Auxh[ni]@@(Pattern[#, _]&/@ivars),
        head@@ivars
    ]//Activate;
    
    (*Implementing the loop:*)
    iiresult = Table[
		(*implelent rule for Derivative:*)
		Derivative[x__][y_Symbol][k__] := $D[{x},y][k]; 
		Derivative[x__][$D[y_List, s_Symbol]][k__] := $D[{x} + y, s][k];
		
		(*alternate between the heads with the gadient definition*)
		Inactive[Set][
			Auxh[i+1]@@(Pattern[#,_]&/@ivars),
			TakeGrad[Auxh[i], dummyvariables, i-1]
		]//Activate;
      
      (*clear the previous head and the Derivative subValues.*)
       Clear[Evaluate[Auxh[i]]]; 
       SubValues[Derivative] = (SubValues[Derivative])[[1]];
       
       (*calculate the numberical value of the gradients in the form of associations*)
       iCalculateGrads[Auxh[i+1], {Orighs}, {uniquehs}, obspoints, i, dim, rule, string],
       {i, ni, order} 
	
	]; (*note that the output format:  { <|{0}->{...}|>,  <|{1}->..., {2}->...|>,  <|{1,1}->...|>, ...}*)
    
    (*Clean vars and hs and return result*)
    ClearAll[Evaluate[Flatten[dummyvariables]//DeleteDuplicates]];
    ClearAll[h1,h2]; ClearAll[Evaluate@ivars];
    
    Remove[Evaluate[Flatten[dummyvariables]//DeleteDuplicates]]; Remove[h1,h2];
    
    (*I want 1 big association in the form <|{0}-> {...}, {1}->{...}, ...|>*)
    Join@@iiresult
]

iGenGrads[x___] := Throw[$Failed, failTag[iGenGrads]]


(* ::Subsubsection::Closed:: *)
(*GenGrads*)


GenGrads//ClearAll


(*iGenGrads[head_, dim_, obspoints_, {Orighs___}, {uniquehs___}, order_, ni_, rule_, string_]*)


GenGrads[head_Symbol, dim_Integer, obspoints_List, {Orighs___Symbol}, {uniquehs___Symbol}, order_Integer, rule_] := (
	(*1 flat <|{0}-> {...}, {1}->{...}, ...|>*)
	iGenGrads[head, dim, obspoints, {Orighs}, {uniquehs},order, 0, rule, "hphc"]
)


GenGrads[head_, {dims__Integer}, {obspoints__List}, {Orighs___Symbol}, {uniquehs___Symbol}, order_Integer, rule_,  "Detectors"]/;(
	Length[{dims}] === Length[{obspoints}] && Length[{dims}]>= 1
) := Table[ (*A list of flat associations, 1 for each detector.*)
	iGenGrads[head, {dims}[[i]], {obspoints}[[i]], {Orighs}, {uniquehs}, order, 0, rule, "FpFc"],
	{i, 1, Length@{dims}}
]

GenGrads[x___] := Throw[$Failed, failTag[GenGrads]]


(* ::Subsubsection::Closed:: *)
(*PartitionLIComponents*)


FillEmpty[{}] := {0}; FillEmpty[x_List]/; (Length@x >= 1) := x;
FillEmpty[x___] := Throw[$Failed, failTag[FillEmpty]]


ClearAll@PartitionLIComponents


PartitionLIComponents::usage="PartitionLIComponents[LIcomponents_List, {varnumbers__}]
LIComponents: list of LIComponents to subdivide
varnumbers: numbers of variables of each function in sequence
Output: List of lists of corresponding derivatives in each coordinate for each LIcomponent.

Example:
>>>With[
	{
		li = {{1,1},{1,2},{1,3},{1,4},{2,2},{2,3},{2,4},{3,3},{3,4},{4,4}}
	},
	PartitionLIComponents[li, {2,2}]
]
	
>>>{
	{{1,1},{1,2},{1},{1},{2,2},{2},{2},{0},{0},{0}},
	{{0},{0},{3},{4},{0},{3},{4},{3,3},{3,4},{4,4}}
}";
	

PartitionLIComponents[LIcomponents_List, {varnumbers__}]/;(MatrixQ[LIcomponents, NumberQ]):= Module[
	{iList, set, cases},
	
	(*create a list of ranges of vars { {1,..., n1}, {1,...,n2}, ... }*)
	iList = Range/@{varnumbers};
	
	(*displace the variables  {0, n1, n1+n2, ...} + { {1,..., n1}, {1,...,n2}, ... } *)
	iList = (FoldList[Plus, 0, {varnumbers}])[[1;;-2]] + iList;
	
	(*split the sets of variables into set[1], set[2], ...*)
	(set[#] = iList[[#]])&/@Range[Length@iList];
	
	(*apply cases to isolate different sets:*)	
	(cases[#] = Cases[Alternatives@@set[#]]/@LIcomponents)&/@Range[Length@iList];
	
	(*Replace {} -> {0}: *)
	(FillEmpty/@cases[#])&/@Range[Length@iList]
]

PartitionLIComponents[x___] := Throw[$Failed]


ClearAll@MakeGradient


MakeGradient::usage="MakeGradient[{gradients__Association}, {varnumbers__Integer}, LIComponents_]

{gradients}: List of associations with all components of all gradients of a function (the Values of the 
association are either numbers or Lists if the function depends on the frequecy)
{varnumbers}: corresponding list of how many vars each function has
LIComponents: LIComponents of the target gradient
	
Output: list with the values of LIComponents, tipically a matrix where each row is the values of the LIComponent
for all frequency points fi";


Clear@idot
idot[v1_List, v2_List] := v1[[1]]*v2[[1]] + v1[[2]] v2[[2]]


(*to fix the problem of zero symbolic derivatives*)
FixZeroDerivative[x_, l_]/; (x === Transpose[0]&&l>0) := ConstantArray[0, {2, l}]
FixZeroDerivative[x_List, l_] := x


MakeGradient[{gradients__Association}, {varnumbers__Integer}, LIComponents_]/;(
	Length[{gradients}] === 2 (*assuming one is the FpFc gradients and the other the hphc*)
) := Module[
	{partitionedList, result, dFpFc, dhphc, Lengthfvec = {gradients}[[1]][[1]][[1]]//Length, igradients2},
	
	igradients2 = FixZeroDerivative[#, Lengthfvec]&/@{gradients}[[2]];
	
	
	partitionedList = PartitionLIComponents[LIComponents, {varnumbers}];
	
	
	(*Create a symbolic structure that pairs the derivatives of FpFc and hphc to get the full gradient*)
	result = idot@@@(MapThread[
		#1/@#2&,
		{{dFpFc, dhphc}, partitionedList}
	]//Transpose);
	
	(*Transpose to get an structure like {
		{dFpFc[{1}], dhphc[{0}]},
		{dFpFc[{2}], dhphc[{0}]},
		{dFpFc[{0}], dhphc[{1}]},
		...
	}
	and then MapApply idot.
	*)
	
	(*change the symbolic constructs by the actual associations and return the result*)
	dFpFc = {gradients}[[1]];
	dhphc = igradients2;
	
	result (*note that Dimensions[result] === {LI, Length[fvec]}*)
]

MakeGradient[x___] := Throw[$Failed, failTag[MakeGradient]]


GradientsList//Clear


GradientsList::usage="GradientsList[{gradients__}, {varnumbers__Integer}, n_Integer]
{gradients}: associations with all gradient components
{varnumbers}: number of variables of each function in the list
n: max gradient order
Output: List of all Full gradients from order 1 to n.

Example:
>>>Block[
	{f, g},
	Head[f] ^= Association; Head[g] ^= Association; 
	GradientsList[{f,g}, {2,2},  2]
]
>>>{
	{
		f[{1}] g[{0}],
		f[{2}] g[{0}],
		f[{0}] g[{3}],
		f[{0}] g[{4}]
	},
	{
		f[{1,1}] g[{0}],
		f[{1,2}] g[{0}],
		f[{1}] g[{3}],
		f[{1}] g[{4}],
		f[{2,2}] g[{0}],
		f[{2}] g[{3}],
		f[{2}] g[{4}],
		f[{0}] g[{3,3}],
		f[{0}] g[{3,4}],
		f[{0}] g[{4,4}]
	}
}";



GradientsList[{gradients__Association}, {varnumbers__Integer}, n_Integer] :=Module[
	
	{dim = Plus@@{varnumbers}, LIComponents, result},
	
	(LIComponents[#] = SymmetrizedIndependentComponents[
		ConstantArray[dim, #], Symmetric[All]
	])&/@Range[n];
	
	(*Note that Dimensions[result] === {n, LI_i, Length[fvec]}, where LI_i -> number of LI comps. of the order i gradient*)
	result = MakeGradient[{gradients}, {varnumbers}, LIComponents[#]]&/@Range[n];
	
	(*we want Dimensions[result] === {n, Length[fvec], LI_i} for the making of DALI tensors*)
	result = Transpose/@result;
	
	
	
	(*one final association: 
		<|"or1"-> m1, "or2"->m2, ...|>;
		where 
		Dimensions[m1] === {Length[fvec], dim}
		Dimensions[m2] === {Length[fvec], dim*(dim+1)/2 },
		...
	*)
	Association@MapThread[
		Rule,
		{("or" <> ToString[#])&/@Range[n], result}
	]
	
	
]

GradientsList[x___] := Throw[$Failed, failTag[GradientsList]]


(* ::Subsection::Closed:: *)
(*Make DALI Tensors*)


(* ::Subsubsection::Closed:: *)
(*GenDaliTerm*)


(*
it is way faster to use 
	Flatten[Transpose[grads1\[ConjugateTranspose]] . grads2]
but this can return non-symmetrical Fisher matrices. I think it is some numerical error.
*)
CiGenDaliTerm = Compile[
	{{gradList1, _Complex, 2}, {gradList2, _Complex, 2}, {PSD, _Real, 1}},
	Sum[
		
		Outer[Times, Conjugate[gradList1[[i]]], gradList2[[i]]]/PSD[[i]],
		
		{i, 1, Length@gradList1}
	],
	
	CompilationTarget->"C",
	RuntimeOptions->"Speed"
]


GenDaliTerm[gradlist1_, gradlist2_, SensitivityVector_,  \[CapitalDelta]f_] := Module[
    {complexSum, inv = SensitivityVector^-1},
    
    (*Likelihood def. eq. 42 of https://arxiv.org/pdf/1809.02293*)
    
    complexSum = CiGenDaliTerm[gradlist1, gradlist2, SensitivityVector]//Flatten;

	4 \[CapitalDelta]f Re[complexSum]
]

GenDaliTerm[x___] := Throw[$Failed, failTag[GenDaliTerm]];


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
(*GenDaliList*)


(*GenDaliList[list_Association, diag_?TrueQ, matrix_List, dim_] := Module[
    {dummy, DaliList, order = Length@list},
    
    DaliList = Do[
        dummy = Divide[list["or" <> ToString[i]], matrix^2]; 
        Do[
            Sow[#,k]&@GenDaliTerm[list["or"<>ToString[k]], dummy,  matrix,  diag],
            {k, i, order}
         ],
       {i, 1, order}
    ]//Reap//Last;
    
    DaliList
]*)


(*GenDaliList[list_Association, order_Integer, diag_/;(diag==False), matrix_List, dim_] := Module[
    {DaliList},
    
    DaliList = Do[
        Sow[#,k]&@GenDaliTerm[list["or" <>ToString[k]], list["or" <> ToString[i]], matrix, diag],
        {i, 1, order}, {k, i, order}
    ]//Reap//Last;
    
    DaliList
    
]*)


GenDaliList[list_Association, Sn_List, \[CapitalDelta]f_] := Module[
    {DaliList, order = Length@list},
    
    DaliList = Do[
        Sow[#,k]&@GenDaliTerm[list["or" <>ToString[k]], list["or" <> ToString[i]],Sn, \[CapitalDelta]f],
        {i, 1, order}, 
        {k, i, order}
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


(* ::Subsubsection:: *)
(**)


(* ::Section::Closed:: *)
(*DALICoefficients*)


(*DALICoefficients[h_, {vars_List, fp_List, n_Integer}, {Cov_, Diag_?BooleanQ}, SymRules_:<||>, NRules_:<||>]/;(
	(MatrixQ[Cov, NumericQ] && Diag===False) || (VectorQ[Cov, NumericQ] && Diag===True) &&
	VectorQ[vars] && Positive[n]
) := Module[
	{SymRulehs, Orighs, Uniquehs, iNRules, gradients, dim = vars[[-1]]},
	
	(*Make iNRules with Unique heads and set up defs from SymRules and NRules:*)
	SymRulehs = Last[#,#]&/@(Flatten[Values@SymRules][[All, 1, 0]]);
	SymRulehs = DeleteDuplicates[SymRulehs];
	Orighs = Last[#,#]&/@(Flatten[Values@NRules][[All, 1, 0]]);
	Orighs = DeleteDuplicates[Orighs];
	
	Uniquehs = Unique[Orighs];
	iNRules = Replace[NRules, Thread@Rule[Orighs, Uniquehs], {4,5}, Heads->True];
	MakeDefs[SymRules, iNRules];  Clear[iNRules];
	
	
	Unprotect[Derivative]; (*Overload of derivative happening in iGenGrads*)
	
	gradients = GenGrads[h, dim, fp, Orighs, Uniquehs, n, {}];
	
	(*Clean definitions:*)
	Remove[Evaluate[Uniquehs]]; GCSymRules[SymRulehs];  Protect[Derivative];
	
	(*This will give you a list of DALIlists:*)
	
	GenDaliList[gradients, n, Diag, Cov, dim] 
]*)


(*Block[
	{h1,h2,l1 = {1,2,Range[5]}, l2 = {1,2,3,Range[5]}},
	
	h1[x_,y_,z_] = Cosh[x] Sin[y] \[ExponentialE]^z;
	h2[x_,y_,z_, w_] = \[ExponentialE]^(w+x) - Tan[y+z];
	
	GenGrads[h1,2, l1, {}, {}, 3, {}]["or3"]//MatrixForm
	
]*)


(*HEADTEST[x_, y_, z_, w_] =  Sin[x y z]*E^(w+y);

points = RandomReal[{0.,1.}, 10];
\[Sigma] = RandomReal[{0.,1.}, 10];

test = DALICoefficients[HEADTEST[x,y,z,w], {{x,y,z}, {1,2,5}, 4}, {\[Sigma], True}, points];
Map[TensorRank, test, {2}]
Clear[HEADTEST, points, \[Sigma], test]*)


(* ::Section::Closed:: *)
(*GWDALICoefficients*)


Clear@iGWDALICoefficients

iGWDALICoefficients[{h__}, detecs_Integer, {{vars__}, {fp__}, n_Integer}, {f0_, f1_, \[CapitalDelta]f_}, {PSD__}, SymRules_, NRules_, returnAll_?BooleanQ]/;(
	Length[{h}] === Length[{vars}] === (Length[{fp}] - (detecs - 1 )) &&
	f1>f0 && MatrixQ[{PSD}, NumericQ] && {SymRules, NRules}[[All,0]] === {Association,Association} && Length[{PSD}] === detecs
) := Module[

	{SymRulehs, Orighs, Uniquehs, iNRules, idetecHs, dims = {vars}[[All, -1]], detectorGradients, fvec, ObsPoints, remainingGradients, result},
	
	(*Make iNRules with Unique heads and set up defs from SymRules and NRules:*)
	EchoTiming[
		SymRulehs = ToExpression/@Keys[SymRules];
		Orighs = ToExpression/@Keys[NRules];
		Uniquehs = Unique[Orighs];,
		"pre work"
	];
	
	EchoTiming[
		(*Orighs = HoldComplete@@Orighs;
		Uniquehs = HoldComplete@@Uniquehs;
		
		Set@@@Transpose[{List@@Orighs, List@@Uniquehs}];
		
		
		iNRules = Association@(Normal[NRules]); (*weird bug: If I do iNRules = NRules the system does not replace the original heads in the association by the unique ones, so I have to convert to normal form and then back to association.*)
		
		List@@(Clear/@Orighs);
		Orighs = List@@Orighs;
		Uniquehs = List@@Uniquehs;*)
		iNRules = NRules/.Thread@Rule[Orighs, Uniquehs];,
		"parsing NRules"
	];
	
	
	EchoTiming[MakeDefs[SymRules, iNRules]; Clear[iNRules];, "MakeDefs"];
	
	
	(*Define some basic quantities:*)
	EchoTiming[
		fvec = Range[f0,f1, \[CapitalDelta]f];
		ObsPoints = Join[#, {fvec}]&/@{fp};
		"Middle work"
	];
	

	Unprotect[Derivative]; (*Overloading of Derivative happening in iGenGrads*)
	
	
	(*Calculate the detector and remaining gradients:*)
	
			
	EchoTiming[
		detectorGradients = GenGrads[
			{h}[[1]], 
			dims[[1]]&/@Range[detecs], 
			ObsPoints[[1;;-2]], 
			Orighs, 
			Uniquehs, 
			n, 
			{Exp[I anything_]-> 1}, 
			"Detectors"
		];,
		"FpFc"
	];
	
	
	EchoTiming[
		remainingGradients = GenGrads[
			{h}[[2]], 
			dims[[2]], 
			ObsPoints[[-1]], 
			Orighs,
			Uniquehs,
			n, 
			{Exp[I anything_]-> 1}
		];,
		"hphc",
		Method->Timing
	];
	
	
	

	(*Clean definitions:*)
	
	Remove[Evaluate[Uniquehs]]; GCSymRules[SymRulehs]; Protect[Derivative];
	
	(*redefine "remainingGradients" indices {1,2} -> {3,4} and so on*)
	EchoTiming[
		remainingGradients = With[
			{d = dims[[1]]},
			KeyMap[Function[{x}, If[x[[1]] > 0, x + d, x]], remainingGradients]
		];, 
		"Redefine labels"
	];
	
	
	
	(*Calculate full gradients for each detector, you have  a list of associations here*)
	EchoTiming[
		detectorGradients = GradientsList[
			{#, remainingGradients}, dims, n
		]&/@detectorGradients;, 
		"Gradient Recombination"
	];
	
	
	EchoTiming[ClearAll[remainingGradients];, "extra cleaning"];
	
	(*Put the matrices in  the form for DALIList*)
	
	(*This will give you a list of DALIlists:*)
	EchoTiming[
		result = MapThread[
			GenDaliList[#1, #2, \[CapitalDelta]f]&,
			{detectorGradients, {PSD}}
		];,
		"DALIList"
	];
	
	If[
		returnAll === True,
		result,
		Total[result]
	]
]


(* ::Section:: *)
(*DALITensors*)


(* ::Subsubsection::Closed:: *)
(*SymRules and NRules*)


SymRules = <||>;
NRules = <||>;


Module[
	{symPv2, symD, symFpFc, nPv2, nD, nFpFc, m},
	m = Quiet[
		DerivativeRulesLoad/@{(*"IMRPhenomPv2",*) "IMRPhenomD", "Detectors"},
		{Part::partw, Part::take}
	];
	
	{(*SymRules["IMRPhenomPv2"],*) SymRules["IMRPhenomD"], SymRules["FpFc"]} = m[[All,1]];
	{(*NRules["IMRPhenomPv2"],*) NRules["IMRPhenomD"], NRules["FpFc"]} = m[[All,2]];
];


NRules[ii\[ScriptA]IMR] = <|
	ii\[ScriptA]IMR -> {
		ii\[ScriptA]IMR[x__] :> Aux\[ScriptA]IMR1[x],
		TagRule[ii\[ScriptA]IMR, $D[{n__}, ii\[ScriptA]IMR], Aux\[ScriptA]IMR2[n]]
	}
|>;


NRules[Aux\[ScriptA]IMR1] = <|
	Aux\[ScriptA]IMR1 -> {
		Aux\[ScriptA]IMR1[x__] :> Module[
			{},
			Transpose[\[ScriptA]IMR[x]]
		]
	}
|>;

NRules[Aux\[ScriptA]IMR2] = <|
	Aux\[ScriptA]IMR2 -> {
		Aux\[ScriptA]IMR2[n__][y__] :> Module[
			{},
			Transpose[$D[{n}, \[ScriptA]IMR][y]]
		]
	}|>


NRules2 = Join[
	NRules["FpFc"], 
	NRules[ii\[ScriptA]IMR], 
	NRules[Aux\[ScriptA]IMR1], 
	NRules[Aux\[ScriptA]IMR2]
];


Protect[SymRules, NRules];


(* ::Subsubsection:: *)
(*ihphc, iFpFc*)


(* ::Text:: *)
(*Mc  = M \[Eta]^(3/5);*)
(*M = Mc \[Eta]^(-3/5);*)


ihphcIMRPhenomD[
	\[ScriptCapitalM]c_, \[Delta]_, \[Chi]s_, \[Chi]a_,
	\[Iota]_, invdL_, tc_, \[Phi]ref_,
	\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_,
	fref_, f_ 
] = Module[
	{G = 4.9254664969309`3.6105383994801805*^-6 (*G/c^3 [s/solarMass]*), \[Omega], \[Omega]ref, M = \[ScriptCapitalM]c ((1-\[Delta]^2)/4)^(-3/5)},
	
	\[Omega] = f M G; 
	\[Omega]ref = fref M G;
	
	ii\[ScriptA]IMR[f,\[ScriptCapitalM]c,\[Delta],\[Chi]s,\[Chi]a,\[Iota]] Exp[
		-I \[CapitalPhi]IMR[\[Omega],\[Omega]ref,\[Delta],\[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4]
	]*Exp[-I (2 \[Pi] f tc - 2 \[Phi]ref)]*invdL
];


ihphcIMRPhenomPv2[
	m1_, m2_, 
	s1x_, s1y_, s1z_,
	s2x_, s2y_, s2z_,
	\[Iota]_, dL_, tc_, \[Phi]ref_,
	\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_,
	fref_, f_
] = ii\[ScriptA]IMR[f, fref, m1, m2, s1x,s1y,s1z,s2x,s2y,s2z,\[Phi]ref,\[Iota]]*Exp[
	-I \[CapitalPhi]IMR[f,fref,m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4]
]*Exp[-I (2 \[Pi] f tc)]/dL;


iFpFc[
	\[Theta]_, \[Phi]_, \[Psi]_,
	pi_, Dij_,
	f_
] = FpFc[f,\[Theta],\[Phi],\[Psi],pi,Dij];


Protect[ihphcIMRPhenomD, ihphcIMRPhenomPv2, iFpFc(*, \[CapitalPhi]IMR, \[ScriptA]IMR, FpFc*)];


(* ::Subsubsection:: *)
(*Utils for Fisher Matrix*)


Clear@order\[Delta]

MapThread[
	(order\[Delta][#1] = #2)&,
	{
		{"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"},
		Range[15]
	}
]; 


vars["Aligned"] = {
	"\[ScriptCapitalM]c", "\[Delta]", "\[Chi]s", "\[Chi]a",
	"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "1/dL", 
	"tc", "\[Phi]ref",
	"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"
};


vars["Precessing"] = {
	"m1", "m2", "s1x","s1y","s1z","s2x", "s2y", "s2z",
	"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "dL",
	"tc", "\[Phi]ref",
	"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"
};


GenErrorMessage["Aligned", True] := Null;

GenErrorMessage["Aligned", False] := Throw[
"\n
PhenomD variables should satisfy the following conditions: \n
1\[LessEqual]M\[LessEqual]\[Infinity] && \!\(\*SuperscriptBox[\(10\), \(-4\)]\)\[LessEqual]\[Eta]\[LessEqual]0.25`&& -1\[LessEqual]s1z\[LessEqual]1 && -1\[LessEqual]s2z\[LessEqual]1 && 0\[LessEqual]\[Iota]\[LessEqual]\[Pi] && 0\[LessEqual]\[Theta]\[LessEqual]\[Pi] && \n
0\[LessEqual]\[Phi]\[LessEqual]2\[Pi] && 0\[LessEqual]\[Psi]\[LessEqual]\[Pi] && \!\(\*SuperscriptBox[\(10\), \(-11\)]\)\[LessEqual]dL\[LessEqual]\[Infinity] && -\[Infinity]\[LessEqual]tc\[LessEqual]\[Infinity] && 0\[LessEqual]\[Phi]ref\[LessEqual]2 \[Pi]
\n"
];

GenErrorMessage["Precessing", True] := Null;

GenErrorMessage["Precessing", False] := Throw[
"
\n
PhenomPv2 variables should satisfy the following conditions: \n
m1 >= m2 && Norm[{s1x, s1y,s1z}]<=1  && Norm[{s2x,s2y, s2z}]<=1 \n
1\[LessEqual] M \[LessEqual]\[Infinity] && \!\(\*SuperscriptBox[\(10\), \(-4\)]\)\[LessEqual]\[Eta]\[LessEqual] 0.25`&& -1\[LessEqual]s1x\[LessEqual]1 && -1\[LessEqual]s1y\[LessEqual]1 && -1\[LessEqual]s1z\[LessEqual]1 &&\n
-1\[LessEqual]s2x\[LessEqual]1 && -1\[LessEqual]s2y\[LessEqual]1 && -1\[LessEqual]s2z\[LessEqual]1 && 0\[LessEqual]\[Iota]\[LessEqual]\[Pi] && 0\[LessEqual]\[Theta]\[LessEqual]\[Pi] &&\n
0\[LessEqual]\[Phi]\[LessEqual]2\[Pi] && 0\[LessEqual]\[Psi]\[LessEqual]\[Pi] && \!\(\*SuperscriptBox[\(10\), \(-11\)]\)\[LessEqual]dL\[LessEqual]\[Infinity] && -\[Infinity]\[LessEqual]tc\[LessEqual]\[Infinity] && 0\[LessEqual]\[Phi]ref\[LessEqual]2 \[Pi] \n"
];


RetrieveFiducial[aligned_String, fp_Association]/;(
	(*there should not be more than 12 variables*)
	Length@Keys[fp] <= 26 && 
	(*the keys must be contained in vars["Aligned"]*)
	ContainsAll[vars["Aligned"],Keys[fp]] &&
	
	aligned == "IMRPhenomD"
) := Module[
	{\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], \[Theta], \[Phi], \[Psi], invdL, tc, \[Phi]ref, \[Delta]p, test, res},
	
	{\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], \[Theta], \[Phi], \[Psi], invdL, tc, \[Phi]ref} = fp/@{
		"\[ScriptCapitalM]c", "\[Delta]", "\[Chi]s", "\[Chi]a",
		"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "1/dL",
		"tc", "\[Phi]ref"
	};
	
	
	(*Test variables:*)
	test = Thread@LessEqual[
		{10^-3, 0, -1, -1, 0, 0, 0,0, 10^-11, -Infinity, 0},
		{\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], \[Theta], \[Phi], \[Psi], invdL, tc, \[Phi]ref},
		{Infinity, 0.9999, 1,1, \[Pi], \[Pi], 2 \[Pi], \[Pi], 10^20,  Infinity, 2 \[Pi]}
	];
	
	(*this collapses to True or False*)
	test = (And@@test)//TrueQ;
	
	GenErrorMessage["Aligned", test];
	
	If[
		Length[Keys[fp]] === 11,
		res = {{\[Theta], \[Phi], \[Psi]}, {\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], invdL, tc, \[Phi]ref}},
		
		\[Delta]p = DeleteElements[Keys[fp], vars["Aligned"][[1;;11]]];(*check for real value*)
		\[Delta]p = fp/@SortBy[\[Delta]p, order\[Delta]];
		
		If[VectorQ[\[Delta]p, RealValuedNumberQ] === False, Throw["\[Delta]pi values should be in the range -\[Infinity]< \[Delta]pi <\[Infinity]"]];
		
		res = {{\[Theta], \[Phi], \[Psi]}, {\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], invdL, tc, \[Phi]ref, Sequence@@\[Delta]p}}
	]
	
]


RetrieveFiducial[precessing_String, fp_Association]/;(
	Length@Keys[fp] <= 30 && 
	
	ContainsAll[vars["Precessing"],Keys[fp]] &&
	
	precessing === "IMRPhenomPv2"
) := Module[
	{m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], \[Theta], \[Phi], \[Psi], dL, tc, \[Phi]ref, \[Delta]p, test, test2},
	
	{
	m1,m2, 
	s1x, s1y, s1z, 
	s2x, s2y, s2z, 
	\[Iota], \[Theta], \[Phi], \[Psi], dL, tc, \[Phi]ref} = fp/@{
		"m1", "m2", 
		"s1x", "s1y","s1z", 
		"s2x", "s2y", "s2z", 
		"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "dL", 
		"tc", "\[Phi]ref"
	};
	
	If[m1==m2, m2 = (0.9999) m1];
	
	test = Thread@LessEqual[
		{1, 1, Sequence@@ConstantArray[-1,6], 0,0,0,0, 10^-11, -Infinity, 0},
		{m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], \[Theta], \[Phi], \[Psi], dL, tc, \[Phi]ref},
		{Infinity, Infinity, Sequence@@ConstantArray[1, 6], \[Pi],\[Pi], 2 \[Pi], \[Pi], Infinity, Infinity, 2 \[Pi]}
	];
	
	test = (And@@test)//TrueQ;
	
	test2 = (m1 >= m2 && Norm[{s1x, s1y, s1z}] <= 1 && Norm[{s2x, s2y, s2z}] <= 1)//TrueQ;
	
	GenErrorMessage["Precessing", TrueQ[test&&test2]];
	
	
	If[
		Length[Keys[fp]] === 15,
		{{\[Theta], \[Phi], \[Psi]}, {m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota],dL, tc, \[Phi]ref}},
		\[Delta]p = DeleteElements[Keys[fp], vars["Precessing"][[1;;15]]]; 
		\[Delta]p = fp/@SortBy[\[Delta]p, order\[Delta]];
		If[VectorQ[\[Delta]p, RealValuedNumberQ]===False, Throw["\[Delta]pi value should be in the range -\[Infinity] < \[Delta]pi <\[Infinity]"]];
		{{\[Theta], \[Phi], \[Psi]}, {m1,m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota],dL, tc, \[Phi]ref, Sequence@@\[Delta]p}}
	]
]

RetrieveFiducial[x___] := Throw[$Failed, failTag[RetrieveFiducial]]


\[Delta]pi = {\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4};


(* ::Text:: *)
(*I don' t know why, but there is a bug in which ```ToExpression``` converts to Global context not the *)
(*private one from the package, so I need to set manually*)


MapThread[
	(iToExpression[#1] =#2)&, 
	{
		{"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"},
		\[Delta]pi
	}
]



\[Delta]pRules = Thread@Rule[
	\[Delta]pi,
	ConstantArray[0, 15]
];


hphcVarNumber["Aligned"] = 8;
hphcVarNumber["Precessing"] = 12;


LIDij = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]];


Clear@headhphc
headhphc["IMRPhenomD"] = ihphcIMRPhenomD;
headhphc["IMRPhenomPv2"] = ihphcIMRPhenomPv2;


iihphcvars["IMRPhenomPv2", fref_] = Join[
	{m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Iota], dL,tc, \[Phi]ref},
	\[Delta]pi,
	{fref, f}
]


iihphcvars["IMRPhenomD", fref_] = Join[
	{\[ScriptCapitalM]c,\[Delta],\[Chi]s,\[Chi]a,\[Iota],invdL,tc,\[Phi]ref},
	\[Delta]pi,
	{fref, f}
]


isAligned["IMRPhenomD"] = True
isAligned["IMRPhenomPv2"] = False


(* ::Subsubsection:: *)
(*Fisher Matrix*)


ClearAll[DALITensors]


make1[0] := 0;
make1[n_]/; n>0 := 1


iconvert["IMRPhenomD"] = "Aligned"
iconvert["IMRPhenomPv2"] = "Precessing"


iDALITensors[
	approximant_String, 
	fiducialPoint_Association, 
	{detectors__Association},
	n_,
	fmin_, fmax_, res_, AllFisherMatrices_
]/;(
	approximant === "IMRPhenomD" ||approximant ===  "IMRPhenomPv2" 
) := Module[
	{
		FiducialFpFc,  Fiducialhphc,
		varsFpFc, ivarshphc, varshphc, iL, \[Delta]pKey, i\[Delta]pRules, ihphc,
		 \[CapitalDelta]f, fvec, PSDs, keysFP = Keys[fiducialPoint], fpDetectors,IFMAX, totalM, alignedOrPrecessing
	},
	(*set the extra keys in order:*)
	
	alignedOrPrecessing = iconvert[approximant];
	
	
	(*get fiducial points*)
	{FiducialFpFc,  Fiducialhphc} = With[
	
		{l = RetrieveFiducial[approximant, fiducialPoint]},
		If[
			Cases[l, Missing, Infinity, Heads->True]==={},
			l,
			Throw["Wrong variables for the fiducial point"]
		]
	];
	
	
	
	varsFpFc = {\[Theta], \[Phi], \[Psi], pi, Dij, f, 3};
	
	ivarshphc = If[
		approximant==="IMRPhenomD",
		iL = DeleteElements[keysFP, vars["Aligned"][[1;;11]]];
		iL = iToExpression/@SortBy[iL, order\[Delta]];
		{\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota], invdL, tc, \[Phi]ref},
		iL = DeleteElements[keysFP, vars["Precessing"][[1;;15]]]; 
		iL = iToExpression/@SortBy[iL, order\[Delta]];
		{m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], dL, tc, \[Phi]ref}
	];
	
	(*Extra \[Delta]p var if it exists*)
	Quiet[\[Delta]pKeys = iL, Set::shape];
	
	(*this will set the variables in case Length[iL] = 0 || Length[iL] =1*)
	
	varshphc[0] = Join[ivarshphc, {f}, {hphcVarNumber[alignedOrPrecessing]}];
	varshphc[1] = Join[ivarshphc, iL, {f}, {hphcVarNumber[alignedOrPrecessing]+Length[iL]}];
	
	
	(*only modifies the rules associated to requested \[Delta]ps*)
	i\[Delta]pRules = DeleteElements[\[Delta]pRules, Thread@Rule[iL, ConstantArray[0, Length[iL]]]];
	
	
	Module[
		{il},
		If[Length[iL]==0, il=0, il=1];
		
		Set@@(set[
			ihphc@@(
				Pattern[#,  Blank[]]&/@varshphc[il][[1;;-2]]
			),
			(headhphc[approximant]@@iihphcvars[approximant, fmin])/.i\[Delta]pRules
		])	
	];
	
	(*Check detector keys*)
	If[
		Cases[{#["Position"], #["DetectorTensor"], #["ASD"]}&/@{detectors}, Missing, Infinity, Heads->True]==={},
		0,
		Throw["Wrong Detector Keys"]
	];
	
	fpDetectors = Join[
		FiducialFpFc, 
		{#["Position"], Extract[#["DetectorTensor"], LIDij]}
	]&/@{detectors};
	
	
	(*this should be an array of real numbers*)
	If[
		VectorQ[Flatten[fpDetectors],RealValuedNumberQ]===False,
		Throw["DetectorTensor and Position should contain only real numbers"]
	];
	
	
	totalM = If[
		isAligned[approximant], 
		Fiducialhphc[[1]] ((1-Fiducialhphc[[2]]^2)/4)^(-3/5),  (*M  = \[ScriptCapitalM]c \[Eta]^(-3/5)*)
		Fiducialhphc[[1]]+Fiducialhphc[[2]]
	];
	
	IFMAX = Min@{fmax, 0.2/(4.9254664969309`3.6105383994801805*^-6 totalM)}//Round;
	
	\[CapitalDelta]f = (IFMAX-fmin)/(res-1);
	
	fvec = Range[fmin, IFMAX, \[CapitalDelta]f];
	
	PSDs = (#["ASD"][fvec])&/@{detectors}; 
	PSDs = PSDs^2;
	
	(*iGWDALICoefficients[{h__}, detecs_, {{vars__}, {fp__}, n_}, {f0_, f1_, \[CapitalDelta]f_}, {PSD__}, SymRules, NRules, returnAll]*)
	
	
	iGWDALICoefficients[
		{iFpFc, ihphc}, 
		Length[{detectors}],
		{
			{varsFpFc, varshphc[make1[Length[iL]]]},
			Join[fpDetectors, {Fiducialhphc}],
			n
		},
		{fmin, IFMAX, \[CapitalDelta]f},
		PSDs,
		SymRules[approximant],
		Join[NRules2, NRules[approximant]],
		AllFisherMatrices
	]
]


Options[DALITensors] = {
	"res" -> 1000,
	"fmin" -> 10,
	"fmax" -> 1024,
	"AllFisherMatrices" -> False
};


DALITensors[x__, OptionsPattern[]] := Module[
	{fmin, fmax, res, AllFisherMatrices},
	{fmin, fmax, res, AllFisherMatrices} = OptionValue[DALITensors, #]&/@{"fmin", "fmax", "res", "AllFisherMatrices"};
	Catch[iDALITensors[x, fmin, fmax, res, AllFisherMatrices]]
]


Protect[DALITensors]


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
