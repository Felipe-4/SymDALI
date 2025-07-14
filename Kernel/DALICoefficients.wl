(* ::Package:: *)

(*SetOptions[EvaluationNotebook[], DefaultNewCellStyle->"Code"]
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
(*Package Header*)


BeginPackage["FelipeBarbosa`SymDALI`DALICoefficients`"];


DALICoefficients::usage="DALICoefficients[h_, {vars_List, fp_List, n_Integer}, {Cov_, Diag_?BooleanQ}, SymRules_:<||>, NRules_:<||>]";


GWDALICoefficients::usage="
Assuming
                      ln(\[ScriptCapitalL]) = - \!\(\*FractionBox[\(1\), \(2\)]\)(g-h, g-h)   and  (u,\[Nu]) \[Congruent] 4 \[ScriptCapitalR] \!\(\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Infinity]\)]\)  \!\(\*FractionBox[\(\*OverscriptBox[\(u\), \(~\)]\[Conjugate]  \((\[ScriptF])\)\\\ \*OverscriptBox[\(\[Nu]\), \(~\)] \((\[ScriptF])\)\), \(\*SubscriptBox[\(S\), \(n\)] \((\[ScriptF])\)\)]\) d\[ScriptF],
the \!\(\*
StyleBox[\"call\", \"Code\"]\)
             GWDALICoefficients[\!\(\*OverscriptBox[\(h\), \(~\)]\), {{\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, {\!\(\*SubscriptBox[\(p0\), \(1\)]\), \!\(\*SubscriptBox[\(p0\), \(2\)]\),...}, n}, Sn, {fmin, fmax, \[CapitalDelta]f}],
generates the list of coefficients for the DALI expansion with respect to {\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...}, around {\!\(\*SubscriptBox[\(p0\), \(1\)]\),\!\(\*SubscriptBox[\(p0\), \(2\)]\),..}
to order n in derivatives, where Sn = {\!\(\*SubscriptBox[\(S\), \(n\)]\)(\!\(\*SubscriptBox[\(f\), \(1\)]\)), \!\(\*SubscriptBox[\(S\), \(n\)]\)(\!\(\*SubscriptBox[\(f\), \(2\)]\)), ...} and \!\(\*OverscriptBox[\(h\), \(~\)]\) = \!\(\*OverscriptBox[\(h\), \(~\)]\)(\!\(\*SubscriptBox[\(p\), \(\[Alpha]\)]\),f).
";

iGWDALICoefficients::usage="iGWDALICoefficients[{h__}, detecs_Integer, {{vars__}, {fp__}, n_Integer}, {f0_, f1_, \[CapitalDelta]f_}, PSD_, SymRules_, NRules_]";


CalculateGrads::usage="iCalculateGrads[head_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List,tensorRank_Integer, dim]

head: function head, such that head[p1,...,pn, X1,...,Xm] gives the array.
Orighs: heads that appear in the explicit expression head[p1,...,pn, X1,...,Xm].
uniquehs: uniqueHeads for which UpValues were defined from NRules. 
OpsPoints: set of values to pass to the function in the form of {pf1,pf2,..., XV1, XV2, ...}
tensorRank: tensorRank of head[p1,...,pn, X1,...,Xm].
dim: dimension of the square array

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
	CalculateGrads[h, {},{}, obs, 2, 2]
]
>>>9
>>>\[LeftAssociation]{1,1}\[Rule]{9,18,27,36,45},{1,2}\[Rule]{9,36,81,144,225},{2,2}\[Rule]{9,9 \!\(\*SqrtBox[\(2\)]\),9 \!\(\*SqrtBox[\(3\)]\),18,9 \!\(\*SqrtBox[\(5\)]\)}\[RightAssociation]

Example2:
>>>Block[
	{h, obs = {2,3,4, Range[0,5]}},
	h[x_,y_,z_,w_] = (x+y+z) w;
	Print[2+3+4];
	CalculateGrads[h, {},{}, obs,0, 3]
]
>>>9
>>>\[LeftAssociation]{0}\[Rule]{0,9,18,27,36,45}\[RightAssociation]";


Begin["`Private`"];


(* ::Section:: *)
(*Definitions*)


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


(* ::Subsection:: *)
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


(*Assuming the lhs is either Condition[$D[{n__}, h], cond] -> something or $D[{n__}, h]-> something*)
GetSymRuleHead[def_Rule]/;def[[1,0]] === Condition := def[[1, 1,-1]]
GetSymRuleHead[def_Rule]/; def[[1,0]] === $D := def[[1,-1]]
GetSymRuleHead[x___] :=Throw[$Failed, failTag[GetSymRuleHead]]


SymUpValues::usage="SymUpValues[{SymRules___Rule}, AuxHead_]
{SymRules}: Flat SymRules List, i.e. Flatten[Values@SymRules]
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

SymUpValues[{}, AuxHead_Symbol] := Null

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


ProcessSymRules::usage="ProcessSymRules[SymRules_Association]
SymRules: Association with SymRules for all heads
Output: list of all heads such that head[x___] = number_i, for some number_i in the SymRules";

ProcessSymRules[SymRules_Association] := Module[
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
	SymUpValues[#1, aux]&/@SymRules;
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
{rules}: list of rules of a particular head. The elements are either Rule[Condition[], rhs] or 
TagRule[tag, Condition[], rhs]
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


(* ::Subsubsection:: *)
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

TakeGrad[functionhead_Symbol,  dummyvariables_List, tensorRank_Integer]/;(
	functionhead[[0]] === Symbol && Positive[tensorRank] && DeleteDuplicates[Flatten[dummyvariables][[All,0]]] === {Symbol}
) := MapIndexed[
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


(* ::Subsubsection:: *)
(*Calculate numerical gradients*)


iCalculateGrads::usage="iCalculateGrads[head_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List,tensorRank_Integer, dim]

head: function head, such that head[p1,...,pn, X1,...,Xm] gives the array.
Orighs: heads that appear in the explicit expression head[p1,...,pn, X1,...,Xm].
uniquehs: uniqueHeads for which UpValues were defined from NRules. 
OpsPoints: set of values to pass to the function in the form of {pf1,pf2,..., XV1, XV2, ...}
tensorRank: tensorRank of head[p1,...,pn, X1,...,Xm].
dim: dimension of the square array

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


iCalculateGrads[head_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List,tensorRank_Integer, dim_, rule_]/;(
	Length[{Orighs}] === Length[{uniquehs}] && Positive@tensorRank && 
	Length[ObsPoints] === Length@FunctionVariables[head]
) := Module[
    {dummyhead, LIComponents, x, dummyvariables, result},
    
    (*Construct pattern variables*)
    dummyvariables = Unique@ConstantArray[x, Length[ObsPoints]];
    
    LIComponents = SymmetrizedIndependentComponents[ConstantArray[dim, tensorRank], Symmetric[All]];
    
    Inactive[Set][(*make dummyhead[x1_, x2_,..., head1_, head2_,...]*)
		dummyhead@@(Pattern[#,_]&/@Join[dummyvariables, {Orighs}]),
		Extract[head@@dummyvariables, LIComponents]//. rule
	]//Activate;
    
    
    result = dummyhead@@Join[ObsPoints, {uniquehs}];
    
    
    (*Clean the dummy variables and return result:*)
    Remove[Evaluate[dummyvariables]];
    
    Association@MapThread[
		Rule,
		{LIComponents, result}
    ]
    
]//Check[#, Throw["TakeGradFail at order" <> ToString[tensorRank]]]&

iCalculateGrads[h_Symbol, {Orighs___Symbol}, {uniquehs___Symbol}, ObsPoints_List, 0, dim_, rule_]/;(
	Length[{Orighs}] === Length[{uniquehs}]
) := Module[
	{dummyhead, dummyvariables, result,x},
    
    (*Construct pattern variables including {Orighs}*)
    dummyvariables = Unique@ConstantArray[x, Length[ObsPoints]];
    
    Inactive[Set][(*make dummyhead[x1_, x2_,..., head1_, head2_,...]*)
		dummyhead@@(Pattern[#,_]&/@Join[dummyvariables, {Orighs}]),
		h@@dummyvariables//. rule
	]//Activate;
    
   result = dummyhead@@Join[ObsPoints, {uniquehs}];
   
   (*Clean the dummy variables and return result:*)
    Remove[Evaluate[dummyvariables]];
   
   <|{0} -> result|>
]//Check[#, Throw["TakeGradFail at order" <> ToString[tensorRank]]]&

iCalculateGrads[x___] := Throw[$Failed, failTag[iCalculateGrads]]


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
	
iGenGrads[head_Symbol, dim_Integer, obspoints_List, {Orighs___Symbol}, {uniquehs___Symbol}, order_Integer, ni_, rule_]/;(
	GenMessage[Positive[order],iGenGrads::negativeOrder ] &&
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
    
    (*Define a function of the form f[p1,...,pn,X1, ...,Xm] to act on:*)
    Inactive[Set][
        Auxh[ni]@@(Pattern[#, _]&/@ivars),
        head@@ivars
    ]//Activate;
    
    (*Implementing the loop:*)
    iiresult = Table[
    (*implelent rule for Derivative:*)
    Derivative[x__][y_Symbol][k__] := $D[{x},y][k]; Derivative[x__][$D[y_List, s_Symbol]][k__] := $D[{x} + y, s][k];
    
        Inactive[Set][
			Auxh[i+1]@@(Pattern[#,_]&/@ivars),
			TakeGrad[Auxh[i], dummyvariables, i-1]
		]//Activate;
      
       Clear[Evaluate[Auxh[i]]]; SubValues[Derivative] = (SubValues[Derivative])[[1]];
       
       iCalculateGrads[Auxh[i+1], {Orighs}, {uniquehs}, obspoints, i, dim, rule],
       {i, ni, order} 
    ];
    
    (*Clean vars and hs and return result*)
    Remove[Evaluate[Flatten[dummyvariables]//DeleteDuplicates]];
    Remove[h1,h2];
    
    iiresult
]

iGenGrads[x___] := Throw[$Failed, failTag[iGenGrads]]


GenGrads::usage="GenGrads[head, dim, obspoints, {Orighs___}, {uniquehs___}, n_]
Output: iGenGrads[head, dim, obspoints, {Orighs___}, {uniquehs___}, n_, 1]


GenGrads[{heads__}, {dims__}, {obspoints__}, {Orighs___}, {uniquehs___}, n_]
Evaluates iGenGrads for the sequence of heads, dimensions and obspoints provided, including the \"0\" order derivative.
Output: List of \"Length[{heads__}]\" associations where each one contains all the gradients of \"{head}[[i]]\".

Example:
>>>Block[
	{h1, h2, l1 = {1,2, Range[5]}, l2 = {1,2,3,Range[5]}},
	
	h1[x_, y_, z_] = Cosh[x]+Sin[y] - \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(z\)]\);
	h2[x_,y_,z_, w_] = \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(w + x\)]\) - Tan[y + z];
	
	GenGrads[{h1,h2}, {2,3}, {l1,l2}, {}, {}, 2]
]

>>>{
	\[LeftAssociation]
		{0}\[Rule]-1+Cosh[1]+Sin[2],
		{1}\[Rule]Sinh[1],
		{2}\[Rule]Cos[2],
		{1,1}\[Rule]Cosh[1],
		{1,2}\[Rule]0,
		{2,2}\[Rule]-Sin[2]
	\[RightAssociation],
	\[LeftAssociation]
		{0}\[Rule]1-Tan[5],
		{1}\[Rule]{\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(2\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(3\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(4\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(5\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(6\)]\)},
		{2}\[Rule]-Sec[5\!\(\*SuperscriptBox[\(]\), \(2\)]\),
		{3}\[Rule]-Sec[5\!\(\*SuperscriptBox[\(]\), \(2\)]\),
		{1,1}\[Rule]{\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(2\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(3\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(4\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(5\)]\),\!\(\*SuperscriptBox[\(\[ExponentialE]\), \(6\)]\)},
		{1,2}\[Rule]0,
		{1,3}\[Rule]0,
		{2,2}\[Rule]-2 Sec[5\!\(\*SuperscriptBox[\(]\), \(2\)]\) Tan[5],
		{2,3}\[Rule]-2 Sec[5\!\(\*SuperscriptBox[\(]\), \(2\)]\) Tan[5],
		{3,3}\[Rule]-2 Sec[5\!\(\*SuperscriptBox[\(]\), \(2\)]\) Tan[5]
	\[RightAssociation]
}";


GenGrads[head_Symbol, dim_Integer, obspoints_List, {Orighs___Symbol}, {uniquehs___Symbol}, order_Integer, rule_] := Module[
	{result, labels},
	
	result = iGenGrads[head, dim, obspoints, {Orighs}, {uniquehs}, order, 1, rule];
	
	(*Transform the list of associations { <|{1} -> ...,{2} -> ...,  ...|>, <|{1,1} -> ..., {1,2} -> ..., ... |>},
	into 1 association of matrices with dimensions {number of obs points, number of LIComponents}*)
	
	result = Transpose[Values[#]]&/@result;
	labels = ("or" <> ToString[#])&/@Range[order];
	
	Association@MapThread[
		Rule,
		{labels, result}
	]
		
]

GenGrads[{heads__Symbol}, {dims__Integer}, {obspoints__List}, {Orighs___Symbol}, {uniquehs___Symbol}, order_Integer, rule_]/;(
	Length[{heads}] === Length[{dims}] === Length[{obspoints}] && Length[{heads}]> 1
) := With[
	{},
	
	Table[
		Association@@iGenGrads[{heads}[[i]], {dims}[[i]], {obspoints}[[i]], {Orighs}, {uniquehs}, order, 0, rule],
		{i, 1, Length@{heads}}
	]
]

GenGrads[x___] := Throw[$Failed, failTag[GenGrads]]


FillEmpty[{}] := {0}; FillEmpty[x_List]/; (Length@x >= 1) := x;
FillEmpty[x___] := Throw[$Failed, failTag[FillEmpty]]


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


MakeGradient::usage="MakeGradient[{gradients__Association}, {varnumbers__Integer}, LIComponents_]

{gradients}: List of associations with all components of all gradients of a function (the Values of the 
association are either numbers or Lists if the function depends on the frequecy)
{varnumbers}: corresponding list of how many vars each function has
LIComponents: LIComponents of the target gradient
	
Output: list with the values of LIComponents, tipically a matrix where each row is the values of the LIComponent
for all frequency points fi";


MakeGradient[{gradients__Association}, {varnumbers__Integer}, LIComponents_] := Module[
	{partitionedList, result},
	
	
	partitionedList = PartitionLIComponents[LIComponents, {varnumbers}];
	
	result = MapThread[
		#1/@#2&,
		{{gradients}, partitionedList}
	]; (*{
			{Ass1comps__}, {Ass2comps__}, ...
		} if Assi came from hi that depend on freq. {Assicomps__} is a matrix. Changing the overall head by Times,
			we can do {Ass1comps__}*{Ass2comps__}*... which will return the list of LI components of the full gradient
			evaluated at all frequency points for each component.
		*)
	
	Times@@result	
]

MakeGradient[x___] := Throw[$Failed, failTag[MakeGradient]]


FillEmpty[{}] := {0}; FillEmpty[x_List]/; (Length@x >= 1) := x;
FillEmpty[x___] := Throw[$Failed, failTag[FillEmpty]]


PartitionLIComponents::usage="
PartitionLIComponents[LIcomponents_List, {varnumbers__}]

LIComponents: list of LIComponents to subdivide
varnumbers: numbers of variables of each function in sequence
	
Output: List of lists of corresponding derivatives in each coordinate for each LIcomponent.

Example:
With[
	{
		li = {{1,1},{1,2},{1,3},{1,4},{2,2},{2,3},{2,4},{3,3},{3,4},{4,4}}
	},
	PartitionLIComponents[li, {2,2}]
]
Returns:	
	{
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


MakeGradient::usage="MakeGradient[{gradients__Association}, {varnumbers__Integer}, LIComponents_]
{gradients}: List of associations with all components of all gradients of a function (the Values of the 
association are either numbers or vectors)
{varnumbers}: corresponding list of how many vars each function has
LIComponents: LIComponents of the target gradient
Output: list with the values of LIComponents.

Example:
>>>Block[
	{h1, h2, LI},
	Head[h1] ^= Association; Head[h2] ^= Association; 
	LI = Echo[SymmetrizedIndependentComponents[{4,4}, Symmetric[All]]];
	MakeGradient[{h1,h2}, {2,2}, LI]
]
>>>{{1,1},{1,2},{1,3},{1,4},{2,2},{2,3},{2,4},{3,3},{3,4},{4,4}}
>>>{
	h1[{1,1}] h2[{0}],
	h1[{1,2}] h2[{0}],
	h1[{1}] h2[{3}],
	h1[{1}] h2[{4}],
	h1[{2,2}] h2[{0}],
	h1[{2}] h2[{3}],
	h1[{2}] h2[{4}],
	h1[{0}] h2[{3,3}],
	h1[{0}] h2[{3,4}],
	h1[{0}] h2[{4,4}]
}";

MakeGradient[{gradients__}, {varnumbers__Integer}, LIComponents_]/;(
	DeleteDuplicates@(Head/@{gradients}) === {Association}
):= Module[
	{partitionedList, result},
	
	
	partitionedList = PartitionLIComponents[LIComponents, {varnumbers}];
	
	result = MapThread[
		#1/@#2&,
		{{gradients}, partitionedList}
	]; (*{
			{Ass1comps__}, {Ass2comps__}, ...
		} if Assi came from hi that depend on freq. {Assicomps__} is a matrix. Changing the overall head by Times,
			we can do {Ass1comps__}*{Ass2comps__}*... which will return the list of LI components of the full gradient
			evaluated at all frequency points for each component.
		*)
	
	Times@@result	
]

MakeGradient[x___] := Throw[$Failed, failTag[MakeGradient]]


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

GradientsList[{gradients__}, {varnumbers__Integer}, n_Integer]/;(
	DeleteDuplicates[Head/@{gradients}] === {Association}
) :=Module[
	{dim = Plus@@{varnumbers}, LIComponents, result},
	
	(LIComponents[#] = SymmetrizedIndependentComponents[
		ConstantArray[dim, #], Symmetric[All]
	])&/@Range[n];
	
	result = MakeGradient[{gradients}, {varnumbers}, LIComponents[#]]&/@Range[n];
	
	(*Transpose the matrices representing gradients and make an association:*)
	result = Transpose/@result;
	
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


GenDaliTerm[grads1_List, grads2_List, matrix_List, True]/;(
	MatrixQ[grads1, NumberQ] &&MatrixQ[grads2, NumberQ]
) := Flatten[Transpose[grads1] . grads2]


GenDaliTerm[grads1_, grads2_, matrix_, False]/;(
	MatrixQ[grads1, NumberQ] &&MatrixQ[grads2, NumberQ]
) := Flatten[Transpose[grads1] . matrix . grads2]


iGenDaliTerm[gradlist1_, gradlist2_, SensitivityVector_] := Sum[
	KroneckerProduct[Conjugate[gradlist1[[i]]],gradlist2[[i]]]/SensitivityVector[[i]],
	{i, 1, Length@gradlist1}
]//Flatten


CiGenDaliTerm = Compile[
	{{gradList1, _Complex, 2}, {gradList2, _Complex, 2}, {PSD, _Real,1}},
	Sum[
	
		Outer[Times, gradList1[[i]], gradList2[[i]]]/PSD[[i]],
		{i, 1, Length@gradList1}
	],
	CompilationTarget->"C",
	RuntimeOptions->"Speed"
]


GenDaliTerm[gradlist1_, gradlist2_, SensitivityVector_,  \[CapitalDelta]f_, "GWs"] := Module[
    {complexSum, inv = SensitivityVector^-1},
    (*Likelihood def. eq. 42 of https://arxiv.org/pdf/1809.02293*)
    (*maybe it would still be faster to do this part with high work precision...*)
    
    (*complexSum = Flatten[gradlist1\[ConjugateTranspose].(gradlist2*inv)];*) (*This dot product is more efficient, but makes Fisher assymetric on tc, \[Phi]c and dL because of numerical errors*)
   
    
    complexSum = CiGenDaliTerm[gradlist1, gradlist2, SensitivityVector]//Flatten;
   
     4 \[CapitalDelta]f Re[complexSum]  (*I think this 4 \[CapitalDelta]f can be just absorbed in the normalization...*)
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


GenDaliList[list_Association, order_Integer, diag_?TrueQ, matrix_List, dim_] := Module[
    {dummy, DaliList},
    
    DaliList = Do[
        dummy = Divide[list["or" <> ToString[i]], matrix^2]; 
        Do[
            Sow[#,k]&@GenDaliTerm[list["or"<>ToString[k]], dummy,  matrix,  diag],
            {k, i, order}
         ],
       {i, 1, order}
    ]//Reap//Last;
    
    DaliList
]


GenDaliList[list_Association, order_Integer, diag_/;(diag==False), matrix_List, dim_] := Module[
    {DaliList},
    
    DaliList = Do[
        Sow[#,k]&@GenDaliTerm[list["or" <>ToString[k]], list["or" <> ToString[i]], matrix, diag],
        {i, 1, order}, {k, i, order}
    ]//Reap//Last;
    
    DaliList
    
]


GenDaliList[list_Association, order_Integer, Sn_List, \[CapitalDelta]f_, "GW"] := Module[
    {DaliList},
    
    DaliList = Do[
        Sow[#,k]&@GenDaliTerm[list["or" <>ToString[k]], list["or" <> ToString[i]],Sn, \[CapitalDelta]f, "GWs"],
        {i, 1, order}, {k, i, order}
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


DALICoefficients[h_, {vars_List, fp_List, n_Integer}, {Cov_, Diag_?BooleanQ}, SymRules_:<||>, NRules_:<||>]/;(
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
]


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

iGWDALICoefficients[{h__}, detecs_Integer, {{vars__}, {fp__}, n_Integer}, {f0_, f1_, \[CapitalDelta]f_}, {PSD__}, SymRules_, NRules_]/;(
	Length[{h}] === Length[{vars}] === (Length[{fp}] - (detecs - 1 )) &&
	f1>f0 && MatrixQ[{PSD}, NumericQ] && {SymRules, NRules}[[All,0]] === {Association,Association} && Length[{PSD}] === detecs
) := Module[
	{SymRulehs, Orighs, Uniquehs, iNRules, idetecHs, dims = {vars}[[All, -1]], detectorGradients, fvec, ObsPoints, remainingGradients, result},
	
	(*Make iNRules with Unique heads and set up defs from SymRules and NRules:*)
	EchoTiming[
	SymRulehs = ToExpression/@Keys[SymRules];
	Orighs =ToExpression/@Keys[NRules];
	
	Uniquehs = Unique[Orighs];,
	"pre work"];
	
	EchoTiming[
		Orighs = HoldComplete@@Orighs;
		Uniquehs = HoldComplete@@Uniquehs;
		
		Set@@@Transpose[{List@@Orighs, List@@Uniquehs}];
		
		
		iNRules = Association@(Normal[NRules]); (*weird bug: If I do iNRules = NRules the system does not replace the original heads in the association by the unique ones, so I have to convert to normal form and then back to association.*)
		
		List@@(Clear/@Orighs);
		Orighs = List@@Orighs;
		Uniquehs = List@@Uniquehs;
		(*iNRules = NRules/.Thread@Rule[Orighs, Uniquehs];*)
		,"parsing NRules"
	];
	
	
	EchoTiming[MakeDefs[SymRules, iNRules]; Clear[iNRules];, "MakeDefs"];
	
	
	(*Define some basic quantities:*)
	EchoTiming[
		idetecHs = ConstantArray[{h}[[1]], detecs];
		fvec = Range[f0,f1, \[CapitalDelta]f];
		ObsPoints = Table[Join[{fp}[[i]], {fvec}], {i, Length@{fp}}];,  
		"Middle work"
	];  (*Join[{#}, {fvec}]&/@{fp};*)
	

	EchoTiming[Unprotect[Derivative];, "pointless"]; (*Overloading of Derivative happening in iGenGrads*)
	(*Calculate the detector and remaining gradients:*)
	
	EchoTiming[detectorGradients = GenGrads[idetecHs, dims[[1]]&/@Range[detecs], ObsPoints[[1;;detecs]], Orighs, Uniquehs, n, {Exp[I anything_]-> 1}];, "detectors"];

	EchoTiming[remainingGradients = GenGrads[{h}[[2;;-1]], dims[[2;;-1]], ObsPoints[[detecs+1;;-1]], Orighs, Uniquehs, n, {Exp[I anything_]-> 1}];, "Core WF", Method->Timing];
	

	(*Clean definitions:*)
	
	EchoTiming[Remove[Evaluate[Uniquehs]]; GCSymRules[SymRulehs]; Protect[Derivative];, "Cleaning"];
	
	(*redefine remaining gradients indices {1,2} -> {3,4} and so on*)
	
	EchoTiming[remainingGradients = MapThread[
		KeyMap[Function[{x}, If[x[[1]] >0, x + #1, x]], #2]&, 
		{FoldList[Plus, 0, dims][[2;;-2]], remainingGradients}
	];, "Redefine labels"];
	
	(*Calculate full gradients for each detector, you have  a list of associations here*)
	EchoTiming[detectorGradients = GradientsList[Join[{#}, remainingGradients], dims, n]&/@detectorGradients;, "Gradient Recombination"];
	EchoTiming[Clear[remainingGradients];, "extra cleaning"];
	
	(*Put the matrices in  the form for DALIList*)
	
	(*This will give you a list of DALIlists:*)
	EchoTiming[result = MapThread[
		GenDaliList[#1, n, #2, \[CapitalDelta]f, "GW"]&,
		{detectorGradients, {PSD}}
	];, "DALIList"];
	
	(*Combine all of them and return:*)
	EchoTiming[Plus@@result, "Combine Detector DALIS"]
]


(* ::Section::Closed:: *)
(*Package Footer*)


End[];
EndPackage[];
