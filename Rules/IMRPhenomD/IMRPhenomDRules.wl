(* ::Package:: *)

Quit


(*Basic stuff:*)
$HistoryLength = 1;
PacletDirectoryLoad[ParentDirectory[NotebookDirectory[], 3]];
<<FelipeBarbosa`SymDALI`

vectorDefs = Module[
	{DALIDir, docDir},
	DALIDir = FindFile["FelipeBarbosa`SymDALI`"]//FileNameDrop//ParentDirectory;
	docDir = FileNameJoin[{DALIDir, "Documentation/English/Tutorials"}];
	Import@FileNameJoin[{docDir, "IMRPhenomDComponentsDefinition.wdx"}]
];


(* ::Section::Closed:: *)
(*Utils*)


PhenomCoeff[\[Eta]_, \[Chi]Pn_, \[Lambda]_?VectorQ] := With[
    {vec = {1, \[Eta], \[Eta]^2}, diff = \[Chi]Pn - 1},
    
    \[Lambda][[1;;2]] . vec[[1;;2]] + diff (\[Lambda][[3;;5]] . vec) + diff^2 (\[Lambda][[6;;8]] . vec) + diff^3 (\[Lambda][[9;;11]] . vec)
];


PhenomDTableV = List[ (*TABLE V in https://arxiv.org/pdf/1508.07253*)
{3931.89794921875`,-17395.7578125`,3132.37548828125`,343965.875`,-1.216256625`*^6,-70698.0078125`,1.383907125`*^6,-3.96627625`*^6,-60017.5234375`,803515.125`,-2.091710375`*^6},
{-40105.4765625`,112253.015625`,23561.6953125`,-3.47618075`*^6,1.1375937`*^7,754313.125`,-1.308476`*^7,3.6444584`*^7,596226.625`,-7.42779`*^6,1.8928978`*^7},
{83208.3515625`,-191237.71875`,-210916.25`,8.717975`*^6,-2.6914942`*^7,-1.988980625`*^6,3.088803`*^7,-8.3908704`*^7,-1.45350325`*^6,1.7063528`*^7,-4.274866`*^7},
{0.8149838447570801`,2.5747554302215576`,1.1610198020935059`,-2.3627772331237793`,6.77103853225708`,0.7570782899856567`,-2.725689649581909`,7.114037990570068`,0.1766934096813202`,-0.797869086265564`,2.116239070892334`},
{0.006927402690052986`,0.030204743146896362`,0.006308024283498526`,-0.12074130773544312`,0.26271599531173706`,0.0034151773434132338`,-0.10779338330030441`,0.27098965644836426`,0.0007374185952357948`,-0.027496211230754852`,0.0733150765299797`},
{1.010344386100769`,0.0008993122028186917`,0.2839491069316864`,-4.049753189086914`,13.207828521728516`,0.10396278649568558`,-7.025059223175049`,24.784893035888672`,0.030932024121284485`,-2.6924023628234863`,9.609374046325684`},
{1.3081616163253784`,-0.0055377297103405`,-0.0678291767835617`,-0.668983519077301`,3.4031479358673096`,-0.05296577513217926`,-0.9923793077468872`,4.820681095123291`,-0.0061341398395597935`,-0.3842925429344177`,1.7561753988265991`},
{2096.552001953125`,1463.749267578125`,1312.54931640625`,18307.330078125`,-43534.14453125`,-833.2889404296875`,32047.3203125`,-108609.453125`,452.2513732910156`,8353.439453125`,-44531.32421875`},
{-10114.056640625`,-44631.01171875`,-6541.30859375`,-266959.21875`,686328.3125`,3405.63720703125`,-437507.71875`,1.631817125`*^6,-7462.6484375`,-114585.25`,674402.5`},
{22933.658203125`,230960.015625`,14961.083984375`,1.194018125`*^6,-3.104224`*^6,-3038.16650390625`,1.87203225`*^6,-7.309145`*^6,42738.23046875`,467502.03125`,-3.0648535`*^6},
{-14621.71484375`,-377812.84375`,-9608.6826171875`,-1.7108925`*^6,4.3329245`*^6,-22366.68359375`,-2.50197175`*^6,1.0274496`*^7,-85360.3046875`,-570025.375`,4.3968445`*^6},
{97.89747619628906`,-42.65972900390625`,153.4842071533203`,-1417.0621337890625`,2752.861328125`,138.7406463623047`,-1433.658447265625`,2857.741943359375`,41.025108337402344`,-423.68072509765625`,850.3594360351562`},
{-3.2827019691467285`,-9.051384925842285`,-12.415450096130371`,55.47164535522461`,-106.05110168457031`,-11.953044891357422`,76.80704498291016`,-155.33172607421875`,-3.412926197052002`,25.572378158569336`,-54.40803527832031`},
{-0.00002515643063816242`,0.000019750257706618868`,-0.000018370670659351163`,0.000021886316972086206`,0.00008250240352936089`,7.157371328503359`*^-6,-0.000055780001275707036`,0.00019142082601319999`,5.447166131489212`*^-6,-0.000032206102332565933`,0.00007974016625666991`},
{43.315147399902344`,638.6332397460938`,-32.857688903808594`,2415.893798828125`,-5766.875`,-61.854591369628906`,2953.9677734375`,-8986.291015625`,-21.571435928344727`,981.2158203125`,-3239.56640625`},
{-0.07020209729671478`,-0.16269798576831818`,-0.18725146353244781`,1.1383136510849`,-2.8334195613861084`,-0.17137955129146576`,1.719754934310913`,-4.539717197418213`,-0.0499834381043911`,0.6062071919441223`,-1.6827696561813354`},
{9.598807334899902`,-397.05438232421875`,16.202125549316406`,-1574.8287353515625`,3600.341064453125`,27.092430114746094`,-1786.4822998046875`,5152.91943359375`,11.17570972442627`,-577.7999267578125`,1808.730712890625`},
{-0.0298948734998703`,1.4022105932235718`,-0.07356049120426178`,0.8337006568908691`,0.22400082647800446`,-0.055202871561050415`,0.5667186379432678`,0.718693196773529`,-0.015507437288761139`,0.15750323235988617`,0.21076816320419312`},
{0.9974408149719238`,-0.007884449325501919`,-0.059046901762485504`,1.3958712816238403`,-4.516631603240967`,-0.055853430181741714`,1.7516579627990723`,-5.990209102630615`,-0.017945336177945137`,0.5965097546577454`,-2.0608880519866943`}
];


Clear@RelativeDiff

Attributes[RelativeDiff] = {Listable};
RelativeDiff[x_,y_]/;x==0 &&y==0 := 0
RelativeDiff[0, y_]/; y!=0 := 1
RelativeDiff[x_, 0]/; x!=0 := 1

RelativeDiff[x_, y_]/; x!=0 &&y!=0 := With[
	{diff = x-y},
	Max[Abs@{diff/x, diff/y}]
]


(* ::Chapter:: *)
(*Phase*)


(* ::Section::Closed:: *)
(*PPN-parameters:*)


(* ::Text:: *)
(*Josiel \[And] Miguel argue that (\[Delta], \[Chi]s, \[Chi]a) is better than (\[Eta], \[Chi]1, \[Chi]2)*)


(* ::Text:: *)
(*\[Delta] -> Sqrt[1 - 4 \[Eta]] *)
(*\[Eta] -> (1 - \[Delta]^2)/4*)


Clear[\[Phi]0, \[Phi]1, \[Phi]2, \[Phi]3, \[Phi]4, \[Phi]5, \[Phi]6, \[Phi]7]


(* Define the coefficients *)
\[Phi]0 = 1;
\[Phi]1 = 0;
\[Phi]2[\[Delta]_] = (3715/756 + 55 \[Eta]/9)//.\[Eta]->(1-\[Delta]^2)/4//Simplify;

\[Phi]3[\[Delta]_, \[Chi]s_, \[Chi]a_] = (-16 \[Pi] + (113 \[Delta] \[Chi]a)/3 + (113/3 - 76 \[Eta]/3) \[Chi]s)//.{\[Eta]->(1-\[Delta]^2)/4}//Simplify;

\[Phi]4[\[Delta]_, \[Chi]s_, \[Chi]a_] = (15293365/508032 + 27145 \[Eta]/504 + 3085 \[Eta]^2/72 + (-405/8 + 200 \[Eta]) \[Chi]a^2 - (405/4) \[Delta] \[Chi]a \[Chi]s + (-405/8 + 5 \[Eta]/2) \[Chi]s^2)//.{
	\[Eta]-> (1-\[Delta]^2)/4
}//Simplify;

\[Phi]5[\[Delta]_, \[Chi]s_, \[Chi]a_] = (1 + Log[\[Pi] M \[Omega]]) * (38645 \[Pi]/756 - 65 \[Pi] \[Eta]/9 + 
    \[Delta] (-732985/2268 - 140 \[Eta]/9) \[Chi]a + 
    (-732985/2268 + 24260 \[Eta]/81 + 340 \[Eta]^2/9) \[Chi]s)//.{
	\[Eta]-> (1-\[Delta]^2)/4, Log[x_] ->0
}//Simplify;

\[Phi]6[\[Delta]_, \[Chi]s_, \[Chi]a_] = (11583231236531/4694215680 - (6848 EulerGamma)/21 - 
   (640 \[Pi]^2)/3 + (-15737765635/3048192 + (2255 \[Pi]^2)/12) \[Eta] + 
   76055 \[Eta]^2/1728 - 127825 \[Eta]^3/1296 - 
   (6848/63) Log[64 \[Pi] M \[Omega]] + (2270/3) \[Pi] \[Delta] \[Chi]a + 
   ((2270 \[Pi])/3 - 520 \[Pi] \[Eta]) \[Chi]s)//.{
	\[Eta]-> (1-\[Delta]^2)/4, Log[x_]-> Log[64]
}//Simplify;

\[Phi]7[\[Delta]_, \[Chi]s_, \[Chi]a_] = (77096675 \[Pi]/254016 + (378515 \[Pi] \[Eta])/1512 - (74045 \[Pi] \[Eta]^2)/756 + 
   \[Delta] (-25150083775/3048192 + (26804935 \[Eta])/6048 - (1985 \[Eta]^2)/48) \[Chi]a + 
   (-25150083775/3048192 + (10566655595 \[Eta])/762048 - 
      (1042165 \[Eta]^2)/3024 + (5345 \[Eta]^3)/36) \[Chi]s)//.{
      \[Eta]-> (1-\[Delta]^2)/4
}//Simplify;


(* ::Section::Closed:: *)
(*Ins-Phase*)


(* ::Text:: *)
(*In  the  notation  of  the  arXiv : 1903.04467  we  want  a  vector*)
(*{\[CurlyPhi]minus2, \[CurlyPhi]0, \[CurlyPhi]1, \[CurlyPhi]2, \[CurlyPhi]3, \[CurlyPhi]4, \[CurlyPhi]5, \[CurlyPhi]5l, \[CurlyPhi]6, \[CurlyPhi]6l, \[CurlyPhi]7} . In  GR  there  is  no  \[CurlyPhi]minus2  or  \[CurlyPhi]1, we  shall  set  their  values  to  1*3/(128 \[Eta])*)
(*and  multiply  by  a  vector  1 + \[Delta]\[CurlyPhi]  where  the  default  value  of  \[Delta]\[CurlyPhi]minus2  and  \[Delta]\[CurlyPhi]1  is - 1 (to recover GR) and -1+\[Delta]\[CurlyPhi] when we want to sample on them:*)


Clear[insVecPhase, \[Omega]InsVecPhase]


Block[
	{},
	insVecPhase[\[Delta]_, \[Chi]s_, \[Chi]a_] = 3*{
		1, (*\[CurlyPhi]minus2*)
		1, (*\[CurlyPhi]0*)
		1, (*\[CurlyPhi]1*)
		\[Phi]2[\[Delta]], (*\[CurlyPhi]2*)
		\[Phi]3[\[Delta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]3*)
		\[Phi]4[\[Delta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]4*)
		\[Phi]5[\[Delta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]5*)
		\[Phi]5[\[Delta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]5l*)
		\[Phi]6[\[Delta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]6*)
		-(6848/63), (*\[CurlyPhi]6l*)
		\[Phi]7[\[Delta], \[Chi]s, \[Chi]a] (*\[CurlyPhi]7*)
	}/(128 \[Eta])//. \[Eta]-> (1-\[Delta]^2)/4//Simplify;
	
	\[Omega]InsVecPhase[\[Omega]_] = Join[
		{(\[Pi] \[Omega])^(-7/3)}, (*-2*)
		(\[Pi] \[Omega])^((#-5)/3)&/@Range[0,5], (*0 to 5*)
		{Log[\[Pi] \[Omega]], (\[Pi] \[Omega])^(1/3), Log[\[Pi] \[Omega]] (\[Pi] \[Omega])^(1/3), (\[Pi] \[Omega])^(2/3)} (*5l, 6, 6l, 7*)
	]
];


(* ::Text:: *)
(*Now, introduce the \[Delta]s:*)
(**)
(*Note that 2010.14529 points out that \[Delta] goes only on the non spinning part: *)
(**)
(*\[CurlyPhi] = \[CurlyPhi]NS + \[CurlyPhi]S -> (1+\[Delta]) \[CurlyPhi]NS + \[CurlyPhi]S*)
(**)
(*note that \[CurlyPhi]NS == \[Phi][\[Delta], 0,0] \[And] \[CurlyPhi]S == \[Phi][\[Delta], \[Chi]s, \[Chi]a] - \[Phi][\[Delta],0,0].*)


Clear@InsVecPhase


Block[
	{v = 1 + {-1 + \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, -1 + \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, 0, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}},
	
	
	Clear@InsVecPhase;
	
	
	InsVecPhase[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_] = (
		v*(insVecPhase[\[Delta], 0, 0]) + (insVecPhase[\[Delta], \[Chi]s, \[Chi]a] - insVecPhase[\[Delta], 0, 0])
	)//Simplify
];


Clear["\[CapitalPhi]*"]


{
	\[CapitalPhi]minus2[\[Delta]_, \[Delta]\[CurlyPhi]minus2_], \[CapitalPhi]0[\[Delta]_, \[Delta]\[CurlyPhi]0_], \[CapitalPhi]1[\[Delta]_, \[Delta]\[CurlyPhi]1_], \[CapitalPhi]2[\[Delta]_, \[Delta]\[CurlyPhi]2_], \[CapitalPhi]3[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]3_], 
	\[CapitalPhi]4[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]4_], \[CapitalPhi]5[\[Delta]_, \[Chi]s_, \[Chi]a_], \[CapitalPhi]5l[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]5l_], \[CapitalPhi]6[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]6_], 
	\[CapitalPhi]6l[\[Delta]_, \[Delta]\[CurlyPhi]6l_], \[CapitalPhi]7[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]7_]
} = InsVecPhase[\[Delta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];


Clear["\[CapitalSigma]*"]


{\[CapitalSigma]1[\[Delta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]2[\[Delta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]3[\[Delta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]4[\[Delta]_, \[Chi]s_, \[Chi]a_]} = \[Eta]^-1*(PhenomCoeff[\[Eta], \[Chi]PN, #]&/@PhenomDTableV[[8;;11]])//.{
	\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]-> (1-\[Delta]^2)/4}//Simplify;


{\[CapitalSigma]1[\[Delta], \[Chi]s, \[Chi]a], \[CapitalSigma]2[\[Delta], \[Chi]s, \[Chi]a],\[CapitalSigma]3[\[Delta], \[Chi]s, \[Chi]a],\[CapitalSigma]4[\[Delta], \[Chi]s, \[Chi]a]} . {\[Omega], 3/4 \[Omega]^(4/3), 3/5 \[Omega]^(5/3), 1/2 \[Omega]^2};


InsExpr = -\[Pi]/4 + {
	U\[CapitalPhi]minus2[\[Delta],\[Delta]\[CurlyPhi]minus2],U\[CapitalPhi]0[\[Delta],\[Delta]\[CurlyPhi]0],U\[CapitalPhi]1[\[Delta],\[Delta]\[CurlyPhi]1],U\[CapitalPhi]2[\[Delta],\[Delta]\[CurlyPhi]2],U\[CapitalPhi]3[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]3],U\[CapitalPhi]4[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]4],U\[CapitalPhi]5[\[Delta], \[Chi]s, \[Chi]a],
	U\[CapitalPhi]5l[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]5l],U\[CapitalPhi]6[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]6],U\[CapitalPhi]6l[\[Delta],\[Delta]\[CurlyPhi]6l],U\[CapitalPhi]7[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]7]
} . \[Omega]InsVecPhase[\[Omega]]  + {U\[CapitalSigma]1[\[Delta], \[Chi]s, \[Chi]a],U\[CapitalSigma]2[\[Delta], \[Chi]s, \[Chi]a],U\[CapitalSigma]3[\[Delta], \[Chi]s, \[Chi]a],U\[CapitalSigma]4[\[Delta], \[Chi]s, \[Chi]a]} . {\[Omega], 3/4 \[Omega]^(4/3), 3/5 \[Omega]^(5/3), 1/2 \[Omega]^2};


AuxiliarDefs = <||>;
Clear@ZeroFunction
ZeroFunction[x__] := 0


AuxiliarDefs["Inspiral"] = HoldForm[{
	Unprotect[Derivative],
	
	Derivative[\[Delta]_, chis_, chia_, d_][\[CapitalPhi]3]/;chis+chia>1 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_, d_][\[CapitalPhi]4]/;chis+chia>2 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_][\[CapitalPhi]5]/;chis+chia>1 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_, d_][\[CapitalPhi]6]/;chis+chia>1 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_, d_][\[CapitalPhi]7]/;chis+chia>1 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_, d_][\[CapitalPhi]5l]/;chis+chia>1 := ZeroFunction,
	
	Derivative[\[Delta]_, chis_, chia_][\[CapitalSigma]1]/; chis+chia>3 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_][\[CapitalSigma]2]/; chis+chia>3 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_][\[CapitalSigma]3]/; chis+chia>3 := ZeroFunction,
	Derivative[\[Delta]_, chis_, chia_][\[CapitalSigma]4]/; chis+chia>3 := ZeroFunction
}];


(* ::Section::Closed:: *)
(*Int-Phase*)


(*Find vectors for Intermediate Phase and their variables:*)
Position[vectorDefs[[All, 1, All, 0]] /. HoldPattern -> Identity, #] & /@ {IntVecPhase, \[Omega]IntVecPhase}

vectorDefs[[3 ;; 4, 1]]

DownValues[IntVecPhase] = {vectorDefs[[3]]};
DownValues[\[Omega]IntVecPhase] = {vectorDefs[[4]]};


(* ::Text:: *)
(*Now  Introduce \[Delta]\[Beta]:*)


Clear["\[Beta]*"]


{\[Beta]1[\[Delta]_, \[Chi]s_, \[Chi]a_], \[Beta]2[\[Delta]_, \[Chi]s_, \[Chi]a_], \[Beta]3[\[Delta]_, \[Chi]s_, \[Chi]a_]} = \[Eta] IntVecPhase[\[Eta], \[Chi]PN]//.{
	\[Chi]PN->\[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]-> (1-\[Delta]^2)/4};


IntVecPhase[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[Beta]2_, \[Delta]\[Beta]3_] = (1 + {0, \[Delta]\[Beta]2, \[Delta]\[Beta]3})*{\[Beta]1[\[Delta], \[Chi]s, \[Chi]a], \[Beta]2[\[Delta], \[Chi]s, \[Chi]a], \[Beta]3[\[Delta], \[Chi]s, \[Chi]a]};


{\[CapitalBeta]1[\[Delta]_, \[Chi]s_, \[Chi]a_], \[CapitalBeta]2[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[Beta]2_], \[CapitalBeta]3[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[Beta]3_]} = (IntVecPhase[\[Delta],\[Chi]s,\[Chi]a,\[Delta]\[Beta]2,\[Delta]\[Beta]3]/\[Eta])//.{
\[Eta]-> (1-\[Delta]^2)/4

}//Simplify;


IntExpr = {U\[CapitalBeta]1[\[Delta], \[Chi]s, \[Chi]a], U\[CapitalBeta]2[\[Delta], \[Chi]s, \[Chi]a, \[Delta]\[Beta]2],U\[CapitalBeta]3[\[Delta], \[Chi]s, \[Chi]a, \[Delta]\[Beta]3]} . \[Omega]IntVecPhase[\[Omega]];


SetOptions[EvaluationNotebook[], Magnification->1.75]


AuxiliarDefs["Intermediate"] = HoldForm[{
	Derivative[\[Delta]_, chis_, chia_][\[CapitalBeta]1]/; chis+chia>3 := ZeroFunction,
	Derivative[\[Delta]_,chis_, chia_, d_][\[CapitalBeta]2]/;chis+chia>3 := ZeroFunction,
	Derivative[\[Delta]_,chis_, chia_, d_][\[CapitalBeta]3]/;chis+chia>3 := ZeroFunction
}];


(* ::Section::Closed:: *)
(*Ringdown and Damping -Phase*)


aeff[\[Delta]_, \[Chi]s_, \[Chi]a_] = Block[
	{S},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S + 2 Sqrt[3.] \[Eta] + (-0.085` S+0.102` S^2-1.355` S^3-0.868` S^4) \[Eta]-4.399` \[Eta]^2+(-5.837` S-2.097` S^2+4.109` S^3+2.064` S^4) \[Eta]^2+9.397` \[Eta]^3-13.181` \[Eta]^4//.\[Eta]-> (1-\[Delta]^2)/4
]//Simplify;

Erad[\[Delta]_, \[Chi]s_, \[Chi]a_] = Block[
	{EradNS=0.0559745 \[Eta]+0.580951 \[Eta]^2-0.960673 \[Eta]^3+3.35241 \[Eta]^4, S},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S = S/(1-2 \[Eta]);
	
	(EradNS (1+(-0.00303023` -2.00661` \[Eta]+7.70506` \[Eta]^2) S))/(1+ (-0.67144` -1.47569` \[Eta] + 7.30468` \[Eta]^2) S)//.\[Eta]-> (1-\[Delta]^2)/4
]//Simplify;

(*re\[Omega] -> interpolation for ringdown and the other is for damping: I think this is frm the PhenomHM paper*)
re\[Omega][\[Chi]_] :=(0.05947169566573468` -0.14989771215394762` \[Chi]+0.09535606290986028` \[Chi]^2+0.02260924869042963` \[Chi]^3-0.02501704155363241` \[Chi]^4-0.005852438240997211` \[Chi]^5+0.0027489038393367993` \[Chi]^6+0.0005821983163192694` \[Chi]^7)/(1-2.8570126619966296` \[Chi]+2.373335413978394` \[Chi]^2-0.6036964688511505` \[Chi]^4+0.0873798215084077` \[Chi]^6);
im\[Omega][\[Chi]_] :=(0.014158792290965177` -0.036989395871554566` \[Chi]+0.026822526296575368` \[Chi]^2+0.0008490933750566702` \[Chi]^3-0.004843996907020524` \[Chi]^4-0.00014745235759327472` \[Chi]^5+0.0001504546201236794` \[Chi]^6)/(1-2.5900842798681376` \[Chi]+1.8952576220623967` \[Chi]^2-0.31416610693042507` \[Chi]^4+0.009002719412204133` \[Chi]^6);

\[Omega]RdDamping[int_, Erad_] = int/(1 - Erad);


Clear[\[Gamma]2, \[Gamma]3]


\[Gamma]2[\[Delta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[6]]//.{\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]-> (1-\[Delta]^2)/4};
\[Gamma]3[\[Delta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[7]]//.{\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]-> (1-\[Delta]^2)/4};


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


AuxiliarDefs["rd \[And] dp"] = HoldForm[{
	Derivative[\[Delta]_, chis_, chia_][\[Gamma]2]/; chis+chia>3 := ZeroFunction,
	Derivative[\[Delta]_,chis_, chia_][\[Gamma]3]/; chis+chia>3 := ZeroFunction,
	Derivative[omRD_, omDamp_, gamma2_, gamma3_][\[Omega]Peak]/;omRD + omDamp>1 := ZeroFunction,
	Derivative[omRD_, omDamp_, gamma2_, gamma3_][\[Omega]Peak]/;gamma3 + omDamp > 2 := ZeroFunction 
}];


(* ::Section::Closed:: *)
(*MR-Phase*)


Position[vectorDefs[[All, 1, All, 0]]/.HoldPattern->Identity, #]&/@{MRVecPhase, \[Omega]MRVecPhase}


vectorDefs[[5;;6,1]]
DownValues[MRVecPhase] = {vectorDefs[[5]]};
DownValues[\[Omega]MRVecPhase] = {vectorDefs[[6]]};


{\[Alpha]1[\[Delta]_, \[Chi]s_, \[Chi]a_], \[Alpha]2[\[Delta]_, \[Chi]s_, \[Chi]a_], \[Alpha]3[\[Delta]_, \[Chi]s_, \[Chi]a_], \[Alpha]4[\[Delta]_, \[Chi]s_, \[Chi]a_]} = \[Eta] MRVecPhase[\[Eta], \[Chi]PN]//.{\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s , \[Eta]-> (1-\[Delta]^2)/4};


\[Alpha]5[\[Delta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[19]]//.{\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]-> (1-\[Delta]^2)/4};


(* ::Text:: *)
(*Introduce  \[Delta]\[Alpha] :*)


MRVecPhase[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[Alpha]2_, \[Delta]\[Alpha]3_, \[Delta]\[Alpha]4_] = (1 + {0, \[Delta]\[Alpha]2, \[Delta]\[Alpha]3, \[Delta]\[Alpha]4})*{\[Alpha]1[\[Delta], \[Chi]s, \[Chi]a], \[Alpha]2[\[Delta], \[Chi]s, \[Chi]a], \[Alpha]3[\[Delta], \[Chi]s, \[Chi]a], \[Alpha]4[\[Delta], \[Chi]s, \[Chi]a]};


{\[CapitalAlpha]1[\[Delta]_, \[Chi]s_, \[Chi]a_], \[CapitalAlpha]2[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[Alpha]2_], \[CapitalAlpha]3[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[Alpha]3_], \[CapitalAlpha]4[\[Delta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[Alpha]4_]} =  MRVecPhase[\[Delta], \[Chi]s, \[Chi]a, \[Delta]\[Alpha]2, \[Delta]\[Alpha]3, \[Delta]\[Alpha]4]/\[Eta]//.\[Eta]-> (1-\[Delta]^2)/4//Simplify;


MRExpr = Block[
	{ringdownFrequency, dampingFrequency, \[Alpha]5},
	\[Alpha]5 = U\[Alpha]5[\[Delta], \[Chi]s, \[Chi]a];
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]];
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]];
	 {U\[CapitalAlpha]1[\[Delta], \[Chi]s, \[Chi]a],U\[CapitalAlpha]2[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[Alpha]2],U\[CapitalAlpha]3[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[Alpha]3],U\[CapitalAlpha]4[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[Alpha]4]} . \[Omega]MRVecPhase[\[Omega], ringdownFrequency, dampingFrequency, \[Alpha]5]
]//Simplify;


AuxiliarDefs["MR"] = HoldForm[{
	Derivative[delta_, chis_, chia_][\[Alpha]5]/; chis+chia>3 := ZeroFunction,
	Derivative[delta_, chis_, chia_][aeff]/; chis+chia>4 := ZeroFunction,
	Derivative[n_, m_][\[Omega]RdDamping]/; n>1 := ZeroFunction,
	
	Derivative[delta_, chis_, chia_][\[CapitalAlpha]1]/; chis+chia>3 := ZeroFunction,
	Derivative[delta_, chis_, chia_, d_][\[CapitalAlpha]2]/; chis+chia>3 := ZeroFunction,
	Derivative[delta_, chis_, chia_, d_][\[CapitalAlpha]3]/; chis+chia>3 := ZeroFunction,
	Derivative[delta_, chis_, chia_, d_][\[CapitalAlpha]4]/; chis+chia>3 := ZeroFunction,
	Protect[Derivative]
}];


(* ::Section::Closed:: *)
(*C(1)*)


C1[beta0_, beta1\[Omega]_, alpha0_, alpha1\[Omega]_, \[Omega]_, ringdownFrequency_] = (
			(beta0 + beta1\[Omega])*us\[Theta][(ringdownFrequency/2 - \[Omega]) (\[Omega]-0.018)] + 
			(alpha0 + alpha1\[Omega])us\[Theta][(\[Omega] - ringdownFrequency/2) (0.2-\[Omega])]
);


Module[{
	ringdownFrequency, dampingFrequency, \[Alpha]5
	},
	\[Alpha]5 = U\[Alpha]5[\[Delta], \[Chi]s, \[Chi]a];
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]];
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]];

	\[Delta]\[Beta]1[d\[Phi]Ins_, d\[Phi]Int_] = d\[Phi]Ins - d\[Phi]Int;
	\[Delta]\[Beta]0[\[CapitalPhi]Ins_, \[CapitalPhi]Int_, \[Delta]\[Beta]1_] = \[CapitalPhi]Ins - \[CapitalPhi]Int - \[Delta]\[Beta]1 0.018;

	\[Delta]\[Alpha]1[d\[CapitalPhi]Int_, d\[CapitalPhi]MR_, \[Delta]\[Beta]1_] = d\[CapitalPhi]Int + \[Delta]\[Beta]1 - d\[CapitalPhi]MR;
	
	\[Delta]\[Alpha]0[\[CapitalPhi]Int_, \[CapitalPhi]MR_, \[Delta]\[Beta]0_, \[Delta]\[Beta]1\[Omega]rd_, \[Delta]\[Alpha]1\[Omega]rd_] := \[CapitalPhi]Int + \[Delta]\[Beta]1\[Omega]rd/2 + \[Delta]\[Beta]0 - \[CapitalPhi]MR - \[Delta]\[Alpha]1\[Omega]rd/2;
	
	Block[
		{dphiIns, dphiInt1,dphiInt2,dphiMR, phiIns, phiInt1, phiInt2,phiMR, beta0, beta1, alpha0, alpha1},
		
		(*derivative functions ###################################################### *)
		dphiIns = D[InsExpr, \[Omega]]//.\[Omega]->0.018;

		dphiInt1 = D[IntExpr, \[Omega]]//.\[Omega]->0.018;
		
		dphiInt2 = D[IntExpr, \[Omega]]//.\[Omega]->ringdownFrequency/2;
		
		dphiMR = D[MRExpr, \[Omega]]//.\[Omega]->ringdownFrequency/2;
		
		(*normal functions ########################################################################*)
		phiIns = InsExpr//.\[Omega]->0.018;

		phiInt1 = IntExpr//.\[Omega]->0.018;
		
		phiInt2 =IntExpr//.\[Omega]-> ringdownFrequency/2;
		
		phiMR = MRExpr//.\[Omega]-> ringdownFrequency/2;
	
		
		beta1 = U\[Delta]\[Beta]1[dphiIns, dphiInt1];
		beta0 = U\[Delta]\[Beta]0[phiIns, phiInt1, beta1];
		
		alpha1 = U\[Delta]\[Alpha]1[dphiInt2, dphiMR, beta1];
		alpha0 = U\[Delta]\[Alpha]0[phiInt2, phiMR, beta0, beta1*ringdownFrequency, alpha1*ringdownFrequency];
		
		
		
		exprC1 = UC1[beta0, beta1*\[Omega], alpha0, alpha1*\[Omega], \[Omega], ringdownFrequency];
	]
];


(* ::Section::Closed:: *)
(*\[Phi]IMR*)


\[Phi]IMR[\[Phi]Ins_, \[Phi]Int_, \[Phi]MR_, c1_, ringdownFrequency_, \[Omega]_] := (
	\[Phi]Ins*us\[Theta][-\[Omega]+0.018] +
	\[Phi]Int*us\[Theta][(ringdownFrequency/2 - \[Omega]) (\[Omega]-0.018)]  +
	\[Phi]MR*us\[Theta][(\[Omega] - ringdownFrequency/2) (0.2-\[Omega])] + 
	c1
)


Module[
	{ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]],
	\[Phi]Ins, \[Phi]Int, \[Phi]MR},
	
	\[Phi]Ins = InsExpr;
	
	\[Phi]Int = IntExpr;
	
	\[Phi]MR = MRExpr;
	
	\[Phi]IMRexpr = U\[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, exprC1, ringdownFrequency, \[Omega]];
	\[Phi]IMRexprRef = U\[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, exprC1, ringdownFrequency, \[Omega]]//.\[Omega]->\[Omega]ref;
]//Simplify;


(* ::Section::Closed:: *)
(*t0*)


\[CapitalDelta]t0[\[Omega]_, d\[Phi]MR_, \[Omega]ref_] = d\[Phi]MR (\[Omega] - \[Omega]ref);

Module[
	{d\[Phi]MR, ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]],
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]], \[Alpha]5, peakFrequency, \[Gamma]2, \[Gamma]3},
	\[Gamma]2 = U\[Gamma]2[\[Delta], \[Chi]s, \[Chi]a];
	\[Gamma]3 = U\[Gamma]3[\[Delta], \[Chi]s, \[Chi]a];
	
	peakFrequency = U\[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2, \[Gamma]3];
	
	\[Alpha]5 = U\[Alpha]5[\[Delta], \[Chi]s, \[Chi]a];
	
	
	d\[Phi]MR =D[MRExpr, \[Omega]]//.\[Omega]->peakFrequency;
		
	\[CapitalDelta]t0expr = U\[CapitalDelta]t0[\[Omega], d\[Phi]MR, \[Omega]ref]
];


(* ::Text:: *)
(*Pos[expri, hj] ={POS1, POS2, ...},  j = 1, ..., N*)


(* ::Section::Closed:: *)
(*Making Block function*)


(* ::Text:: *)
(*Now, all function definitions, starting at Primitives and ending on \[Phi]IMR:*)


FunctionDefs = <||>;


FunctionDefs["Ins+Int"]=With[{
	rules =Thread@Rule[
		{-2,0,1,2,3,4,5,5l,6,6l,7},
		{
			\[CapitalPhi]minus2[\[Delta], \[Delta]\[CurlyPhi]minus2],\[CapitalPhi]0[\[Delta],\[Delta]\[CurlyPhi]0],\[CapitalPhi]1[\[Delta],\[Delta]\[CurlyPhi]1],\[CapitalPhi]2[\[Delta],\[Delta]\[CurlyPhi]2],\[CapitalPhi]3[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]3],\[CapitalPhi]4[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]4],
			\[CapitalPhi]5[\[Delta], \[Chi]s, \[Chi]a],\[CapitalPhi]5l[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]5l],\[CapitalPhi]6[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]6],\[CapitalPhi]6l[\[Delta],\[Delta]\[CurlyPhi]6l],\[CapitalPhi]7[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]7]
		}
	]},

 
   HoldForm[{
	\[CapitalPhi]minus2[\[Delta]_, \[Delta]\[CurlyPhi]minus2_]=-2,
	\[CapitalPhi]0[\[Delta]_,\[Delta]\[CurlyPhi]0_]=0,
	\[CapitalPhi]1[\[Delta]_,\[Delta]\[CurlyPhi]1_]=1,
	\[CapitalPhi]2[\[Delta]_,\[Delta]\[CurlyPhi]2_]=2,
	\[CapitalPhi]3[\[Delta]_,\[Chi]s_,\[Chi]a_,\[Delta]\[CurlyPhi]3_]=3,
	\[CapitalPhi]4[\[Delta]_,\[Chi]s_,\[Chi]a_,\[Delta]\[CurlyPhi]4_]=4,
	\[CapitalPhi]5[\[Delta]_,\[Chi]s_,\[Chi]a_]=5,
	\[CapitalPhi]5l[\[Delta]_,\[Chi]s_,\[Chi]a_,\[Delta]\[CurlyPhi]5l_]=5l,
	\[CapitalPhi]6[\[Delta]_,\[Chi]s_,\[Chi]a_,\[Delta]\[CurlyPhi]6_]=6,
	\[CapitalPhi]6l[\[Delta]_,\[Delta]\[CurlyPhi]6l_]=6l,
	\[CapitalPhi]7[\[Delta]_,\[Chi]s_,\[Chi]a_,\[Delta]\[CurlyPhi]7_]=7,
	\[CapitalSigma]1[\[Delta]_,\[Chi]s_,\[Chi]a_] = A, 
	\[CapitalSigma]2[\[Delta]_,\[Chi]s_,\[Chi]a_] = B, 
	\[CapitalSigma]3[\[Delta]_,\[Chi]s_,\[Chi]a_] = c, 
	\[CapitalSigma]4[\[Delta]_,\[Chi]s_,\[Chi]a_] = d,
	\[CapitalBeta]1[\[Delta]_,\[Chi]s_,\[Chi]a_] = x,
	\[CapitalBeta]2[\[Delta]_,\[Chi]s_,\[Chi]a_, \[Delta]\[Beta]2_] = y,
	\[CapitalBeta]3[\[Delta]_,\[Chi]s_,\[Chi]a_, \[Delta]\[Beta]3_] =z
	
}]/.Join[
	rules,
	{x-> \[CapitalBeta]1[\[Delta],\[Chi]s, \[Chi]a], y-> \[CapitalBeta]2[\[Delta],\[Chi]s, \[Chi]a, \[Delta]\[Beta]2],z-> \[CapitalBeta]3[\[Delta],\[Chi]s, \[Chi]a, \[Delta]\[Beta]3]},
	{A ->\[CapitalSigma]1[\[Delta],\[Chi]s, \[Chi]a], B -> \[CapitalSigma]2[\[Delta],\[Chi]s, \[Chi]a], c -> \[CapitalSigma]3[\[Delta],\[Chi]s, \[Chi]a], d -> \[CapitalSigma]4[\[Delta],\[Chi]s, \[Chi]a] }
]

];


FunctionDefs["MR"] = HoldForm[{
	aeff[\[Delta]_,\[Chi]s_, \[Chi]a_] := x,
	Erad[\[Delta]_,\[Chi]s_, \[Chi]a_] := y,
	
	re\[Omega][\[Chi]_] := z,
	im\[Omega][\[Chi]_] := w,
	\[Omega]RdDamping[int_, Erad_] := \[ScriptX],
	\[Gamma]2[\[Delta]_,\[Chi]s_, \[Chi]a_] := \[ScriptY],
	\[Gamma]3[\[Delta]_,\[Chi]s_, \[Chi]a_] := \[ScriptZ],
	
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = \[ScriptW],
	\[Alpha]5[\[Delta]_,\[Chi]s_, \[Chi]a_] = \[ScriptA],
	
	\[CapitalAlpha]1[\[Delta]_, \[Chi]s_, \[Chi]a_] = 1,
	\[CapitalAlpha]2[\[Delta]_,\[Chi]s_, \[Chi]a_,\[Delta]\[Alpha]2_] =2 ,
	\[CapitalAlpha]3[\[Delta]_,\[Chi]s_, \[Chi]a_,\[Delta]\[Alpha]3_] = 3,
	\[CapitalAlpha]4[\[Delta]_,\[Chi]s_, \[Chi]a_,\[Delta]\[Alpha]4_] = 4
	
}]/.{
	x -> aeff[\[Delta],\[Chi]s, \[Chi]a], y-> Erad[\[Delta],\[Chi]s, \[Chi]a], 
	z-> re\[Omega][\[Chi]], w-> im\[Omega][\[Chi]], \[ScriptX]-> \[Omega]RdDamping[int, Erad], \[ScriptY]-> \[Gamma]2[\[Delta],\[Chi]s, \[Chi]a], \[ScriptZ]-> \[Gamma]3[\[Delta],\[Chi]s, \[Chi]a], 
	\[ScriptW]-> \[Omega]Peak[\[Omega]RD, \[Omega]DAMP, \[Gamma]2, \[Gamma]3], \[ScriptA]-> \[Alpha]5[\[Delta],\[Chi]s, \[Chi]a], 
	1 -> \[CapitalAlpha]1[\[Delta], \[Chi]s, \[Chi]a], 2-> \[CapitalAlpha]2[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[Alpha]2], 3 ->\[CapitalAlpha]3[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[Alpha]3], 4 -> \[CapitalAlpha]4[\[Delta], \[Chi]s, \[Chi]a,\[Delta]\[Alpha]4]
};


FunctionDefs["C1"] = HoldForm[{
	
	C1[beta0_, beta1\[Omega]_, alpha0_, alpha1\[Omega]_, \[Omega]_,ringdownFrequency_]  := w,
	
	\[Delta]\[Beta]1[d\[Phi]Ins_, d\[Phi]Int_] := d\[Phi]Ins - d\[Phi]Int,
	
	\[Delta]\[Beta]0[\[CapitalPhi]Ins_, \[CapitalPhi]Int_, \[Delta]\[Beta]1_] := \[CapitalPhi]Ins - \[CapitalPhi]Int - \[Delta]\[Beta]1 0.018,

	\[Delta]\[Alpha]1[d\[CapitalPhi]Int_, d\[CapitalPhi]MR_, \[Delta]\[Beta]1_] := d\[CapitalPhi]Int + \[Delta]\[Beta]1 - d\[CapitalPhi]MR,
	
	\[Delta]\[Alpha]0[\[CapitalPhi]Int_, \[CapitalPhi]MR_, \[Delta]\[Beta]0_, \[Delta]\[Beta]1\[Omega]rd_, \[Delta]\[Alpha]1\[Omega]rd_] := \[CapitalPhi]Int + \[Delta]\[Beta]1\[Omega]rd/2 + \[Delta]\[Beta]0 - \[CapitalPhi]MR - \[Delta]\[Alpha]1\[Omega]rd/2
}]//.{
	w -> C1[beta0, beta1\[Omega], alpha0, alpha1\[Omega], \[Omega],ringdownFrequency]
};


FunctionDefs["\[Phi]IMR+t0"] = HoldForm[{
	\[Phi]IMR[\[Phi]Ins_, \[Phi]Int_, \[Phi]MR_, c1_, ringdownFrequency_, \[Omega]_] := x,
	\[CapitalDelta]t0[\[Omega]_, d\[Phi]MR_, \[Omega]ref_] := d\[Phi]MR (\[Omega] - \[Omega]ref)
}]//.{
	x-> \[Phi]IMR[\[Phi]Ins, \[Phi]Int, \[Phi]MR, c1, ringdownFrequency, \[Omega]]
};


FunctionDefs["Commute D"] = HoldForm[{
	(*distribute derivatives*)
	Unprotect[D],
	DownValues[D] = {},
	D[\[Phi]IMR[x__] - \[Phi]IMR[y__] - \[CapitalDelta]t0[z__], n___] := D[\[Phi]IMR[x],n] - D[\[Phi]IMR[y],n] - D[\[CapitalDelta]t0[z], n],
	Protect[D],

	(*Commute D with \[Phi]IMR*)
	\[Phi]IMR/: D[\[Phi]IMR[x__], n___] := With[
	{args = D[#, n]&/@({x}[[1;;4]])},
	\[Phi]IMR@@(Join[args, {x}[[5;;6]]])
	],
	
	(*Commute D with C1*)
	C1/: D[C1[x__], n___] := With[
		{args = D[#, n]&/@({x}[[1;;4]])},
		C1@@(Join[args, {x}[[5;;6]]])
	],
	(*Commute with \[Delta]\[Beta]0, \[Delta]\[Beta]1, \[Delta]\[Alpha]0*)
	\[Delta]\[Beta]0/: D[\[Delta]\[Beta]0[x__], n___] := \[Delta]\[Beta]0@@(D[#, n]&/@{x}),
	\[Delta]\[Beta]1/: D[\[Delta]\[Beta]1[x__], n___] := \[Delta]\[Beta]1@@(D[#,n]&/@{x}),
	\[Delta]\[Alpha]0/: D[\[Delta]\[Alpha]0[x__], n___] := \[Delta]\[Alpha]0@@(D[#,n]&/@{x})
}];


(* ::Text:: *)
(*Make rules to eliminate the U from expressions:*)


(*Make rules to eliminate U*)
AllFunctionHeads = (FunctionDefs//Values)[[All, 1, All, 1,0]]//Flatten//DeleteDuplicates;
AllFunctionHeads = DeleteElements[AllFunctionHeads, {D, DownValues, Symbol}];


Module[{undefined =  ("U" <> #)&/@(ToString/@AllFunctionHeads)},
	
	undefined = ToExpression[undefined,InputForm];
	Headrules = Thread@Rule[undefined, AllFunctionHeads]
];


(* ::Text:: *)
(*Make  the  final  expression :*)


FinalExpr = Module[
	{dummy},
	dummy = HoldForm[
		{phiIMR - phiIMRref - deltat}
	]/.{phiIMR -> \[Phi]IMRexpr, phiIMRref -> \[Phi]IMRexprRef, deltat -> \[CapitalDelta]t0expr};

	dummy//.Headrules
];


(* ::Text:: *)
(*The $Block  structure : *)
(*$Block[*)
(*	{{heads__}, localvars__},*)
(*	localvardefs__;*)
(*	functiondefs__;*)
(*	expr*)
(*] ,*)
(*following that:*)


defs = Module[
	{functionDefs = FunctionDefs//Values, orderedDefs},
	
	functionDefs = Join[functionDefs, Values[AuxiliarDefs]];
	
	functionDefs = HoldForm@@@#&/@functionDefs; (*{HoldForm[HoldForm[]], HoldForm[HoldForm[]],...}*)
	orderedDefs = {functionDefs, HoldForm@@@FinalExpr}//Flatten; (*{HoldForm[HoldForm[]], HoldForm[HoldForm[]],...}*)
	orderedDefs = Flatten[HoldForm@@orderedDefs]; (*HoldForm[allDefs]*)
	orderedDefs = CompoundExpression@@@(HoldForm[Evaluate@orderedDefs]) 
];


declarations = HoldForm[Evaluate@{AllFunctionHeads}];


$blockexpr = $Block@@@HoldForm[Join[declarations, defs]//Evaluate];


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {\[Omega], \[Delta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4};


(*All derivatives up to order 3*)

derivatives = Combinations[vars, 3];

derivatives//Length


(* ::Text:: *)
(*Bcs \[Delta]pi appears linearly, all terms with more than 1 derivative in \[Delta]pi are 0.*)
(*Also, there is no correlation between \[Delta]\[CurlyPhi]i and \[Chi]s, \[Chi]a bcs the \[Delta]\[Beta]1 term cancels in the C(1) contribution (\[Delta]\[Alpha]1 \[Proportional] \[Delta]\[Beta]1) and there is no correlation between \[Delta]\[CurlyPhi]i \[And] {\[Chi]s, \[Chi]a} originally in \[Phi]Ins*)


numberOf\[Delta]p[elem_List] := Module[
	{\[Delta]ps = {\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}},
	
	Count[
		elem, 
		x_/; Or@@(Equal[x, #]&/@\[Delta]ps)
	]
]

numberOf\[Delta]\[CurlyPhi]\[Chi]s\[Chi]a[elem_List] := Module[
	{l, dummy1,dummy2, \[Delta]ps = {\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}},
	
	l = Join[\[Delta]ps, {\[Chi]s, \[Chi]a}];
	
	(*number of \[Delta]ps and/or \[Chi]s, \[Chi]a*)
	dummy1 = Count[
		elem, 
		x_/; Or@@(Equal[x, #]&/@l)
	];
	
	(*number of \[Delta]ps in the element*)
	dummy2 = Count[
		elem, 
		x_/; Or@@(Equal[x, #]&/@\[Delta]ps)
	];
	
	(*I want the number to be > 0 only if there is at least 1 \[Delta]p in elem*)
	If[
		dummy2 == 0,
		0, 
		dummy1
	]
	
]

Filter\[Delta]ps[derivatives_List]:= Module[
	{iL, pos, iList},
	
	iL = numberOf\[Delta]p/@derivatives;
	
	pos = Position[iL, x_/; x<=1];
	
	iList = Extract[derivatives, pos]; (*exclude \[Delta]pi^2*)
	
	iL = numberOf\[Delta]\[CurlyPhi]\[Chi]s\[Chi]a/@iList;
	
	pos = Position[iL, x_/; x<=1];
	
	Extract[iList, pos] (*exclude \[Delta]\[CurlyPhi]i*spins*)
	
]

derivatives = Filter\[Delta]ps[derivatives];

derivatives//Length

PrependTo[derivatives, {}];

derivatives//Length


Clear@expr;

With[{
	vars = {\[Omega], \[Omega]ref, \[Delta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}
	},
	expr =  Hold[
	{\[CapitalPhi]IMR,  vars, derivatives},
	Evaluate@$blockexpr,
	"IncludeZeroDerivative"->False
]//.HoldForm[X_] :> X;
]


<<FelipeBarbosa`SymDALI`


res1 = EchoTiming[DerivativeRules@@expr];


Unprotect[Derivative];
SubValues[Derivative] = Drop[SubValues[Derivative], {2, -1}]
Protect[Derivative];


test[i_] := Block[
	{testF,\[Omega], \[Eta],\[Delta],\[Chi]s, \[Chi]a,
	\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4,
	\[Omega]ref, M, G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
	 f},
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = ConstantArray[0, 16];
	
	testF[\[Omega]_, \[Delta]_, \[Chi]s_, \[Chi]a_] = res1[[i,2]]; DownValues[testF] = DownValues[testF]//. HoldForm[x_]:> x;
	
	\[Eta] = RandomReal[{0.1, 0.24}];
	\[Delta] = Sqrt[1-4 \[Eta]];
	{\[Chi]s, \[Chi]a} = RandomReal[{-1,1}, 2];
	M = RandomReal[{20,100}];
	f = RandomReal[{10, 0.2/(M G)}];
	\[Omega] = f M G;
	\[Omega]ref = 20 M G;
	
	testF[\[Omega], \[Delta], \[Chi]s, \[Chi]a]
]//Quiet


SetDirectory[NotebookDirectory[]];


Export["Phase_Ds_order_0_to_3.mx", res1]


(*Compiler`$CCompilerOptions = {
    "ShellCommandFunction" -> Print, 
    "SystemCompileOptions" -> "-march=native -finline-functions -funroll-loops -flto -fPIC -O3 -fno-fast-math -ftree-vectorize -fvect-cost-model=dynamic"
};
*)

(*Compiler`$CCompilerOptions={
	"ShellCommandFunction"->Print, 
	"SystemCompileOptions"->"-march=native -finline-functions -funroll-loops -flto -fPIC -O3"};*)

(*Compiler`$CCompilerOptions =Compiler`$CCompilerOptions={
	"ShellCommandFunction"->Print, 
	"SystemCompileOptions"->" -fPIC -O2"};*)


(* ::Section::Closed:: *)
(*Testing the function*)


SetDirectory[NotebookDirectory[]];
res = Import["Phase_Ds_order_0_to_3.mx"];


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	ClearAll[Test\[CapitalPsi]];
	Test\[CapitalPsi][
		\[Omega]ref_, \[Omega]_, \[Delta]_, \[Chi]s_, \[Chi]a_, 
		\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_
	] = res[[1,2]]//.Join[{G -> g}];

]
DownValues[Test\[CapitalPsi]] = DownValues[Test\[CapitalPsi]]//.HoldForm[x_]:> x;


(* ::Subsection::Closed:: *)
(*Testing the Phase against Ripple*)


DeleteObject/@ExternalSessions[]

Clear@python

python = StartExternalSession["Python"];
ExternalEvaluate[python,"
import numpy as np

from ripplegw.waveforms import IMRPhenomD as IMRD
from ripplegw.waveforms import IMRPhenomD_utils as IMRD_utils
"]


Clear[helperCoeffs, RippleCoeffs]
helperCoeffs = ExternalFunction[python, "def Coeffs(theta):
	a = IMRD_utils.get_coeffs(theta)
	return np.array(a)"
];

RippleCoeffs[\[Theta]_] := helperCoeffs[\[Theta]]//Normal


Clear[helperTransitionFrequencies, RippleTransitionFrequencies]
helperTransitionFrequencies = ExternalFunction[python, "def transitionfrequencies(theta, gamma2, gamma3):
	a = IMRD_utils.get_transition_frequencies(theta, gamma2, gamma3)
	return np.array(a)
"];

RippleTransitionFrequencies[\[Theta]_, \[Gamma]2_, \[Gamma]3_] := helperTransitionFrequencies[\[Theta], \[Gamma]2, \[Gamma]3]//Normal


Clear[helperArg, Argument]
helperArg = ExternalFunction[python, "def h0(f, \[Theta]in, \[Theta]extr, coeffs,  fref):
	b = np.array(f)
	h0, Psi = IMRD._gen_IMRPhenomD(b, \[Theta]in, \[Theta]extr, coeffs,  fref)
	return np.array(Psi)"
];

(*the function takes the arguments: Mc, eta, chi1, chi2, dist_mpc, tc, phic, inclination*)
Argument[f_, \[Theta]in_, \[Theta]ex_, coeffs_, fref_] := helperArg[f, \[Theta]in, \[Theta]ex, coeffs, fref]//Normal


Test := Module[
	{f, M, \[Chi]1, \[Chi]2,pos, m1, m2,\[Eta],tc, \[Phi]c, Ripple, MMA, \[Omega], \[Theta]in, \[Theta]ex, coeffs,  \[Delta], \[Chi]s, \[Chi]a,
	G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],diff, transition,
	ringdown,\[Omega]ref,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4},

	f = Range[20, 2048, 1.];
	{m1, m2} = ReverseSort@RandomReal[{10,120},2];
	
	M = (m1+m2);
	
	\[Eta] = (m1 m2)/M^2;
	tc =0;  
	\[Phi]c = 0; 
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2];
	
	\[Delta] = Sqrt[1-4 \[Eta]];
	\[Chi]s = (\[Chi]1+\[Chi]2)/2; \[Chi]a = (\[Chi]1-\[Chi]2)/2;
	\[Omega] = G M f;
	\[Omega]ref = 20 G M;
	
	

	\[Theta]ex = {1, tc, \[Phi]c}//N;
	\[Theta]in = {m1, m2, \[Chi]1, \[Chi]2};
	coeffs = RippleCoeffs[\[Theta]in];
	transition = RippleTransitionFrequencies[\[Theta]in,  Sequence@@coeffs[[6;;7]] ];

	Ripple = Argument[f, \[Theta]in, \[Theta]ex, coeffs, 20.];
	
	
	MMA = Test\[CapitalPsi][\[Omega]ref, \[Omega], \[Delta], \[Chi]s, \[Chi]a, Sequence@@ConstantArray[0, 15]];
	
	diff = RelativeDiff@@{MMA, Ripple};
	ringdown = transition[[-2]] M G;
	
	
	Ripple = Riffle[\[Omega], Ripple]//Partition[#,2]&;
	MMA = Riffle[\[Omega],MMA]//Partition[#,2]&;

	pos = FirstPosition[\[Omega], x_/; x>=0.19]//Last; (*0.2 is the upper cutoff for IMRPhenomD. *)
	diff = Riffle[\[Omega], diff]//Partition[#,2]&;
	{
		ListLinePlot[
			Take[diff, pos], 
			GridLines->{{{0.018,Red}, {ringdown/2, Red}}, None}, 
			PlotRange->All, ImageSize->Medium, Background->White
		],
		ListLinePlot[
			{Take[Ripple, pos], Take[MMA, pos]}, GridLines->{{{0.018,Red}, {ringdown/2, Red}}, None},PlotRange->All,
			PlotLegends->{"Python", "MMA"}, ImageSize->Medium, Background->White
		]
	}
	

]


Test


(* ::Subsection::Closed:: *)
(*Comparing numerical x Symbolic derivatives*)


NGrad//Clear

NGrad[f_, vars_, n_] := Module[
	{h = 1. 10^-6, dummy, Point1, Point2, denominator},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = If[vars[[i]]==0, h, vars[[i]] + h vars[[i]]];
		dummy,
		{i, n+1, Length@vars}
	];
	
	
	Point2 = ConstantArray[vars, (Length[vars] - n)];
	
	denominator = Table[If[vars[[i]]==0, h, h vars[[i]]], {i, n+1, Length@vars}];

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



Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	ClearAll[TestGrad\[CapitalPsi]];
	
	TestGrad\[CapitalPsi][\[Omega]ref_, \[Omega]_, \[Delta]_, \[Chi]s_, \[Chi]a_, 
		\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_
		] = res[[2;;20, 2]]//.G -> g;

]
DownValues[TestGrad\[CapitalPsi]] = DownValues[TestGrad\[CapitalPsi]]//.HoldForm[x_]:> x;


Clear@Test
Test := Module[
	{
		\[Delta], \[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100}, 2]];
	\[Eta] = (m1 m2)/(m1+m2)^2;
	\[Delta] = Sqrt[1 - 4 \[Eta]];
	{\[Chi]s,\[Chi]a} = RandomReal[{-1,1}, 2];
	
	
	f = RandomReal[{10., 0.2/(G (m1+m2))}];
	\[Omega] = (m1+m2) G f;
	\[Omega]ref = 10 (m1+m2) G;
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4} = RandomReal[{-10, 10}, 15];
	
	Symbolic = TestGrad\[CapitalPsi][\[Omega]ref, \[Omega], \[Delta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4];
	vars = {\[Omega]ref, \[Omega], \[Delta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4};
	Numeric = NGrad[Test\[CapitalPsi], vars, 1];
	
	RelativeDiff@@{Symbolic, Numeric}
	

	
]


Test//ScientificForm


Table[Test, {100}]//MinMax


(* ::Section::Closed:: *)
(*Compiling*)


SetDirectory[NotebookDirectory[]];
Ds = Import["Phase_Ds_order_0_to_3.mx"];


(* ::Subsection::Closed:: *)
(*Compiling Phase terms:*)


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis[x_HoldForm] := Module[
	{dummy, vars},
	vars = {{\[Omega],  _Real,  1}, \[Omega]ref, \[Delta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4};
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
	]//.{HoldForm[y_] :> y, us\[Theta]->UnitStep, Rational->Divide, Piecewise[{{val_, cond_}}, else_]-> If[cond, val, else]};
	
	Compile@@dummy
]


<<CompiledFunctionTools`


compiledDs = MapAt[
	compileThis,
	Ds, 
	{All, 2}
];


<<CCompilerDriver`
<<CCodeGenerator`


$CCompilerDefaultDirectory = FileNameJoin[{
	ParentDirectory[NotebookDirectory[], 2],
	"/LibraryResources/",
	$SystemID,
	"DerivativeRules/IMRPhenomD/NRules"
}]


list = MapIndexed[
	LibraryGenerate[#1[[2]], "phi" <> ToString[#2//First]]&,
	compiledDs
];


(* ::Subsection::Closed:: *)
(*Defining phase terms for RosettaStone*)


Module[
	{d, direc = ParentDirectory[NotebookDirectory[], 2], dummy, system=$SystemID},
	dummy = FileNameJoin[{direc, "/LibraryResources", system,"DerivativeRules/IMRPhenomD/NRules"}];
	d =FileNames["phi*", {
		dummy
	}]; 
	
	list = SortBy[(StringReplace[FileBaseName[#], "phi"->""]//ToExpression)&]@d;
	
	(*THIS IS IMPORTANT IT TAKES FROM THE FILE NAME EVERETHING BEFORE "SymDALI/...":*)
	list = FileNameDrop[#, 5]&/@list
];


SetDirectory[NotebookDirectory[]]


Dterms =Module[{Ds = Import["Phase_Ds_order_0_to_3.mx"]}, Ds[[All,1]] ];


phis = Block[
	{Cvariables},
	Cvariables= ConstantArray[{Real, 0},19];
	Cvariables = Join[{{Real,1}}, Cvariables];

MapThread[
	(#1 -> LF[
		#2, 
		FileBaseName[#2], 
		Cvariables,
		{Real, 1}
	])&,
	{
		Dterms, 
		list
	}
]
];


(*put the rules in the correct format*)
phis2 = phis//.{($D[{n__}, \[CapitalPhi]IMR][x__] -> LF[y__]) :>  TagRule[\[CapitalPhi]IMR, $D[{n}, \[CapitalPhi]IMR][x],  LF[y]]};


phis2[[1,1]]


(*more than 1 derivative in \[Delta]pi is zero because each appears linearly*)
sym1 = $D[{x__}, \[CapitalPhi]IMR]/; Total[{x}[[6;;-1]]] > 1 -> 0

(*more than 1 derivative in \[Delta]\[CurlyPhi]i \[And] \[Chi]s/\[Chi]a is also zero*)
sym2 = $D[{x__}, \[CapitalPhi]IMR]/; (Total[{x}[[6;;15]]] ==1 && Total[{x}[[4;;5]]]> 0) -> 0

(*non-zero derivatives in \[Omega] and \[Omega]ref are automatically 0*)
sym3 = $D[{x__}, \[CapitalPhi]IMR]/;({x}[[1]] > 0 && {x}[[2]] > 0)  ->  0


(*
	Consider that \[CapitalPsi] = \[CapitalPhi]IMR[\[Omega]] - \[CapitalPhi][\[Omega]ref] - t0 (\[Omega]-\[Omega]ref) and that \[Omega] variables has to be a vector.
	The derivative with respect to \[Omega]ref is minus the derivative of \[Omega] evaluated at the value of \[Omega]ref passed as a vector
*)


					(*\[Omega] is at position 1 and \[Omega]ref is at position 2*)
$D[{x__}, \[CapitalPhi]IMR][y__]/; {x}[[1]] === 0 && {x}[[2]] > 0  -> Aux\[CapitalPhi]IMR2;

AuxRule2 = Aux\[CapitalPhi]IMR2[x__][y__] :>  Module[
	{ds, args},
	
	ds = Join[
		{x}[[{2}]], 
		{0},
		{x}[[3;;-1]]
	];
	
	args = Join[
		{{y}[[{2}]]},
		{0.001}, (*note that this number is irrelevant as the derivative in \[Omega] kills the contribution from \[Omega]ref*)
		{y}[[3;;-1]]
	];
	
	-Last[
		$D[ds, \[CapitalPhi]IMR]@@args
	]
]


phis3 = Join[
	phis2,
	{
		TagRule[\[CapitalPhi]IMR, $D[{x__}, \[CapitalPhi]IMR]/; {x}[[1]] === 0 && {x}[[2]] > 0 , Aux\[CapitalPhi]IMR2[x]]
	}
];


(* ::Chapter:: *)
(*Amplitude*)


(* ::Section::Closed:: *)
(*PPN*)


(* Coefficients *)
A0 = 1;

A1 = 0;

A2 = -323/224 + (451 \[Eta])/168//.{\[Eta] -> (1-\[Delta]^2)/4}//Simplify;

A3 = (27 \[Delta] \[Chi]a)/8 + (27/8 - (11 \[Eta])/6) \[Chi]s//.{\[Eta] -> (1-\[Delta]^2)/4}//Simplify;

A4 = (-27312085/8128512 - (1975055 \[Eta])/338688 + (105271 \[Eta]^2)/24192 +
     (-81/32 + 8 \[Eta]) \[Chi]a^2 - (81/16) \[Delta] \[Chi]a \[Chi]s +
     (-81/32 + (17 \[Eta])/8) \[Chi]s^2)//.{\[Eta] -> (1-\[Delta]^2)/4}//Simplify;
     
A5 = (-85 \[Pi]/64 + (85 \[Pi] \[Eta])/16 +
     \[Delta] (285197/16128 - (1579 \[Eta])/4032) \[Chi]a +
     (285197/16128 - (15317 \[Eta])/672 - (2227 \[Eta]^2)/1008) \[Chi]s)//.{\[Eta] -> (1-\[Delta]^2)/4}//Simplify;
     
     
A6 = (-177520268561/8583708672 + 
     ((545384828789/5007163392) - (205 \[Pi]^2)/48) \[Eta] - 
     (3248849057 \[Eta]^2)/178827264 + 
     (34473079 \[Eta]^3)/6386688 +
     (1614569/64512 - (1873643 \[Eta])/16128 + (2167 \[Eta]^2)/42) \[Chi]a^2 +
     (31 \[Pi]/12 - (7 \[Pi] \[Eta])/3) \[Chi]s +
     (1614569/64512 - (61391 \[Eta])/1344 + (57451 \[Eta]^2)/4032) \[Chi]s^2 +
     \[Delta] \[Chi]a (31 \[Pi]/12 + ((1614569/32256) - (165961 \[Eta])/2688) \[Chi]s))//.{\[Eta] -> (1-\[Delta]^2)/4}//Simplify;



(* ::Section::Closed:: *)
(*Inspiral*)


insVecAmplitude[\[Delta]_, \[Chi]s_,\[Chi]a_] = Sqrt[\[Eta]] {A0, A1, A2, A3, A4, A5, A6}//.{\[Eta]->(1-\[Delta]^2)/4}//Simplify;
\[Omega]insVecAmplitude[\[Omega]_] = (\[Pi] \[Omega])^(#/3)&/@Range[0,6];
(*Include \[Omega]^(-7/6)*)
\[Omega]insVecAmplitude[\[Omega]_] = \[Omega]^(-7/6) \[Omega]insVecAmplitude[\[Omega]];


(* ::Text:: *)
(*Light is something*)


{\[ScriptCapitalA]0[\[Delta]_, \[ScriptCapitalM]c_], \[ScriptCapitalA]1[\[Delta]_, \[ScriptCapitalM]c_], \[ScriptCapitalA]2[\[Delta]_, \[ScriptCapitalM]c_], \[ScriptCapitalA]3[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_], \[ScriptCapitalA]4[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_], \[ScriptCapitalA]5[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_], \[ScriptCapitalA]6[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_]} = Block[
	{G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]List, M },
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5)//.\[Eta]->(1-\[Delta]^2)/4;
	
	\[Omega]List = \[Omega]insVecAmplitude[\[Omega]]//. \[Omega] -> M G;
	
	insVecAmplitude[\[Delta], \[Chi]s, \[Chi]a]*M^2*\[Omega]List
];


insvec  = Sqrt[\[Eta]] (PhenomCoeff[\[Eta], \[Chi]PN, #]&/@PhenomDTableV[[1;;3]])//.{\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]->(1-\[Delta]^2)/4};
\[Omega]insvec = \[Omega]^(-7/6)*(\[Omega]^((#+6)/3)&/@Range[3]);


{\[Rho]1[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_],\[Rho]2[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_], \[Rho]3[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_]}= Block[
	{G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]List, M},
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5)//.\[Eta]->(1-\[Delta]^2)/4;
	
	\[Omega]List = \[Omega]insvec//.\[Omega]-> G M;
	
	insvec*M^2*\[Omega]List


];


InsAmpExpr  = {
	U\[ScriptCapitalA]0[\[Delta], \[ScriptCapitalM]c], 0, U\[ScriptCapitalA]2[\[Delta], \[ScriptCapitalM]c],
	U\[ScriptCapitalA]3[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],U\[ScriptCapitalA]4[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],
	U\[ScriptCapitalA]5[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],U\[ScriptCapitalA]6[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a]} . (\[Omega]insVecAmplitude[\[Omega]]//. {\[Omega] -> f, \[Pi]->1}) + {U\[Rho]1[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],U\[Rho]2[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],U\[Rho]3[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a]} . (\[Omega]insvec//.\[Omega]->f);


AuxiliarDefs = <||>;

Clear@ZeroFunction

ZeroFunction[x__] := 0


AuxiliarDefs["Inspiral"] = HoldForm[{
	Unprotect[Derivative],
	
	Derivative[delta_, Mc_, chis_, chia_][\[ScriptCapitalA]3]/;chis+chia>1 := ZeroFunction,
	Derivative[delta_, Mc_, chis_, chia_][\[ScriptCapitalA]4]/;chis+chia>2 := ZeroFunction,
	Derivative[delta_, Mc_, chis_, chia_][\[ScriptCapitalA]5]/;chis+chia>1 := ZeroFunction,
	Derivative[delta_, Mc_, chis_, chia_][\[ScriptCapitalA]6]/;chis+chia>2 := ZeroFunction,
	
	Derivative[delta_, Mc_, chis_, chia_][\[Rho]1]/;chis+chia>3 := ZeroFunction,
	Derivative[delta_, Mc_, chis_, chia_][\[Rho]2]/;chis+chia>3 := ZeroFunction,
	Derivative[delta_, Mc_, chis_, chia_][\[Rho]3]/;chis+chia>3 := ZeroFunction
}];


(* ::Section::Closed:: *)
(*Ringdown and Damping -Phase*)


aeff[\[Delta]_, \[Chi]s_, \[Chi]a_] = Block[
	{S},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S + 2 Sqrt[3.] \[Eta] + (-0.085` S+0.102` S^2-1.355` S^3-0.868` S^4) \[Eta]-4.399` \[Eta]^2+(-5.837` S-2.097` S^2+4.109` S^3+2.064` S^4) \[Eta]^2+9.397` \[Eta]^3-13.181` \[Eta]^4//.\[Eta]-> (1-\[Delta]^2)/4
]//Simplify;

Erad[\[Delta]_, \[Chi]s_, \[Chi]a_] = Block[
	{EradNS=0.0559745 \[Eta]+0.580951 \[Eta]^2-0.960673 \[Eta]^3+3.35241 \[Eta]^4, S},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S = S/(1-2 \[Eta]);
	
	(EradNS (1+(-0.00303023` -2.00661` \[Eta]+7.70506` \[Eta]^2) S))/(1+ (-0.67144` -1.47569` \[Eta] + 7.30468` \[Eta]^2) S)//.\[Eta]-> (1-\[Delta]^2)/4
]//Simplify;

(*re\[Omega] -> interpolation for ringdown and the other is for damping: I think this is frm the PhenomHM paper*)
re\[Omega][\[Chi]_] :=(0.05947169566573468` -0.14989771215394762` \[Chi]+0.09535606290986028` \[Chi]^2+0.02260924869042963` \[Chi]^3-0.02501704155363241` \[Chi]^4-0.005852438240997211` \[Chi]^5+0.0027489038393367993` \[Chi]^6+0.0005821983163192694` \[Chi]^7)/(1-2.8570126619966296` \[Chi]+2.373335413978394` \[Chi]^2-0.6036964688511505` \[Chi]^4+0.0873798215084077` \[Chi]^6);
im\[Omega][\[Chi]_] :=(0.014158792290965177` -0.036989395871554566` \[Chi]+0.026822526296575368` \[Chi]^2+0.0008490933750566702` \[Chi]^3-0.004843996907020524` \[Chi]^4-0.00014745235759327472` \[Chi]^5+0.0001504546201236794` \[Chi]^6)/(1-2.5900842798681376` \[Chi]+1.8952576220623967` \[Chi]^2-0.31416610693042507` \[Chi]^4+0.009002719412204133` \[Chi]^6);

\[Omega]RdDamping[int_, Erad_] = int/(1 - Erad);


\[Gamma]2[\[Delta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[6]]//.{\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]-> (1-\[Delta]^2)/4};
\[Gamma]3[\[Delta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[7]]//.{\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s, \[Eta]-> (1-\[Delta]^2)/4};


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


AuxiliarDefs["rd \[And] dp"] = HoldForm[{
	Derivative[\[Delta]_, chis_, chia_][\[Gamma]2]/; chis+chia>3 := ZeroFunction,
	Derivative[\[Delta]_,chis_, chia_][\[Gamma]3]/; chis+chia>3 := ZeroFunction,
	Derivative[omRD_, omDamp_, gamma2_, gamma3_][\[Omega]Peak]/;omRD + omDamp>1 := ZeroFunction,
	Derivative[omRD_, omDamp_, gamma2_, gamma3_][\[Omega]Peak]/;gamma3 + omDamp > 2 := ZeroFunction 
}];


(* ::Section::Closed:: *)
(*MR Amplitude*)


\[Gamma]1[\[Delta]_, \[Chi]s_, \[Chi]a_] = Block[
	{dummy = PhenomDTableV[[5]], res},
	
	res = PhenomCoeff[\[Eta], \[Chi]PN, dummy]//.{\[Chi]PN -> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};
	
	res Sqrt[\[Eta]]//.\[Eta]->(1 - \[Delta]^2)/4
]; 


Block[
	{
		ringdownFrequency, dampingFrequency, \[Gamma]1 = U\[Gamma]1[\[Delta], \[Chi]s, \[Chi]a], \[Gamma]2 = U\[Gamma]2[\[Delta], \[Chi]s, \[Chi]a], \[Gamma]3 = U\[Gamma]3[\[Delta], \[Chi]s, \[Chi]a],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega],M
	},
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5)//.\[Eta]->(1-\[Delta]^2)/4;
	
	\[Omega] = f M G;
	
	
	ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]];
	
	
	dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]];
	
	
	MRAmpExpr = \[Omega]^(-7/6) (M^2 \[Gamma]1 \[Gamma]3 dampingFrequency)/((\[Omega] - ringdownFrequency)^2 + (\[Gamma]3 dampingFrequency)^2) Exp[-((\[Gamma]2 (\[Omega] - ringdownFrequency))/(\[Gamma]3 dampingFrequency))];


]


AuxiliarDefs["MR"] = HoldForm[{
	Derivative[delta_, chis_, chia_][\[Gamma]1]/; chis+chia>3 := ZeroFunction,
	Derivative[delta_, chis_, chia_][aeff]/; chis+chia>4 := ZeroFunction,
	Derivative[n_, m_][\[Omega]RdDamping]/; n>1 := ZeroFunction
	(*Protect[Derivative]*)
}];


(* ::Section::Closed:: *)
(*Intermediate*)


Block[
	{\[Omega],M},
	
	(*M = \[ScriptCapitalM]c \[Eta]^(-3/5)//.\[Eta]->(1-\[Delta]^2)/4;*)
	
	\[Omega] = G M f;
	
	
	IntAmpfrequency = Sqrt[\[Eta]] M^2 \[Omega]^(-7/6) {1, \[Omega], \[Omega]^2, \[Omega]^3, \[Omega]^4}//Simplify;
]


With[{
	point1 = IntAmpfrequency//.f->f1,
	point2 = IntAmpfrequency//.f->f2,
	point3 = IntAmpfrequency//.f->f3,
	point1D = D[IntAmpfrequency, f]//.f->f1,
	point3D = D[IntAmpfrequency, f]//.f->f3(*,
	M = \[ScriptCapitalM]c ((1-\[Delta]^2)/4)^(-3/5)*)
	},
	
	sol = LinearSolve[{point1, point2, point3, point1D, point3D}, {v1, f2^(-7/6) v2,v3,d1,d3}]//Simplify;
	sol = (sol//.{f1 -> 0.014/(M G)})//Simplify;
]


Clear["\[Delta]*"]


{
	\[Delta]0[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]1[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]2[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]3[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]4[\[Eta]_, M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_]
} = sol;


v2[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = With[
	{
		dummy =(PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[4]])//.{\[Chi]PN-> \[Delta] \[Chi]a +(1-76 \[Eta]/113) \[Chi]s, \[Eta] -> (1-\[Delta]^2)/4},
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
		M = \[ScriptCapitalM]c \[Eta]^(-3/5)
		
	},
	dummy Sqrt[\[Eta]] M^2 (G M)^(-7/6)//.{ \[Eta] -> (1-\[Delta]^2)/4}//Simplify
];


ClearAll[\[CapitalDelta]0, \[CapitalDelta]1, \[CapitalDelta]2, \[CapitalDelta]3, \[CapitalDelta]4]

{
		\[CapitalDelta]0[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]1[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]2[M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]3[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]4[M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_]
	} = { 
		\[Delta]0[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[Delta]1[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[Delta]2[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[Delta]3[\[Eta],M,f2, f3,v1,v2,v3,d1,d3],
		\[Delta]4[\[Eta],M,f2, f3,v1,v2,v3,d1,d3]
		(*\[Delta]i ~ 1/Sqrt[\[Eta]] so \[CapitalDelta]i is not dependent on \[Eta]*)
		}*((Sqrt[\[Eta]] M^2 \[Omega]^(-7/6))*{1, \[Omega], \[Omega]^2, \[Omega]^3, \[Omega]^4}//. \[Omega]-> G M)//Simplify;
	


{
		\[CapitalDelta]0[ M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]1[ M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]2[ M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]3[ M_, f2_,  f3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]4[ M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_]
	} = {
		\[CapitalDelta]0[M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]1[M,f2, f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]2[M,f2, f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]3[M,f2, f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]4[M,f2, f3,v1,v2,v3,d1,d3]
	}//.G-> UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]//Simplify;


Block[
	{ G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
		 
		ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Delta],\[Chi]s, \[Chi]a]], UErad[\[Delta],\[Chi]s, \[Chi]a]],
		dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Delta],\[Chi]s, \[Chi]a]], UErad[\[Delta],\[Chi]s, \[Chi]a]],  peakFrequency, \[Gamma]2, \[Gamma]3,
		f1, f2,f3, v1,v2,v3, d1, d3, \[Omega], M = \[ScriptCapitalM]c \[Eta]^(-3/5)
	},
	
	\[Omega] = M G f;

	\[Gamma]2 = U\[Gamma]2[\[Delta],\[Chi]s, \[Chi]a];
	\[Gamma]3 = U\[Gamma]3[\[Delta],\[Chi]s, \[Chi]a];
	
	peakFrequency = U\[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2, \[Gamma]3];
	f1 = 0.014/(G M);  f3 = peakFrequency/(M G);  f2 =(0.014+peakFrequency)/(2 G M);
	
	v1 = InsAmpExpr//.f-> f1; v3 = MRAmpExpr//.f->f3; v2 = Uv2[\[Delta], \[ScriptCapitalM]c, \[Chi]s, \[Chi]a];
	d1 = D[InsAmpExpr, f]//.f-> f1; d3 = D[MRAmpExpr, f]//.f->f3;
	
	
	IntAmpExpr = {
		U\[CapitalDelta]0[M,f2,f3,v1,v2,v3,d1,d3],
		U\[CapitalDelta]1[M,f2,f3,v1,v2,v3,d1 ,d3],
		U\[CapitalDelta]2[M,f2,f3,v1,v2,v3,d1,d3],
		U\[CapitalDelta]3[M,f2,f3,v1,v2,v3,d1,d3],
		U\[CapitalDelta]4[M,f2,f3,v1,v2,v3,d1,d3]
	} . (f^(-7/6)*{1, f, f^2, f^3, f^4})//.\[Eta]->(1-\[Delta]^2)/4//Simplify;
]


AuxiliarDefs["intermediate"] = HoldForm[{
	Derivative[delta_, Mc_, chis_, chia_][v2]/; chis+chia>3 := ZeroFunction,
	
	
	Derivative[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_][\[CapitalDelta]0]/;v1+v2+v3+d1+d3>1 := ZeroFunction,
	Derivative[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_][\[CapitalDelta]1]/;v1+v2+v3+d1+d3>1 := ZeroFunction,
	Derivative[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_][\[CapitalDelta]2]/;v1+v2+v3+d1+d3>1 := ZeroFunction,
	Derivative[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_][\[CapitalDelta]3]/;v1+v2+v3+d1+d3>1 := ZeroFunction,
	Derivative[M_, f2_, f3_, v1_, v2_, v3_, d1_, d3_][\[CapitalDelta]4]/;v1+v2+v3+d1+d3>1 := ZeroFunction,
	
	Protect[Derivative]
}];


(* ::Section::Closed:: *)
(*\[ScriptCapitalA]IMR*)


Clear@\[ScriptCapitalA]IMR


\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins_, \[ScriptCapitalA]Int_, \[ScriptCapitalA]MR_, f_, M_, \[Omega]Peak_] = Block[
	{c = UnitConvert["SpeedOfLight", ("Gigaparsecs")/("Seconds")][[1]], 
	G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], \[Omega]},
	
	\[Omega] = f M G;
	
	
	Sqrt[5/6] c/(2 \[Pi]^(2/3)) G^2 (
		\[ScriptCapitalA]Ins us\[Theta][0.014 - \[Omega]] + \[ScriptCapitalA]Int us\[Theta][(\[Omega] - 0.014) (\[Omega]Peak - \[Omega])] + \[ScriptCapitalA]MR us\[Theta][(0.2 - \[Omega]) (\[Omega]- \[Omega]Peak)]
	)
];


\[ScriptCapitalA]IMRExpr = Block[
	{ 
		G= UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], 
		ringdownFrequency = U\[Omega]RdDamping[Ure\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]],
		dampingFrequency = U\[Omega]RdDamping[Uim\[Omega][Uaeff[\[Delta], \[Chi]s, \[Chi]a]], UErad[\[Delta], \[Chi]s, \[Chi]a]],  peakFrequency, \[Gamma]2, \[Gamma]3,
		f1, f2,f3, v1,v2,v3, d1, d3, \[Omega], M = \[ScriptCapitalM]c \[Eta]^(-3/5)
	},
	
	\[Omega] = M G f;

	\[Gamma]2 = U\[Gamma]2[\[Delta], \[Chi]s, \[Chi]a];
	\[Gamma]3 = U\[Gamma]3[\[Delta], \[Chi]s, \[Chi]a];
	
	peakFrequency = U\[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2, \[Gamma]3];
	
	U\[ScriptCapitalA]IMR[InsAmpExpr, IntAmpExpr, MRAmpExpr, f, M, peakFrequency]//.\[Eta]->(1-\[Delta]^2)/4
];


(* ::Section::Closed:: *)
(*Making Block function*)


AmpDefs = <||>;


AmpDefs["Ins"] = Block[
	{rule, sym, defs},
	sym = {a0,a2,a3,a4,a5,a6,r1,r2,r3};
	defs = {
		\[ScriptCapitalA]0[\[Delta], \[ScriptCapitalM]c],\[ScriptCapitalA]2[\[Delta], \[ScriptCapitalM]c],\[ScriptCapitalA]3[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],\[ScriptCapitalA]4[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],\[ScriptCapitalA]5[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],\[ScriptCapitalA]6[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],
		\[Rho]1[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],\[Rho]2[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a],\[Rho]3[\[Delta], \[ScriptCapitalM]c,\[Chi]s,\[Chi]a]
	};
	
	rule = MapThread[
		Rule,
		{sym, defs}
	];
	
	
	
	HoldForm[{
		\[ScriptCapitalA]0[\[Delta]_,  \[ScriptCapitalM]c_] = a0, 
		\[ScriptCapitalA]2[\[Delta]_,  \[ScriptCapitalM]c_] = a2, 
		\[ScriptCapitalA]3[\[Delta]_,  \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = a3, 
		\[ScriptCapitalA]4[\[Delta]_,  \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = a4, 
		\[ScriptCapitalA]5[\[Delta]_,  \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = a5, 
		\[ScriptCapitalA]6[\[Delta]_,  \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = a6,
		\[Rho]1[\[Delta]_,  \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = r1,
		\[Rho]2[\[Delta]_,  \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = r2, 
		\[Rho]3[\[Delta]_,  \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = r3
	}]/.rule
];


AmpDefs["MR"] = HoldForm[{
	\[Gamma]1[\[Delta]_, \[Chi]s_, \[Chi]a_] = g1,
	\[Gamma]2[\[Delta]_, \[Chi]s_, \[Chi]a_] = g2,
	\[Gamma]3[\[Delta]_, \[Chi]s_, \[Chi]a_] = g3,
	aeff[\[Delta]_, \[Chi]s_, \[Chi]a_] = a ,
	re\[Omega][\[Chi]_] = r,
	im\[Omega][\[Chi]_] = i,
	Erad[\[Delta]_, \[Chi]s_, \[Chi]a_] = erad,
	\[Omega]RdDamping[int_, Erad_] = s,
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ]
}]/.{
	g1 -> \[Gamma]1[\[Delta], \[Chi]s, \[Chi]a], g2 -> \[Gamma]2[\[Delta], \[Chi]s, \[Chi]a], g3 -> \[Gamma]3[\[Delta], \[Chi]s, \[Chi]a],
	a-> aeff[\[Delta], \[Chi]s, \[Chi]a], r -> re\[Omega][\[Chi]], i-> im\[Omega][\[Chi]], erad-> Erad[\[Delta], \[Chi]s, \[Chi]a], s-> \[Omega]RdDamping[int, Erad]
};


AmpDefs["Int"] = Block[
	{rule, defs},
	
	defs = {
		\[CapitalDelta]0[M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]1[M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]2[M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]3[M,f2,f3,v1,v2,v3,d1,d3],
		\[CapitalDelta]4[M,f2,f3,v1,v2,v3,d1,d3],
		v2[\[Delta], \[ScriptCapitalM]c, \[Chi]s, \[Chi]a]
	};
	
	rule = MapThread[
		Rule,
		{{D0,D1,D2,D3,D4,V2},defs}
	];

	HoldForm[{
		\[CapitalDelta]0[M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] = D0,
		\[CapitalDelta]1[M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] = D1,
		\[CapitalDelta]2[M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] = D2,
		\[CapitalDelta]3[M_,f2_, f3_, v1_, v2_, v3_, d1_, d3_] = D3,
		\[CapitalDelta]4[M_,f2_,f3_, v1_, v2_, v3_, d1_, d3_] = D4,
		v2[\[Delta]_, \[ScriptCapitalM]c_, \[Chi]s_, \[Chi]a_] = V2
	}]/.rule

];


AmpDefs["\[ScriptCapitalA]IMR"] = HoldForm[{
	\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins_, \[ScriptCapitalA]Int_, \[ScriptCapitalA]MR_, f_, M_, \[Omega]Peak_] = x
}]/.x->\[ScriptCapitalA]IMR[\[ScriptCapitalA]Ins, \[ScriptCapitalA]Int, \[ScriptCapitalA]MR, f, M, \[Omega]Peak];


AmpDefs["Commute D"] = HoldForm[{
	(*Commute D with \[Phi]IMR*)
	\[ScriptCapitalA]IMR/: D[\[ScriptCapitalA]IMR[x__], n___] := With[
		{args = D[#, n]&/@({x}[[1;;3]])},
		\[ScriptCapitalA]IMR@@(Join[args, {x}[[4;;6]]])
	]
}];


(*Make rules to eliminate U*)
AllFunctionHeads = (AmpDefs//Values)[[All, 1, All, 1,0]]//Flatten//DeleteDuplicates;
AllFunctionHeads = DeleteElements[AllFunctionHeads, {Symbol}];


Module[{undefined =  ("U" <> #)&/@(ToString/@AllFunctionHeads)},
	
	undefined = ToExpression[undefined,InputForm];
	Headrules = Thread@Rule[undefined, AllFunctionHeads]
];


FinalExpr = Module[
	{dummy},
	dummy = HoldForm[
		{aIMR}
	]/.aIMR-> \[ScriptCapitalA]IMRExpr*{(1+Cos[\[Iota]]^2)/2, I Cos[\[Iota]]};

	dummy//.Headrules
];


defs = Module[
	{functionDefs = AmpDefs//Values, orderedDefs},
	
	functionDefs = Join[functionDefs,Values@AuxiliarDefs ];
	

	functionDefs = HoldForm@@@#&/@functionDefs; (*{HoldForm[HoldForm[]], HoldForm[HoldForm[]],...}*)
	orderedDefs = {functionDefs, HoldForm@@@FinalExpr}//Flatten; (*{HoldForm[HoldForm[]], HoldForm[HoldForm[]],...}*)
	orderedDefs = Flatten[HoldForm@@orderedDefs]; (*HoldForm[allDefs]*)
	orderedDefs = CompoundExpression@@@(HoldForm[Evaluate@orderedDefs]) 
];


declarations = HoldForm[Evaluate@{AllFunctionHeads}];


$blockexpr = $Block@@@HoldForm[Join[declarations, defs]//Evaluate];


expr =  Hold[
	{test, {\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota]}, 3},
	Evaluate@$blockexpr
]//.HoldForm[X_] :> X;


<<FelipeBarbosa`SymDALI`


res1 = EchoTiming[DerivativeRules@@expr];


Unprotect[Derivative];
SubValues[Derivative] = Drop[SubValues[Derivative], {2, -1}]
Protect[Derivative];


SetDirectory[NotebookDirectory[]]
Export["Amp_Ds_order_0_to_3.wdx", res1]


Export["Amp_Ds_order_0_to_3.mx", res1]


itest[i_] := Block[
	{testF, \[Phi]1, \[Phi]2, ND\[Omega],\[Omega], \[Eta], \[Chi]s, \[Chi]a, DM, \[Delta], \[ScriptCapitalM]c,
	M, G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],
    grad, ngrad, f},
	
	testF[\[ScriptCapitalM]c_, \[Delta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = res1[[i, 2]]; DownValues[testF] = DownValues[testF]//. HoldForm[x_]:> x;

	
	\[Eta] = RandomReal[{0.1, 0.24}];
	\[Delta] = Sqrt[1-4 \[Eta]];
	{\[Chi]s, \[Chi]a} = RandomReal[{-1,1}, 2];
	M = RandomReal[{20, 100}];
	\[ScriptCapitalM]c = M \[Eta]^(3/5);
	f = RandomReal[{20, 0.2/(G M)}];

	
	testF[\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, 2 \[Iota]]
]//Quiet


(*Compiler`$CCompilerOptions =Compiler`$CCompilerOptions={
	"ShellCommandFunction"->Print, 
	"SystemCompileOptions"->" -fPIC -O2", 
	"CleanIntermediate"->True};*)


(* ::Section::Closed:: *)
(*Testing the function*)


SetDirectory[NotebookDirectory[]];
res = Import["Amp_Ds_order_0_to_3.mx"];


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	ClearAll[Test\[ScriptCapitalA]];
	Test\[ScriptCapitalA][f_, \[ScriptCapitalM]c_, \[Delta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = res[[1,2]]//.Join[{G -> g}];

]
DownValues[Test\[ScriptCapitalA]] = DownValues[Test\[ScriptCapitalA]]//.HoldForm[x_]:> x;


(* ::Subsection::Closed:: *)
(*Testing against Ripple*)


Clear@python
DeleteObject/@ExternalSessions[]


python = StartExternalSession[
	{

		"System" -> "Python"
    }
];


ExternalEvaluate[python, {
	   "import numpy as np",
            "from ripplegw.waveforms import IMRPhenomD as IMRD",
            "from ripplegw.waveforms import IMRPhenomD_utils as IMRD_utils"
       }]


helperCoeffs = ExternalFunction[python, "def Coeffs(theta):
	a = IMRD_utils.get_coeffs(theta)
	return np.array(a)"
];

RippleCoeffs[\[Theta]_] := helperCoeffs[\[Theta]]//Normal

helperArg = ExternalFunction[python, "def h0(f, \[Theta]in, \[Theta]extr, coeffs,  fref):
	b = np.array(f)
	a, psi = IMRD._gen_IMRPhenomD(b, \[Theta]in, \[Theta]extr, coeffs,  fref)
	return np.array(a)"
];


Amplitude[f_, \[Theta]in_, \[Theta]ex_, coeffs_, fref_] := Abs[helperArg[f, \[Theta]in, \[Theta]ex, coeffs, fref]//Normal]


helperTransitionFrequencies = ExternalFunction[python, "def transitionfrequencies(theta, gamma2, gamma3):
	a = IMRD_utils.get_transition_frequencies(theta, gamma2, gamma3)
	return np.array(a)
"];

RippleTransitionFrequencies[\[Theta]_, \[Gamma]2_, \[Gamma]3_] := helperTransitionFrequencies[\[Theta], \[Gamma]2, \[Gamma]3]//Normal


test\[ScriptCapitalA] := Module[
	{f, M, \[Chi]1, \[Chi]2,pos, m1, m2,\[Eta], \[Theta], Ripple, MMA, \[Omega], \[Theta]in, \[Theta]ex, coeffs,  ripplefreqs,\[Omega]p,posPeak,diff,
	\[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a,
	G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]
	},

	{m1, m2} = ReverseSort@RandomReal[{10,120},2];
	M = (m1+m2); f = Range[20, 0.2/(M G), 0.25];
	\[Eta] = (m1 m2)/M^2;
	\[ScriptCapitalM]c = M \[Eta]^(3/5);
	\[Delta] = Sqrt[1-4 \[Eta]];
	\[Chi]s = (\[Chi]1+\[Chi]2)/2; \[Chi]a = (\[Chi]1-\[Chi]2)/2;
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2];
	\[Omega] = G M f;


	\[Theta]ex = {1, 0, 0};
	\[Theta]in = {m1, m2, \[Chi]1, \[Chi]2};
	coeffs = RippleCoeffs[\[Theta]in];
	ripplefreqs = RippleTransitionFrequencies[\[Theta]in, coeffs[[6]], coeffs[[7]]];
	\[Omega]p = ripplefreqs[[4]] G M;
	

	Ripple = Amplitude[f, \[Theta]in, \[Theta]ex, coeffs,  20];
	MMA = 10^3 Test\[ScriptCapitalA][f, \[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, 0][[1]];
	diff = RelativeDiff[Ripple, MMA];
	
	
	Ripple = Riffle[\[Omega], Ripple]//Partition[#,2]&;
		(*\[ScriptCapitalA] ~ 1/dL, MMA unities are GPC and dL is set to 1 Mpc*)
		
	MMA = Riffle[\[Omega],MMA]//Partition[#,2]&;
	
	diff = Riffle[\[Omega],diff]//Partition[#,2]&;

	pos = FirstPosition[\[Omega], x_/; x>=0.2]//Last; (*0.2 is the upper cutoff for IMRPhenomD. *)
	posPeak = FirstPosition[\[Omega], x_/; x>=\[Omega]p]//Last;

	{
		ListLinePlot[
			diff, 
			GridLines->{{{0.014,Red}, {\[Omega]p, Red}}, None}, 
			PlotRange->All, ImageSize->Medium, Background->White, ScalingFunctions->"Log10"
		],
		
		ListLinePlot[
			{MMA, Ripple},  PlotLegends->{"MMA", "Python"}, ImageSize->Medium,PlotRange->All,
			GridLines -> {{{0.014, Red}, {\[Omega]p, Red}}, None},
			Background->White
		]
	}
]



test\[ScriptCapitalA]


(* ::Subsection:: *)
(*Comparing numerical x Symbolic derivatives*)


NGrad//Clear

NGrad[f_, vars_, n_] := Module[
	{h = 1. 10^-6, dummy, Point1, Point2, denominator},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = If[vars[[i]]==0, h, vars[[i]] + h vars[[i]]];
		dummy,
		{i, n+1, Length@vars}
	];
	
	
	Point2 = ConstantArray[vars, (Length[vars] - n)];
	
	denominator = Table[If[vars[[i]]==0, h, h vars[[i]]], {i, n+1, Length@vars}];

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



Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	ClearAll[TestGrad\[ScriptCapitalA]];
	
	TestGrad\[ScriptCapitalA][f_, \[ScriptCapitalM]c_, \[Delta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = res[[2;;6, 2]]//.G->g;

]
DownValues[TestGrad\[ScriptCapitalA]] = DownValues[TestGrad\[ScriptCapitalA]]//.HoldForm[x_]:> x;


Clear@Test

Test := Module[
	{
		 M, m1, m2, \[ScriptCapitalM]c, \[Delta], f, \[Eta],\[Iota], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100}, 2]]; M = m1+m2;
	\[Eta] = (m1 m2)/(m1+m2)^2;
	M = m1+m2;
	\[ScriptCapitalM]c = M \[Eta]^(3/5); \[Delta] = Sqrt[1-4 \[Eta]];
	{\[Chi]s,\[Chi]a} = RandomReal[{-1,1}, 2];
	
	f = RandomReal[{10., 0.2/(G (m1+m2))}];
	\[Iota] = RandomReal[{0,\[Pi]}];
	
	Symbolic = TestGrad\[ScriptCapitalA][f, \[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota]][[All,1]];
	vars = {f, \[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota]};
	Numeric = NGrad[Test\[ScriptCapitalA], vars, 1][[All,1]];
	
	RelativeDiff@@{Symbolic, Numeric}
]


Test//ScientificForm


Table[Round@Test, {100}]//DeleteDuplicates


(* ::Section:: *)
(*Compiling*)


SetDirectory[NotebookDirectory[]]


Ds = Import["Amp_Ds_order_0_to_3.mx"];


(* ::Subsection::Closed:: *)
(*Compiling Amp terms:*)


Ds[[1,1]]


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{f, _Real, 1}, \[ScriptCapitalM]c, \[Delta], \[Chi]s, \[Chi]a, \[Iota]}, 
		Evaluate[x],
		RuntimeOptions->{
			"CatchMachineOverflow"->False,
			"CatchMachineIntegerOverflow"->False,
			"EvaluateSymbolically"->False,
			"RuntimeErrorHandler"->None,
			"WarningMessages"->True
		}
	]//.{HoldForm[y_] :> y, us\[Theta]->UnitStep, Rational[any_, any2_] :> Divide[any, any2], Complex[x1_, x2_] :> x1+I x2};
	
	Compile@@dummy
]


<<CompiledFunctionTools`


compiledDs = MapAt[
	compileThis,
	Ds, 
	{All, 2}
];


CompilePrint[compiledDs[[20,2]]]


<<CCompilerDriver`
<<CCodeGenerator`


$CCompilerDefaultDirectory = FileNameJoin[{
	ParentDirectory[NotebookDirectory[], 2],
	"/LibraryResources/",
	$SystemID,
	"/DerivativeRules/IMRPhenomD/NRules"
}]


list = MapIndexed[
	LibraryGenerate[#1[[2]], "a" <> ToString[#2//First]]&,
	compiledDs
];


(* ::Subsection::Closed:: *)
(*Defining amplitude for Rosetta stone*)


Module[
	{d, direc = ParentDirectory[NotebookDirectory[], 2], system = $SystemID, dummy},
	dummy = FileNameJoin[{direc,"/LibraryResources", system,"/DerivativeRules/IMRPhenomD/NRules"}];
	d =FileNames["a*", {
		dummy
	}];
	
	list\[ScriptCapitalA] = SortBy[(StringReplace[FileBaseName[#], "a"->""]//ToExpression)&]@d;
	
	list\[ScriptCapitalA] = FileNameDrop[#,5]&/@list\[ScriptCapitalA];
];


Dterms\[ScriptCapitalA] =Module[{Ds = Import["Amp_Ds_order_0_to_3.mx"]}, Ds[[All,1]]];


(* ::Text:: *)
(*Add f_ to the patterns : *)


Dterms\[ScriptCapitalA] = (Join[#[[0]][f_], #]&/@Dterms\[ScriptCapitalA]);


(* ::Text:: *)
(*There are no derivatives in f, so we have to correct for that: *)


Dterms\[ScriptCapitalA] = Dterms\[ScriptCapitalA]/.$D[{n__}, test][y__] -> $D[{0, n}, test][y];


\[ScriptCapitalA]s = MapThread[
	(#1 -> LF[
		#2, 
		FileBaseName[#2], 
		{{Real, 1}, {Real, 0}, {Real, 0}, {Real, 0}, {Real, 0}, {Real, 0}},
		{Complex, 2}
	])&,
	{Dterms\[ScriptCapitalA], list\[ScriptCapitalA]}
];


\[ScriptCapitalA]s = \[ScriptCapitalA]s//.test->\[ScriptA]IMR;


\[ScriptCapitalA]s2 = \[ScriptCapitalA]s//.{($D[{n__}, \[ScriptA]IMR][x__] -> LF[y__]) :>  TagRule[\[ScriptA]IMR, $D[{n}, \[ScriptA]IMR][x],  LF[y]]};


iNRules = <||>;


iNRules[\[CapitalPhi]IMR] = phis3;
iNRules[\[ScriptA]IMR] = \[ScriptCapitalA]s2;
iNRules[Aux\[CapitalPhi]IMR2] = {AuxRule2};


ParentDirectory[NotebookDirectory[],2]


Module[
	{name = ParentDirectory[NotebookDirectory[],2], system = $SystemID},
	
	name = FileNameJoin[{name,"LibraryResources/", system,  "/DerivativeRules/IMRPhenomD/NRules/RosettaStone.wdx"}];
	
	Export[name, iNRules]
]


iSymRules = <||>;

iSymRules[\[CapitalPhi]IMR] = {sym1, sym2, sym3};


Module[
	{name = ParentDirectory[NotebookDirectory[],2], system = $SystemID},
	
	name = FileNameJoin[{name,"LibraryResources/",  system,  "/DerivativeRules/IMRPhenomD/SymRules/file.wdx"}];
	
	Export[name, iSymRules]
]
