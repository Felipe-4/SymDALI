(* ::Package:: *)

Quit


SetOptions[EvaluationNotebook[], LightDark->"Light"]


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


DZeroFunction//ClearAll

DZeroFunction[\[Omega]_?VectorQ, x__] := ConstantArray[0, Length@\[Omega]]
DZeroFunction[\[Omega]_?NumberQ, x__] := 0


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


(*function to compare complex vectors and generate plots*)
Clear@CompareComplexVectors

CompareComplexVectors[{v1_, v2_}, {name1_, name2_}, dataRange_, plotRange_]/;(
	Length[v1] == Length[v2] && Length[v1] == Length[dataRange]
) := Module[
	{Rev1, Rev2, Imv1, Imv2, Rediff, Imdiff, rePlot, imPlot, reDiffPlot, imDiffPlot},
	
	(*real parts*)
	Rev1 = Re[v1];
	Rev2 = Re[v2];
	
	(*im parts*)
	Imv1 = Im[v1];
	Imv2 = Im[v2];
	
	Rediff = RelativeDiff@@{Rev1, Rev2};
	Imdiff = RelativeDiff@@{Imv2, Imv1};
	
	
	
	rePlot = ListLinePlot[
		{Rev1, Rev2},
		PlotLegends->Placed[{name1, name2}, {Right, Top}],
		ImageSize->Medium,
		
		
		PlotLabel->Style["Re", Black],
		PlotRange->{plotRange, All},
		DataRange->MinMax[dataRange], 
		Frame->True,
		Background->White
	];
	imPlot = ListLinePlot[
		{Imv1, Imv2},
		PlotLegends->Placed[{name1, name2}, {Right, Top}],
		ImageSize->Medium,
		PlotLabel->Style["Im", Black],
		PlotRange->{plotRange, All},
		DataRange->MinMax[dataRange],
		Frame->True,
		Background->White
		];
		
	reDiffPlot = ListLinePlot[
		Rediff,
		ImageSize->Medium,
		PlotLabel->Style["RelativeDiff Re part", Black],
		PlotRange->{plotRange, All},
		DataRange->MinMax[dataRange],
		Frame->True,
		ScalingFunctions->"Log10",
		Background->White
	];
		
	imDiffPlot = ListLinePlot[
		Imdiff,
		ImageSize->Medium,
		PlotLabel->Style["RelativeDiff Im part", Black],
		PlotRange->{plotRange, All},
		DataRange->MinMax[dataRange],
		GridLinesStyle->Directive[Red, 13, Dashed],
		Frame->True,
		ScalingFunctions->"Log10",
		Background->White
	];
	
	Grid[{{rePlot, imPlot}, {reDiffPlot, imDiffPlot}}]


]


(* ::Chapter::Closed:: *)
(*Phase*)


(* ::Section::Closed:: *)
(*Inspiral *)


(* ::Subsection::Closed:: *)
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
\[Phi]2[\[Eta]_] = (3715/756 + 55 \[Eta]/9);

\[Phi]3[\[Eta]_, \[Chi]s_, \[Chi]a_] = (-16 \[Pi] + (113 \[Delta] \[Chi]a)/3 + (113/3 - 76 \[Eta]/3) \[Chi]s)//.{\[Delta]->Sqrt[1-4 \[Eta]]}//Simplify;

\[Phi]4[\[Eta]_, \[Chi]s_, \[Chi]a_] = (15293365/508032 + 27145 \[Eta]/504 + 3085 \[Eta]^2/72 + (-405/8 + 200 \[Eta]) \[Chi]a^2 - (405/4) \[Delta] \[Chi]a \[Chi]s + (-405/8 + 5 \[Eta]/2) \[Chi]s^2)//.{
	\[Delta]-> Sqrt[1 - 4 \[Eta]]
}//Simplify;

\[Phi]5[\[Eta]_, \[Chi]s_, \[Chi]a_] = (1 + Log[\[Pi] M \[Omega]]) * (38645 \[Pi]/756 - 65 \[Pi] \[Eta]/9 + 
    \[Delta] (-732985/2268 - 140 \[Eta]/9) \[Chi]a + 
    (-732985/2268 + 24260 \[Eta]/81 + 340 \[Eta]^2/9) \[Chi]s)//.{
	\[Delta]->Sqrt[1- 4 \[Eta]],
	Log[x_] ->0
}//Simplify;

\[Phi]6[\[Eta]_, \[Chi]s_, \[Chi]a_] = (11583231236531/4694215680 - (6848 EulerGamma)/21 - 
   (640 \[Pi]^2)/3 + (-15737765635/3048192 + (2255 \[Pi]^2)/12) \[Eta] + 
   76055 \[Eta]^2/1728 - 127825 \[Eta]^3/1296 - 
   (6848/63) Log[64 \[Pi] M \[Omega]] + (2270/3) \[Pi] \[Delta] \[Chi]a + 
   ((2270 \[Pi])/3 - 520 \[Pi] \[Eta]) \[Chi]s)//.{
	\[Delta]->Sqrt[1 - 4 \[Eta]],
	Log[x_]-> Log[64]
}//Simplify;

\[Phi]7[\[Eta]_, \[Chi]s_, \[Chi]a_] = (77096675 \[Pi]/254016 + (378515 \[Pi] \[Eta])/1512 - (74045 \[Pi] \[Eta]^2)/756 + 
   \[Delta] (-25150083775/3048192 + (26804935 \[Eta])/6048 - (1985 \[Eta]^2)/48) \[Chi]a + 
   (-25150083775/3048192 + (10566655595 \[Eta])/762048 - 
      (1042165 \[Eta]^2)/3024 + (5345 \[Eta]^3)/36) \[Chi]s)//.{
      \[Delta]->Sqrt[1- 4 \[Eta]]
}//Simplify;


(* ::Subsection::Closed:: *)
(*Ins-Phase*)


(* ::Text:: *)
(*In  the  notation  of  the  arXiv : 1903.04467  we  want  a  vector*)
(*{\[CurlyPhi]minus2, \[CurlyPhi]0, \[CurlyPhi]1, \[CurlyPhi]2, \[CurlyPhi]3, \[CurlyPhi]4, \[CurlyPhi]5, \[CurlyPhi]5l, \[CurlyPhi]6, \[CurlyPhi]6l, \[CurlyPhi]7} . In  GR  there  is  no  \[CurlyPhi]minus2  or  \[CurlyPhi]1, we  shall  set  their  values  to  1*3/(128 \[Eta])*)
(*and  multiply  by  a  vector  1 + \[Delta]\[CurlyPhi]  where  the  default  value  of  \[Delta]\[CurlyPhi]minus2  and  \[Delta]\[CurlyPhi]1  is - 1 (to recover GR) and -1+\[Delta]\[CurlyPhi] when we want to sample on them:*)


Clear[insVecPhase, \[Omega]InsVecPhase]


Block[
	{},
	insVecPhase[\[Eta]_, \[Chi]s_, \[Chi]a_] = 3*{
		1, (*\[CurlyPhi]minus2*)
		1, (*\[CurlyPhi]0*)
		1, (*\[CurlyPhi]1*)
		\[Phi]2[\[Eta]], (*\[CurlyPhi]2*)
		\[Phi]3[\[Eta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]3*)
		\[Phi]4[\[Eta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]4*)
		\[Phi]5[\[Eta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]5*)
		\[Phi]5[\[Eta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]5l*)
		\[Phi]6[\[Eta], \[Chi]s, \[Chi]a], (*\[CurlyPhi]6*)
		-(6848/63), (*\[CurlyPhi]6l*)
		\[Phi]7[\[Eta], \[Chi]s, \[Chi]a] (*\[CurlyPhi]7*)
	}/(128 \[Eta]);
	
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
	{v = 1 + {-1 + \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, (-1 + \[Delta]\[CurlyPhi]1), \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, 0, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}},
	
	
	Clear@InsVecPhase;
	
	
	InsVecPhase[\[Eta]_, \[Chi]s_, \[Chi]a_,\[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_] = (
		v*(insVecPhase[\[Eta], 0, 0]) + (insVecPhase[\[Eta], \[Chi]s, \[Chi]a] - insVecPhase[\[Eta], 0, 0])
	)//Simplify
];


Clear["\[CapitalPhi]*"]


{
	\[CapitalPhi]minus2[\[Eta]_, \[Delta]\[CurlyPhi]minus2_], \[CapitalPhi]0[\[Eta]_, \[Delta]\[CurlyPhi]0_], \[CapitalPhi]1[\[Eta]_, \[Delta]\[CurlyPhi]1_], \[CapitalPhi]2[\[Eta]_, \[Delta]\[CurlyPhi]2_], \[CapitalPhi]3[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]3_], 
	\[CapitalPhi]4[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]4_], \[CapitalPhi]5[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalPhi]5l[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]5l_], \[CapitalPhi]6[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]6_], 
	\[CapitalPhi]6l[\[Eta]_, \[Delta]\[CurlyPhi]6l_], \[CapitalPhi]7[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]7_]
} = InsVecPhase[\[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];


Clear["\[CapitalSigma]*"]

{\[CapitalSigma]1[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]2[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]3[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]4[\[Eta]_, \[Chi]s_, \[Chi]a_]} = \[Eta]^-1*(PhenomCoeff[\[Eta], \[Chi]PN, #]&/@PhenomDTableV[[8;;11]])//.{
	\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s}//Simplify;


{\[CapitalSigma]1[\[Delta], \[Chi]s, \[Chi]a], \[CapitalSigma]2[\[Delta], \[Chi]s, \[Chi]a],\[CapitalSigma]3[\[Delta], \[Chi]s, \[Chi]a],\[CapitalSigma]4[\[Delta], \[Chi]s, \[Chi]a]} . {\[Omega], 3/4 \[Omega]^(4/3), 3/5 \[Omega]^(5/3), 1/2 \[Omega]^2};


InsExpr = -\[Pi]/4 + {
	U\[CapitalPhi]minus2[\[Eta],\[Delta]\[CurlyPhi]minus2],U\[CapitalPhi]0[\[Eta],\[Delta]\[CurlyPhi]0],U\[CapitalPhi]1[\[Eta],\[Delta]\[CurlyPhi]1],U\[CapitalPhi]2[\[Eta],\[Delta]\[CurlyPhi]2],U\[CapitalPhi]3[\[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]3],U\[CapitalPhi]4[\[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]4],U\[CapitalPhi]5[\[Eta], \[Chi]s, \[Chi]a],
	U\[CapitalPhi]5l[\[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]5l],U\[CapitalPhi]6[\[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]6],U\[CapitalPhi]6l[\[Eta],\[Delta]\[CurlyPhi]6l],U\[CapitalPhi]7[\[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]7]
} . \[Omega]InsVecPhase[\[Omega]]  + {U\[CapitalSigma]1[\[Eta], \[Chi]s, \[Chi]a],U\[CapitalSigma]2[\[Eta], \[Chi]s, \[Chi]a],U\[CapitalSigma]3[\[Eta], \[Chi]s, \[Chi]a],U\[CapitalSigma]4[\[Eta], \[Chi]s, \[Chi]a]} . {\[Omega], 3/4 \[Omega]^(4/3), 3/5 \[Omega]^(5/3), 1/2 \[Omega]^2};


InsExpr


\[CapitalPhi]0[\[Eta], \[Delta]\[CurlyPhi]0]


InsExpr;


Expr = HoldComplete[
	{{
		\[CapitalPhi]minus2, \[CapitalPhi]0, \[CapitalPhi]1, \[CapitalPhi]2, \[CapitalPhi]3, \[CapitalPhi]4, \[CapitalPhi]5, \[CapitalPhi]5l, \[CapitalPhi]6, \[CapitalPhi]6l, \[CapitalPhi]7,
		\[CapitalSigma]1, \[CapitalSigma]2, \[CapitalSigma]3, \[CapitalSigma]4
	}},
	
	\[CapitalPhi]minus2[\[Eta]_, \[Delta]\[CurlyPhi]minus2_] := pminus2;
	\[CapitalPhi]0[\[Eta]_, \[Delta]\[CurlyPhi]0_] := p0;
	\[CapitalPhi]1[\[Eta]_, \[Delta]\[CurlyPhi]1_] := p1;
	\[CapitalPhi]2[\[Eta]_, \[Delta]\[CurlyPhi]2_] := p2;
	\[CapitalPhi]3[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]3_] := p3;
	\[CapitalPhi]4[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]4_] := p4;
	\[CapitalPhi]5[\[Eta]_, \[Chi]s_, \[Chi]a_] :=p5;
	\[CapitalPhi]5l[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]5l_] :=p5l;
	\[CapitalPhi]6[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]6_] := p6;
	\[CapitalPhi]6l[\[Eta]_, \[Delta]\[CurlyPhi]6l_] := p6l;
	\[CapitalPhi]7[\[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]7_] := p7;
	
	\[CapitalSigma]1[\[Eta]_, \[Chi]s_, \[Chi]a_] :=  S1;
	\[CapitalSigma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := S2;
	\[CapitalSigma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := S3;
	\[CapitalSigma]4[\[Eta]_, \[Chi]s_, \[Chi]a_] := S4;
	
	-(\[Pi]/4)+
	(
		\[Omega] \[CapitalSigma]1[\[Eta],\[Chi]s,\[Chi]a] + 3/4 \[Omega]^(4/3) \[CapitalSigma]2[\[Eta],\[Chi]s,\[Chi]a]+3/5 \[Omega]^(5/3) \[CapitalSigma]3[\[Eta],\[Chi]s,\[Chi]a]+1/2 \[Omega]^2 \[CapitalSigma]4[\[Eta],\[Chi]s,\[Chi]a]
	) +
	(
		\[CapitalPhi]0[\[Eta],\[Delta]\[CurlyPhi]0]/(\[Pi]^(5/3) \[Omega]^(5/3))+\[CapitalPhi]1[\[Eta],\[Delta]\[CurlyPhi]1]/(\[Pi]^(4/3) \[Omega]^(4/3))+\[CapitalPhi]2[\[Eta],\[Delta]\[CurlyPhi]2]/(\[Pi] \[Omega])+\[CapitalPhi]3[\[Eta],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]3]/(\[Pi]^(2/3) \[Omega]^(2/3))+
		\[CapitalPhi]4[\[Eta],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]4]/(\[Pi]^(1/3) \[Omega]^(1/3))+\[CapitalPhi]5[\[Eta],\[Chi]s,\[Chi]a]+Log[\[Pi] \[Omega]] \[CapitalPhi]5l[\[Eta],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]5l]+
		\[Pi]^(1/3) \[Omega]^(1/3) \[CapitalPhi]6[\[Eta],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]6]+\[Pi]^(1/3) \[Omega]^(1/3) Log[\[Pi] \[Omega]] \[CapitalPhi]6l[\[Eta],\[Delta]\[CurlyPhi]6l]+
		\[Pi]^(2/3) \[Omega]^(2/3) \[CapitalPhi]7[\[Eta],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]7] + \[CapitalPhi]minus2[\[Eta],\[Delta]\[CurlyPhi]minus2]/(\[Pi]^(7/3) \[Omega]^(7/3))
	)
]//.{
	pminus2 ->  \[CapitalPhi]minus2[\[Eta], \[Delta]\[CurlyPhi]minus2], 
	p1 -> \[CapitalPhi]1[\[Eta], \[Delta]\[CurlyPhi]1],
	p0-> \[CapitalPhi]0[\[Eta], \[Delta]\[CurlyPhi]0], p2-> \[CapitalPhi]2[\[Eta],  \[Delta]\[CurlyPhi]2], p3-> \[CapitalPhi]3[\[Eta], \[Chi]s, \[Chi]a,  \[Delta]\[CurlyPhi]3], p4-> \[CapitalPhi]4[\[Eta], \[Chi]s, \[Chi]a,  \[Delta]\[CurlyPhi]4], p5-> \[CapitalPhi]5[\[Eta], \[Chi]s, \[Chi]a],
	p5l-> \[CapitalPhi]5l[\[Eta], \[Chi]s, \[Chi]a,  \[Delta]\[CurlyPhi]5l], p6-> \[CapitalPhi]6[\[Eta], \[Chi]s, \[Chi]a,  \[Delta]\[CurlyPhi]6], p6l-> \[CapitalPhi]6l[\[Eta],  \[Delta]\[CurlyPhi]6l],p7-> \[CapitalPhi]7[\[Eta], \[Chi]s, \[Chi]a,  \[Delta]\[CurlyPhi]7],
	
	S1 ->\[CapitalSigma]1[\[Eta], \[Chi]s, \[Chi]a], S2 ->\[CapitalSigma]2[\[Eta], \[Chi]s, \[Chi]a], S3->\[CapitalSigma]3[\[Eta], \[Chi]s, \[Chi]a], S4->\[CapitalSigma]4[\[Eta], \[Chi]s, \[Chi]a]
};


Expr = $Block@@@(HoldForm[Evaluate@Expr]);


{\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7};


(*All derivatives up to order 6*)

derivatives = Combinations[vars, 5];

derivatives//Length


(* ::Text:: *)
(*Bcs \[Delta]pi appears linearly, all terms with more than 1 derivative in \[Delta]pi are 0.*)


numberOf\[Delta]p[elem_List] := Module[
	{\[Delta]ps = {\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}},
	
	Count[
		elem, 
		x_/; Or@@(Equal[x, #]&/@\[Delta]ps)
	]
]


Filter\[Delta]ps[derivatives_List]:= Module[
	{iL, pos, iList},
	
	iL = numberOf\[Delta]p/@derivatives;
	
	pos = Position[iL, x_/; x<=1];
	
	iList = Extract[derivatives, pos]; (*exclude \[Delta]pi^2*)
	
	iList
	
]

derivatives = Filter\[Delta]ps[derivatives];

derivatives//Length

PrependTo[derivatives, {}];

derivatives//Length


D\[CapitalPhi]Ins//ClearAll

expr = Hold[
	{D\[CapitalPhi]Ins, {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7},
	 derivatives
	},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalPhi]resInspiral = EchoTiming[DerivativeRules@@expr];


\[CapitalPhi]resInspiral[[-1]]


Clear@T

T[i_] := Block[
	{test},
	
	test[
		\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_,
		\[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_
	] = \[CapitalPhi]resInspiral[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, 
		\[Delta]\[CurlyPhi]minus2$->\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0$->\[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1$->\[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2$->\[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3$->\[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4$->\[Delta]\[CurlyPhi]4,
		\[Delta]\[CurlyPhi]5$->\[Delta]\[CurlyPhi]5, \[Delta]\[CurlyPhi]5l$->\[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6$->\[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l$->\[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7$->\[Delta]\[CurlyPhi]7
	};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, 0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`]
]


Table[T[i], {i, Length@\[CapitalPhi]resInspiral}];//AbsoluteTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, 
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

compileThis[x_Integer] := DZeroFunction


c\[CapitalPhi]resInspiral = MapAt[
	compileThis,
	\[CapitalPhi]resInspiral,
	{All,2}
];


Do[
	c\[CapitalPhi]resInspiral[[i]][[2, -2]] = Function[{\[Omega],\[Eta],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, X]//.{X -> ToString[c\[CapitalPhi]resInspiral[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPhi]resInspiral}
] (*The error comes from derivatives that are zero.*)


Clear@T

T[i_] := Block[
	{test},
	
	test[
		\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_,
		\[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_
	] := c\[CapitalPhi]resInspiral[[i,2]][\[Omega],\[Eta],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];
	
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, 0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`,0.1`]
]


Length@\[CapitalPhi]resInspiral


Table[T[i], {i, Length@\[CapitalPhi]resInspiral}];//AbsoluteTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = c\[CapitalPhi]resInspiral[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@D\[CapitalPhi]Ins = list;
]


(* ::Text:: *)
(*Set DownValues:*)


D\[CapitalPhi]Ins[
	\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, 
	\[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_
] = c\[CapitalPhi]resInspiral[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];


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
{\[Beta]1[\[Eta]_, \[Chi]s_, \[Chi]a_], \[Beta]2[\[Eta]_, \[Chi]s_, \[Chi]a_],  \[Beta]3[\[Eta]_, \[Chi]s_, \[Chi]a_]} = IntVecPhase[\[Eta], \[Chi]PN]//.{
	\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s
};


Expr = HoldComplete[
	{{\[Beta]1, \[Beta]2, \[Beta]3}},
	
	\[Beta]1[\[Eta]_, \[Chi]s_, \[Chi]a_] := b1;
	\[Beta]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := b2;
	\[Beta]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := b3;
	
	\[Omega] \[Beta]1[\[Eta],\[Chi]s,\[Chi]a]+ Log[\[Omega]] \[Beta]2[\[Eta],\[Chi]s,\[Chi]a]- 1/3 \[Beta]3[\[Eta],\[Chi]s,\[Chi]a] \[Omega]^-3
]//.{
	b1->\[Beta]1[\[Eta], \[Chi]s, \[Chi]a],
	b2->\[Beta]2[\[Eta], \[Chi]s, \[Chi]a],
	b3->\[Beta]3[\[Eta], \[Chi]s, \[Chi]a]
};


Expr = $Block@@@(HoldForm[Evaluate@Expr]);


D\[CapitalPhi]Int//ClearAll
expr = Hold[
	{D\[CapitalPhi]Int, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 5},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalPhi]resInt = EchoTiming[DerivativeRules@@expr];


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = \[CapitalPhi]resInt[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@\[CapitalPhi]resInt}]//RepeatedTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a}, 
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

compileThis[x_Integer] := DZeroFunction


c\[CapitalPhi]resInt = MapAt[
	compileThis,
	\[CapitalPhi]resInt,
	{All,2}
];


Do[
	c\[CapitalPhi]resInt[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[c\[CapitalPhi]resInt[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPhi]resInt}
] (*The error comes from derivatives that are zero.*)


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := c\[CapitalPhi]resInt[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@c\[CapitalPhi]resInt}]//RepeatedTiming


Module[
	{list = c\[CapitalPhi]resInt[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@D\[CapitalPhi]Int = list;
]


D\[CapitalPhi]Int[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = c\[CapitalPhi]resInt[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


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


\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[6]]//.{
	\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s,
	\[Delta]->Sqrt[1-4 \[Eta]]
};
\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[7]]//.{
	\[Chi]PN-> \[Delta] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s,
	\[Delta]->Sqrt[1-4 \[Eta]]
};


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


(* ::Section::Closed:: *)
(*MR-Phase*)


Clear@Erad
Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] = Block[
	{EradNS=0.0559745 \[Eta]+0.580951 \[Eta]^2-0.960673 \[Eta]^3+3.35241 \[Eta]^4, S, \[Delta]=Sqrt[1 - 4 \[Eta]]},
	
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S = S/(1-2 \[Eta]);
	
	(EradNS (1+(-0.00303023` -2.00661` \[Eta]+7.70506` \[Eta]^2) S))/(1+ (-0.67144` -1.47569` \[Eta] + 7.30468` \[Eta]^2) S)(*//.\[Eta]-> (1-\[Delta]^2)/4*)
]//Simplify;


Clear@aeff
aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] = Block[
	{S, \[Delta]=Sqrt[1-4 \[Eta]]},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S + 2 Sqrt[3.] \[Eta] + (-0.085` S+0.102` S^2-1.355` S^3-0.868` S^4) \[Eta]-4.399` \[Eta]^2+(-5.837` S-2.097` S^2+4.109` S^3+2.064` S^4) \[Eta]^2+9.397` \[Eta]^3-13.181` \[Eta]^4(*//.\[Eta]-> (1-\[Delta]^2)/4*)
]//Simplify;


{\[Alpha]1[\[Eta]_, \[Chi]s_, \[Chi]a_], \[Alpha]2[\[Eta]_, \[Chi]s_, \[Chi]a_], \[Alpha]3[\[Eta]_, \[Chi]s_, \[Chi]a_], \[Alpha]4[\[Eta]_, \[Chi]s_, \[Chi]a_]} = \[Eta]^-1*(PhenomCoeff[\[Eta], \[Chi]PN, #]&/@PhenomDTableV[[-5;;-2]])//.{
	\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s}//Simplify;


\[Alpha]5[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[19]]//.{\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};


Expr = HoldComplete[
	{{\[Omega]RdDampinglm, re\[Omega], im\[Omega], Erad, aeff, \[Alpha]1, \[Alpha]2, \[Alpha]3, \[Alpha]4, \[Alpha]5}},
	
	Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
	aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
	re\[Omega][\[Kappa]_] := re;
	im\[Omega][\[Kappa]_] := im;
	\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
	
	\[Alpha]1[\[Eta]_, \[Chi]s_, \[Chi]a_] := a1;
	\[Alpha]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := a2;
	\[Alpha]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := a3;
	\[Alpha]4[\[Eta]_, \[Chi]s_, \[Chi]a_] := a4;
	\[Alpha]5[\[Eta]_, \[Chi]s_, \[Chi]a_] := a5;
	
	(
		\[Omega] \[Alpha]1[\[Eta],\[Chi]s,\[Chi]a]-
		\[Alpha]2[\[Eta],\[Chi]s,\[Chi]a]/\[Omega]+
		4/3 \[Omega]^(3/4) \[Alpha]3[\[Eta],\[Chi]s,\[Chi]a]+
		( ArcTan[((\[Omega]-ringdownFrequency \[Alpha]5[\[Eta],\[Chi]s,\[Chi]a]))/(dampingFrequency)] \[Alpha]4[\[Eta],\[Chi]s,\[Chi]a])
	)
]//.{
	a1->\[Alpha]1[\[Eta], \[Chi]s, \[Chi]a], a2->\[Alpha]2[\[Eta], \[Chi]s, \[Chi]a], a3->\[Alpha]3[\[Eta], \[Chi]s, \[Chi]a],a4->\[Alpha]4[\[Eta], \[Chi]s, \[Chi]a], a5->\[Alpha]5[\[Eta], \[Chi]s, \[Chi]a],
	
	erad->Erad[\[Eta], \[Chi]s, \[Chi]a], ae->aeff[\[Eta], \[Chi]s, \[Chi]a],
	re->re\[Omega][\[Kappa]], im->im\[Omega][\[Kappa]],
	ringdownFrequency :> \[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	dampingFrequency :> \[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]]
};

Expr = $Block@@@(HoldForm[Evaluate@Expr]);


ClearAll@D\[CapitalPhi]MR
expr = Hold[
	{D\[CapitalPhi]MR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 5},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalPhi]resMR = EchoTiming[DerivativeRules@@expr];


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = \[CapitalPhi]resMR[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a
	};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@\[CapitalPhi]resMR}]//AbsoluteTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a}, 
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


c\[CapitalPhi]resMR = MapAt[
	compileThis,
	\[CapitalPhi]resMR,
	{All,2}
];


Do[
	c\[CapitalPhi]resMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[c\[CapitalPhi]resMR[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPhi]resMR}
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := c\[CapitalPhi]resMR[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@c\[CapitalPhi]resMR}]//AbsoluteTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = c\[CapitalPhi]resMR[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@D\[CapitalPhi]MR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


D\[CapitalPhi]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = c\[CapitalPhi]resMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


(* ::Section::Closed:: *)
(*\[CapitalPhi]IMR*)


Block[
	{\[Beta]0, \[Beta]1, \[Alpha]0, \[Alpha]1},
	
	\[Beta]1 = Derivative[1,0,0,0, Sequence@@ConstantArray[0, 10]][\[Phi]Ins][0.018`, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7] - Derivative[1,0,0,0][\[Phi]Int][0.018`, \[Eta], \[Chi]s, \[Chi]a];
	\[Beta]0 = \[Phi]Ins[0.018, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7] - \[Phi]Int[0.018`, \[Eta], \[Chi]s, \[Chi]a] - \[Beta]1 0.018`;
	
	\[Alpha]1 = (
		Derivative[1,0,0,0][\[Phi]Int][ringdownFrequency/2, \[Eta], \[Chi]s, \[Chi]a] +
		\[Beta]1 - 
		Derivative[1,0,0,0][\[Phi]MR][ringdownFrequency/2,  \[Eta], \[Chi]s, \[Chi]a]
	)//Simplify;
	
	\[Alpha]0 = (
		\[Phi]Int[ringdownFrequency/2, \[Eta], \[Chi]s, \[Chi]a] + 
		\[Beta]0 + \[Beta]1 ringdownFrequency/2 - 
		\[Phi]MR[ringdownFrequency/2, \[Eta], \[Chi]s, \[Chi]a] - 
		\[Alpha]1 ringdownFrequency/2
	)//Simplify;
	
	Expr = HoldComplete[
		{{\[Omega]RdDamping, re\[Omega], Erad, aeff, \[Phi]Ins, \[Phi]Int, \[Phi]MR}},
	
		Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
		aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
		re\[Omega][\[Kappa]_] := re;
		\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
		
		\[Phi]Ins[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] := U\[CapitalPhi]Ins[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7];
		\[Phi]Int[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[CapitalPhi]Int[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
		\[Phi]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[CapitalPhi]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	
		\[Phi]Ins[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7]*us\[Theta][0.018`- \[Omega]] + 
		(\[Phi]Int[\[Omega], \[Eta], \[Chi]s, \[Chi]a] + b0 + b1 \[Omega])*us\[Theta][(\[Omega]-0.018`) (ringdownFrequency/2-\[Omega])] + 
		(\[Phi]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a] + a0 + a1 \[Omega])*us\[Theta][(\[Omega]-ringdownFrequency/2) (0.2-\[Omega])]
		
	]//.{
		b0->\[Beta]0, b1->\[Beta]1, 
		a0->\[Alpha]0, a1->\[Alpha]1,
		erad->Erad[\[Eta], \[Chi]s, \[Chi]a], ae->aeff[\[Eta], \[Chi]s, \[Chi]a], kappa->\[Kappa][aeff, l, m],
		re->re\[Omega][\[Kappa]], 
		ringdownFrequency :> \[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]]
	};
]

Expr = $Block@@@(HoldForm[Evaluate@Expr]);


vars = {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7};


derivatives = Combinations[vars, 4];

derivatives//Length


(* ::Text:: *)
(*Bcs \[Delta]pi appears linearly, all terms with more than 1 derivative in \[Delta]pi are 0.*)


numberOf\[Delta]p[elem_List] := Module[
	{\[Delta]ps = {\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}},
	
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


D\[CapitalPhi]IMR//ClearAll
expr = Hold[
	{D\[CapitalPhi]IMR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}, derivatives},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalPhi]resIMR = EchoTiming[DerivativeRules@@expr];


Module[
	{temp, rule},
	
	rule = {
		Derivative[n__][x_] -> $D[{n}, x], 
		
		$D[{n__}, x_][y1_, y2__]/;Length[y1]==0 && y1[[0]] === Real :> Last[$D[{n}, x][{y1}, y2]],
		$D[{n__}, x_][y1_, y2__]/;Length[y1]>0 && y1[[0]] === Times :> Last[$D[{n}, x][{y1}, y2]],
		
		
		U\[CapitalPhi]Ins[y1_, y2__]/; y1[[0]] === Real :> Last[U\[CapitalPhi]Ins[{y1}, y2]],
		
		
		U\[CapitalPhi]Int[y1_, y2__]/;y1[[0]] === Real :> Last[U\[CapitalPhi]Int[{y1}, y2]],
		U\[CapitalPhi]Int[y1_, y2__]/;Length[y1]>0 && y1[[0]] === Times :> Last[U\[CapitalPhi]Int[{y1}, y2]],
		
		
		U\[CapitalPhi]MR[y1_, y2__]/;Length[y1]>0 && y1[[0]] ===Times :> Last[U\[CapitalPhi]MR[{y1}, y2]]
	
	
	(*, U\[ScriptCapitalA]Ins->\[ScriptCapitalA]Ins, U\[ScriptCapitalA]MR->\[ScriptCapitalA]MR*)};
	
	
	temp = \[CapitalPhi]resIMR//.rule;
	
	\[CapitalPhi]resIMR2 = temp//.{U\[CapitalPhi]Ins->D\[CapitalPhi]Ins, U\[CapitalPhi]Int->D\[CapitalPhi]Int, U\[CapitalPhi]MR->D\[CapitalPhi]MR};
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] = \[CapitalPhi]resIMR2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a,
		\[Delta]\[CurlyPhi]minus2$->\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0$->\[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1$->\[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2$->\[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3$->\[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4$->\[Delta]\[CurlyPhi]4,
		\[Delta]\[CurlyPhi]5$->\[Delta]\[CurlyPhi]5, \[Delta]\[CurlyPhi]5l$->\[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6$->\[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l$->\[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7$->\[Delta]\[CurlyPhi]7
	};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2,Sequence@@ConstantArray[1., 10]]
]


Table[T[i], {i, Length@\[CapitalPhi]resIMR2}]//AbsoluteTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, 
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


c\[CapitalPhi]resIMR = MapAt[
	compileThis,
	\[CapitalPhi]resIMR2,
	{All,2}
];


Do[
	c\[CapitalPhi]resIMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, X]//.{X -> ToString[c\[CapitalPhi]resIMR[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPhi]resIMR}
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_,  \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] := c\[CapitalPhi]resIMR[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, Sequence@@ConstantArray[1., 10]]
]


Table[T[i], {i, Length@c\[CapitalPhi]resIMR}]//RepeatedTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = c\[CapitalPhi]resIMR[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	UpValues@D\[CapitalPhi]IMR = {};
	UpValues@D\[CapitalPhi]IMR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


D\[CapitalPhi]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] = c\[CapitalPhi]resIMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];


(* ::Section::Closed:: *)
(*\[CapitalPsi]*)


Clear@\[Phi]0


Expr = HoldComplete[
	{{
		 Erad, aeff, \[Omega]RdDamping, \[Omega]Peak, \[Phi]MR, re\[Omega], im\[Omega], \[Gamma]2, \[Gamma]3, \[Phi]IMR
	}},
	
	\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := 1.010344386100769` +0.0008993122028186917` \[Eta]+(0.2839491069316864` -4.049753189086914` \[Eta]+13.207828521728516` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)+(0.10396278649568558` -7.025059223175049` \[Eta]+24.784893035888672` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^2+(0.030932024121284485` -2.6924023628234863` \[Eta]+9.609374046325684` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^3;
	\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := 1.3081616163253784` -0.0055377297103405` \[Eta]+(-0.0678291767835617` -0.668983519077301` \[Eta]+3.4031479358673096` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)+(-0.05296577513217926` -0.9923793077468872` \[Eta]+4.820681095123291` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^2+(-0.0061341398395597935` -0.3842925429344177` \[Eta]+1.7561753988265991` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^3;
	
	\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
	\[Omega]Peak[\[Omega]RD_,\[Omega]DAMP_,\[Gamma]2_,\[Gamma]3_] := If[\[Gamma]2<=1,\[Omega]RD+(\[Omega]DAMP \[Gamma]3 (Sqrt[1-\[Gamma]2^2]-1))/\[Gamma]2,\[Omega]RD-(\[Omega]DAMP \[Gamma]3)/\[Gamma]2];
	
	\[Phi]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[CapitalPhi]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	\[Phi]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_] := U\[CapitalPhi]IMR[
		\[Omega], \[Eta], \[Chi]s, \[Chi]a,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7
	];
	
	re\[Omega][\[Kappa]_] := (0.05947169566573468` -0.14989771215394762` \[Kappa]+0.09535606290986028` \[Kappa]^2+0.02260924869042963` \[Kappa]^3-0.02501704155363241` \[Kappa]^4-0.005852438240997211` \[Kappa]^5+0.0027489038393367993` \[Kappa]^6+0.0005821983163192694` \[Kappa]^7)/(1-2.8570126619966296` \[Kappa]+2.373335413978394` \[Kappa]^2-0.6036964688511505` \[Kappa]^4+0.0873798215084077` \[Kappa]^6);
	im\[Omega][\[Kappa]_] := (0.014158792290965177` -0.036989395871554566` \[Kappa]+0.026822526296575368` \[Kappa]^2+0.0008490933750566702` \[Kappa]^3-0.004843996907020524` \[Kappa]^4-0.00014745235759327472` \[Kappa]^5+0.0001504546201236794` \[Kappa]^6)/(1-2.5900842798681376` \[Kappa]+1.8952576220623967` \[Kappa]^2-0.31416610693042507` \[Kappa]^4+0.009002719412204133` \[Kappa]^6);
	
	Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
	aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
	
	\[Phi]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7] - \[Delta]t0 - \[Phi]0
	
	
]//.{
	erad -> Erad[\[Eta], \[Chi]s, \[Chi]a], 
	ae -> aeff[\[Eta], \[Chi]s, \[Chi]a],
	
	Peak :> \[Omega]Peak[
		\[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]], 
		\[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]], 
		\[Gamma]2[\[Eta], \[Chi]s, \[Chi]a], 
		\[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]
	],
	
	
	\[Delta]t0 -> Derivative[1, 0, 0, 0][\[Phi]MR][
		Peak, \[Eta], \[Chi]s, \[Chi]a
	]*(\[Omega]-\[Omega]ref),
	
	\[Phi]0 :>  \[Phi]IMR[\[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7]

};

Expr = $Block@@@(HoldForm[Evaluate@Expr]);


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7};


(*All derivatives up to order 3*)

derivatives = Combinations[vars, 4];

derivatives//Length


(* ::Text:: *)
(*Bcs \[Delta]pi appears linearly, all terms with more than 1 derivative in \[Delta]pi are 0.*)
(*Also, there is no correlation between \[Delta]\[CurlyPhi]i and \[Chi]s, \[Chi]a bcs the \[Delta]\[Beta]1 term cancels in the C(1) contribution (\[Delta]\[Alpha]1 \[Proportional] \[Delta]\[Beta]1) and there is no correlation between \[Delta]\[CurlyPhi]i \[And] {\[Chi]s, \[Chi]a} originally in \[Phi]Ins*)


numberOf\[Delta]p[elem_List] := Module[
	{\[Delta]ps = {\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}},
	
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


D\[CapitalPsi]//ClearAll
expr = Hold[
	{D\[CapitalPsi], {\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4, \[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7}, derivatives},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalPsi]res = EchoTiming[DerivativeRules@@expr];


Module[
	{rule, exprs, temp},
	
	(*exprs ={
		(0.036` \[Omega]rdlm)/(m $x2568), 
		$x2568, 
		(\[Omega]rdlm (-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568))/(\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568)+(0.018` \[Omega]rdlm (2/m-(-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568)/(\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568)))/$x2568, 
		(0.018` \[Omega]rdlm (-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568))/((\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568) $x2568)+(0.018` \[Omega]rdlm (2/m-(-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568)/(\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568)))/$x2568
	}//.{$x2568->$x4476};*)
	
	rule = {
		Derivative[n__][U\[CapitalPhi]IMR][x__] -> $D[{n}, U\[CapitalPhi]IMR][x],
		Derivative[n__][U\[CapitalPhi]MR][x__] -> $D[{n}, U\[CapitalPhi]MR][x],
		
		$D[{n__}, U\[CapitalPhi]IMR][\[Omega]ref, y2__] :> Last[$D[{n}, U\[CapitalPhi]IMR][{\[Omega]ref}, y2]],
		U\[CapitalPhi]IMR[\[Omega]ref, y2__] :>  Last[U\[CapitalPhi]IMR[{\[Omega]ref}, y2]],
		
		$D[{n__}, U\[CapitalPhi]MR][y1_, y2__]/; y1[[0]] =!= List :> Last[$D[{n}, U\[CapitalPhi]MR][{y1}, y2]]
		
	};
	
	
	
	
	temp =  \[CapitalPsi]res//.rule;
	
	\[CapitalPsi]res2 = temp//.{U\[CapitalPhi]IMR->D\[CapitalPhi]IMR, U\[CapitalPhi]MR-> D\[CapitalPhi]MR};
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] = \[CapitalPsi]res2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Omega]ref$->\[Omega]ref, \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a,
		\[Delta]\[CurlyPhi]minus2$->\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0$->\[Delta]\[CurlyPhi]0, \[Delta]\[CurlyPhi]1$->\[Delta]\[CurlyPhi]1, \[Delta]\[CurlyPhi]2$->\[Delta]\[CurlyPhi]2, \[Delta]\[CurlyPhi]3$->\[Delta]\[CurlyPhi]3, \[Delta]\[CurlyPhi]4$->\[Delta]\[CurlyPhi]4,
		\[Delta]\[CurlyPhi]5$->\[Delta]\[CurlyPhi]5, \[Delta]\[CurlyPhi]5l$->\[Delta]\[CurlyPhi]5l, \[Delta]\[CurlyPhi]6$->\[Delta]\[CurlyPhi]6, \[Delta]\[CurlyPhi]6l$->\[Delta]\[CurlyPhi]6l, \[Delta]\[CurlyPhi]7$->\[Delta]\[CurlyPhi]7
	};
	
	test[Range[10^-4, 0.124, 10^-4], 10^-4, 0.2, 0.1, 0.2,Sequence@@ConstantArray[1., 10]]
]


Table[T[i], {i, Length@\[CapitalPsi]res2}]//AbsoluteTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, 
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


c\[CapitalPsi]2 = MapAt[
	compileThis,
	\[CapitalPsi]res2,
	{All,2}
];


Do[
	c\[CapitalPsi]2[[i]][[2, -2]] = Function[{\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, X]//.{X -> ToString[c\[CapitalPsi]2[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPsi]2}
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_,  \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] := c\[CapitalPsi]2[[i,2]][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];
	
	test[Range[10^-4, 0.124, 10^-4], 10^-4, 0.2, 0.1, 0.2, Sequence@@ConstantArray[1., 10]]
]


Table[T[i], {i, Length@c\[CapitalPsi]2}]//RepeatedTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = c\[CapitalPsi]2[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	UpValues@D\[CapitalPsi] = {};
	
	UpValues@D\[CapitalPsi] = list;
]


(* ::Text:: *)
(*Set DownValues:*)


D\[CapitalPsi][\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_] = c\[CapitalPsi]2[[1, 2]][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];


(* ::Section::Closed:: *)
(*Testing *)


(* ::Subsection::Closed:: *)
(*Testing phase accuracy*)


(* ::Subsubsection::Closed:: *)
(*My def*)


Clear@MMA

MMA[\[ScriptCapitalM]c_, \[Eta]_, \[Chi]1_, \[Chi]2_] := Module[
	{
		M, \[Omega], \[Chi]s, \[Chi]a, freq, \[Omega]ref, res,  G =  4.925490947641267`*^-6
	},
	
	\[Chi]s = (\[Chi]1+\[Chi]2)/2; 
	\[Chi]a = (\[Chi]1-\[Chi]2)/2;
	
	
	freq =Range[5., 1024., 1];
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5);
	
	\[Omega] = M freq G;
	
	\[Omega]ref = Min[\[Omega]];
	
	
	res = D\[CapitalPsi][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, Sequence@@ConstantArray[0,10]];
	
	res
]


(* ::Subsubsection::Closed:: *)
(*Testing the Phase against Ripple*)


DeleteObject/@ExternalSessions[]

Clear@python

python = StartExternalSession[{
	"Python",
	"Evaluator" -> "/Users/felipe/anaconda3/envs/RIPPLE/bin/python"
}];

ExternalEvaluate[python,"
import numpy as np

from ripplegw.waveforms import IMRPhenomD as IMRD
from ripplegw.waveforms import IMRPhenomD_utils as IMRD_utils
"]


ExternalEvaluate[python,"
import ripplegw
print(ripplegw.__file__)
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
	{
		f, M, \[Chi]1, \[Chi]2,pos, m1, m2,\[Eta],tc, \[Phi]c, Ripple, mMA, \[Omega], \[Theta]in, \[Theta]ex, coeffs,  \[Delta], \[Chi]s, \[Chi]a,t0,
		diff, transition,
		ringdown,\[Omega]ref, G =  4.925490947641267`*^-6
	},
	
	
	M = RandomReal[{5, 200}];
	\[Eta] = RandomReal[{0.125, 0.25}];
	{m1, m2} = M/2 {1 + Sqrt[1- 4 \[Eta]], 1 - Sqrt[1- 4 \[Eta]]};
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2];
	
	\[Omega] = G M Range[5, 1024, 1.];
	
	tc =0; 
	\[Phi]c = 0;
	
	\[Theta]ex = {1, tc, \[Phi]c}//N;
	\[Theta]in = {m1, m2, \[Chi]1, \[Chi]2};
	
	coeffs = RippleCoeffs[\[Theta]in];
	
	transition = RippleTransitionFrequencies[\[Theta]in,  Sequence@@coeffs[[6;;7]] ];

	Ripple = Argument[Range[5, 1024, 1.], \[Theta]in, \[Theta]ex, coeffs, 5.];
	
	mMA = MMA[M \[Eta]^(3/5), \[Eta], \[Chi]1, \[Chi]2];
	
	diff = RelativeDiff@@{mMA, Ripple};
	ringdown = transition[[-2]] M G;
	
	
	Ripple = Riffle[\[Omega], Ripple]//Partition[#,2]&;
	mMA = Riffle[\[Omega],mMA]//Partition[#,2]&;

	pos = FirstPosition[\[Omega], x_/; x>=0.19]//Last; (*0.2 is the upper cutoff for IMRPhenomD. *)
	diff = Riffle[\[Omega], diff]//Partition[#,2]&;
	
	{
		ListLinePlot[
			Take[diff, pos], 
			GridLines->{{{0.018,Red}, {ringdown/2, Red}}, None}, 
			PlotRange->All, ImageSize->Medium, Background->White,
			ScalingFunctions->"Log10"
		],
		ListLinePlot[
			{Take[Ripple, pos], Take[mMA, pos]}, GridLines->{{{0.018,Red}, {ringdown/2, Red}}, None},PlotRange->All,
			PlotLegends->{"Python", "MMA"}, ImageSize->Medium, Background->White
		]
	}
	

]


Test


(* ::Subsubsection::Closed:: *)
(*lal hphc function*)


DeleteObject/@ExternalSessions[]

Clear@python

python = StartExternalSession[{
	"Python",
	"Evaluator" -> "/Users/felipe/anaconda3/envs/TIGER/bin/python"
}];

ExternalEvaluate[python,"
import lalsimulation
import bilby

import warnings
import numpy as np
warnings.filterwarnings(\"ignore\", \"Wswiglal-redir-stdio\")
import bilby_tgr
"]


ExternalEvaluate[python,"

# Set the duration and sampling frequency of the data segment that we're
# going to inject the signal into
duration = 4.0
sampling_frequency = 2048.0
minimum_frequency = 5


# Fixed arguments passed into the source model
waveform_arguments = dict(
    waveform_approximant=\"IMRPhenomD\",
    reference_frequency=5,
    minimum_frequency=5,
)

WG = bilby.gw.WaveformGenerator(
    duration=4,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby_tgr.tiger.source.lal_binary_black_hole_TIGER_PhenomP,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments,
)
"]


lalPhase = ExternalFunction[python,"
def lal_phase(m1, m2, chi1, chi2, delta_pi):

    f = np.linspace(5, 1024, 4077)

    params = {
        \"mass_1\": m1,
        \"mass_2\": m2,
        \"chi_1\": chi1,
        \"chi_2\": chi2,
        \"luminosity_distance\": 1.0,
        \"theta_jn\": 1.0,
        \"phase\": 0,

        \"dchi_0\":0,
        \"dchi_1\":0,
        \"dchi_2\": 0,
        \"dchi_3\": 0,
        \"dchi_4\": 0,
        \"dchi_5l\": 0,
        \"dchi_6\": 0,
        \"dchi_6l\": 0,
        \"dchi_7\": delta_pi,

        \"dbeta_2\": 0,
        \"dbeta_3\": 0,

        \"dalpha_2\": 0,
        \"dalpha_3\": 0,
        \"dalpha_4\": 0,
        \"dalpha_5\": 0,
    }

    hs = WG.frequency_domain_strain(params)

    hp_bilby = hs[\"plus\"][20:]
    hc_bilby = hs[\"cross\"][20:]

    phase = hp_bilby / np.abs(hp_bilby)

    return phase
"]


(* ::Subsubsection::Closed:: *)
(*Test*)


(* ::Text:: *)
(*I believe Josiel is using frd and fdamp as calculated in the original PhenomD implementation. I am not convinced this is right. I find it more likely that those should be replaced by the frd and fdamp 22 formulas attached to Z22.*)


Clear@Test

Test[M_] := Module[
	{
		m1, m2, mc, eta, s1z, s2z, \[Iota], \[Phi]Ref, lalRes, fmax, fmin, myDef, amp, t0, chis, chia, psi,  
		G =  4.925490947641267`*^-6, \[Delta]p=RandomReal[{-10,10}]
	},
	
	eta = RandomReal[{0.1, 0.25}];
	m1 = M/2 (1 + Sqrt[1-4 eta]);
	m2 = M/2 (1 - Sqrt[1-4 eta]);
	
	
	mc = (m1+m2) eta^(3/5);
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	
	fmax =  1024;
	fmin=5;
	
	lalRes = lalPhase[m1, m2, s1z, s2z, \[Delta]p]//Normal;
	
	myDef = D\[CapitalPsi][
		M Range[5., 1024, 0.25] G, 
		5 M G,
		eta,
		(s1z+s2z)/2, (s1z-s2z)/2, 
		0, 0,0, 0, 0, 0, 0, 0, 0, \[Delta]p
	];
		
	{lalRes, Exp[-I myDef]}
	
]


Module[
	{M = RandomReal[{50,200}], list, G =  4.925490947641267`*^-6},
	
	list = Test[M];
	
	ListLinePlot[
		RelativeDiff@@list,
		ScalingFunctions->"Log10"
	]//Echo;
	
	CompareComplexVectors[list, {"lal", "MMA"}, Range[5,1024,0.25], {5, 0.2/(G M)}]
	
]


(* ::Subsection::Closed:: *)
(*Comparing numerical x Symbolic derivatives*)


(* ::Subsubsection::Closed:: *)
(*First Order*)


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


Test\[CapitalPsi][\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_,\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_] := D\[CapitalPsi][
	{\[Omega]}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7
]


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	ClearAll[TestGrad\[CapitalPsi]];
	
	TestGrad\[CapitalPsi][
		\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, 
		\[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_
		] = \[CapitalPsi]res2[[2;;15, 2]];

]
DownValues[TestGrad\[CapitalPsi]] = DownValues[TestGrad\[CapitalPsi]]//.HoldForm[x_]:> x;


Clear@Test
Test := Module[
	{
		\[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4, M
	},
	
	M = RandomReal[{50,200}];
	
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2]; 
	
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	{\[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7} = RandomReal[{-1, 1}, 10];
	
	Symbolic = TestGrad\[CapitalPsi][{\[Omega]}, \[Omega]ref,  \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7];
	vars = {\[Omega]ref, \[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2, \[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7};
	Numeric = NGrad[Test\[CapitalPsi], vars, 1];
	
	RelativeDiff@@{Symbolic, Numeric}//Flatten
	

	
]


Test//ScientificForm


Table[Test, {1000}]//MinMax


(* ::Subsubsection::Closed:: *)
(*Second Order*)


(* ::Text:: *)
(*Since First order checks out we can Use the symbolic gradient of first order and compare its numerical derivatives against the symbolic *)
(*gradients of order 2.*)


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	rule = Join[
		{G -> g}, Thread@Rule[
			{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, ConstantArray[0, 10]
		]
	];
	
	ClearAll[iTestGrad\[CapitalPsi]O1];
	
	iTestGrad\[CapitalPsi]O1[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = \[CapitalPsi]res2[[2;;5, 2]]//.rule;

]
DownValues[iTestGrad\[CapitalPsi]O1] = DownValues[iTestGrad\[CapitalPsi]O1]//.HoldForm[x_]:> x;

TestGrad\[CapitalPsi]O1[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := iTestGrad\[CapitalPsi]O1[\[Omega]ref, {\[Omega]}, \[Eta], \[Chi]s, \[Chi]a]


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	rule = Join[
		{G -> g}, Thread@Rule[
			{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, ConstantArray[0, 10]
		]
	];
	
	ClearAll[TestGrad\[CapitalPsi]O2];
	
	TestGrad\[CapitalPsi]O2[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_
		] = Join[ (*these give the order: 
				{{\[Omega],\[Omega]},{\[Eta],\[Omega]},{\[Chi]s,\[Omega]},{\[Chi]a,\[Omega]},{\[Eta],\[Eta]},{\[Eta],\[Chi]s},{\[Eta],\[Chi]a},{\[Chi]s,\[Chi]s},{\[Chi]a,\[Chi]s},{\[Chi]a,\[Chi]a}}
			*)
			\[CapitalPsi]res2[[16;;19, 2]],
			\[CapitalPsi]res2[[30;;32, 2]],
			\[CapitalPsi]res2[[43;;45, 2]]
		]//.rule;

]
DownValues[TestGrad\[CapitalPsi]O2] = DownValues[TestGrad\[CapitalPsi]O2]//.HoldForm[x_]:> x;


Clear@Test

Test := Module[
	{
		\[Delta], \[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		li = SymmetrizedIndependentComponents[{4,4}, Symmetric[All]],M
	},
	
	M = RandomReal[{50, 200}];
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1},2];
	{\[Chi]s,\[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	Symbolic = TestGrad\[CapitalPsi]O2[\[Omega]ref, {\[Omega]}, \[Eta], \[Chi]s, \[Chi]a]//Flatten;
	
	vars = {\[Omega]ref, \[Omega], \[Eta], \[Chi]s, \[Chi]a};
	
	Numeric = Extract[
		ArrayReshape[NGrad[TestGrad\[CapitalPsi]O1, vars, 1], {4,4}],
		li
	];
	
	
	RelativeDiff@@{Symbolic, Numeric}
	

	
]


Test//ScientificForm


Table[Test, {1000}]//MinMax


(* ::Subsubsection::Closed:: *)
(*Third Order*)


(* ::Text:: *)
(*Second order checks out so we can use the symbolic gradient of second order and compare its numerical derivatives against the symbolic gradients of order 3.*)


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	rule = Join[
		{G -> g}, Thread@Rule[
			{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, ConstantArray[0, 10]
		]
	];
	
	ClearAll[iTestGrad\[CapitalPsi]O2];
	
	iTestGrad\[CapitalPsi]O2[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_
		] = Join[ (*these give the order: 
				{{\[Omega],\[Omega]},{\[Eta],\[Omega]},{\[Chi]s,\[Omega]},{\[Chi]a,\[Omega]},{\[Eta],\[Eta]},{\[Eta],\[Chi]s},{\[Eta],\[Chi]a},{\[Chi]s,\[Chi]s},{\[Chi]a,\[Chi]s},{\[Chi]a,\[Chi]a}}
			*)
			\[CapitalPsi]res2[[16;;19, 2]],
			\[CapitalPsi]res2[[30;;32, 2]],
			\[CapitalPsi]res2[[43;;45, 2]]
		]//.rule;

]
DownValues[iTestGrad\[CapitalPsi]O2] = DownValues[iTestGrad\[CapitalPsi]O2]//.HoldForm[x_]:> x;

TestGrad\[CapitalPsi]O2//Clear

TestGrad\[CapitalPsi]O2[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := iTestGrad\[CapitalPsi]O2[\[Omega]ref, {\[Omega]}, \[Eta], \[Chi]s, \[Chi]a]//Flatten


SymmetricTestGradOrder2[x__] := Module[
	{li = SymmetrizedIndependentComponents[{4, 4}, Symmetric[All]],rule, values = TestGrad\[CapitalPsi]O2[x], rules},
	
	rules = MapThread[Rule, {li, values}];
	
	SymmetrizedArray[rules, {4,4},  Symmetric[All]]
]


(* ::Text:: *)
(*Just make sure derivatives are in the same order of \[CapitalPsi]*)


positions = Position[derivatives, x_/; Length[x]===3&&numberOf\[Delta]p[x]===0];


2


Clear@Test

Test := Module[
	{
		\[Delta], \[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		li = SymmetrizedIndependentComponents[{4, 4, 4}, Symmetric[All]], M
	},
	
	M = RandomReal[{50, 200}];
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	
	{\[Chi]1,\[Chi]2} = RandomReal[{-1,1}, 2];
	{\[Chi]s,\[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	Symbolic = TestGrad\[CapitalPsi]O3[\[Omega]ref, \[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	vars = {\[Omega]ref, \[Omega], \[Eta], \[Chi]s, \[Chi]a};
	
	Numeric = Extract[NGrad[SymmetricTestGradOrder2, vars, 1], li];
	
	RelativeDiff@@{Symbolic, Numeric}
]


Test//ScientificForm


Table[Test, {10^3}]//MinMax//EchoTiming


(* ::Subsubsection::Closed:: *)
(*Fourth Order*)


(* ::Text:: *)
(*Second order checks out so we can use the symbolic gradient of second order and compare its numerical derivatives against the symbolic gradients of order 3.*)


positions = Position[derivatives, x_/; Length[x]===3&&numberOf\[Delta]p[x]===0];


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	rule = Join[
		{G -> g}, Thread@Rule[
			{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, ConstantArray[0, 10]
		]
	];
	
	ClearAll[iTestGrad\[CapitalPsi]O3];
	
	iTestGrad\[CapitalPsi]O3[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = \[CapitalPsi]res2[[positions//Flatten, 2]]//.rule;

]
DownValues[iTestGrad\[CapitalPsi]O3] = DownValues[iTestGrad\[CapitalPsi]O3]//.HoldForm[x_]:> x;

TestGrad\[CapitalPsi]O3//Clear

TestGrad\[CapitalPsi]O3[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := iTestGrad\[CapitalPsi]O3[\[Omega]ref, {\[Omega]}, \[Eta], \[Chi]s, \[Chi]a]//Flatten


SymmetricTestGradOrder3[x__] := Module[
	{li = SymmetrizedIndependentComponents[{4, 4, 4}, Symmetric[All]],rule, values = TestGrad\[CapitalPsi]O3[x], rules},
	
	rules = MapThread[Rule, {li, values}];
	
	SymmetrizedArray[rules, {4,4, 4},  Symmetric[All]]
]


positions2 = Position[derivatives, x_/; Length[x]===4&&numberOf\[Delta]p[x]===0];


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	rule = Join[
		{G -> g}, Thread@Rule[
			{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7}, ConstantArray[0, 10]
		]
	];
	
	ClearAll[iTestGrad\[CapitalPsi]O4];
	
	iTestGrad\[CapitalPsi]O4[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = \[CapitalPsi]res2[[positions2//Flatten, 2]]//.rule;

]
DownValues[iTestGrad\[CapitalPsi]O4] = DownValues[iTestGrad\[CapitalPsi]O4]//.HoldForm[x_]:> x;

TestGrad\[CapitalPsi]O4//Clear

TestGrad\[CapitalPsi]O4[\[Omega]ref_, \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := iTestGrad\[CapitalPsi]O4[\[Omega]ref, {\[Omega]}, \[Eta], \[Chi]s, \[Chi]a]//Flatten


Clear@Test

Test := Module[
	{
		\[Delta], \[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		li = SymmetrizedIndependentComponents[{4, 4, 4, 4}, Symmetric[All]], M
	},
	
	M = RandomReal[{50, 200}];
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	
	{\[Chi]1,\[Chi]2} = RandomReal[{-1,1}, 2];
	{\[Chi]s,\[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	Symbolic = TestGrad\[CapitalPsi]O4[\[Omega]ref, \[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	vars = {\[Omega]ref, \[Omega], \[Eta], \[Chi]s, \[Chi]a};
	
	Numeric = Extract[NGrad[SymmetricTestGradOrder3, vars, 1], li];
	
	RelativeDiff@@{Symbolic, Numeric}
]


Test//ScientificForm


Table[Test, {10^3}]//MinMax//EchoTiming


(* ::Subsection::Closed:: *)
(*Testing against the old implementation*)


(* ::Subsubsection::Closed:: *)
(*First order derivatives*)


test = FelipeBarbosa`SymDALI`DALICoefficients`Private`test;


list = list//.{
	FelipeBarbosa`SymDALI`DALICoefficients`Private`\[CapitalPhi]IMR -> \[CapitalPhi]IMR
};


TagSetDelayed@@@list[[2;;-1]];


list = test[FelipeBarbosa`SymDALI`DALICoefficients`Private`\[CapitalPhi]IMR];


\[Phi]IMR[
	\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, 
	\[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_, \[Delta]\[CurlyPhi]1_, \[Delta]\[CurlyPhi]2_, \[Delta]\[CurlyPhi]3_, \[Delta]\[CurlyPhi]4_, \[Delta]\[CurlyPhi]5l_, \[Delta]\[CurlyPhi]6_, \[Delta]\[CurlyPhi]6l_, \[Delta]\[CurlyPhi]7_,
	\[Delta]\[Beta]2_, \[Delta]\[Beta]3_, \[Delta]\[Alpha]2_, \[Delta]\[Alpha]3_, \[Delta]\[Alpha]4_
] = \[CapitalPhi]IMR[
	\[Omega],\[Omega]ref, Sqrt[1 - 4 \[Eta]],\[Chi]s,\[Chi]a,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4
]


Unprotect@Derivative;
Derivative[n__][\[CapitalPhi]IMR] := $D[{n}, \[CapitalPhi]IMR];
Derivative[n__][D\[CapitalPsi]] := $D[{n}, D\[CapitalPsi]];
Derivative[x__][$D[y_List, s_Symbol]][k__] := $D[{x} + y, s][k];
Protect@Derivative;


derivatives[[2;;15]]


DownValues@D\[CapitalPsi]={};


T := Block[
	{
		headList1, headList2, 
		\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,
		G = 4.925490947641267`*^-6, M, \[Chi]1, \[Chi]2
	},
	
	
	headList1 = D[D\[CapitalPsi][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7], #]&/@(derivatives[[2;;15]]);
	headList2 = D[\[Phi]IMR[\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,0,0,0,0,0], #]&/@(derivatives[[2;;15]]);
	
	
	M = RandomReal[{50, 200}];
	
	\[Omega] = RandomReal[{6, 0.15/(G M)}] G M//List;
	\[Omega]ref = 5 G M;
	\[Eta] = RandomReal[{0.1, 0.24999}];
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1},2];
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7} = RandomReal[{-10,10},10];
	
	RelativeDiff@@{
		headList1//Flatten,
		headList2//Flatten
	}
	
]


tt = Table[T, {10^6}];


tt//Max


(* ::Subsubsection::Closed:: *)
(*second order derivatives:*)


Clear@T

T := Block[
	{
		headList1, headList2, 
		\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,
		G = 4.925490947641267`*^-6, M, \[Chi]1, \[Chi]2
	},
	
	Block[
		{D\[CapitalPsi], \[CapitalPhi]IMR},
		headList1 = D[D\[CapitalPsi][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7], Sequence@@#]&/@(derivatives[[16;;45]]);
		headList2 = D[\[Phi]IMR[\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,0,0,0,0,0], Sequence@@#]&/@(derivatives[[16;;45]]);
	
	];
	
	M = RandomReal[{50, 200}];
	
	\[Omega] = RandomReal[{6, 0.15/(G M)}] G M//List;
	\[Omega]ref = 5 G M;
	\[Eta] = RandomReal[{0.1, 0.24999}];
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1},2];
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7} = RandomReal[{-10,10},10];
	
	RelativeDiff@@{
		headList1//Flatten,
		headList2//Flatten
	}
	
	
	
]


tt = Table[T, {10^5}]//EchoTiming;


tt//Max


Flatten@tt//Dimensions


Count[Flatten@tt, x_/; x< 10^-10]


ListPlot[Flatten@tt, ScalingFunctions->"Log10"]


(* ::Subsubsection::Closed:: *)
(*Third order: *)


Clear@T

T := Block[
	{
		headList1, headList2, 
		\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,
		G = 4.925490947641267`*^-6, M, \[Chi]1, \[Chi]2
	},
	
	Block[
		{D\[CapitalPsi], \[CapitalPhi]IMR},
		headList1 = D[D\[CapitalPsi][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7], Sequence@@#]&/@(derivatives[[46;;95]]);
		headList2 = D[\[Phi]IMR[\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,0,0,0,0,0], Sequence@@#]&/@(derivatives[[46;;95]]);
	
	];
	
	M = RandomReal[{50, 200}];
	
	\[Omega] = RandomReal[{6, 0.15/(G M)}] G M//List;
	\[Omega]ref = 5 G M;
	\[Eta] = RandomReal[{0.1, 0.24999}];
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1},2];
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7} = RandomReal[{-10,10},10];
	
	RelativeDiff@@{
		headList1//Flatten,
		headList2//Flatten
	}
	
]


tt = Table[T, 3 10^4]//EchoTiming;


tt//Max


Count[Flatten@tt, x_/; x<10^-10]/Length[Flatten[tt]]//N


ListPlot[Flatten@tt, ScalingFunctions->"Log10"]


(* ::Chapter::Closed:: *)
(*Amplitude*)


(* ::Section::Closed:: *)
(*Inspiral*)


(* ::Subsection::Closed:: *)
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



(* ::Subsection::Closed:: *)
(*Inspiral*)


insVecAmplitude[\[Delta]_, \[Chi]s_,\[Chi]a_] = Sqrt[\[Eta]] {A0, A1, A2, A3, A4, A5, A6}//.{\[Eta]->(1-\[Delta]^2)/4}//Simplify;
\[Omega]insVecAmplitude[\[Omega]_] = (\[Pi] \[Omega])^(#/3)&/@Range[0,6];
(*Include \[Omega]^(-7/6)*)
\[Omega]insVecAmplitude[\[Omega]_] = \[Omega]^(-7/6) \[Omega]insVecAmplitude[\[Omega]];


Clear["\[ScriptCapitalA]*"]

{\[ScriptCapitalA]0[\[Eta]_], \[ScriptCapitalA]1[\[Eta]_], \[ScriptCapitalA]2[\[Eta]_], \[ScriptCapitalA]3[\[Eta]_, \[Chi]s_, \[Chi]a_], \[ScriptCapitalA]4[\[Eta]_, \[Chi]s_, \[Chi]a_], \[ScriptCapitalA]5[\[Eta]_, \[Chi]s_, \[Chi]a_], \[ScriptCapitalA]6[\[Eta]_, \[Chi]s_, \[Chi]a_]} = Block[
	{},
	
	(insVecAmplitude[\[Delta], \[Chi]s, \[Chi]a]//.\[Delta]->Sqrt[1 - 4 \[Eta]])//Simplify
];


insvec  = Sqrt[\[Eta]] (PhenomCoeff[\[Eta], \[Chi]PN, #]&/@PhenomDTableV[[1;;3]])//.{\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};
\[Omega]insvec =(\[Omega]^((#+6)/3)&/@Range[3])*\[Omega]^(-7/6);


Clear["\[Rho]*"]

{\[Rho]1[\[Eta]_, \[Chi]s_, \[Chi]a_], \[Rho]2[\[Eta]_, \[Chi]s_, \[Chi]a_], \[Rho]3[\[Eta]_, \[Chi]s_, \[Chi]a_]}=  insvec//Simplify;


Expr = HoldComplete[
	{{\[ScriptCapitalA]0, \[ScriptCapitalA]2, \[ScriptCapitalA]3, \[ScriptCapitalA]4, \[ScriptCapitalA]5, \[ScriptCapitalA]6, \[Rho]1, \[Rho]2, \[Rho]3}}, 
	
	\[ScriptCapitalA]0[\[Eta]_] := a0;
	\[ScriptCapitalA]2[\[Eta]_] := a2;
	\[ScriptCapitalA]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := a3;
	\[ScriptCapitalA]4[\[Eta]_, \[Chi]s_, \[Chi]a_] := a4;
	\[ScriptCapitalA]5[\[Eta]_, \[Chi]s_, \[Chi]a_] :=a5;
	\[ScriptCapitalA]6[\[Eta]_, \[Chi]s_, \[Chi]a_] := a6;
	
	\[Rho]1[\[Eta]_, \[Chi]s_, \[Chi]a_] :=  r1;
	\[Rho]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := r2;
	\[Rho]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := r3;
	
	(
		\[ScriptCapitalA]0[\[Eta]]/\[Omega]^(7/6)+(\[Pi]^(2/3) \[ScriptCapitalA]2[\[Eta]])/Sqrt[\[Omega]]+(\[Pi] \[ScriptCapitalA]3[\[Eta],\[Chi]s,\[Chi]a])/\[Omega]^(1/6)+\[Pi]^(4/3) \[Omega]^(1/6) \[ScriptCapitalA]4[\[Eta],\[Chi]s,\[Chi]a]+
		\[Pi]^(5/3) Sqrt[\[Omega]] \[ScriptCapitalA]5[\[Eta],\[Chi]s,\[Chi]a]+\[Pi]^2 \[Omega]^(5/6) \[ScriptCapitalA]6[\[Eta],\[Chi]s,\[Chi]a]
	) + (\[Omega]^(7/6) \[Rho]1[\[Eta],\[Chi]s,\[Chi]a]+\[Omega]^(3/2) \[Rho]2[\[Eta],\[Chi]s,\[Chi]a]+\[Omega]^(11/6) \[Rho]3[\[Eta],\[Chi]s,\[Chi]a])
]//.{
	a0-> \[ScriptCapitalA]0[\[Eta]], a2-> \[ScriptCapitalA]2[\[Eta]], a3-> \[ScriptCapitalA]3[\[Eta], \[Chi]s, \[Chi]a], a4-> \[ScriptCapitalA]4[\[Eta], \[Chi]s, \[Chi]a], a5-> \[ScriptCapitalA]5[\[Eta], \[Chi]s, \[Chi]a],
	a6-> \[ScriptCapitalA]6[\[Eta], \[Chi]s, \[Chi]a], r1 ->\[Rho]1[\[Eta], \[Chi]s, \[Chi]a], r2 ->\[Rho]2[\[Eta], \[Chi]s, \[Chi]a], r3->\[Rho]3[\[Eta], \[Chi]s, \[Chi]a]
};



Expr = $Block@@@(HoldForm[Evaluate@Expr]);


D\[ScriptCapitalA]Ins//ClearAll

expr = Hold[
	{D\[ScriptCapitalA]Ins, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 5},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


resInspiral = EchoTiming[DerivativeRules@@expr];


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = resInspiral[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[0.1, 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@resInspiral}]//RepeatedTiming


(* ::Text:: *)
(*NOW we COMPILE TO MMA Virtual machine*)


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a}, 
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

compileThis[x_Integer] := DZeroFunction


cresInspiral = MapAt[
	compileThis,
	resInspiral,
	{All,2}
];


Do[
	cresInspiral[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[cresInspiral[[i, 1]][[All,1]]]},
	{i, Length@cresInspiral}
] (*The error comes from derivatives that are zero.*)


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := cresInspiral[[i,2]][{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[0.1, 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@resInspiral}]//RepeatedTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = cresInspiral[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@D\[ScriptCapitalA]Ins = list;
]


(* ::Text:: *)
(*Set DownValues:*)


D\[ScriptCapitalA]Ins[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresInspiral[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


(* ::Section::Closed:: *)
(*Ringdown and Damping -Phase*)


(* ::Subsection::Closed:: *)
(*Commom to PhenomD*)


aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] = Block[
	{S, \[Delta]=Sqrt[1- 4 \[Eta]]},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S + 2 Sqrt[3.] \[Eta] + (-0.085` S+0.102` S^2-1.355` S^3-0.868` S^4) \[Eta]-4.399` \[Eta]^2+(-5.837` S-2.097` S^2+4.109` S^3+2.064` S^4) \[Eta]^2+9.397` \[Eta]^3-13.181` \[Eta]^4
]//Simplify;

Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] = Block[
	{EradNS=0.0559745 \[Eta]+0.580951 \[Eta]^2-0.960673 \[Eta]^3+3.35241 \[Eta]^4, S, \[Delta]=Sqrt[1- 4 \[Eta]]},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S = S/(1-2 \[Eta]);
	
	(EradNS (1+(-0.00303023` -2.00661` \[Eta]+7.70506` \[Eta]^2) S))/(1+ (-0.67144` -1.47569` \[Eta] + 7.30468` \[Eta]^2) S)//.\[Eta]-> (1-\[Delta]^2)/4
]//Simplify;

(*re\[Omega] -> interpolation for ringdown and the other is for damping: I think this is frm the PhenomHM paper*)
re\[Omega][\[Chi]_] :=(0.05947169566573468` -0.14989771215394762` \[Chi]+0.09535606290986028` \[Chi]^2+0.02260924869042963` \[Chi]^3-0.02501704155363241` \[Chi]^4-0.005852438240997211` \[Chi]^5+0.0027489038393367993` \[Chi]^6+0.0005821983163192694` \[Chi]^7)/(1-2.8570126619966296` \[Chi]+2.373335413978394` \[Chi]^2-0.6036964688511505` \[Chi]^4+0.0873798215084077` \[Chi]^6);
im\[Omega][\[Chi]_] :=(0.014158792290965177` -0.036989395871554566` \[Chi]+0.026822526296575368` \[Chi]^2+0.0008490933750566702` \[Chi]^3-0.004843996907020524` \[Chi]^4-0.00014745235759327472` \[Chi]^5+0.0001504546201236794` \[Chi]^6)/(1-2.5900842798681376` \[Chi]+1.8952576220623967` \[Chi]^2-0.31416610693042507` \[Chi]^4+0.009002719412204133` \[Chi]^6);

\[Omega]RdDamping[int_, Erad_] = int/(1 - Erad);


\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[6]]//.{\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};
\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[7]]//.{\[Chi]PN-> Sqrt[1 -4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};


\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] = If[\[Gamma]2<=1, \[Omega]RD +( \[Omega]DAMP \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2, \[Omega]RD +( - \[Omega]DAMP \[Gamma]3 )/\[Gamma]2 ];


(* ::Section::Closed:: *)
(*MR Amplitude*)


Clear[\[Gamma]1, \[Gamma]2, \[Gamma]3]
\[Gamma]1[\[Eta]_, \[Chi]s_, \[Chi]a_] = Block[
	{dummy = PhenomDTableV[[5]], res},
	
	res = PhenomCoeff[\[Eta], \[Chi]PN, dummy]//.{\[Chi]PN -> Sqrt[1 - 4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};
	
	res Sqrt[\[Eta]](*//.\[Eta]->(1 - \[Delta]^2)/4*)
]; 


\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[6]]//.{\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};
\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[7]]//.{\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};


Clear@Erad
Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] = Block[
	{EradNS=0.0559745 \[Eta]+0.580951 \[Eta]^2-0.960673 \[Eta]^3+3.35241 \[Eta]^4, S, \[Delta]=Sqrt[1 - 4 \[Eta]]},
	
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S = S/(1-2 \[Eta]);
	
	(EradNS (1+(-0.00303023` -2.00661` \[Eta]+7.70506` \[Eta]^2) S))/(1+ (-0.67144` -1.47569` \[Eta] + 7.30468` \[Eta]^2) S)(*//.\[Eta]-> (1-\[Delta]^2)/4*)
]//Simplify;


Clear@aeff
aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] = Block[
	{S, \[Delta]=Sqrt[1-4 \[Eta]]},
	(*S = 1/4 (1+\[Delta])^2 \[Chi]1 + 1/4 (1-\[Delta])^2 \[Chi]2;*)
	S = \[Chi]s/2 + \[Chi]a \[Delta] + \[Delta]^2 \[Chi]s/2;
	S + 2 Sqrt[3.] \[Eta] + (-0.085` S+0.102` S^2-1.355` S^3-0.868` S^4) \[Eta]-4.399` \[Eta]^2+(-5.837` S-2.097` S^2+4.109` S^3+2.064` S^4) \[Eta]^2+9.397` \[Eta]^3-13.181` \[Eta]^4(*//.\[Eta]-> (1-\[Delta]^2)/4*)
]//Simplify;


Expr = HoldComplete[
	{{\[Omega]RdDamping, re\[Omega], im\[Omega],  Erad, \[Gamma]1, \[Gamma]2, \[Gamma]3, aeff}},

	\[Gamma]1[\[Eta]_, \[Chi]s_, \[Chi]a_] := g1;
	\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := g2;
	\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := g3;
	
	Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
	aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
	
	re\[Omega][\[Kappa]_] := re;
	im\[Omega][\[Kappa]_] := im;
	\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
	
	
	(dampingFrequency E^(-((\[Gamma]2[\[Eta], \[Chi]s, \[Chi]a] (-ringdownFrequency+\[Omega]))/(dampingFrequency \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]))) \[Gamma]1[\[Eta], \[Chi]s, \[Chi]a] \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a])/(\[Omega]^(7/6) (dampingFrequency^2 \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]^2+(-ringdownFrequency+\[Omega])^2))
]//.{
	g1->\[Gamma]1[\[Eta], \[Chi]s, \[Chi]a], g2->\[Gamma]2[\[Eta], \[Chi]s, \[Chi]a], g3->\[Gamma]3[\[Eta], \[Chi]s, \[Chi]a],
	erad->Erad[\[Eta], \[Chi]s, \[Chi]a], ae->aeff[\[Eta], \[Chi]s, \[Chi]a], kappa->\[Kappa][aeff, l, m],
	re->re\[Omega][\[Kappa]], im->im\[Omega][\[Kappa]],
	ringdownFrequency :> \[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	dampingFrequency :> \[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]]
};


Expr = $Block@@@(HoldForm[Evaluate@Expr]);


D\[ScriptCapitalA]MR//ClearAll

expr = Hold[
	{D\[ScriptCapitalA]MR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 5},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


resMR = EchoTiming[DerivativeRules@@expr];


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = resMR[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[0.1, 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@resMR}]//RepeatedTiming


(* ::Text:: *)
(*Now you compile:*)


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a}, 
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


cresMR = MapAt[
	compileThis,
	resMR,
	{All,2}
];


Do[
	cresMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[cresMR[[i, 1]][[All,1]]]},
	{i, Length@cresMR}
] 


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := cresMR[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[0.0001, 0.0001 1024, 0.0001], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@cresMR}]//RepeatedTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = cresMR[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	UpValues@D\[ScriptCapitalA]MR = {};
	UpValues@D\[ScriptCapitalA]MR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


D\[ScriptCapitalA]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


(* ::Section::Closed:: *)
(*Intermediate*)


Block[
	{\[Omega]},
	IntAmpfrequency = \[Omega]^(-7/6) Sqrt[\[Eta]] {1, \[Omega], \[Omega]^2, \[Omega]^3, \[Omega]^4}//Simplify;
]


With[{
	point1 = IntAmpfrequency//.\[Omega]->\[Omega]1,
	point2 = IntAmpfrequency//.\[Omega]->\[Omega]2,
	point3 = IntAmpfrequency//.\[Omega]->\[Omega]3,
	point1D = D[IntAmpfrequency, \[Omega]]//.\[Omega]->\[Omega]1,
	point3D = D[IntAmpfrequency, \[Omega]]//.\[Omega]->\[Omega]3
	},
	
	sol = LinearSolve[{point1, point2, point3, point1D, point3D}, {v1, \[Omega]2^(-7/6) v2, v3,d1,d3}]//Simplify;
	sol = (sol//.{\[Omega]1 -> 0.014})//Simplify;
]


Clear["\[Delta]*"]


{
	\[Delta]0[\[Eta]_,  f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]1[\[Eta]_,  f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]2[\[Eta]_,  f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]3[\[Eta]_,  f2_, f3_, v1_, v2_, v3_, d1_, d3_],
	\[Delta]4[\[Eta]_,  f2_, f3_, v1_, v2_, v3_, d1_, d3_]
} = sol;


Clear@v2

v2[\[Eta]_, \[Chi]s_, \[Chi]a_] = With[
	{
		dummy =(PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[4]])//.{\[Chi]PN-> Sqrt[1 - 4 \[Eta]] \[Chi]a +(1-76 \[Eta]/113) \[Chi]s}
	},
	dummy Sqrt[\[Eta]]//Simplify
];


ClearAll[\[CapitalDelta]0, \[CapitalDelta]1, \[CapitalDelta]2, \[CapitalDelta]3, \[CapitalDelta]4]

{
		\[CapitalDelta]0[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]1[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]2[ \[Omega]2_,  \[Omega]3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]3[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_],
		\[CapitalDelta]4[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_]
	} = { 
		\[Delta]0[\[Eta], \[Omega]2,\[Omega]3,v1,v2,v3,d1,d3],
		\[Delta]1[\[Eta], \[Omega]2,\[Omega]3,v1,v2,v3,d1,d3],
		\[Delta]2[\[Eta], \[Omega]2,\[Omega]3,v1,v2,v3,d1,d3],
		\[Delta]3[\[Eta], \[Omega]2,\[Omega]3,v1,v2,v3,d1,d3],
		\[Delta]4[\[Eta], \[Omega]2,\[Omega]3,v1,v2,v3,d1,d3]
		(*\[Delta]i ~ 1/Sqrt[\[Eta]] so \[CapitalDelta]i is not dependent on \[Eta]*)
		}*Sqrt[\[Eta]];
	


Expr = HoldComplete[
	{{
		\[Omega]RdDamping, re\[Omega],  aeff, Erad, im\[Omega], \[Gamma]2, \[Gamma]3, v2, \[Omega]Peak, \[ScriptA]Ins, \[ScriptA]MR, 
		\[CapitalDelta]0, \[CapitalDelta]1, \[CapitalDelta]2, \[CapitalDelta]3, \[CapitalDelta]4
	}},
	
	\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
	re\[Omega][\[Kappa]_] := (0.05947169566573468` -0.14989771215394762` \[Kappa]+0.09535606290986028` \[Kappa]^2+0.02260924869042963` \[Kappa]^3-0.02501704155363241` \[Kappa]^4-0.005852438240997211` \[Kappa]^5+0.0027489038393367993` \[Kappa]^6+0.0005821983163192694` \[Kappa]^7)/(1-2.8570126619966296` \[Kappa]+2.373335413978394` \[Kappa]^2-0.6036964688511505` \[Kappa]^4+0.0873798215084077` \[Kappa]^6);
	im\[Omega][\[Kappa]_] := (0.014158792290965177` -0.036989395871554566` \[Kappa]+0.026822526296575368` \[Kappa]^2+0.0008490933750566702` \[Kappa]^3-0.004843996907020524` \[Kappa]^4-0.00014745235759327472` \[Kappa]^5+0.0001504546201236794` \[Kappa]^6)/(1-2.5900842798681376` \[Kappa]+1.8952576220623967` \[Kappa]^2-0.31416610693042507` \[Kappa]^4+0.009002719412204133` \[Kappa]^6);
	
	
	aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := aef;
	Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
	
	\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := g2;
	\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := g3;
	v2[\[Eta]_, \[Chi]s_, \[Chi]a_] := ve2base;
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] := omegaPbase;
	
	\[ScriptA]Ins[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[ScriptCapitalA]Ins[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	\[ScriptA]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[ScriptCapitalA]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	\[CapitalDelta]0[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_] := delta0;
	\[CapitalDelta]1[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_] := delta1;
	\[CapitalDelta]2[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_] := delta2;
	\[CapitalDelta]3[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_] := delta3;
	\[CapitalDelta]4[ \[Omega]2_, \[Omega]3_, v1_, v2_, v3_, d1_, d3_] := delta4;
	
	
	
	(
		\[CapitalDelta]0[\[CapitalOmega]2,\[CapitalOmega]3,V1,V2,V3,D1,D3]/\[Omega]^(7/6)+\[CapitalDelta]1[\[CapitalOmega]2,\[CapitalOmega]3,V1,V2,V3,D1,D3]/\[Omega]^(1/6)+ \[Omega]^(5/6) \[CapitalDelta]2[\[CapitalOmega]2,\[CapitalOmega]3,V1,V2,V3,D1,D3]+ 
		\[Omega]^(11/6) \[CapitalDelta]3[\[CapitalOmega]2,\[CapitalOmega]3,V1,V2,V3,D1,D3]+\[Omega]^(17/6) \[CapitalDelta]4[\[CapitalOmega]2,\[CapitalOmega]3,V1,V2,V3,D1,D3]
	)
]//.{
	\[CapitalOmega]2 :> (omegaP+0.014)/2,
	\[CapitalOmega]3:> omegaP,
	V1-> \[ScriptA]Ins[0.014, \[Eta], \[Chi]s, \[Chi]a],
	V2->ve2,
	V3->\[ScriptA]MR[omegaP, \[Eta], \[Chi]s, \[Chi]a],
	D1->Derivative[1, 0, 0, 0][\[ScriptA]Ins][0.014, \[Eta], \[Chi]s, \[Chi]a],
	D3->Derivative[1, 0, 0, 0][\[ScriptA]MR][omegaP, \[Eta], \[Chi]s, \[Chi]a]
};


Expr = Expr//.{
	 delta0 -> \[CapitalDelta]0[\[Omega]2, \[Omega]3, v1, v2, v3, d1, d3],
	 delta1 -> \[CapitalDelta]1[\[Omega]2, \[Omega]3, v1, v2, v3, d1, d3],
	 delta2 -> \[CapitalDelta]2[\[Omega]2, \[Omega]3, v1, v2, v3, d1, d3],
	 delta3 -> \[CapitalDelta]3[\[Omega]2, \[Omega]3, v1, v2, v3, d1, d3],
	 delta4 -> \[CapitalDelta]4[\[Omega]2, \[Omega]3, v1, v2, v3, d1, d3],
	 
	 aef ->aeff[\[Eta], \[Chi]s, \[Chi]a],
	 erad ->Erad[\[Eta],\[Chi]s, \[Chi]a],
	 g2-> \[Gamma]2[\[Eta], \[Chi]s, \[Chi]a],
	 g3-> \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a],
	 ve2base-> v2[\[Eta], \[Chi]s, \[Chi]a],
	 ve2 :> v2[\[Eta], \[Chi]s, \[Chi]a],
	 omegaPbase->\[Omega]Peak[\[Omega]RD, \[Omega]DAMP, \[Gamma]2, \[Gamma]3],
	 omegaP :> \[Omega]Peak[
		\[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
		\[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]], 
		\[Gamma]2[\[Eta], \[Chi]s, \[Chi]a],  \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]
	]
};


Expr = $Block@@@(HoldForm[Evaluate@Expr]);


D\[ScriptCapitalA]Int//ClearAll

expr = Hold[
	{D\[ScriptCapitalA]Int, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 4},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


resInt = EchoTiming[DerivativeRules@@expr];


Module[
	{rule, temp},
	
	
	rule = {
		Derivative[n__][x_] -> $D[{n}, x], 
		$D[{n__}, x_][y1_, y2__]/;Length[y1]==0 :> Last[$D[{n}, x][{y1}, y2]],
		U\[ScriptCapitalA]Ins[y1_, y2__]/;Length[y1]==0 :> Last[U\[ScriptCapitalA]Ins[{y1}, y2]],
		U\[ScriptCapitalA]MR[y1_, y2__]/;Length[y1]==0 :> Last[U\[ScriptCapitalA]MR[{y1}, y2]]
	
	
	(*, U\[ScriptCapitalA]Ins->\[ScriptCapitalA]Ins, U\[ScriptCapitalA]MR->\[ScriptCapitalA]MR*)};
	
	
	(*resInt2 = resInt/.Join[listIns, listMR];*)
	temp = resInt//.rule;
	
	resInt2 = temp//.{U\[ScriptCapitalA]Ins->D\[ScriptCapitalA]Ins,U\[ScriptCapitalA]MR->D\[ScriptCapitalA]MR};
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = resInt2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[0.0001, 0.0001 1024, 0.0001], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@resInt2}]//AbsoluteTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a}, 
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


cresInt = MapAt[
	compileThis,
	resInt2,
	{All,2}
];


Do[
	cresInt[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[cresInt[[i, 1]][[All,1]]]},
	{i, Length@cresInt}
] 


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := cresInt[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[1. 10^-4, 1. 10^-4 1024, 1. 10^-4], 0.2, 0.1, 0.2]
]


Table[T[i], {i,Length@cresInt}]//RepeatedTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = cresInt[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@D\[ScriptCapitalA]Int = list;
]


(* ::Text:: *)
(*Set DownValues:*)


D\[ScriptCapitalA]Int[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresInt[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


(* ::Section::Closed:: *)
(*\[ScriptCapitalA]IMR*)


Expr = HoldComplete[
	{
		{\[Omega]RdDamping, re\[Omega], aeff, Erad, im\[Omega], \[Gamma]2, \[Gamma]3, \[Omega]Peak, \[ScriptA]Ins, \[ScriptA]MR, \[ScriptA]Int}
	},
	
	\[Omega]RdDamping[int_, Erad_] := int/((1-Erad));
	re\[Omega][\[Kappa]_] := (0.05947169566573468` -0.14989771215394762` \[Kappa]+0.09535606290986028` \[Kappa]^2+0.02260924869042963` \[Kappa]^3-0.02501704155363241` \[Kappa]^4-0.005852438240997211` \[Kappa]^5+0.0027489038393367993` \[Kappa]^6+0.0005821983163192694` \[Kappa]^7)/(1-2.8570126619966296` \[Kappa]+2.373335413978394` \[Kappa]^2-0.6036964688511505` \[Kappa]^4+0.0873798215084077` \[Kappa]^6);
	im\[Omega][\[Kappa]_] := (0.014158792290965177` -0.036989395871554566` \[Kappa]+0.026822526296575368` \[Kappa]^2+0.0008490933750566702` \[Kappa]^3-0.004843996907020524` \[Kappa]^4-0.00014745235759327472` \[Kappa]^5+0.0001504546201236794` \[Kappa]^6)/(1-2.5900842798681376` \[Kappa]+1.8952576220623967` \[Kappa]^2-0.31416610693042507` \[Kappa]^4+0.009002719412204133` \[Kappa]^6);
	
	
	aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := aef;
	Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
	
	\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := g2;
	\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := g3;
	\[Omega]Peak[\[Omega]RD_, \[Omega]DAMP_, \[Gamma]2_, \[Gamma]3_] := omegaPbase;
	
	\[ScriptA]Ins[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[ScriptCapitalA]Ins[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	\[ScriptA]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[ScriptCapitalA]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	\[ScriptA]Int[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[ScriptCapitalA]Int[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	Sqrt[5/6] c/(2 \[Pi]^(2/3)) G^2 (
			\[ScriptA]Ins[\[Omega], \[Eta], \[Chi]s, \[Chi]a] us\[Theta][0.014 - \[Omega]] + 
			\[ScriptA]Int[\[Omega], \[Eta], \[Chi]s, \[Chi]a] us\[Theta][(\[Omega] - 0.014) (omegaP - \[Omega])] +
			\[ScriptA]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a] us\[Theta][(0.2 - \[Omega]) (\[Omega]- omegaP)]
	)*{(1+Cos[\[Iota]]^2)/2, -I Cos[\[Iota]]}
]//.{
	omegaP :> \[Omega]Peak[
		\[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
		\[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]], 
		\[Gamma]2[\[Eta], \[Chi]s, \[Chi]a],  \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]
	],
	aef ->aeff[\[Eta], \[Chi]s, \[Chi]a],
	erad ->Erad[\[Eta],\[Chi]s, \[Chi]a],
	g2-> \[Gamma]2[\[Eta], \[Chi]s, \[Chi]a],
	g3-> \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a],
	omegaPbase->\[Omega]Peak[\[Omega]RD, \[Omega]DAMP, \[Gamma]2, \[Gamma]3],
	c-> UnitConvert["SpeedOfLight", ("Gigaparsecs")/("Seconds")][[1]],
	G->4.925490947641267`*^-6
};


Expr = $Block@@@(HoldForm[Evaluate@Expr]);


D\[ScriptCapitalA]IMR//Clear

expr = Hold[
	{D\[ScriptCapitalA]IMR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]}, 4},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


resIMR = EchoTiming[DerivativeRules@@expr];


resIMR2 = resIMR;

Do[
	
	resIMR2[[i,2]]  = resIMR[[i,2]]/.List[x1_, x2_] :> Transpose[List[x1, x2]],
	
	{i, Length@resIMR}
]


Module[
	{rule, temp},
	
	rule = {
		Derivative[n__][x_] -> $D[{n}, x], 
		U\[ScriptCapitalA]Ins->D\[ScriptCapitalA]Ins,U\[ScriptCapitalA]MR->D\[ScriptCapitalA]MR, U\[ScriptCapitalA]Int->D\[ScriptCapitalA]Int,
		us\[Theta]->UnitStep
	};
	
	
	
	resIMR2 = resIMR2//.rule;
]


T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = resIMR2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, \[Iota]$->\[Iota]};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1,0.2, 1.]
]


Table[T[i], {i, Length@resIMR}]//RepeatedTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]}, 
		Evaluate[x],
		RuntimeOptions->{
			"CatchMachineOverflow"->False,
			"CatchMachineIntegerOverflow"->False,
			"EvaluateSymbolically"->False,
			"RuntimeErrorHandler"->None,
			"WarningMessages"->True
		},
		RuntimeAttributes->{Listable}
	]//.{HoldForm[y_] :> y, us\[Theta]->UnitStep, Rational[any_, any2_] :> Divide[any, any2], Complex[x1_, x2_] :> x1+I x2};
	
	Compile@@dummy
]


cresIMR = MapAt[
	compileThis,
	resIMR2,
	{All,2}
];


Do[
	cresIMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]}, X]//.{X -> ToString[cresIMR[[i, 1]][[All,1]]]},
	{i, Length@cresInt}
] 


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] := cresIMR[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, \[Iota]$->\[Iota]};
	
	test[Range[1. 10^-4, 1. 10^-4 1024, 1. 10^-4], 0.2, 0.1, 0.2, 0]
]


Table[T[i], {i, Length@resIMR}]//RepeatedTiming


Module[
	{list = cresIMR[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	UpValues@D\[ScriptCapitalA]IMR = {};
	UpValues@D\[ScriptCapitalA]IMR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


DownValues@D\[ScriptCapitalA]IMR = {};
D\[ScriptCapitalA]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = cresIMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]];


(* ::Section::Closed:: *)
(*Testing*)


(* ::Subsubsection::Closed:: *)
(*lal hphc function*)


DeleteObject/@ExternalSessions[]

Clear@python

python = StartExternalSession[{
	"Python",
	"Evaluator" -> "/Users/felipe/anaconda3/envs/TIGER/bin/python"
}];

ExternalEvaluate[python,"
import lalsimulation
import bilby

import warnings
import numpy as np
warnings.filterwarnings(\"ignore\", \"Wswiglal-redir-stdio\")
import bilby_tgr
"]


ExternalEvaluate[python,"

# Set the duration and sampling frequency of the data segment that we're
# going to inject the signal into
duration = 4.0
sampling_frequency = 2048.0
minimum_frequency = 5


# Fixed arguments passed into the source model
waveform_arguments = dict(
    waveform_approximant=\"IMRPhenomD\",
    reference_frequency=5,
    minimum_frequency=5,
)

WG = bilby.gw.WaveformGenerator(
    duration=4,
    sampling_frequency=sampling_frequency,
    frequency_domain_source_model=bilby_tgr.tiger.source.lal_binary_black_hole_TIGER_PhenomP,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments,
)
"]


lalAmp = ExternalFunction[python,"
def lal_phase(m1, m2, chi1, chi2, iota):

    f = np.linspace(5, 1024, 4077)

    params = {
        \"mass_1\": m1,
        \"mass_2\": m2,
        \"chi_1\": chi1,
        \"chi_2\": chi2,
        \"luminosity_distance\": 1.0,
        \"theta_jn\": iota,
        \"phase\": 0,

        \"dchi_0\":0,
        \"dchi_1\":0,
        \"dchi_2\": 0,
        \"dchi_3\": 0,
        \"dchi_4\": 0,
        \"dchi_5l\": 0,
        \"dchi_6\": 0,
        \"dchi_6l\": 0,
        \"dchi_7\":0,

        \"dbeta_2\": 0,
        \"dbeta_3\": 0,

        \"dalpha_2\": 0,
        \"dalpha_3\": 0,
        \"dalpha_4\": 0,
        \"dalpha_5\": 0,
    }

    hs = WG.frequency_domain_strain(params)

    hp_bilby = hs[\"plus\"][20:]
    hc_bilby = hs[\"cross\"][20:]

    amp = np.abs(hc_bilby)

    return amp
"]


(* ::Subsubsection::Closed:: *)
(*Test*)


(* ::Text:: *)
(*I believe Josiel is using frd and fdamp as calculated in the original PhenomD implementation. I am not convinced this is right. I find it more likely that those should be replaced by the frd and fdamp 22 formulas attached to Z22.*)


DownValues@D\[ScriptCapitalA]IMR


Clear@Test

Test[M_] := Module[
	{
		m1, m2, mc, eta, s1z, s2z, \[Iota], \[Phi]Ref, lalRes, fmax, fmin, myDef, amp, t0, chis, chia, psi,  
		G =  4.925490947641267`*^-6
	},
	
	\[Iota] = RandomReal[{0, \[Pi]}];
	eta = RandomReal[{0.1, 0.25}];
	
	m1 = M/2 (1 + Sqrt[1-4 eta]);
	m2 = M/2 (1 - Sqrt[1-4 eta]);
	
	
	mc = (m1+m2) eta^(3/5);
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	
	fmax =  1024;
	fmin=5;
	
	lalRes = lalAmp[m1, m2, s1z, s2z, \[Iota]]//Normal;
	
	myDef = M^2/10^-3 D\[ScriptCapitalA]IMR[
		M Range[5., 1024, 0.25] G, 
		eta,
		(s1z+s2z)/2, (s1z-s2z)/2,
		\[Iota]
	][[All,2]]//Abs;
		
	{lalRes, myDef}
	
]


Module[
	{M = RandomReal[{50,200}], list, G =  4.925490947641267`*^-6},
	
	list = Test[M];
	
	ListLinePlot[
		RelativeDiff@@list,
		ScalingFunctions->"Log10",
		PlotRange->{{5, 0.2/(G M)}, All},
		DataRange->{5, 1024},
		PlotLabel->"hp"
	]
	
]


(* ::Subsection::Closed:: *)
(*Comparing numerical x Symbolic derivatives*)


(* ::Subsubsection::Closed:: *)
(*First Order*)


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


Test\[ScriptCapitalA][\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_,\[Iota]_] := D\[ScriptCapitalA]IMR[{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]]


Block[
	{}, 
	
	ClearAll[TestGrad\[ScriptCapitalA]];
	
	TestGrad\[ScriptCapitalA][\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = resIMR2[[2;;6, 2]];

]
DownValues[TestGrad\[ScriptCapitalA]] = DownValues[TestGrad\[ScriptCapitalA]]//.HoldForm[x_]:> x;


Clear@Test
Test := Module[
	{
		\[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		M, \[Iota]=RandomReal[{0, \[Pi]}]
	},
	
	M = RandomReal[{50,200}];
	
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1}, 2]; 
	
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	Symbolic = TestGrad\[ScriptCapitalA][{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]];
	vars = { \[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]};
	Numeric = NGrad[Test\[ScriptCapitalA], vars, 0];
	
	RelativeDiff@@{Symbolic, Numeric}
	

	
]


Test//ScientificForm


Table[Test, {1000}]//MinMax


(* ::Subsubsection::Closed:: *)
(*Second Order*)


(* ::Text:: *)
(*Since First order checks out we can Use the symbolic gradient of first order and compare its numerical derivatives against the symbolic *)
(*gradients of order 2.*)


Block[
	{},
	
	ClearAll[iTestGrad\[CapitalPsi]O1];
	
	iTestGrad\[CapitalPsi]O1[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = resIMR2[[2;;6, 2]];

]
DownValues[iTestGrad\[CapitalPsi]O1] = DownValues[iTestGrad\[CapitalPsi]O1]//.HoldForm[x_]:> x;

TestGrad\[CapitalPsi]O1//Clear
TestGrad\[CapitalPsi]O1[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] := iTestGrad\[CapitalPsi]O1[{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]]


Range[7,22]//Length


Block[{rule, g = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]]}, 
	
	ClearAll[TestGrad\[CapitalPsi]O2];
	
	TestGrad\[CapitalPsi]O2[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_
		] = resIMR2[[7;;21, 2]];

]
DownValues[TestGrad\[CapitalPsi]O2] = DownValues[TestGrad\[CapitalPsi]O2]//.HoldForm[x_]:> x;


Clear@Test

Test := Module[
	{
		\[Delta], \[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		li = SymmetrizedIndependentComponents[{5,5}, Symmetric[All]],M, \[Iota] = RandomReal[{0, \[Pi]}]
	},
	
	M = RandomReal[{50, 200}];
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1},2];
	{\[Chi]s,\[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	Symbolic = ArrayReshape[TestGrad\[CapitalPsi]O2[{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]], {15,2}];
	
	vars = {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]};
	
	Numeric = Extract[
		ArrayReshape[NGrad[TestGrad\[CapitalPsi]O1, vars, 0], {5, 5, 2}],
		li
	];
	
	RelativeDiff@@{Symbolic, Numeric}
	

	
]


Test//ScientificForm


Table[Test, {10^3}]//MinMax


(* ::Subsubsection::Closed:: *)
(*Third Order*)


(* ::Text:: *)
(*Second order checks out so we can use the symbolic gradient of second order and compare its numerical derivatives against the symbolic gradients of order 3.*)


Block[
	{}, 
	
	ClearAll[iTestGrad\[ScriptCapitalA]O2];
	
	iTestGrad\[ScriptCapitalA]O2[ \[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_
		] = resIMR2[[7;;21, 2]];

]
DownValues[iTestGrad\[ScriptCapitalA]O2] = DownValues[iTestGrad\[ScriptCapitalA]O2]//.HoldForm[x_]:> x;

TestGrad\[ScriptCapitalA]O2//Clear

TestGrad\[ScriptCapitalA]O2[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] := (
	ArrayReshape[iTestGrad\[ScriptCapitalA]O2[{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]], {15, 2}]
);


SymmetrizedIndependentComponents[{5, 5, 2}, Symmetric[{1,2}]]//Length


SymmetricTestGradOrder2[x__] := Module[
	{
		li = SymmetrizedIndependentComponents[{5, 5, 2}, Symmetric[{1,2}]],
		rule, values = TestGrad\[ScriptCapitalA]O2[x]//Flatten, rules
	},
	
	rules = MapThread[Rule, {li, values}];
	
	SymmetrizedArray[rules, {5, 5, 2},  Symmetric[{1,2}]]
]


Block[
	{},
	
	ClearAll[iTestGrad\[ScriptCapitalA]O3];
	
	iTestGrad\[ScriptCapitalA]O3[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = resIMR2[[22;;56, 2]];

]
DownValues[iTestGrad\[ScriptCapitalA]O3] = DownValues[iTestGrad\[ScriptCapitalA]O3]//.HoldForm[x_]:> x;

TestGrad\[ScriptCapitalA]O3//Clear
TestGrad\[ScriptCapitalA]O3[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] := iTestGrad\[ScriptCapitalA]O3[{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]]


(* ::Text:: *)
(*Just make sure derivatives are in the same order of \[CapitalPsi]*)


Length@Range[22,56]


Clear@Test

Test := Module[
	{
		\[Delta], \[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		li = SymmetrizedIndependentComponents[{5, 5, 5, 2}, Symmetric[{1,2,3}]], M, \[Iota] = RandomReal[{0,\[Pi]}]
	},
	
	
	M = RandomReal[{50, 200}];
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	
	{\[Chi]1,\[Chi]2} = RandomReal[{-1,1}, 2];
	{\[Chi]s,\[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	Symbolic = ArrayReshape[TestGrad\[ScriptCapitalA]O3[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]], {35,2}]//Flatten;
	
	vars = {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]};
	
	Numeric = Extract[NGrad[SymmetricTestGradOrder2, vars, 0],li];
	
	RelativeDiff@@{Symbolic, Numeric}
]


Test//ScientificForm


Table[Test, {10^3}]//MinMax//EchoTiming


(* ::Subsubsection::Closed:: *)
(*Fourth Order*)


(* ::Text:: *)
(*Second order checks out so we can use the symbolic gradient of second order and compare its numerical derivatives against the symbolic gradients of order 3.*)


Block[
	{}, 
	
	ClearAll[iTestGrad\[ScriptCapitalA]O3];
	
	iTestGrad\[ScriptCapitalA]O3[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = resIMR2[[22;;56, 2]];

]
DownValues[iTestGrad\[ScriptCapitalA]O3] = DownValues[iTestGrad\[ScriptCapitalA]O3]//.HoldForm[x_]:> x;

TestGrad\[ScriptCapitalA]O3//Clear

TestGrad\[ScriptCapitalA]O3[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] := iTestGrad\[ScriptCapitalA]O3[{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]]//Flatten


SymmetricTestGradOrder3[x__] := Module[
	{
		li = SymmetrizedIndependentComponents[{5, 5, 5, 2}, Symmetric[{1,2,3}]],
		rule, 
		values = TestGrad\[ScriptCapitalA]O3[x]//Flatten, rules
	},
	
	rules = MapThread[Rule, {li, values}];
	
	SymmetrizedArray[rules, {5, 5, 5,2},  Symmetric[{1,2,3}]]
]


Block[
	{},
	
	ClearAll[iTestGrad\[ScriptCapitalA]O4];
	
	iTestGrad\[ScriptCapitalA]O4[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = resIMR2[[57;;-1, 2]];

]
DownValues[iTestGrad\[ScriptCapitalA]O4] = DownValues[iTestGrad\[ScriptCapitalA]O4]//.HoldForm[x_]:> x;

TestGrad\[ScriptCapitalA]O4//Clear

TestGrad\[ScriptCapitalA]O4[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] := iTestGrad\[ScriptCapitalA]O4[{\[Omega]}, \[Eta], \[Chi]s, \[Chi]a, \[Iota]]//Flatten


Clear@Test

Test := Module[
	{
		\[Delta], \[Omega]ref, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		li = SymmetrizedIndependentComponents[{5, 5, 5, 5, 2}, Symmetric[{1,2,3,4}]], M, \[Iota]=RandomReal[{0, \[Pi]}]
	},
	
	M = RandomReal[{50, 200}];
	\[Eta] = RandomReal[{0.125, 0.24999}];
	
	
	{\[Chi]1,\[Chi]2} = RandomReal[{-1,1}, 2];
	{\[Chi]s,\[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	
	
	f = RandomReal[{10., 0.2/(G M)}];
	\[Omega] = M G f;
	\[Omega]ref = 10 M G;
	
	Symbolic = ArrayReshape[TestGrad\[ScriptCapitalA]O4[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]], {70,2}];
	
	vars = {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]};
	
	Numeric = ArrayReshape[
		Extract[NGrad[SymmetricTestGradOrder3, vars, 0], li],
		{70,2}
	];

	RelativeDiff@@{Symbolic, Numeric}
]


Test//ScientificForm


Table[Test, {10^3}]//MinMax//EchoTiming


(* ::Subsection::Closed:: *)
(*Comparing to old implementation*)


test = FelipeBarbosa`SymDALI`DALICoefficients`Private`test;


list = test[FelipeBarbosa`SymDALI`DALICoefficients`Private`\[ScriptA]IMR];


list = list//.{
	FelipeBarbosa`SymDALI`DALICoefficients`Private`\[ScriptA]IMR -> \[ScriptA]IMR
};


TagSetDelayed@@@list[[2;;-1]];


list[[1]]


aIMR[f_, \[ScriptCapitalM]c_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = \[ScriptA]IMR[f, \[ScriptCapitalM]c, Sqrt[1 - 4 \[Eta]], \[Chi]s, \[Chi]a,\[Iota]]


DownValues@d\[ScriptCapitalA]IMR


DownValues@d\[ScriptCapitalA]IMR = {HoldPattern[d\[ScriptCapitalA]IMR[f_,\[ScriptCapitalM]c_,\[Eta]_,\[Chi]s_,\[Chi]a_,\[Iota]_]]:>(\[ScriptCapitalM]c^2 D\[ScriptCapitalA]IMR[(4.9254664969309`3.6105383994801805*^-6 f \[ScriptCapitalM]c)/\[Eta]^(3/5),\[Eta],\[Chi]s,\[Chi]a,\[Iota]])/\[Eta]^(6/5)};


(*d\[ScriptCapitalA]IMR//Clear
d\[ScriptCapitalA]IMR[f_, \[ScriptCapitalM]c_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_] = Module[
	{G = 4.9254664969309`3.6105383994801805*^-6, M, \[Omega]},
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5);
	
	\[Omega] = f M G;
	
	M^2 D\[ScriptCapitalA]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota]]
]*)


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


vars = {\[ScriptCapitalM]c, \[Eta], \[Chi]s, \[Chi]a, \[Iota]};


(*All derivatives up to order 3*)

derivatives = Combinations[vars, 3];

derivatives//Length


Unprotect@Derivative;
Derivative[n__][s_Symbol] := $D[{n}, s];
(*Derivative[n__][D\[ScriptCapitalA]IMR] := $D[{n}, D\[ScriptCapitalA]IMR];*)
Derivative[x__][$D[y_List, s_Symbol]][k__] := $D[{x} + y, s][k];
Protect@Derivative;


DownValues@D\[ScriptCapitalA]IMR = {};


derivatives[[1;;5]]


(* ::Subsubsection::Closed:: *)
(*First order*)


T := Block[
	{
		headList1, headList2, 
		 \[Eta], \[Chi]s, \[Chi]a,
		G = 4.9254664969309`3.6105383994801805*^-6, M, \[Chi]1, \[Chi]2, f, \[Iota], \[ScriptCapitalM]c
	},
	
	
	Block[
		{D\[ScriptCapitalA]IMR, \[ScriptA]IMR},
		headList1 = D[d\[ScriptCapitalA]IMR[f, \[ScriptCapitalM]c, \[Eta], \[Chi]s, \[Chi]a, \[Iota]], #]&/@(derivatives[[1;;5]]);
		headList2 = D[aIMR[f, \[ScriptCapitalM]c, \[Eta],  \[Chi]s, \[Chi]a, \[Iota]], #]&/@(derivatives[[1;;5]]);
	];
	
	
	M = RandomReal[{50, 200}];
	
	f = RandomReal[{6, 0.15/(G M)}]//List;
	\[Eta] = RandomReal[{0.1, 0.24999}];
	\[ScriptCapitalM]c = M \[Eta]^(3/5);
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-1,1},2];
	
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	
	RelativeDiff@@{
		headList1//Flatten,
		(Transpose/@headList2)//Flatten
	}
]


tt = Table[T, 10^5]//EchoTiming;


tt//Max//ScientificForm


ListPlot[Flatten@tt, ScalingFunctions->"Log10"]


(* ::Subsubsection::Closed:: *)
(*Second order*)


derivatives[[6;;20]]


T := Block[
	{
		headList1, headList2, 
		 \[Eta], \[Chi]s, \[Chi]a,
		G = 4.9254664969309`3.6105383994801805*^-6, M, \[Chi]1, \[Chi]2, f, \[Iota], \[ScriptCapitalM]c
	},
	
	
	Block[
		{D\[ScriptCapitalA]IMR, \[ScriptA]IMR},
		headList1 = D[d\[ScriptCapitalA]IMR[f, \[ScriptCapitalM]c, \[Eta], \[Chi]s, \[Chi]a, \[Iota]], Sequence@@#]&/@(derivatives[[6;;20]]);
		headList2 = D[aIMR[f, \[ScriptCapitalM]c, \[Eta],  \[Chi]s, \[Chi]a, \[Iota]], Sequence@@#]&/@(derivatives[[6;;20]]);
	];
	
	
	M = RandomReal[{50, 200}];
	
	f = RandomReal[{6, 0.15/(G M)}]//List;
	\[Eta] = RandomReal[{0.1, 0.24999}];
	\[ScriptCapitalM]c = M \[Eta]^(3/5);
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-0.5,0.5},2];
	
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	
	RelativeDiff@@{
		headList1//Flatten,
		(Transpose/@headList2)//Flatten
	}
	
]


T


tt = Table[T, 10^5]//EchoTiming;


tt//Max//ScientificForm


Count[Flatten@tt, x_/; x<10^-5]/Length[Flatten[tt]]//N


ListPlot[Flatten@tt, ScalingFunctions->"Log10"]


(* ::Subsubsection::Closed:: *)
(*Third order*)


derivatives[[21;;-1]]


T := Block[
	{
		headList1, headList2, 
		 \[Eta], \[Chi]s, \[Chi]a,
		G = 4.9254664969309`3.6105383994801805*^-6, M, \[Chi]1, \[Chi]2, f, \[Iota], \[ScriptCapitalM]c
	},
	
	
	Block[
		{D\[ScriptCapitalA]IMR, \[ScriptA]IMR},
		headList1 = D[d\[ScriptCapitalA]IMR[f, \[ScriptCapitalM]c, \[Eta], \[Chi]s, \[Chi]a, \[Iota]], Sequence@@#]&/@(derivatives[[21;;]]);
		headList2 = D[aIMR[f, \[ScriptCapitalM]c, \[Eta],  \[Chi]s, \[Chi]a, \[Iota]], Sequence@@#]&/@(derivatives[[21;;]]);
	];
	
	
	M = RandomReal[{50, 200}];
	
	f = RandomReal[{6, 0.15/(G M)}]//List;
	
	\[Eta] = RandomReal[{0.1, 0.24999}];
	
	\[ScriptCapitalM]c = M \[Eta]^(3/5);
	
	{\[Chi]1, \[Chi]2} = RandomReal[{-0.5,0.5},2];
	
	{\[Chi]s, \[Chi]a} = {(\[Chi]1+\[Chi]2)/2, (\[Chi]1-\[Chi]2)/2};
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	
	RelativeDiff@@{
		headList1//Flatten,
		(Transpose/@headList2)//Flatten
	}
	
]


T


tt = Table[T, 10^4]//EchoTiming;


tt//Max//ScientificForm


Count[Flatten@tt, x_/; x<10^-5]/Length[Flatten[tt]]//N


ListPlot[Flatten@tt, ScalingFunctions->"Log10"]


(* ::Chapter:: *)
(*Saving defs*)


(* ::Subsection::Closed:: *)
(*New*)


(* ::Text:: *)
(*THIS ADDS SIDE EFFECTS:*)


DownValues[HMhphc]


newVars = {
	\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota],\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7
}//.x_Symbol/; Context[x] ==="Global`" :> ToExpression[
	"FelipeBarbosa`SymDALI`DALICoefficients`Private`"<>ToString[x]
]


newHeads = {
			DZeroFunction, 
			D\[ScriptCapitalA]Ins, D\[ScriptCapitalA]Int, D\[ScriptCapitalA]MR, D\[ScriptCapitalA]IMR,
			D\[CapitalPhi]Ins, D\[CapitalPhi]Int, D\[CapitalPhi]MR, D\[CapitalPhi]IMR, D\[CapitalPsi]
}//.x_Symbol/; Context[x] ==="Global`" :> ToExpression[
	"FelipeBarbosa`SymDALI`DALICoefficients`Private`"<>ToString[x]
]


varRules = MapThread[
	Rule, 
	{
		{\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota],\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7},
		newVars
	}
];


HeadRules = MapThread[
	Rule,
	{{
			DZeroFunction, 
			D\[ScriptCapitalA]Ins, D\[ScriptCapitalA]Int, D\[ScriptCapitalA]MR, D\[ScriptCapitalA]IMR,
			D\[CapitalPhi]Ins, D\[CapitalPhi]Int, D\[CapitalPhi]MR, D\[CapitalPhi]IMR, D\[CapitalPsi]
	}, 
	newHeads}
]


ChangeContext//Clear


ChangeContext[x_Symbol] := Module[
	{
		newSymbol, newContext = "FelipeBarbosa`SymDALI`DALICoefficients`Private`",
		newDownValues, newUpValues
	},
	
	newSymbol = newContext <> ToString[x]//ToExpression;
	
	newDownValues = DownValues[x]//.Join[
		varRules, HeadRules
	];
	
	(Hold[DownValues[ss] = newDownValues]//.ss->newSymbol)//ReleaseHold;
	
	newUpValues = UpValues[x]//.Join[
		varRules, HeadRules
	];
	
	(Hold[UpValues[ss] = newUpValues]//.ss->newSymbol)//ReleaseHold;
	
]


ChangeContext[D\[ScriptCapitalA]Ins]


ChangeContext/@{
			DZeroFunction, 
			D\[ScriptCapitalA]Ins, D\[ScriptCapitalA]Int, D\[ScriptCapitalA]MR, D\[ScriptCapitalA]IMR,
			D\[CapitalPhi]Ins, D\[CapitalPhi]Int, D\[CapitalPhi]MR, D\[CapitalPhi]IMR, D\[CapitalPsi]
}


Module[
	{SymDALIDir = NotebookDirectory[]//ParentDirectory[#, 2]&, fileName, file},
	
	
	fileName = FileNameJoin[{
		SymDALIDir, 
		"LibraryResources",
		$SystemID,
		"DerivativeRules/IMRPhenomD/MMA/Defs.mx"
	}];
	
	file = CreateFile[fileName];
	
	DumpSave[
		file,
		{
			FelipeBarbosa`SymDALI`DALICoefficients`Private`DZeroFunction,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[ScriptCapitalA]Ins,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[ScriptCapitalA]Int,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[ScriptCapitalA]MR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[ScriptCapitalA]IMR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[CapitalPhi]Ins,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[CapitalPhi]Int,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[CapitalPhi]MR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[CapitalPhi]IMR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`D\[CapitalPsi]
		}
	]
	
]


DownValues@FelipeBarbosa`SymDALI`DALICoefficients`Private`DZeroFunction


FelipeBarbosa`SymDALI`DALICoefficients`Private`DZeroFunction//ClearAll


Get["/Users/felipe/Documents/GitHub/SymDALI/LibraryResources/MacOSX-ARM64/DerivativeRules/IMRPhenomD/MMA/Defs.mx"]


DownValues@FelipeBarbosa`SymDALI`DALICoefficients`Private`DZeroFunction


Import["/Users/felipe/Documents/GitHub/SymDALI/LibraryResources/MacOSX-ARM64/DerivativeRules/IMRPhenomD/SymRules/file.wdx"]
