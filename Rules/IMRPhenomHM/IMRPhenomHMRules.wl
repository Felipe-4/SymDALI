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


HMZeroFunction//ClearAll

HMZeroFunction[\[Omega]_?VectorQ, x__] := ConstantArray[0, Length@\[Omega]]
HMZeroFunction[\[Omega]_?NumberQ, x__] := 0


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
		PlotRange->All,
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


HM\[ScriptCapitalA]Ins//ClearAll

expr = Hold[
	{HM\[ScriptCapitalA]Ins, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 4},
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

compileThis[x_Integer] := HMZeroFunction


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
	
	UpValues@HM\[ScriptCapitalA]Ins = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[ScriptCapitalA]Ins[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresInspiral[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


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


(* ::Subsection::Closed:: *)
(*Exclusive to PhenomHM*)


\[Kappa][aeff_, l_, m_] := (Log[2-aeff]/Log[3])^(1/(2+ l- Abs[m]))


Z22[\[Kappa]22_] := (
	1 + 
	\[Kappa]22*1.557847*Exp[2.903124 I] + 
	\[Kappa]22^2*1.95097051*Exp[5.920970 I] + 
	\[Kappa]22^3*2.09971716 * Exp[2.760585 I] + 
	\[Kappa]22^4*1.41094660 * Exp[5.914340 I]+ 
	\[Kappa]22^5*0.41063923*Exp[2.795235 I]
)

Z32[\[Kappa]32_] := (
	1.022464*Exp[0.004870 I] + 
	\[Kappa]32*0.24731213*Exp[0.665292 I] + 
	\[Kappa]32^2*1.70468239*Exp[3.138283 I] + 
	\[Kappa]32^3*0.94604882*Exp[0.163247 I] + 
	\[Kappa]32^4 * 1.53189884*Exp[5.703573 I] + 
	\[Kappa]32^5 * 2.28052668*Exp[2.685231 I] + 
	\[Kappa]32^6 * 0.92150314*Exp[5.841704 I]
)


Z44[\[Kappa]44_] := (
	2 + 
	\[Kappa]44*2.658908*Exp[3.002787 I] + 
	\[Kappa]44^2 * 2.97825567 * Exp[6.050955 I] + 
	\[Kappa]44^3 * 3.21842350 * Exp[2.877514 I] + 
	\[Kappa]44^4 * 2.12764967 * Exp[5.989669 I] + 
	\[Kappa]44^5 * 0.60338186 * Exp[2.830031 I]
)

Z21[\[Kappa]21_] := (
	0.589113 * Exp[0.043525 I] + 
	\[Kappa]21*0.18896353  * Exp[2.289868 I] + 
	\[Kappa]21^2*1.15012965  * Exp[5.810057 I] + 
	\[Kappa]21^3*6.04585476  * Exp[2.741967 I] + 
	\[Kappa]21^4*11.12627777 * Exp[5.844130 I] + 
	\[Kappa]21^5*9.34711461  * Exp[2.669372 I] + 
	\[Kappa]21^6*3.03838318  * Exp[5.791518 I]
)


Z33[\[Kappa]33_] := (
	1.5 +
	\[Kappa]33 * 2.095657   * Exp[2.964973 I] +
	\[Kappa]33^2 * 2.46964352 * Exp[5.996734 I] + 
	\[Kappa]33^3 * 2.66552551 * Exp[2.817591 I] + 
	\[Kappa]33^4 * 1.75836443 * Exp[5.932693 I] + 
	\[Kappa]33^5 * 0.49905688 * Exp[2.781658 I]
)



Z43[\[Kappa]43_] := (
	1.5 + 
	\[Kappa]43 * 0.205046   * Exp[0.595328 I] + 
	\[Kappa]43^2* 3.10333396 * Exp[3.016200 I] + 
	\[Kappa]43^3 * 4.23612166 * Exp[6.038842 I] + 
	\[Kappa]43^4* 3.02890198 * Exp[2.826239 I] + 
	\[Kappa]43^5 * 0.90843949 * Exp[5.915164 I]
)


(* ::Text:: *)
(*This folllow the reasoning of re\[Omega] and im\[Omega] of last subsection:*)


re\[Omega]22[\[Kappa]22_] = Re[Z22[\[Kappa]22]]//Simplify[#, Assumptions-> \[Kappa]22 \[Element] Reals]&;
im\[Omega]22[\[Kappa]22_] = Im[Z22[\[Kappa]22]]//Simplify[#, Assumptions-> \[Kappa]22 \[Element] Reals]&;

re\[Omega]32[\[Kappa]32_] = Re[Z32[\[Kappa]32]]//Simplify[#, Assumptions-> \[Kappa]32 \[Element] Reals]&;
im\[Omega]32[\[Kappa]32_] = Im[Z32[\[Kappa]32]]//Simplify[#, Assumptions-> \[Kappa]32 \[Element] Reals]&;


re\[Omega]44[\[Kappa]44_] = Re[Z44[\[Kappa]44]]//Simplify[#, Assumptions-> \[Kappa]44 \[Element] Reals]&;
im\[Omega]44[\[Kappa]44_] = Im[Z44[\[Kappa]44]]//Simplify[#, Assumptions-> \[Kappa]44 \[Element] Reals]&;

re\[Omega]21[\[Kappa]21_] = Re[Z21[\[Kappa]21]]//Simplify[#, Assumptions-> \[Kappa]21 \[Element] Reals]&;
im\[Omega]21[\[Kappa]21_] = Im[Z21[\[Kappa]21]]//Simplify[#, Assumptions-> \[Kappa]21 \[Element] Reals]&;


re\[Omega]33[\[Kappa]33_] = Re[Z33[\[Kappa]33]]//Simplify[#, Assumptions-> \[Kappa]33 \[Element] Reals]&;
im\[Omega]33[\[Kappa]33_] = Im[Z33[\[Kappa]33]]//Simplify[#, Assumptions-> \[Kappa]33 \[Element] Reals]&;

re\[Omega]43[\[Kappa]43_] = Re[Z43[\[Kappa]43]]//Simplify[#, Assumptions-> \[Kappa]43 \[Element] Reals]&;
im\[Omega]43[\[Kappa]43_] = Im[Z43[\[Kappa]43]]//Simplify[#, Assumptions-> \[Kappa]43 \[Element] Reals]&;


\[Omega]RdDampinglm[int_, Erad_] := int/(2 \[Pi] (1-Erad))


ClearAll@f22\[ScriptCapitalA]

Attributes[f22\[ScriptCapitalA]] = {Listable};

f22\[ScriptCapitalA][\[Omega]_, \[Omega]rd22_, \[Omega]rdlm_, m_] := Module[
	{f0, case1, case2, case3},
	
	f0 = 0.014 \[Omega]rdlm/\[Omega]rd22;
	
	
	case1 =2/m \[Omega];
	
	case2 = (\[Omega]rd22 - 2 f0/m)/(\[Omega]rdlm-f0) (\[Omega]-f0) + 2 f0/m;
	
	case3 = \[Omega]  - (\[Omega]rdlm - \[Omega]rd22);
	
	Which[
		\[Omega] <= f0, case1,
		f0 < \[Omega] <= \[Omega]rdlm, case2,
		\[Omega]rdlm<\[Omega], case3
	]
]


ClearAll@f22\[ScriptCapitalA]

f22\[ScriptCapitalA][\[Omega]_, \[Omega]rd22_, \[Omega]rdlm_, m_] := Module[
	{f0, case1, case2, case3},
	
	f0 = 0.014 \[Omega]rdlm/\[Omega]rd22;
	
	
	case1 =2/m \[Omega];
	
	case2 = (\[Omega]rd22 - 2 f0/m)/(\[Omega]rdlm-f0) (\[Omega]-f0) + 2 f0/m;
	
	case3 = \[Omega]  - (\[Omega]rdlm - \[Omega]rd22);
	
	us\[Theta][f0-\[Omega]] case1 + us\[Theta][(\[Omega]-f0) (\[Omega]rdlm-\[Omega])] case2 + us\[Theta][\[Omega]-\[Omega]rdlm] case3//Simplify
]


f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rdlm, m]


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


HM\[ScriptCapitalA]MR//ClearAll

expr = Hold[
	{HM\[ScriptCapitalA]MR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 4},
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
	
	UpValues@HM\[ScriptCapitalA]MR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[ScriptCapitalA]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


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


HM\[ScriptCapitalA]Int//ClearAll

expr = Hold[
	{HM\[ScriptCapitalA]Int, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 3},
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
	
	
	(*resInt2 = resInt/.Join[listIns, listMR];	*)
	temp = resInt//.rule;
	
	resInt2 = temp//.{U\[ScriptCapitalA]Ins->HM\[ScriptCapitalA]Ins,U\[ScriptCapitalA]MR->HM\[ScriptCapitalA]MR};
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = resInt2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[0.0001, 0.0001 1024, 0.0001], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@resInt2}]//RepeatedTiming


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
	
	UpValues@HM\[ScriptCapitalA]Int = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[ScriptCapitalA]Int[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresInt[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


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
	)
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


HM\[ScriptCapitalA]IMR//Clear

expr = Hold[
	{HM\[ScriptCapitalA]IMR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 3},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


resIMR = EchoTiming[DerivativeRules@@expr];


Module[
	{rule, temp},
	
	rule = {
		Derivative[n__][x_] -> $D[{n}, x], 
		U\[ScriptCapitalA]Ins->HM\[ScriptCapitalA]Ins,U\[ScriptCapitalA]MR->HM\[ScriptCapitalA]MR, U\[ScriptCapitalA]Int->HM\[ScriptCapitalA]Int,
		us\[Theta]->UnitStep
	};
	
	
	
	resIMR2 = resIMR//.rule;
]


T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = resIMR2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1,0.2]
]


Table[T[i], {i, Length@resIMR}]//RepeatedTiming


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
	cresIMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[cresIMR[[i, 1]][[All,1]]]},
	{i, Length@cresInt}
] 


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := cresIMR[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[1. 10^-4, 1. 10^-4 1024, 1. 10^-4], 0.2, 0.1, 0.2]
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
	
	UpValues@HM\[ScriptCapitalA]IMR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[ScriptCapitalA]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresIMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


(* ::Section::Closed:: *)
(*Hlm*)


H21[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, m_] := Module[
	
	{v = (2 \[Pi] \[Omega]/m)^(1/3), h, \[Delta] = Sqrt[1 - 4 \[Eta]], hstar},
	
	h = Sqrt[2]/3 ( 
		v \[Delta] - 3/2 v^2 (\[Chi]a + \[Delta] \[Chi]s) + v^3 \[Delta] (335/672 + 117/56 \[Eta]) + 
		v^4 (
			\[Chi]a (3427/1344 - 2101/336 \[Eta]) + \[Delta] \[Chi]s (3427/1344 - 956/336 \[Eta]) + \[Delta] (-I/2 -\[Pi] - 2 I 0.69314718056)
		)
		
	);
	
	hstar = Sqrt[2]/3 ( 
		v \[Delta] - 3/2 v^2 (\[Chi]a + \[Delta] \[Chi]s) + v^3 \[Delta] (335/672 + 117/56 \[Eta]) + 
		v^4 (
			\[Chi]a (3427/1344 - 2101/336 \[Eta]) + \[Delta] \[Chi]s (3427/1344 - 956/336 \[Eta]) + \[Delta] (I/2 -\[Pi] + 2 I 0.69314718056)
		)
		
	);
	
	Simplify[Sqrt[hstar h], Assumptions->\[Delta]\[Element]Reals&&\[Omega]\[Element]Reals&&\[Chi]s\[Element]Reals && \[Chi]a\[Element]Reals]
];


H21[\[Omega], \[Eta], \[Chi]s, \[Chi]a, m]


Clear@H33
H33[\[Omega]_, \[Eta]_, m_] := Module[
	{v = (2 \[Pi] \[Omega]/m)^(1/3), h, \[Delta]=Sqrt[1  - 4  \[Eta]]},
	
	3/4 Sqrt[5/7] v \[Delta]
]


H33[\[Omega], \[Eta], 3]/H33[2 \[Omega]/3, \[Eta], 3]


Clear@H32
H32[\[Omega]_, \[Eta]_, m_] := Module[
	{v = (2 \[Pi] \[Omega]/m)^(1/3), \[Delta] = Sqrt[1-4 \[Eta]], h},
	
	1/3 Sqrt[5/7] v^2 (1-3 \[Eta])
]


H32[\[Omega], \[Eta], 2]


Clear@H44
H44[\[Omega]_, \[Eta]_,  m_] := Module[
	{v = (2 \[Pi] \[Omega]/m)^(1/3), h},
	
	4/9 Sqrt[10/7] v^2 (1-3 \[Eta])
]


H44[\[Omega], \[Eta], 4]/H44[\[Omega]/2, \[Eta], 4]


Clear@H43
H43[\[Omega]_, \[Eta]_, m_] := Module[
	{v = (2 \[Pi] \[Omega]/m)^(1/3), \[Delta]=Sqrt[1-4 \[Eta]], h},
	
	3/4 Sqrt[3/35] v^3 \[Delta] (1-2 \[Eta])
]


H43[\[Omega], \[Eta], 3]/H43[2 \[Omega]/3, \[Eta], 3]


(* ::Section::Closed:: *)
(*g[\[Iota]]*)


(* ::Subsection::Closed:: *)
(*Y *)


ClearAll[SphericalHarm];

(* Define specific (l, m) branches *)
SphericalHarm[theta_, phi_, 2, 2]  := Sqrt[5/(64 Pi)] * (1 + Cos[theta])^2 * Exp[I * 2 * phi];
SphericalHarm[theta_, phi_, 2, 1]  := Sqrt[5/(16 Pi)] * Sin[theta] * (1 + Cos[theta]) * Exp[I * 1 * phi];
SphericalHarm[theta_, phi_, 2, -1] := Sqrt[5/(16 Pi)] * Sin[theta] * (1 - Cos[theta]) * Exp[I * (-1) * phi];
SphericalHarm[theta_, phi_, 2, -2] := Sqrt[5/(64 Pi)] * (1 - Cos[theta])^2 * Exp[I * (-2) * phi];

SphericalHarm[theta_, phi_, 3, 3]  := -Sqrt[21/(2 Pi)] * Cos[theta/2]^5 * Sin[theta/2] * Exp[I * 3 * phi];
SphericalHarm[theta_, phi_, 3, 2]  := Sqrt[7/Pi] * Cos[theta/2]^4 * (3 Cos[theta] - 2)/2 * Exp[I * 2 * phi];
SphericalHarm[theta_, phi_, 3, -2] := Sqrt[7/Pi] * Sin[theta/2]^4 * (3 Cos[theta] + 2)/2 * Exp[I * (-2) * phi];
SphericalHarm[theta_, phi_, 3, -3] := Sqrt[21/(2 Pi)] * Sin[theta/2]^5 * Cos[theta/2] * Exp[I * (-3) * phi];

SphericalHarm[theta_, phi_, 4, 4]  := 3 Sqrt[7/Pi] * Cos[theta/2]^6 * Sin[theta/2]^2 * Exp[I * 4 * phi];
SphericalHarm[theta_, phi_, 4, 3]  := -3 Sqrt[7/(2 Pi)] * Cos[theta/2]^5 * (2 Cos[theta] - 1) * Sin[theta/2] * Exp[I * 3 * phi];
SphericalHarm[theta_, phi_, 4, -3] := 3 Sqrt[7/(2 Pi)] * Sin[theta/2]^5 * (2 Cos[theta] + 1) * Cos[theta/2] * Exp[I * (-3) * phi];
SphericalHarm[theta_, phi_, 4, -4] := 3 Sqrt[7/Pi] * Sin[theta/2]^6 * Cos[theta/2]^2 * Exp[I * (-4) * phi];

(* Default branch for any other l, m combinations *)
SphericalHarm[theta_, phi_, l_, m_] := 0;


(* ::Subsection::Closed:: *)
(*g*)


g22Plus[\[Iota]_] = 1/2 (SphericalHarm[\[Iota], 0, 2, 2] + (-1)^2 SphericalHarm[\[Iota], 0, 2, -2])  1/2 Sqrt[(64 \[Pi])/5]//Simplify
g21Plus[\[Iota]_] = 1/2 (SphericalHarm[\[Iota], 0, 2, 1] + (-1)^2 SphericalHarm[\[Iota], 0, 2, -1])  1/2 Sqrt[(64 \[Pi])/5]//Simplify

g33Plus[\[Iota]_] = 1/2 (SphericalHarm[\[Iota], 0, 3, 3] + (-1)^3 SphericalHarm[\[Iota], 0, 3, -3])  1/2 Sqrt[(64 \[Pi])/5]//Simplify
g32Plus[\[Iota]_] = 1/2 (SphericalHarm[\[Iota], 0, 3, 2] + (-1)^3 SphericalHarm[\[Iota], 0, 3, -2])  1/2 Sqrt[(64 \[Pi])/5]//Simplify


g44Plus[\[Iota]_] = 1/2 (SphericalHarm[\[Iota], 0, 4, 4] + (-1)^4 SphericalHarm[\[Iota], 0, 4, -4])  1/2 Sqrt[(64 \[Pi])/5]//Simplify
g43Plus[\[Iota]_] = 1/2 (SphericalHarm[\[Iota], 0, 4, 3] + (-1)^4 SphericalHarm[\[Iota], 0, 4, -3])  1/2 Sqrt[(64 \[Pi])/5]//Simplify


g22Cross[\[Iota]_] = -(I/2) (SphericalHarm[\[Iota], 0, 2, 2] - (-1)^2 SphericalHarm[\[Iota], 0, 2, -2])  1/2 Sqrt[(64 \[Pi])/5]//Simplify
g21Cross[\[Iota]_] = -(I/2) (SphericalHarm[\[Iota], 0, 2, 1] - (-1)^2 SphericalHarm[\[Iota], 0, 2, -1])  1/2 Sqrt[(64 \[Pi])/5]//Simplify


g33Cross[\[Iota]_] = -(I/2) (SphericalHarm[\[Iota], 0, 3, 3] - (-1)^3 SphericalHarm[\[Iota], 0, 3, -3])  1/2 Sqrt[(64 \[Pi])/5]//Simplify
g32Cross[\[Iota]_] = -(I/2) (SphericalHarm[\[Iota], 0, 3, 2] - (-1)^3 SphericalHarm[\[Iota], 0, 3, -2])  1/2 Sqrt[(64 \[Pi])/5]//Simplify


g44Cross[\[Iota]_] = -(I/2) (SphericalHarm[\[Iota], 0, 4, 4] - (-1)^4 SphericalHarm[\[Iota], 0, 4, -4])  1/2 Sqrt[(64 \[Pi])/5]//Simplify
g43Cross[\[Iota]_] = -(I/2) (SphericalHarm[\[Iota], 0, 4, 3] - (-1)^4 SphericalHarm[\[Iota], 0, 4, -3])  1/2 Sqrt[(64 \[Pi])/5]//Simplify


(* ::Section::Closed:: *)
(*Testing against Josiel*)


(* ::Subsection::Closed:: *)
(*MyDef*)


PhenomD\[ScriptCapitalA][\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = cresIMR[[1,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


Clear@MMA
MMA[\[ScriptCapitalM]c_, \[Eta]_, \[Chi]1_, \[Chi]2_] := Block[
	{
		\[Omega]rd22, \[Omega]rd21, \[Omega]rd33, \[Omega]rd32, \[Omega]rd44, \[Omega]rd43,\[Chi]s, \[Chi]a,  freq, \[Omega], M, \[Omega]ref,
		\[ScriptCapitalA]22, \[ScriptCapitalA]21, \[ScriptCapitalA]33, \[ScriptCapitalA]32, \[ScriptCapitalA]44,\[ScriptCapitalA]43,
		G = 4.925490947641267`*^-6
	},
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5);
	
	\[Chi]s = (\[Chi]1+\[Chi]2)/2; \[Chi]a = (\[Chi]1-\[Chi]2)/2;
	
	\[Omega]rd22 = \[Omega]RdDampinglm[re\[Omega]22[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]rd21 = \[Omega]RdDampinglm[re\[Omega]21[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 1]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	\[Omega]rd33 = \[Omega]RdDampinglm[re\[Omega]33[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]rd32 = \[Omega]RdDampinglm[re\[Omega]32[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	\[Omega]rd44 = \[Omega]RdDampinglm[re\[Omega]44[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 4]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]rd43 = \[Omega]RdDampinglm[re\[Omega]43[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4,3]],  Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	freq =Range[1, 1024.];
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5);
	
	\[Omega] = M freq G;
	
	\[ScriptCapitalA]22 = PhenomD\[ScriptCapitalA][\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	\[ScriptCapitalA]21 = PhenomD\[ScriptCapitalA][f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd21, 1], \[Eta], \[Chi]s, \[Chi]a]*(H21[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd21, 1], \[Eta], \[Chi]s, \[Chi]a, 1]*H21[\[Omega], \[Eta], \[Chi]s, \[Chi]a, 1]/H21[2 \[Omega], \[Eta], \[Chi]s, \[Chi]a,1]);
	
	\[ScriptCapitalA]33 = PhenomD\[ScriptCapitalA][f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd33, 3],  \[Eta], \[Chi]s, \[Chi]a]*(H33[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd33, 3], \[Eta], 3]*(3/2)^(1/3));
	\[ScriptCapitalA]32 = PhenomD\[ScriptCapitalA][f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd32, 2], \[Eta], \[Chi]s, \[Chi]a]*(H32[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd32, 2], \[Eta], 2]);
	
	\[ScriptCapitalA]44 = PhenomD\[ScriptCapitalA][f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd44, 4], \[Eta], \[Chi]s, \[Chi]a]*(H44[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd44, 4], \[Eta], 4]*2^(2/3));
	\[ScriptCapitalA]43 = PhenomD\[ScriptCapitalA][f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd43, 3], \[Eta], \[Chi]s, \[Chi]a]*(H43[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd43, 3], \[Eta], 3]*3/2);
	
	M^2 {
		\[ScriptCapitalA]22, Re/@\[ScriptCapitalA]21, \[ScriptCapitalA]33, \[ScriptCapitalA]32, \[ScriptCapitalA]44, \[ScriptCapitalA]43
	
	}
	
	
	
]


(* ::Subsection::Closed:: *)
(*Josiel def*)


DeleteObject/@ExternalSessions[]
Clear@python
python = StartExternalSession[{
	"Python", (*you should change the evaluator to your installation:*)
	(*"Evaluator"-> "/home/cosmo-ufes/anaconda3/envs/gwfast_env/bin/python"*)
	"Evaluator" -> "/Users/felipe/anaconda3/envs/josiel/bin/python"
}];


ExternalEvaluate[python,"
import numpy as np
from GWDALI_v1.lib.IMRPhenomHM import josiel_Amp

freq = np.linspace(1, 1024,1024)
"]


josiel = ExternalFunction[python,"
def gwDALI_Amp(Mc, eta, chi1, chi2):
    phi0=0
    t0=0
    sx1=0
    sy1=0
    sz1=chi1

    sx2=0
    sy2=0
    sz2=chi2

    dL=1
    iota=0

    all_modes = josiel_Amp(
        dL,iota,phi0,t0,
        Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,
        freq
    )

    phi_22 = all_modes[2][2]
    phi_21 = all_modes[2][1]
    phi_33 = all_modes[3][3]
    phi_32 = all_modes[3][2]
    phi_44 = all_modes[4][4]
    phi_43 = all_modes[4][3]
    

    
    return np.array([
        phi_22, phi_21, 
        phi_33, phi_32, 
        phi_44, phi_43        
    ])
"]


(* ::Subsection::Closed:: *)
(*Test*)


(* ::Text:: *)
(*I believe Josiel is using frd and fdamp as calculated in the original PhenomD implementation. I am not convinced this is right. I find it more likely that those should be replaced by the frd and fdamp 22 formulas attached to Z22.*)


Clear@Test

Test[M_] := Module[
	
	{mc, eta, chi1, chi2},
	
	eta = RandomReal[{0.01, 0.25}];
	
	mc= M eta^(3/5);
	
	{chi1, chi2} = RandomReal[{-1,1}, 2];
	
	{
		(josiel[mc, eta, chi1, chi2]//Normal)//Re,
		(MMA[mc, eta, chi1, chi2])
	}
	
]


G =  4.925490947641267`*^-6;


plot[a_, M_] := {ListPlot[
	RelativeDiff@@(Extract[a, {{1,1}, {2,1}}]), 
	PlotRange->All,
	GridLines->{{0.014/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"22"
],



ListPlot[
	RelativeDiff@@(Extract[a, {{1,2}, {2,2}}]), 
	PlotRange->All,
	GridLines->{{0.014/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"21"
],



ListPlot[
	RelativeDiff@@(Extract[a, {{1,3}, {2,3}}]), 
	PlotRange->All,
	GridLines->{{0.014/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"33"
],


ListPlot[
	RelativeDiff@@(Extract[a, {{1,4}, {2,4}}]), 
	PlotRange->All,
	GridLines->{{0.014/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"32"
],




ListPlot[
	RelativeDiff@@(Extract[a, {{1,5}, {2,5}}]), 
	PlotRange->All,
	GridLines->{{0.014/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"44"
],



ListPlot[
	RelativeDiff@@(Extract[a, {{1,6}, {2,6}}]), 
	PlotRange->All,
	GridLines->{{0.014/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"43"
]}


Make:= Module[
	{M = RandomReal[{5,100}],a},
	Echo[0.2/(G M)];
	a = Test[M];
	plot[a, M]
]


Make


(* ::Chapter::Closed:: *)
(*Phase*)


(* ::Section::Closed:: *)
(*Inspiral*)


(* ::Subsection:: *)
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
	\[CapitalPhi]minus2[\[Eta]_], \[CapitalPhi]0[\[Eta]_], \[CapitalPhi]1[\[Eta]_], \[CapitalPhi]2[\[Eta]_], \[CapitalPhi]3[\[Eta]_, \[Chi]s_, \[Chi]a_], 
	\[CapitalPhi]4[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalPhi]5[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalPhi]5l[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalPhi]6[\[Eta]_, \[Chi]s_, \[Chi]a_], 
	\[CapitalPhi]6l[\[Eta]_], \[CapitalPhi]7[\[Eta]_, \[Chi]s_, \[Chi]a_]
} = InsVecPhase[\[Eta], \[Chi]s, \[Chi]a, 0,0,0,0,0,0,0,0,0,0];


Clear["\[CapitalSigma]*"]

{\[CapitalSigma]1[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]2[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]3[\[Eta]_, \[Chi]s_, \[Chi]a_], \[CapitalSigma]4[\[Eta]_, \[Chi]s_, \[Chi]a_]} = \[Eta]^-1*(PhenomCoeff[\[Eta], \[Chi]PN, #]&/@PhenomDTableV[[8;;11]])//.{
	\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s}//Simplify;


Expr = HoldComplete[
	{{\[CapitalPhi]0, \[CapitalPhi]2, \[CapitalPhi]3, \[CapitalPhi]4, \[CapitalPhi]5, \[CapitalPhi]5l, \[CapitalPhi]6, \[CapitalSigma]1, \[CapitalSigma]2, \[CapitalSigma]3, \[CapitalSigma]4}},
	
	\[CapitalPhi]0[\[Eta]_] := p0;
	\[CapitalPhi]2[\[Eta]_] := p2;
	\[CapitalPhi]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := p3;
	\[CapitalPhi]4[\[Eta]_, \[Chi]s_, \[Chi]a_] := p4;
	\[CapitalPhi]5[\[Eta]_, \[Chi]s_, \[Chi]a_] :=p5;
	\[CapitalPhi]5l[\[Eta]_, \[Chi]s_, \[Chi]a_] :=p5l;
	\[CapitalPhi]6[\[Eta]_, \[Chi]s_, \[Chi]a_] := p6;
	\[CapitalPhi]6l[\[Eta]_] := p6l;
	\[CapitalPhi]7[\[Eta]_, \[Chi]s_, \[Chi]a_] := p7;
	
	\[CapitalSigma]1[\[Eta]_, \[Chi]s_, \[Chi]a_] :=  S1;
	\[CapitalSigma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := S2;
	\[CapitalSigma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := S3;
	\[CapitalSigma]4[\[Eta]_, \[Chi]s_, \[Chi]a_] := S4;
	
	-\[Pi]/4 + 
	(
		\[CapitalPhi]0[\[Eta]]/(\[Pi]^(5/3) \[Omega]^(5/3))+\[CapitalPhi]2[\[Eta]]/(\[Pi] \[Omega])+\[CapitalPhi]3[\[Eta],\[Chi]s,\[Chi]a]/(\[Pi]^(2/3) \[Omega]^(2/3))+\[CapitalPhi]4[\[Eta],\[Chi]s,\[Chi]a]/(\[Pi]^(1/3) \[Omega]^(1/3))+\[CapitalPhi]5[\[Eta],\[Chi]s,\[Chi]a]+
		Log[\[Pi] \[Omega]] \[CapitalPhi]5l[\[Eta],\[Chi]s,\[Chi]a]+\[Pi]^(1/3) \[Omega]^(1/3) \[CapitalPhi]6[\[Eta],\[Chi]s,\[Chi]a]+\[Pi]^(1/3) \[Omega]^(1/3) Log[\[Pi] \[Omega]] \[CapitalPhi]6l[\[Eta]]+ \[Pi]^(2/3) \[Omega]^(2/3) \[CapitalPhi]7[\[Eta],\[Chi]s,\[Chi]a]
	) + (
		\[Omega] \[CapitalSigma]1[\[Eta],\[Chi]s,\[Chi]a]+3/4 \[Omega]^(4/3) \[CapitalSigma]2[\[Eta],\[Chi]s,\[Chi]a]+3/5 \[Omega]^(5/3) \[CapitalSigma]3[\[Eta],\[Chi]s,\[Chi]a]+1/2 \[Omega]^2 \[CapitalSigma]4[\[Eta],\[Chi]s,\[Chi]a]
	)
]//.{
	p0-> \[CapitalPhi]0[\[Eta]], p2-> \[CapitalPhi]2[\[Eta]], p3-> \[CapitalPhi]3[\[Eta], \[Chi]s, \[Chi]a], p4-> \[CapitalPhi]4[\[Eta], \[Chi]s, \[Chi]a], p5-> \[CapitalPhi]5[\[Eta], \[Chi]s, \[Chi]a],
	p5l-> \[CapitalPhi]5l[\[Eta], \[Chi]s, \[Chi]a], p6-> \[CapitalPhi]6[\[Eta], \[Chi]s, \[Chi]a], p6l-> \[CapitalPhi]6l[\[Eta]],p7-> \[CapitalPhi]7[\[Eta], \[Chi]s, \[Chi]a],
	
	S1 ->\[CapitalSigma]1[\[Eta], \[Chi]s, \[Chi]a], S2 ->\[CapitalSigma]2[\[Eta], \[Chi]s, \[Chi]a], S3->\[CapitalSigma]3[\[Eta], \[Chi]s, \[Chi]a], S4->\[CapitalSigma]4[\[Eta], \[Chi]s, \[Chi]a]
};


Expr = $Block@@@(HoldForm[Evaluate@Expr]);


HM\[CapitalPhi]Ins//ClearAll

expr = Hold[
	{HM\[CapitalPhi]Ins, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 4},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalPhi]resInspiral = EchoTiming[DerivativeRules@@expr];


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = \[CapitalPhi]resInspiral[[i,2]];
	
	DownValues[test] = DownValues[test]//.{HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@\[CapitalPhi]resInspiral}]


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

compileThis[x_Integer] := ZeroFunction


c\[CapitalPhi]resInspiral = MapAt[
	compileThis,
	\[CapitalPhi]resInspiral,
	{All,2}
];


Do[
	c\[CapitalPhi]resInspiral[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[cresInspiral[[i, 1]][[All,1]]]},
	{i, Length@cresInspiral}
] (*The error comes from derivatives that are zero.*)


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := c\[CapitalPhi]resInspiral[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	
	
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@resInspiral}];//RepeatedTiming


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
	
	UpValues@HM\[CapitalPhi]Ins = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[CapitalPhi]Ins[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = c\[CapitalPhi]resInspiral[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


(* ::Section::Closed:: *)
(*Int-Phase*)


(*Find vectors for Intermediate Phase and their variables:*)
Position[vectorDefs[[All, 1, All, 0]] /. HoldPattern -> Identity, #] & /@ {IntVecPhase, \[Omega]IntVecPhase}

vectorDefs[[3 ;; 4, 1]]

DownValues[IntVecPhase] = {vectorDefs[[3]]};
DownValues[\[Omega]IntVecPhase] = {vectorDefs[[4]]};


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


HM\[CapitalPhi]Int//ClearAll
expr = Hold[
	{HM\[CapitalPhi]Int, {\[Omega], \[Eta], \[Chi]s, \[Chi]a}, 4},
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

compileThis[x_Integer] := HMZeroFunction


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


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = c\[CapitalPhi]resInt[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@HM\[CapitalPhi]Int = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[CapitalPhi]Int[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] = c\[CapitalPhi]resInt[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a];


(* ::Section::Closed:: *)
(*Ringdown and Damping frequency*)


(* ::Subsection::Closed:: *)
(*Commom to PhenomD*)


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


(* ::Subsection::Closed:: *)
(*Exclusive to PhenomHM*)


\[Kappa][aeff_, l_, m_] := (Log[2-aeff]/Log[3])^(1/(2+ l- Abs[m]))


Z22[\[Kappa]22_] := (
	1 + 
	\[Kappa]22*1.557847*Exp[2.903124 I] + 
	\[Kappa]22^2*1.95097051*Exp[5.920970 I] + 
	\[Kappa]22^3*2.09971716 * Exp[2.760585 I] + 
	\[Kappa]22^4*1.41094660 * Exp[5.914340 I]+ 
	\[Kappa]22^5*0.41063923*Exp[2.795235 I]
)

Z32[\[Kappa]32_] := (
	1.022464*Exp[0.004870 I] + 
	\[Kappa]32*0.24731213*Exp[0.665292 I] + 
	\[Kappa]32^2*1.70468239*Exp[3.138283 I] + 
	\[Kappa]32^3*0.94604882*Exp[0.163247 I] + 
	\[Kappa]32^4 * 1.53189884*Exp[5.703573 I] + 
	\[Kappa]32^5 * 2.28052668*Exp[2.685231 I] + 
	\[Kappa]32^6 * 0.92150314*Exp[5.841704 I]
)


Z44[\[Kappa]44_] := (
	2 + 
	\[Kappa]44*2.658908*Exp[3.002787 I] + 
	\[Kappa]44^2 * 2.97825567 * Exp[6.050955 I] + 
	\[Kappa]44^3 * 3.21842350 * Exp[2.877514 I] + 
	\[Kappa]44^4 * 2.12764967 * Exp[5.989669 I] + 
	\[Kappa]44^5 * 0.60338186 * Exp[2.830031 I]
)

Z21[\[Kappa]21_] := (
	0.589113 * Exp[0.043525 I] + 
	\[Kappa]21*0.18896353  * Exp[2.289868 I] + 
	\[Kappa]21^2*1.15012965  * Exp[5.810057 I] + 
	\[Kappa]21^3*6.04585476  * Exp[2.741967 I] + 
	\[Kappa]21^4*11.12627777 * Exp[5.844130 I] + 
	\[Kappa]21^5*9.34711461  * Exp[2.669372 I] + 
	\[Kappa]21^6*3.03838318  * Exp[5.791518 I]
)


Z33[\[Kappa]33_] := (
	1.5 +
	\[Kappa]33 * 2.095657   * Exp[2.964973 I] +
	\[Kappa]33^2 * 2.46964352 * Exp[5.996734 I] + 
	\[Kappa]33^3 * 2.66552551 * Exp[2.817591 I] + 
	\[Kappa]33^4 * 1.75836443 * Exp[5.932693 I] + 
	\[Kappa]33^5 * 0.49905688 * Exp[2.781658 I]
)



Z43[\[Kappa]43_] := (
	1.5 + 
	\[Kappa]43 * 0.205046   * Exp[0.595328 I] + 
	\[Kappa]43^2* 3.10333396 * Exp[3.016200 I] + 
	\[Kappa]43^3 * 4.23612166 * Exp[6.038842 I] + 
	\[Kappa]43^4* 3.02890198 * Exp[2.826239 I] + 
	\[Kappa]43^5 * 0.90843949 * Exp[5.915164 I]
)


(* ::Text:: *)
(*This folllow the reasoning of re\[Omega] and im\[Omega] of last subsection:*)


re\[Omega]22[\[Kappa]22_] = Re[Z22[\[Kappa]22]]//Simplify[#, Assumptions-> \[Kappa]22 \[Element] Reals]&;
im\[Omega]22[\[Kappa]22_] = Im[Z22[\[Kappa]22]]//Simplify[#, Assumptions-> \[Kappa]22 \[Element] Reals]&;

re\[Omega]32[\[Kappa]32_] = Re[Z32[\[Kappa]32]]//Simplify[#, Assumptions-> \[Kappa]32 \[Element] Reals]&;
im\[Omega]32[\[Kappa]32_] = Im[Z32[\[Kappa]32]]//Simplify[#, Assumptions-> \[Kappa]32 \[Element] Reals]&;


re\[Omega]44[\[Kappa]44_] = Re[Z44[\[Kappa]44]]//Simplify[#, Assumptions-> \[Kappa]44 \[Element] Reals]&;
im\[Omega]44[\[Kappa]44_] = Im[Z44[\[Kappa]44]]//Simplify[#, Assumptions-> \[Kappa]44 \[Element] Reals]&;

re\[Omega]21[\[Kappa]21_] = Re[Z21[\[Kappa]21]]//Simplify[#, Assumptions-> \[Kappa]21 \[Element] Reals]&;
im\[Omega]21[\[Kappa]21_] = Im[Z21[\[Kappa]21]]//Simplify[#, Assumptions-> \[Kappa]21 \[Element] Reals]&;


re\[Omega]33[\[Kappa]33_] = Re[Z33[\[Kappa]33]]//Simplify[#, Assumptions-> \[Kappa]33 \[Element] Reals]&;
im\[Omega]33[\[Kappa]33_] = Im[Z33[\[Kappa]33]]//Simplify[#, Assumptions-> \[Kappa]33 \[Element] Reals]&;

re\[Omega]43[\[Kappa]43_] = Re[Z43[\[Kappa]43]]//Simplify[#, Assumptions-> \[Kappa]43 \[Element] Reals]&;
im\[Omega]43[\[Kappa]43_] = Im[Z43[\[Kappa]43]]//Simplify[#, Assumptions-> \[Kappa]43 \[Element] Reals]&;


\[Omega]RdDampinglm[int_, Erad_] := int/(2 \[Pi] (1-Erad))


f22[\[Omega]_, \[Omega]rd22_, \[Omega]rdlm_, m_] := Module[
	{f0, case1, case2, case3},
	
	f0 = 0.018 \[Omega]rdlm/\[Omega]rd22;
	
	
	case1 =2/m \[Omega];
	
	case2 = (\[Omega]rd22-2 f0/m)/(\[Omega]rdlm-f0) (\[Omega]-f0) + 2 f0/m;
	
	case3 = \[Omega]rd22/\[Omega]rdlm \[Omega];
	
	Which[
		\[Omega] <= f0, case1,
		f0 < \[Omega] <= \[Omega]rdlm, case2,
		\[Omega]rdlm<\[Omega], case3
	]
]


(* ::Section::Closed:: *)
(*MR phase*)


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
		(\[Rho]lm ArcTan[((\[Omega]-ringdownFrequency \[Alpha]5[\[Eta],\[Chi]s,\[Chi]a]))/(dampingFrequency \[Rho]lm \[Tau]lm)] \[Alpha]4[\[Eta],\[Chi]s,\[Chi]a])
	)
]//.{
	a1->\[Alpha]1[\[Eta], \[Chi]s, \[Chi]a], a2->\[Alpha]2[\[Eta], \[Chi]s, \[Chi]a], a3->\[Alpha]3[\[Eta], \[Chi]s, \[Chi]a],a4->\[Alpha]4[\[Eta], \[Chi]s, \[Chi]a], a5->\[Alpha]5[\[Eta], \[Chi]s, \[Chi]a],
	
	erad->Erad[\[Eta], \[Chi]s, \[Chi]a], ae->aeff[\[Eta], \[Chi]s, \[Chi]a],
	re->re\[Omega][\[Kappa]], im->im\[Omega][\[Kappa]],
	ringdownFrequency :> \[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	dampingFrequency :> \[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]]
};

Expr = $Block@@@(HoldForm[Evaluate@Expr]);


ClearAll@HM\[CapitalPhi]MR
expr = Hold[
	{HM\[CapitalPhi]MR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm}, 4},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalPhi]resMR = EchoTiming[DerivativeRules@@expr];


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] = \[CapitalPhi]resMR[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, \[Rho]lm$->\[Rho]lm, \[Tau]lm$->\[Tau]lm
	};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, 1.2,1.1]
]


Table[T[i], {i, Length@\[CapitalPhi]resMR}]//RepeatedTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm}, 
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
	c\[CapitalPhi]resMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm}, X]//.{X -> ToString[c\[CapitalPhi]resMR[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPhi]resMR}
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] := c\[CapitalPhi]resMR[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, 1., 1.]
]


Table[T[i], {i, Length@c\[CapitalPhi]resMR}]//RepeatedTiming


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
	
	UpValues@HM\[CapitalPhi]MR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[CapitalPhi]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] = c\[CapitalPhi]resMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];


(* ::Section::Closed:: *)
(*\[CapitalPhi]IMR*)


Block[
	{\[Beta]0, \[Beta]1, \[Alpha]0, \[Alpha]1},
	
	\[Beta]1 = Derivative[1,0,0,0][\[Phi]Ins][0.018`, \[Eta], \[Chi]s, \[Chi]a] - Derivative[1,0,0,0][\[Phi]Int][0.018`, \[Eta], \[Chi]s, \[Chi]a];
	\[Beta]0 = \[Phi]Ins[0.018, \[Eta], \[Chi]s, \[Chi]a] - \[Phi]Int[0.018`, \[Eta], \[Chi]s, \[Chi]a] - \[Beta]1 0.018`;
	
	\[Alpha]1 = (
		Derivative[1,0,0,0][\[Phi]Int][ringdownFrequency/2, \[Eta], \[Chi]s, \[Chi]a] +
		\[Beta]1 - 
		Derivative[1,0,0,0,0,0][\[Phi]MR][ringdownFrequency/2,  \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm]
	)//Simplify;
	
	\[Alpha]0 = (
		\[Phi]Int[ringdownFrequency/2, \[Eta], \[Chi]s, \[Chi]a] + 
		\[Beta]0 + \[Beta]1 ringdownFrequency/2 - 
		\[Phi]MR[ringdownFrequency/2, \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm] - 
		\[Alpha]1 ringdownFrequency/2
	)//Simplify;
	
	Expr = HoldComplete[
		{{\[Omega]RdDamping, re\[Omega], Erad, aeff, \[Phi]Ins, \[Phi]Int, \[Phi]MR}},
	
		Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
		aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
		re\[Omega][\[Kappa]_] := re;
		\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
		
		\[Phi]Ins[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[CapitalPhi]Ins[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
		\[Phi]Int[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[CapitalPhi]Int[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
		\[Phi]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] := U\[CapitalPhi]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];
	
	
		\[Phi]Ins[\[Omega], \[Eta], \[Chi]s, \[Chi]a]*us\[Theta][0.018`- \[Omega]] + 
		(\[Phi]Int[\[Omega], \[Eta], \[Chi]s, \[Chi]a] + b0 + b1 \[Omega])*us\[Theta][(\[Omega]-0.018`) (ringdownFrequency/2-\[Omega])] + 
		(\[Phi]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm] + a0 + a1 \[Omega])*us\[Theta][(\[Omega]-ringdownFrequency/2) (0.2-\[Omega])]
		
	]//.{
		b0->\[Beta]0, b1->\[Beta]1, 
		a0->\[Alpha]0, a1->\[Alpha]1,
		erad->Erad[\[Eta], \[Chi]s, \[Chi]a], ae->aeff[\[Eta], \[Chi]s, \[Chi]a], kappa->\[Kappa][aeff, l, m],
		re->re\[Omega][\[Kappa]], 
		ringdownFrequency :> \[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]]
	};
]

Expr = $Block@@@(HoldForm[Evaluate@Expr]);


HM\[CapitalPhi]IMR//ClearAll
expr = Hold[
	{HM\[CapitalPhi]IMR, {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm}, 3},
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
	
	\[CapitalPhi]resIMR2 = temp//.{U\[CapitalPhi]Ins->HM\[CapitalPhi]Ins, U\[CapitalPhi]Int->HM\[CapitalPhi]Int, U\[CapitalPhi]MR->HM\[CapitalPhi]MR};
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] = \[CapitalPhi]resIMR2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, \[Rho]lm$->\[Rho]lm, \[Tau]lm$->\[Tau]lm
	};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2,1.,1.]
]


Table[T[i], {i, Length@\[CapitalPhi]resIMR2}]//RepeatedTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm}, 
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
	c\[CapitalPhi]resIMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[c\[CapitalPhi]resIMR[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPhi]resIMR}
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] := c\[CapitalPhi]resIMR[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, 1., 1.]
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
	
	UpValues@HM\[CapitalPhi]IMR = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[CapitalPhi]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] = c\[CapitalPhi]resIMR[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];


(* ::Section::Closed:: *)
(*\[CapitalPsi]Dlm*)


Block[
	{C01, C02, \[Omega]rd22, f0, am, bm, ar, ai=2/m, \[Rho]lm, \[Tau]lm},
	
	
	\[Omega]rd22 = ringdownFrequency;
	
	f0 = 0.018 \[Omega]rdlm/\[Omega]rd22;
	
	am = (\[Omega]rd22-2 f0/m)/(\[Omega]rdlm-f0);
	bm = f0 (2/m - am);
	ar = \[Omega]rd22/\[Omega]rdlm;
	
	\[Rho]lm = \[Omega]rd22/\[Omega]rdlm;
	\[Tau]lm = \[Omega]damplm/DampingFrequency;
	
	
	C01 = 1/ai \[Phi]IMR[ai f0, \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm] - 1/am \[Phi]IMR[am f0 + bm, \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];
	C02 = 1/am \[Phi]IMR[am \[Omega]rdlm + bm, \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm] - 1/ar \[Phi]IMR[ar \[Omega]rdlm, \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm] + C01;
	
	
	Expr = HoldComplete[
		{{\[Omega]RdDampinglm, re\[Omega]22, im\[Omega]22, \[Kappa], Erad, aeff, \[Phi]IMR, shift}},
	
		Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
		aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
		\[Kappa][aeff_, l_, m_] := kappa;
		re\[Omega]22[\[Kappa]_] := re22;
		im\[Omega]22[\[Kappa]_] := im22;
		\[Omega]RdDampinglm[int_, Erad_] := int/(2 \[Pi] (1-Erad));
		
		shift[m_] := Arg[Exp[I \[Pi]/2 (2-m)]];
		
		\[Phi]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, rholm_, taulm_] := U\[CapitalPhi]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, rholm, taulm];
	
		shift[m] + 
		
		\[Phi]IMR[Ai \[Omega], \[Eta], \[Chi]s, \[Chi]a, rlm, tlm]/Ai*us\[Theta][F0 - \[Omega]] +
		
		(\[Phi]IMR[Am \[Omega] + Bm, \[Eta], \[Chi]s, \[Chi]a, rlm, tlm]/Am + c01)*us\[Theta][(\[Omega]-F0) (\[Omega]rdlm-\[Omega])] + 
		
		(\[Phi]IMR[Ar \[Omega], \[Eta], \[Chi]s, \[Chi]a, rlm, tlm]/Ar + c02)*us\[Theta][(\[Omega] - \[Omega]rdlm) (0.2-\[Omega])] 
		
	]//.{
		Ai->ai, Bm->bm, Am->am, Ar->ar,c01->C01, c02->C02,F0->f0,
		
		erad->Erad[\[Eta], \[Chi]s, \[Chi]a], ae->aeff[\[Eta], \[Chi]s, \[Chi]a], kappa->\[Kappa][aeff, l, m],
		
		re22->re\[Omega]22[\[Kappa]], im22->im\[Omega]22[\[Kappa]],
		
		ringdownFrequency :> \[Omega]RdDampinglm[re\[Omega]22[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
		DampingFrequency :> \[Omega]RdDampinglm[im\[Omega]22[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
		
		rlm-> ringdownFrequency/\[Omega]rdlm,
		tlm -> \[Omega]damplm/DampingFrequency
	};
]

Expr = $Block@@@(HoldForm[Evaluate@Expr]);


Combinations[vars_List, n_Integer]/;n>0 := Module[
	{result},
	result  = Table[
			(Sort/@Tuples[vars, i])//DeleteDuplicates,
			{i,n}
	]//Flatten[#,1]&
]


derivatives = Combinations[{\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rdlm, \[Omega]damplm}, 3];
PrependTo[derivatives, {}];


ClearAll@HM\[CapitalPsi]lm

expr = Hold[
	{HM\[CapitalPsi]lm, {\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rdlm, \[Omega]damplm, m}, derivatives},
	Evaluate@Expr,
	"IncludeZeroDerivative"->False
]//.HoldForm[X_] :> X;


\[CapitalPsi]reslm = EchoTiming[DerivativeRules@@expr];


Module[
	{rule, exprs, temp},
	
	exprs ={
		(0.036` \[Omega]rdlm)/(m $x2568), 
		$x2568, 
		(\[Omega]rdlm (-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568))/(\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568)+(0.018` \[Omega]rdlm (2/m-(-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568)/(\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568)))/$x2568, 
		(0.018` \[Omega]rdlm (-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568))/((\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568) $x2568)+(0.018` \[Omega]rdlm (2/m-(-((0.036` \[Omega]rdlm)/(m $x2568))+$x2568)/(\[Omega]rdlm-(0.018` \[Omega]rdlm)/$x2568)))/$x2568
	}//.{$x2568->$x4476};
	
	rule = {
		Derivative[n__][U\[CapitalPhi]IMR][x__] -> $D[{n}, U\[CapitalPhi]IMR][x],
		$D[{n__}, U\[CapitalPhi]IMR][y1_, y2__]/;(
			y1[[0]]=!=List && Length[Cases[y1-exprs, x_/;x==0]]>0
		) :> Last[$D[{n}, U\[CapitalPhi]IMR][{y1}, y2]],
		
		U\[CapitalPhi]IMR[y1_, y2__]/;(
			y1[[0]]=!=List && Length[Cases[y1-exprs, x_/;x==0]]>0
		)  :>  Last[U\[CapitalPhi]IMR[{y1}, y2]]
	};
	
	
	
	
	temp =  \[CapitalPsi]reslm//.rule;
	
	\[CapitalPsi]reslm2 = temp//.U\[CapitalPhi]IMR->HM\[CapitalPhi]IMR;
]


Clear@T

T[i_] := Block[
	{test, ringdownFrequency,DampingFrequency, eta, chis, chia},
	{eta, chis, chia} = {0.2, 0.1, 0.2};
	
	ringdownFrequency = \[Omega]RdDampinglm[re\[Omega]22[\[Kappa][aeff[eta, chis, chia], 2, 2]], Erad[eta, chis, chia]];
	DampingFrequency = \[Omega]RdDampinglm[im\[Omega]22[\[Kappa][aeff[eta, chis, chia], 2, 2]], Erad[eta, chis, chia]];
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Omega]rdlm_, \[Omega]damplm_, m_] = \[CapitalPsi]reslm2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, \[Omega]rdlm$->\[Omega]rdlm, \[Omega]damplm$->\[Omega]damplm, m$->m
	};
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, ringdownFrequency, DampingFrequency, 2]
]


Table[T[i], {i, Length@\[CapitalPsi]reslm2}]//RepeatedTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Eta], \[Chi]s, \[Chi]a, \[Omega]rdlm, \[Omega]damplm, m}, 
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


c\[CapitalPsi]reslm = MapAt[
	compileThis,
	\[CapitalPsi]reslm2,
	{All,2}
];


Do[
	c\[CapitalPhi]resIMR[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rdlm, \[Omega]damplm, m}, X]//.{X -> ToString[c\[CapitalPhi]resIMR[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalPsi]reslm}
]


Clear@T

T[i_] := Block[
	{test, ringdownFrequency,DampingFrequency, eta, chis, chia},
	{eta, chis, chia} = {0.2, 0.1, 0.2};
	
	ringdownFrequency = \[Omega]RdDampinglm[re\[Omega]22[\[Kappa][aeff[eta, chis, chia], 2, 2]], Erad[eta, chis, chia]];
	DampingFrequency = \[Omega]RdDampinglm[im\[Omega]22[\[Kappa][aeff[eta, chis, chia], 2, 2]], Erad[eta, chis, chia]];
	
	test[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_, m_] := c\[CapitalPsi]reslm[[i,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm, m];
	
	
	test[Range[10^-4, 0.124, 10^-4], 0.2, 0.1, 0.2, ringdownFrequency, DampingFrequency, 2]
]


Table[T[i], {i, Length@c\[CapitalPsi]reslm}]//RepeatedTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = c\[CapitalPsi]reslm[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@HM\[CapitalPsi]lm = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HM\[CapitalPsi]lm[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Omega]rdlm_, \[Omega]damplm_, m_] = c\[CapitalPsi]reslm[[1, 2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rdlm, \[Omega]damplm, m];


(* ::Section::Closed:: *)
(*t0*)


Clear[\[Gamma]2, \[Gamma]3]
\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[6]]//.{\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};
\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] = PhenomCoeff[\[Eta], \[Chi]PN, #]&@PhenomDTableV[[7]]//.{\[Chi]PN-> Sqrt[1-4 \[Eta]] \[Chi]a + (1-76 \[Eta]/113) \[Chi]s};


Expr = HoldComplete[
		{{\[Gamma]2, \[Gamma]3, \[Omega]Peak, \[Omega]RdDampinglm, re\[Omega], im\[Omega],  Erad, aeff, \[Phi]MR}},
	
		Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
		aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
		re\[Omega][\[Kappa]_] := re;
		im\[Omega][\[Kappa]_] := im;
		\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
		\[Omega]Peak[\[Omega]RD_,\[Omega]DAMP_,\[Gamma]2_,\[Gamma]3_] := If[\[Gamma]2<=1,\[Omega]RD+(\[Omega]DAMP \[Gamma]3 (Sqrt[1-\[Gamma]2^2]-1))/\[Gamma]2,\[Omega]RD-(\[Omega]DAMP \[Gamma]3)/\[Gamma]2];
		\[Phi]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] := U\[CapitalPhi]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];
		\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := g2;
		\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := g3;
		
	
		Derivative[1, 0,0,0,0,0][\[Phi]MR][Peak, \[Eta], \[Chi]s, \[Chi]a, 1, 1]*(\[Omega]-\[Omega]ref)
		
	]//.{
		erad->Erad[\[Eta], \[Chi]s, \[Chi]a], ae->aeff[\[Eta], \[Chi]s, \[Chi]a], kappa->\[Kappa][aeff, l, m],
		re->re\[Omega][\[Kappa]], im->im\[Omega][\[Kappa]],
		Peak :> \[Omega]Peak[ringdownFrequency, dampingFrequency, \[Gamma]2[\[Eta], \[Chi]s, \[Chi]a], \[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]],
		ringdownFrequency :> \[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
		dampingFrequency :> \[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
		g2->\[Gamma]2[\[Eta], \[Chi]s, \[Chi]a],g3->\[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]
	};


Expr = $Block@@@(HoldForm[Evaluate@Expr]);


\[CapitalDelta]t0//ClearAll

expr = Hold[
	{\[CapitalDelta]t0, {\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a}, 3},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


\[CapitalDelta]t0res = EchoTiming[DerivativeRules@@expr];


Module[
	{rule, temp},
	
	rule = {
		Derivative[n__][U\[CapitalPhi]MR][y1_, y2__]/;y1[[0]]===Symbol :> Last[$D[{n}, U\[CapitalPhi]MR][{y1}, y2]]
	};
	
	
	
	temp = \[CapitalDelta]t0res//.rule;
	
	\[CapitalDelta]t0res2 = temp//.U\[CapitalPhi]MR->HM\[CapitalPhi]MR;
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_] = \[CapitalDelta]t0res2[[i,2]];
	
	DownValues[test] = DownValues[test]//.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, \[Omega]ref$->\[Omega]ref
	};
	
	test[0.1, 0.02, 0.1, 0.2, 0.12]
]


Table[T[i], {i, Length@\[CapitalDelta]t0res}]


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll
compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a}, 
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

compileThis[x_Integer] := HMZeroFunction


c\[CapitalDelta]t0res = MapAt[
	compileThis,
	\[CapitalDelta]t0res2,
	{All,2}
];


Do[
	c\[CapitalDelta]t0res[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a}, X]//.{X -> ToString[c\[CapitalDelta]t0res[[i, 1]][[All,1]]]},
	{i, Length@c\[CapitalDelta]t0res}
]


Clear@T

T[i_] := Block[
	{test},
	
	
	test[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_] = c\[CapitalDelta]t0res[[i,2]][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a];
	
	
	test[Range[10^-4, 0.124, 10^-4], 10^-4, 0.2, 0.1, 0.2]
]


Table[T[i], {i, Length@c\[CapitalDelta]t0res}]//RepeatedTiming


Module[
	{list = c\[CapitalDelta]t0res[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@\[CapitalDelta]t0 = list;
]


(* ::Text:: *)
(*Set DownValues:*)


\[CapitalDelta]t0[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_] = c\[CapitalDelta]t0res[[1, 2]][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a];


$D[{1, 0, 0, 0, 0}, \[CapitalDelta]t0][{0.001, 0.002}, 0.001, 0.2, 0.1, 0.3]


(* ::Section::Closed:: *)
(*Testing against Josiel*)


(* ::Subsection::Closed:: *)
(*My def*)


PhenomD\[CapitalPsi]lm[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Omega]rdlm_,\[Omega]damplm_,m_] = c\[CapitalPsi]reslm[[1,2]][\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rdlm,\[Omega]damplm,m];


PhenomD\[CapitalDelta]t0[\[Omega]_, \[Omega]ref_,\[Eta]_,\[Chi]s_,\[Chi]a_] = c\[CapitalDelta]t0res[[1,2]][\[Omega], \[Omega]ref,\[Eta],\[Chi]s,\[Chi]a];


MMAPhase[\[ScriptCapitalM]c_, \[Eta]_, \[Chi]1_, \[Chi]2_] := Module[
	{
		\[Omega]rd22, \[Omega]rd21, \[Omega]rd33, \[Omega]rd32, \[Omega]rd44, \[Omega]rd43,
		\[Omega]d22, \[Omega]d21, \[Omega]d33, \[Omega]d32, \[Omega]d44, \[Omega]d43, \[Chi]s, \[Chi]a, freq, \[Omega], M, \[Omega]ref,
		\[Psi]22, \[Psi]21, \[Psi]33, \[Psi]32, \[Psi]44,\[Psi]43, \[CapitalDelta]t0,\[Phi]0,
		G = 4.925490947641267`*^-6
	},
	
	
	\[Chi]s = (\[Chi]1+\[Chi]2)/2; \[Chi]a = (\[Chi]1-\[Chi]2)/2;
	
	
	
	\[Omega]rd22 = \[Omega]RdDampinglm[re\[Omega]22[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]d22 = \[Omega]RdDampinglm[im\[Omega]22[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	\[Omega]rd21 = \[Omega]RdDampinglm[re\[Omega]21[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 1]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]d21 = \[Omega]RdDampinglm[im\[Omega]21[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 1]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	\[Omega]rd33 = \[Omega]RdDampinglm[re\[Omega]33[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]d33 = \[Omega]RdDampinglm[im\[Omega]33[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	\[Omega]rd32 = \[Omega]RdDampinglm[re\[Omega]32[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]d32 = \[Omega]RdDampinglm[im\[Omega]32[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	\[Omega]rd44 = \[Omega]RdDampinglm[re\[Omega]44[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 4]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]d44 =\[Omega]RdDampinglm[im\[Omega]44[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 4]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	
	\[Omega]rd43 = \[Omega]RdDampinglm[re\[Omega]43[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	\[Omega]d43 = \[Omega]RdDampinglm[im\[Omega]43[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]];
	
	
	freq = Range[1,1024];
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5);
	
	\[Omega] = M freq G;
	
	\[Omega]ref = M freq[[1]] G;
	
	\[CapitalDelta]t0 = PhenomD\[CapitalDelta]t0[\[Omega], \[Omega]ref, \[Eta], \[Chi]s,\[Chi]a];
	
	\[Phi]0 =  Last[HM\[CapitalPhi]IMR[{\[Omega]ref}, \[Eta], \[Chi]s, \[Chi]a, 1, 1]/2];
	
	
	\[Psi]22 = HM\[CapitalPhi]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, 1, 1] - \[CapitalDelta]t0 - 2 \[Phi]0;
	\[Psi]21 = PhenomD\[CapitalPsi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd21, \[Omega]d21, 1]- \[CapitalDelta]t0 -  \[Phi]0;
	
	\[Psi]33 = PhenomD\[CapitalPsi]lm[\[Omega], \[Eta],\[Chi]s, \[Chi]a, \[Omega]rd33, \[Omega]d33, 3]- \[CapitalDelta]t0 - 3 \[Phi]0;
	\[Psi]32 = PhenomD\[CapitalPsi]lm[\[Omega], \[Eta],\[Chi]s, \[Chi]a, \[Omega]rd32, \[Omega]d32, 2] - \[CapitalDelta]t0 - 2 \[Phi]0;
	
	\[Psi]44 = PhenomD\[CapitalPsi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd44, \[Omega]d44, 4]- \[CapitalDelta]t0 - 4 \[Phi]0;
	\[Psi]43 = PhenomD\[CapitalPsi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd43, \[Omega]d43, 3]- \[CapitalDelta]t0 - 3 \[Phi]0;
	
	{
		\[Psi]22, \[Psi]21, \[Psi]33, \[Psi]32, \[Psi]44, \[Psi]43	
	
	}
	
	
	
]


(* ::Subsection::Closed:: *)
(*Josiel def*)


DeleteObject/@ExternalSessions[]
Clear@python
python = StartExternalSession[{
	"Python", (*you should change the evaluator to your installation:*)
	(*"Evaluator"-> "/home/cosmo-ufes/anaconda3/envs/gwfast_env/bin/python"*)
	"Evaluator" -> "/Users/felipe/anaconda3/envs/josiel/bin/python"
}];


ExternalEvaluate[python,"
import numpy as np
from GWDALI_v1.lib.IMRPhenomHM import josiel_phase

freq = np.linspace(1, 1024, 1024)
"]


Range[10,1024]//Length


josiel = ExternalFunction[python,"
def gwDALI_phase(Mc, eta, chi1, chi2):
    phi0=0
    t0=0
    sx1=0
    sy1=0
    sz1=chi1

    sx2=0
    sy2=0
    sz2=chi2

    dL=1
    iota=0

    all_modes = josiel_phase(
        dL,iota,phi0,t0,
        Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,
        freq
    )

    phi_22 = all_modes[2][2]
    phi_21 = all_modes[2][1]
    phi_33 = all_modes[3][3]
    phi_32 = all_modes[3][2]
    phi_44 = all_modes[4][4]
    phi_43 = all_modes[4][3]
    

    
    return np.array([
        phi_22, phi_21, 
        phi_33, phi_32, 
        phi_44, phi_43        
    ])
"]


(* ::Subsection::Closed:: *)
(*Test*)


(* ::Text:: *)
(*I believe Josiel is using frd and fdamp as calculated in the original PhenomD implementation. I am not convinced this is right. I find it more likely that those should be replaced by the frd and fdamp 22 formulas attached to Z22.*)


Clear@Test
Test[M_] := Module[
	
	{mc, eta, chi1, chi2, \[Chi]s,\[Chi]a},
	
	eta = RandomReal[{0.1, 0.25}]//Echo;
	
	mc= M eta^(3/5);
	
	{chi1, chi2} = RandomReal[{-0.9,0.9}, 2]//Echo;
	
	\[Chi]s=(chi1+chi2)/2; \[Chi]a=(chi1-chi2)/2;
	

	
	(*\[Omega]RdDamping[re\[Omega][aeff[eta, \[Chi]s, \[Chi]a]], Erad[eta, \[Chi]s, \[Chi]a]]//Echo;
	PhenomD\[CapitalDelta]t0[{1}, 0, eta,\[Chi]s,\[Chi]a]//Echo;*)
	
	{
		(josiel[mc, eta, chi1, chi2]//Normal),
		MMAPhase[mc, eta, chi1, chi2]
	}
	
]


plot[a_, M_] := {ListPlot[
	RelativeDiff@@(Extract[a, {{1,1}, {2,1}}]), 
	PlotRange->All,
	GridLines->{{0.018/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"22"
],



ListPlot[
	RelativeDiff@@(Extract[a, {{1,2}, {2,2}}]), 
	PlotRange->All,
	GridLines->{{0.018/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"21"
],



ListPlot[
	RelativeDiff@@(Extract[a, {{1,3}, {2,3}}]), 
	PlotRange->All,
	GridLines->{{0.018/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"33"
],


ListPlot[
	RelativeDiff@@(Extract[a, {{1,4}, {2,4}}]), 
	PlotRange->All,
	GridLines->{{0.018/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"32"
],




ListPlot[
	RelativeDiff@@(Extract[a, {{1,5}, {2,5}}]), 
	PlotRange->All,
	GridLines->{{0.018/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"44"
],



ListPlot[
	RelativeDiff@@(Extract[a, {{1,6}, {2,6}}]), 
	PlotRange->All,
	GridLines->{{0.018/(G M), 0.2/(G M)}},
	GridLinesStyle->Directive[Black, Thick],
	ScalingFunctions->"Log10",
	PlotLabel->"43"
]}


Make:= Module[
	{M = RandomReal[{5,100}],a},
	Echo[M];
	a = Test[M];
	plot[a, M]
]


Make


(* ::Chapter:: *)
(*Making Final function*)


(* ::Section:: *)
(*Symbolic expression*)


\[Phi]0//Clear


Expr = HoldComplete[
	{{
		\[Omega]RdDampinglm, \[Kappa], Erad, aeff, re\[Omega]22, im\[Omega]22, re\[Omega]21,im\[Omega]21, re\[Omega]33,im\[Omega]33, re\[Omega]32,im\[Omega]32,
		re\[Omega]44, im\[Omega]44, re\[Omega]43, im\[Omega]43, H21, H33, H32, H44, H43, f\[ScriptCapitalA], \[Psi]lm, \[ScriptA]IMR, deltat0,
		\[Omega]RdDamping, \[Omega]Peak, \[Phi]MR, re\[Omega], im\[Omega], \[Gamma]2, \[Gamma]3, \[Phi]IMR
	}},
	
	\[Gamma]2[\[Eta]_, \[Chi]s_, \[Chi]a_] := 1.010344386100769` +0.0008993122028186917` \[Eta]+(0.2839491069316864` -4.049753189086914` \[Eta]+13.207828521728516` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)+(0.10396278649568558` -7.025059223175049` \[Eta]+24.784893035888672` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^2+(0.030932024121284485` -2.6924023628234863` \[Eta]+9.609374046325684` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^3;
	\[Gamma]3[\[Eta]_, \[Chi]s_, \[Chi]a_] := 1.3081616163253784` -0.0055377297103405` \[Eta]+(-0.0678291767835617` -0.668983519077301` \[Eta]+3.4031479358673096` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)+(-0.05296577513217926` -0.9923793077468872` \[Eta]+4.820681095123291` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^2+(-0.0061341398395597935` -0.3842925429344177` \[Eta]+1.7561753988265991` \[Eta]^2) (-1+Sqrt[1-4 \[Eta]] \[Chi]a+(1-(76 \[Eta])/113) \[Chi]s)^3;
	
	\[Omega]RdDamping[int_, Erad_] := int/(1-Erad);
	\[Omega]Peak[\[Omega]RD_,\[Omega]DAMP_,\[Gamma]2_,\[Gamma]3_] := If[\[Gamma]2<=1,\[Omega]RD+(\[Omega]DAMP \[Gamma]3 (Sqrt[1-\[Gamma]2^2]-1))/\[Gamma]2,\[Omega]RD-(\[Omega]DAMP \[Gamma]3)/\[Gamma]2];
	
	\[Phi]MR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] := U\[CapitalPhi]MR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];
	\[Phi]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Rho]lm_, \[Tau]lm_] := U\[CapitalPhi]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Rho]lm, \[Tau]lm];
	
	re\[Omega][\[Kappa]_] := (0.05947169566573468` -0.14989771215394762` \[Kappa]+0.09535606290986028` \[Kappa]^2+0.02260924869042963` \[Kappa]^3-0.02501704155363241` \[Kappa]^4-0.005852438240997211` \[Kappa]^5+0.0027489038393367993` \[Kappa]^6+0.0005821983163192694` \[Kappa]^7)/(1-2.8570126619966296` \[Kappa]+2.373335413978394` \[Kappa]^2-0.6036964688511505` \[Kappa]^4+0.0873798215084077` \[Kappa]^6);
	im\[Omega][\[Kappa]_] := (0.014158792290965177` -0.036989395871554566` \[Kappa]+0.026822526296575368` \[Kappa]^2+0.0008490933750566702` \[Kappa]^3-0.004843996907020524` \[Kappa]^4-0.00014745235759327472` \[Kappa]^5+0.0001504546201236794` \[Kappa]^6)/(1-2.5900842798681376` \[Kappa]+1.8952576220623967` \[Kappa]^2-0.31416610693042507` \[Kappa]^4+0.009002719412204133` \[Kappa]^6);
	
	\[Omega]RdDampinglm[int_,Erad_] := int/(2 \[Pi] (1-Erad));
	\[Kappa][aeff_, l_, m_] := Log[3]^(-(1/(2+l-Abs[m]))) Log[2-aeff]^(1/(2+l-Abs[m]));
	Erad[\[Eta]_, \[Chi]s_, \[Chi]a_] := erad;
	aeff[\[Eta]_, \[Chi]s_, \[Chi]a_] := ae;
	
	re\[Omega]22[\[Kappa]_] := r22;
	im\[Omega]22[\[Kappa]_] := i22;
	
	re\[Omega]21[\[Kappa]_] := r21;
	im\[Omega]21[\[Kappa]_] := i21;
	
	re\[Omega]33[\[Kappa]_] := r33;
	im\[Omega]33[\[Kappa]_] := i33;
	
	re\[Omega]32[\[Kappa]_] := r32;
	im\[Omega]32[\[Kappa]_] := i32;
	
	re\[Omega]44[\[Kappa]_] := r44;
	im\[Omega]44[\[Kappa]_] := i44;
	
	re\[Omega]43[\[Kappa]_] := r43;
	im\[Omega]43[\[Kappa]_] := i43;
	
	H21[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := h21;
	H33[\[Omega]_, \[Eta]_] := 1/2 Sqrt[5/7] (3/2)^(2/3) \[Pi]^(1/3) Sqrt[1-4 \[Eta]] \[Omega]^(1/3);
	H32[\[Omega]_, \[Eta]_] := 1/3 Sqrt[5/7] \[Pi]^(2/3) (1-3 \[Eta]) \[Omega]^(2/3);
	H44[\[Omega]_, \[Eta]_] := 2/9 Sqrt[5/7] 2^(5/6) \[Pi]^(2/3) (1-3 \[Eta]) \[Omega]^(2/3);
	H43[\[Omega]_, \[Eta]_] := 1/2 Sqrt[3/35] \[Pi] Sqrt[1-4 \[Eta]] (1-2 \[Eta]) \[Omega];
	
	f22\[ScriptCapitalA][\[Omega]_, ringdownFrequency_, \[Omega]rdlm_, m_] := f22a;
	
	\[Psi]lm[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Omega]rdlm_, \[Omega]damplm_, m_] := U\[CapitalPsi]lm[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rdlm, \[Omega]damplm, m];
	\[ScriptA]IMR[\[Omega]_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[ScriptCapitalA]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a];
	deltat0[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_] := U\[CapitalDelta]t0[\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a];
	
	
	
	
	
	{
	(
		gp22 \[ScriptA]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Phi]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, 1, 1] - \[Delta]t0 - 2 \[Phi]0)] + 
		gp21 \[ScriptA]IMR[f\[ScriptCapitalA]21, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd21, \[Omega]d21, 1]- \[Delta]t0 -  \[Phi]0)]*ratio21 + 
		gp33 \[ScriptA]IMR[f\[ScriptCapitalA]33, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega], \[Eta],\[Chi]s, \[Chi]a, \[Omega]rd33, \[Omega]d33, 3]- \[Delta]t0 - 3 \[Phi]0)]*ratio33 + 
		gp32 \[ScriptA]IMR[f\[ScriptCapitalA]32, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rd32, \[Omega]d32, 2] - \[Delta]t0 - 2 \[Phi]0)]*ratio32 + 
		gp44 \[ScriptA]IMR[f\[ScriptCapitalA]44, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd44, \[Omega]d44, 4]- \[Delta]t0 - 4 \[Phi]0)]*ratio44 + 
		gp43 \[ScriptA]IMR[f\[ScriptCapitalA]43, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd43, \[Omega]d43, 3]- \[Delta]t0 - 3 \[Phi]0)]*ratio43
	), 
	(
		gc22 \[ScriptA]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Phi]IMR[\[Omega], \[Eta], \[Chi]s, \[Chi]a, 1, 1] - \[Delta]t0 - 2 \[Phi]0)] + 
		gc21 \[ScriptA]IMR[f\[ScriptCapitalA]21, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rd21, \[Omega]d21, 1]- \[Delta]t0 -  \[Phi]0)]*ratio21 + 
		gc33 \[ScriptA]IMR[f\[ScriptCapitalA]33, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega], \[Eta],\[Chi]s, \[Chi]a, \[Omega]rd33, \[Omega]d33, 3]- \[Delta]t0 - 3 \[Phi]0)]*ratio33 + 
		gc32 \[ScriptA]IMR[f\[ScriptCapitalA]32, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Omega]rd32, \[Omega]d32, 2] - \[Delta]t0 - 2 \[Phi]0)]*ratio32 + 
		gc44 \[ScriptA]IMR[f\[ScriptCapitalA]44, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd44, \[Omega]d44, 4]- \[Delta]t0 - 4 \[Phi]0)]*ratio44 + 
		gc43 \[ScriptA]IMR[f\[ScriptCapitalA]43, \[Eta], \[Chi]s, \[Chi]a] Exp[-I (\[Psi]lm[\[Omega],\[Eta], \[Chi]s, \[Chi]a, \[Omega]rd43, \[Omega]d43, 3]- \[Delta]t0 - 3 \[Phi]0)]*ratio43
	)
	}
]//.{
	erad ->Erad[\[Eta], \[Chi]s, \[Chi]a], ae -> aeff[\[Eta], \[Chi]s, \[Chi]a],
	r22->re\[Omega]22[\[Kappa]], i22->im\[Omega]22[\[Kappa]], r21->re\[Omega]21[\[Kappa]], i21->im\[Omega]21[\[Kappa]],
	r33->re\[Omega]33[\[Kappa]], i33->im\[Omega]33[\[Kappa]], r32->re\[Omega]32[\[Kappa]], i32->im\[Omega]32[\[Kappa]],
	r44->re\[Omega]44[\[Kappa]], i44->im\[Omega]44[\[Kappa]], r43->re\[Omega]43[\[Kappa]], i43->im\[Omega]43[\[Kappa]],
	h21-> H21[\[Omega],\[Eta],\[Chi]s,\[Chi]a, 1],
	f22a -> f22\[ScriptCapitalA][\[Omega], ringdownFrequency, \[Omega]rdlm, m]

};


NormalExpressionRules = {
	gp22 -> g22Plus[\[Iota]], gp21->g21Plus[\[Iota]], gp33 -> g33Plus[\[Iota]], 
	gp32->g32Plus[\[Iota]], gp44 -> g44Plus[\[Iota]], gp43->g43Plus[\[Iota]],
	
	gc22 -> g22Cross[\[Iota]], gc21->g21Cross[\[Iota]], gc33 -> g33Cross[\[Iota]], 
	gc32->g32Cross[\[Iota]], gc44 -> g44Cross[\[Iota]], gc43->g43Cross[\[Iota]],
	
	f\[ScriptCapitalA]21 :>  f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd21, 1],
	f\[ScriptCapitalA]33 :>  f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd33, 3],
	f\[ScriptCapitalA]32 :>  f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd32, 2],
	f\[ScriptCapitalA]44 :>  f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd44, 4],
	f\[ScriptCapitalA]43 :>  f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd43, 3],
	
	
	\[Omega]rd22 :> \[Omega]RdDampinglm[re\[Omega]22[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	\[Omega]d22 :> \[Omega]RdDampinglm[im\[Omega]22[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	
	\[Omega]rd21 :> \[Omega]RdDampinglm[re\[Omega]21[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 1]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	\[Omega]d21 :> \[Omega]RdDampinglm[im\[Omega]21[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 2, 1]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	
	\[Omega]rd33 :> \[Omega]RdDampinglm[re\[Omega]33[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	\[Omega]d33 :> \[Omega]RdDampinglm[im\[Omega]33[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	
	
	\[Omega]rd32 :> \[Omega]RdDampinglm[re\[Omega]32[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	\[Omega]d32 :> \[Omega]RdDampinglm[im\[Omega]32[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 3, 2]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	
	\[Omega]rd44 :> \[Omega]RdDampinglm[re\[Omega]44[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 4]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	\[Omega]d44 :> \[Omega]RdDampinglm[im\[Omega]44[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 4]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	
	\[Omega]rd43 :> \[Omega]RdDampinglm[re\[Omega]43[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	\[Omega]d43 :> \[Omega]RdDampinglm[im\[Omega]43[\[Kappa][aeff[\[Eta], \[Chi]s, \[Chi]a], 4, 3]], Erad[\[Eta], \[Chi]s, \[Chi]a]],
	
	Peak :> \[Omega]Peak[
		\[Omega]RdDamping[re\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]], 
		\[Omega]RdDamping[im\[Omega][aeff[\[Eta], \[Chi]s, \[Chi]a]], Erad[\[Eta], \[Chi]s, \[Chi]a]], 
		\[Gamma]2[\[Eta], \[Chi]s, \[Chi]a], 
		\[Gamma]3[\[Eta], \[Chi]s, \[Chi]a]
	],
	
	
	\[Delta]t0 -> Derivative[1, 0, 0,0,0,0][\[Phi]MR][
		Peak, \[Eta], \[Chi]s, \[Chi]a, 1, 1
	]*(\[Omega]-\[Omega]ref),
	
	\[Phi]0 :>  \[Phi]IMR[\[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, 1, 1]/2 + \[Phi]Ref,
	
	ratio21 :> H21[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd21, 1], \[Eta], \[Chi]s, \[Chi]a]*H21[\[Omega], \[Eta], \[Chi]s, \[Chi]a]/H21[2 \[Omega], \[Eta], \[Chi]s, \[Chi]a],
	ratio33 :>H33[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22,  \[Omega]rd33, 3], \[Eta]]*(3/2)^(1/3),
	ratio32 :> (H32[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22,\[Omega]rd32, 2], \[Eta]]),
	ratio44 :> H44[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22, \[Omega]rd44, 4], \[Eta]]*2^(2/3),
	ratio43:> H43[f22\[ScriptCapitalA][\[Omega], \[Omega]rd22,  \[Omega]rd43, 3], \[Eta]]*3/2
};


Expr = Expr//.NormalExpressionRules;
Expr = $Block@@@(HoldForm[Evaluate@Expr]);
ClearAll[hphcres3]


HMhphc//ClearAll

expr = Hold[
	{HMhphc, {\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref}, 3},
	Evaluate@Expr
]//.HoldForm[X_] :> X;


hphcres = EchoTiming[DerivativeRules@@expr];


hphcres2 = hphcres;
Do[
	
	hphcres2[[i,2]]  = hphcres[[i,2]]/.List[x1_, x2_] :> Transpose[List[x1, x2]],
	{i, Length@hphcres}
]


Module[
	{rule, temp},
	
	rule = {Derivative[n__][s_][x__] :> $D[{n}, s][x]};
	
	
	temp = hphcres2//.rule;
	
	rule =  {
		U\[CapitalPhi]IMR[\[Omega]ref, y__] :> Last[U\[CapitalPhi]IMR[{\[Omega]ref}, y]],
		$D[{n__}, U\[CapitalPhi]IMR][\[Omega]ref, y__] :> Last[$D[{n},  U\[CapitalPhi]IMR][{\[Omega]ref},y]],
		
		$D[{n__}, U\[CapitalPhi]MR][y1_, y__]/; y1[[0]]=!=List :> Last[$D[{n},  U\[CapitalPhi]MR][{y1},y]]  
	};
	
	temp = temp//.rule;
	
	
	hphcres3 = temp//.{U\[ScriptCapitalA]IMR->HM\[ScriptCapitalA]IMR, U\[CapitalPsi]lm->HM\[CapitalPsi]lm, U\[CapitalDelta]t0->\[CapitalDelta]t0, U\[CapitalPhi]MR->HM\[CapitalPhi]MR, U\[CapitalPhi]IMR->HM\[CapitalPhi]IMR};
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] = hphcres3[[i,2]];
	
	DownValues[test] = DownValues[test]/.{
		HoldForm[X_] :> X, \[Omega]$->\[Omega], \[Eta]$->\[Eta], \[Chi]s$->\[Chi]s, \[Chi]a$->\[Chi]a, \[Omega]ref$->\[Omega]ref, \[Iota]$->\[Iota], \[Phi]Ref$->\[Phi]Ref
	};
	
	test[Range[10^-4, 0.124, 10^-4], 0.02, 0.1, 0.2, 0.12, 1., 1.]
]


Table[T[i], {i, Length@hphcres}]//AbsoluteTiming


(*Function to take HoldForm[Block[...]] and make library functions*)
compileThis//ClearAll

compileThis[x_HoldForm] := Module[
	{dummy},
	dummy = Hold[
		{{\[Omega], _Real, 1}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref}, 
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


chphcres = MapAt[
	compileThis,
	hphcres3,
	{All,2}
];


Do[
	chphcres[[i]][[2, -2]] = Function[{\[Omega], \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref}, X]//.{X -> ToString[chphcres[[i, 1]][[All,1]]]},
	{i, Length@chphcres}
]


Clear@T

T[i_] := Block[
	{test},
	
	test[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] := chphcres[[i, 2]][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref];
	
	test[Range[10^-4, 0.124, 10^-4], 10^-4, 0.2, 0.1, 0.2, 1., 1.]
]


(*All or3 gradients included*)
Table[T[i], {i, Length@chphcres}]//AbsoluteTiming


(* ::Text:: *)
(*Set UpValues:*)


Module[
	{list = chphcres[[2;;-1]]},
	
	list[[All,1]] = list[[All,1, 0]];
	
	list = list//.Rule->RuleDelayed;
	list = MapAt[
		HoldPattern, 
		list, 
		{All,1}
	];
	
	UpValues@HMhphc = list;
]


(* ::Text:: *)
(*Set DownValues:*)


HMhphc[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] = chphcres[[1, 2]][\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref];


(* ::Section:: *)
(*Testing*)


ClearAll[Testhphc];
Testhphc[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] = HMhphc[\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref];


(* ::Subsection::Closed:: *)
(*Testing hphc*)


(* ::Subsubsection::Closed:: *)
(*My def*)


Clear@MMA

MMA[\[ScriptCapitalM]c_, \[Eta]_, \[Chi]1_, \[Chi]2_, \[Iota]_, \[Phi]Ref_] := Module[
	{
		M, \[Omega], \[Chi]s, \[Chi]a, freq, \[Omega]ref, res, G =  4.925490947641267`*^-6
	},
	
	\[Chi]s = (\[Chi]1+\[Chi]2)/2; 
	\[Chi]a = (\[Chi]1-\[Chi]2)/2;
	
	
	freq =Range[1., 1024., 1];
	
	M = \[ScriptCapitalM]c \[Eta]^(-3/5);
	
	\[Omega] = M freq G;
	
	\[Omega]ref = Min[\[Omega]];
	
	
	res =  M^2 Testhphc[\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref]//Transpose;
	
	res
]


(* ::Subsubsection::Closed:: *)
(*lal hphc function*)


Clear[python]
DeleteObject/@ExternalSessions[];

python = StartExternalSession[{
	"Python", (*you should change the evaluator to your installation:*)
	(*"Evaluator"-> "/home/cosmo-ufes/anaconda3/envs/gwfast_env/bin/python"*)
	"Evaluator" -> "/Users/felipe/anaconda3/envs/GWFAST/bin/python"
}];
ExternalEvaluate[python,"
import numpy as np

import jax
import jax.numpy as jnp


from jax import grad, vmap
from functools import partial


import lalsimulation as lalsim
import lal

jax.config.update(\"jax_enable_x64\", True)
"]


lal = ExternalFunction[python, "
def lal_hp(m1, m2, s1z, s2z, dL, iota, phiRef, deltaF,
           f_min, f_max, f_ref, appr):

    approximant = lalsim.SimInspiralGetApproximantFromString(appr)
    
    dL_m = dL * 3.08568 * (10**22)  # convert megaparsecs to meters

    # solar mass to Kg
    m1_kg = m1 * 1.9884 * (10**30)
    m2_kg = m2 * 1.9884 * (10**30)
    
    hp, hc = lalsim.SimInspiralChooseFDWaveform(
        m1_kg,
        m2_kg,
        0,#s1x,
        0,#s1y,
        s1z,
        0,#s2x,
        0,#s2y,
        s2z,
        dL_m,
        iota,
        phiRef,
        0,
        0,
        0,
        deltaF,
        f_min,
        f_max,
        f_ref,
        None,
        approximant,
    )

    hP = hp.data.data
    hC = hc.data.data    

    return [hP.tolist(), hC.tolist()]

"]


(* ::Subsubsection::Closed:: *)
(*Josiel def*)


(*DeleteObject/@ExternalSessions[]*)
Clear@python2

python2 = StartExternalSession[{
	"Python", (*you should change the evaluator to your installation:*)
	(*"Evaluator"-> "/home/cosmo-ufes/anaconda3/envs/gwfast_env/bin/python"*)
	"Evaluator" -> "/Users/felipe/anaconda3/envs/josiel/bin/python"
}];


Range[1, 1024]//Length


ExternalEvaluate[python2,"
import numpy as np
from GWDALI_v1.lib.IMRPhenomHM import hphx_IMRPhenomHM

freq = np.linspace(1, 1024, 1024)
"]


josiel = ExternalFunction[python2,"
def gwDALI_phase(Mc, eta, chi1, chi2, iota):
    phi0=0
    t0=0
    sx1=0
    sy1=0
    sz1=chi1

    sx2=0
    sy2=0
    sz2=chi2

    dL=1

    hp, hx = hphx_IMRPhenomHM(
        dL,iota,phi0,t0,
        Mc,eta,sx1,sy1,sz1,sx2,sy2,sz2,
        freq
    )

    
    

    
    return np.array([
        hp, hx      
    ])
"]


(* ::Subsubsection::Closed:: *)
(*Test*)


(* ::Text:: *)
(*I believe Josiel is using frd and fdamp as calculated in the original PhenomD implementation. I am not convinced this is right. I find it more likely that those should be replaced by the frd and fdamp 22 formulas attached to Z22.*)


Clear@Test
Test[M_] := Module[
	
	{m1, m2, mc, eta, s1z, s2z, \[Iota], \[Phi]Ref, lalRes, fmax, fmin, myDef},
	
	
	
	eta = RandomReal[{0.1, 0.25}];
	m1 = M/2 (1 + Sqrt[1-4 eta]);
	m2 = M/2 (1 - Sqrt[1-4 eta]);
	
	
	mc = (m1+m2) eta^(3/5);
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	
	\[Phi]Ref = 0; 
	
	mc = M eta^(3/5);
	
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	fmax =  1024;
	fmin=1;
	
	lalRes = lal[m1, m2, s1z, s2z, 10^3, \[Iota], \[Phi]Ref, 1, fmin, fmax, fmin, "IMRPhenomHM"]//Normal;
	myDef = hphcIMRPhenomHM[
		Range[1, 1024], 
		mc, 
		Sqrt[1-4 eta], (s1z+s2z)/2, (s1z-s2z)/2, \[Iota], 0, \[Phi]Ref, 
		1
	];
	
	
	{
		
		josiel[mc, eta, s1z, s2z, \[Iota]]//Normal,
		{lalRes[[1, 2;;1025]], lalRes[[2, 2;;1025]]},
		(*MMA[mc, eta, s1z, s2z, \[Iota], \[Phi]Ref]*)
		myDef
	}
	
]


plot[a_, M_] := Module[
	{fcutoff},
	fcutoff = 0.1/(G M);
(*all hp comparisons:*)
	{
(*josielvslal*)CompareComplexVectors[Extract[a, {{1,1}, {2,1}}], {"josiel", "lal"}, Range[1, 1024, 1], {1, fcutoff}],
(*josielvsmine*)CompareComplexVectors[Extract[a, {{1,1}, {3,1}}], {"josiel", "MMA"}, Range[1, 1024, 1], {1, fcutoff}]
(*lalvsmine*)CompareComplexVectors[Extract[a, {{2,1}, {3,1}}], {"lal", "MMA"}, Range[1, 1024, 1], {1, fcutoff}]
	}
]


Make:= Module[
	{M = RandomReal[{50,100(*450*)}],a},
	Echo[M];
	Echo[0.014/(G M)];
	a = Test[M];
	plot[a, M]
]


(* ::Text:: *)
(*Josiel's code does not agree perfectly with lal suite either.*)


Make


Clear@Test
Test[M_] := Module[
	
	{m1, m2, mc, eta, s1z, s2z, \[Iota], \[Phi]Ref, lalRes, fmax, fmin, myDef},
	
	
	
	eta = RandomReal[{0.1, 0.25}];
	m1 = M/2 (1 + Sqrt[1-4 eta]);
	m2 = M/2 (1 - Sqrt[1-4 eta]);
	
	
	mc = (m1+m2) eta^(3/5);
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}]; 
	
	mc = M eta^(3/5);
	
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	fmax =  1024;
	fmin=1;
	
	lalRes = lal[m1, m2, s1z, s2z, 10^3, \[Iota], \[Phi]Ref, 1, fmin, fmax, fmin, "IMRPhenomHM"]//Normal;
	
	
	{
		{lalRes[[1, 2;;1025]], lalRes[[2, 2;;1025]]},
		MMA[mc, eta, s1z, s2z, \[Iota], \[Phi]Ref]
	}
	
]


plot[a_, M_] := Module[
	{fcutoff},
	fcutoff = 0.2/(G M);
(*all hp comparisons:*)
	{
	CompareComplexVectors[Extract[a, {{1,1}, {2,1}}], {"hp_lal", "MMA"}, Range[1, 1024, 1], {1, fcutoff}],
	CompareComplexVectors[Extract[a, {{1,2}, {2,2}}], {"hc_lal", "MMA"}, Range[1, 1024, 1], {1, fcutoff}]
	}
]


(* ::Text:: *)
(*Check agreement on the mass range associated to the fiducial point of event 12:*)


{225.9445206399997, 211.44219146000015}//Total


Make:= Module[
	{M = RandomReal[{300,500}],a},
	Echo[M];
	Echo[0.014/(G M)];
	a = Test[M];
	plot[a, M]
]


(* ::Text:: *)
(*Josiel's code does not agree perfectly with lal suite either.*)


Extract[Make, {{1,1,1}, {2,1,1}}]


(* ::Subsection:: *)
(*Testing derivatives*)


(* ::Subsubsection:: *)
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


ClearAll[TestGradhphc];

TestGradhphc[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] = hphcres3[[2;;8, 2]];


DownValues[TestGradhphc] = DownValues[TestGradhphc]//.HoldForm[x_]:> x;


iTesthphc[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] := Testhphc[{\[Omega]}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref];


Clear@Test

Test := Module[
	{
		 M, \[Omega],\[Omega]ref, \[Eta],\[Iota], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		 hpSym, hpN, hcSym, hcN,\[Chi]1, \[Chi]2,\[Phi]Ref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars
	},
	M = RandomReal[{5,100}];
	\[Eta] = RandomReal[{0.1, 0.249}];
	
	
	
	{\[Chi]1,\[Chi]2} = RandomReal[{-0.9,0.9}, 2];
	\[Chi]s = (\[Chi]1+\[Chi]2)/2; \[Chi]a = (\[Chi]1-\[Chi]2)/2;
	
	\[Omega] = RandomReal[{1 G M, 0.1}];
	\[Omega]ref = 1 G M;
	\[Iota] = RandomReal[{0,\[Pi]}];
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	
	Symbolic = TestGradhphc[{\[Omega]}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref];
	Symbolic = ArrayReshape[Symbolic, {7,2}];
	vars = {\[Omega], \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref};
	
	Numeric = NGrad[iTesthphc, vars, 0];
	Numeric = ArrayReshape[Numeric, {7,2}];
	
	hpSym = Symbolic[[All,1]];
	hpN = Numeric[[All,1]];
	
	hcSym = Symbolic[[All,2]];
	hcN = Numeric[[All,2]];
	
	{
		RelativeDiff@@{hpSym, hpN},
		RelativeDiff@@{hcSym, hcN}
	}
]


Test


Table[Round@Test, {100}]


(* ::Subsubsection:: *)
(*Second Order*)


(* ::Text:: *)
(*Since First order checks out we can Use the symbolic gradient of first order and compare its numerical derivatives against the symbolic *)
(*gradients of order 2.*)


ClearAll[TestGrad\[ScriptCapitalA]O2];

TestGrad\[ScriptCapitalA]O2[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] = hphcres3[[9;;36, 2]];

DownValues[TestGrad\[ScriptCapitalA]O2] = DownValues[TestGrad\[ScriptCapitalA]O2]//.HoldForm[x_]:> x;


iTestGradhphc[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_]  := TestGradhphc[{\[Omega]}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref] 


Clear@Test

Test := Module[
	{
		 M, chi1, chi2, \[Eta], \[Iota], \[Omega], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,fref,
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], vars,
		li = SymmetrizedIndependentComponents[{7,7}, Symmetric[All]],
		hpS, hcS, hpN, hcN, \[Phi]Ref
	},
	
	M = RandomReal[{5,100}];
	\[Eta] = RandomReal[{0.1, 0.24}];
	
	{chi1,chi2} = RandomReal[{-1,1}, 2];
	\[Chi]s = (chi1+chi2)/2; \[Chi]a = (chi1-chi2)/2;
	
	
	\[Iota] = RandomReal[{0,\[Pi]}];
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	
	\[Omega] = RandomReal[{1, 1024}] G M;
	
	Symbolic = ArrayReshape[TestGrad\[ScriptCapitalA]O2[{\[Omega]}, G M, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref], {28,2}];
	
	hpS = Symbolic[[All,1]];
	hcS = Symbolic[[All,2]];
	
	vars = {\[Omega], G M, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref};
	
	Numeric = ArrayReshape[NGrad[iTestGradhphc, vars, 0], {7, 7, 2}];
	
	hpN = Extract[
		Numeric[[All,All,1]], 
		li
	];
	
	hcN = Extract[
		Numeric[[All,All,2]], 
		li
	];
	
	
	{
		RelativeDiff@@{hpN, hpS},
		RelativeDiff@@{hcN, hcS}
	}
	
	
]


Test//ScientificForm


Table[Round@Test, {100}]


(* ::Subsubsection:: *)
(*Third Order*)


(* ::Text:: *)
(*Second order checks out so we can use the symbolic gradient of second order and compare its numerical derivatives against the symbolic gradients of order 3.*)


SymmetricTestGradOrder2[x__] := Module[
	{
		li = SymmetrizedIndependentComponents[{7, 7, 2}, Symmetric[{1,2}]], 
		rule, symbolicO2,
		(*values = TestGrad\[ScriptCapitalA]O2[x][[All, 1]], *)
		rules
	},
	symbolicO2 = TestGrad\[ScriptCapitalA]O2[x];
	
	symbolicO2 = ArrayReshape[symbolicO2, {28, 2}];
	
	rules = MapThread[Rule, {li, Flatten[symbolicO2]}];
	
	SymmetrizedArray[rules, {7, 7, 2}, Symmetric[{1,2}]]
]


iSymmetricTestGradOrder2[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] := SymmetricTestGradOrder2[{\[Omega]}, \[Omega]ref, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref]


ClearAll[TestGrad\[ScriptCapitalA]O3];

TestGrad\[ScriptCapitalA]O3[\[Omega]_, \[Omega]ref_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, \[Phi]Ref_] = hphcres3[[37;;-1, 2]];

DownValues[TestGrad\[ScriptCapitalA]O3] = DownValues[TestGrad\[ScriptCapitalA]O3]//.HoldForm[x_]:> x;


hphcres3[[37;;-1, 1]]//Length


Clear@Test

Test := Module[
	{
		chi1, chi2, m1, m2, \[Chi]1, \[Chi]2, f,\[Omega], \[Eta], Symbolic,Numeric, RRe, MMARe, RIm, MMAIm, \[Chi]s, \[Chi]a,diff,
		li = SymmetrizedIndependentComponents[{7,7,7}, Symmetric[All]], M, \[ScriptCapitalM]c, \[Iota],
		vars, \[Phi]Ref
	},
	
	\[Eta] = RandomReal[{0.1, 0.249}];
	
	M  = RandomReal[{5, 100}];
	
	
	{chi1, chi2} = RandomReal[{-1,1}, 2];
	\[Chi]s = (chi1+chi2)/2; \[Chi]a = (chi1-chi2)/2;
	
	\[Iota] = RandomReal[{0,\[Pi]}];
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	
	
	
	\[Omega] = RandomReal[{1, 1024}] G M;
	
	Symbolic = TestGrad\[ScriptCapitalA]O3[{\[Omega]}, G M, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref];
	
	Symbolic= ArrayReshape[Symbolic, {84, 2}];
	
	
	
	vars = {\[Omega], G M, \[Eta], \[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref};
	
	

	
	Numeric = Extract[NGrad[iSymmetricTestGradOrder2, vars, 0], li];
	
	{
		RelativeDiff@@{Symbolic[[All,1]], Numeric[[All,1]]}, (*hp*)
		RelativeDiff@@{Symbolic[[All,2]], Numeric[[All,2]]} (*hc*)
	}
]


Round@Test//ScientificForm


Table[Round@Test, {10}]//EchoTiming


list = %;


(* ::Section:: *)
(*Saving defs*)


(* ::Subsection:: *)
(*New*)


(* ::Text:: *)
(*THIS ADDS SIDE EFFECTS:*)


"FelipeBarbosa`SymDALI`DALICoefficients`Private`" <>ToString[HM\[ScriptCapitalA]Ins]//ToExpression


DownValues[HMhphc]


newVars = {\[Omega], \[Omega]ref, \[Eta],\[Chi]s, \[Chi]a, \[Iota], \[Phi]Ref, \[Rho]lm,\[Tau]lm, \[Omega]rdlm,\[Omega]damplm,m}//.x_Symbol/; Context[x] ==="Global`" :> ToExpression[
	"FelipeBarbosa`SymDALI`DALICoefficients`Private`"<>ToString[x]
]


newHeads = {
			HMZeroFunction, 
			HM\[ScriptCapitalA]Ins, HM\[ScriptCapitalA]Int, HM\[ScriptCapitalA]MR, HM\[ScriptCapitalA]IMR,
			HM\[CapitalPhi]Ins, HM\[CapitalPhi]Int, HM\[CapitalPhi]MR, HM\[CapitalPhi]IMR, HM\[CapitalPsi]lm, 
			HMhphc
}//.x_Symbol/; Context[x] ==="Global`" :> ToExpression[
	"FelipeBarbosa`SymDALI`DALICoefficients`Private`"<>ToString[x]
]


varRules = MapThread[
	Rule, 
	{{\[Omega], \[Omega]ref, \[Eta],\[Chi]s,\[Chi]a, \[Iota], \[Phi]Ref, \[Rho]lm,\[Tau]lm, \[Omega]rdlm,\[Omega]damplm,m}, newVars}
]


HeadRules = MapThread[
	Rule,
	{{
			HMZeroFunction, 
			HM\[ScriptCapitalA]Ins, HM\[ScriptCapitalA]Int, HM\[ScriptCapitalA]MR, HM\[ScriptCapitalA]IMR,
			HM\[CapitalPhi]Ins, HM\[CapitalPhi]Int, HM\[CapitalPhi]MR, HM\[CapitalPhi]IMR, HM\[CapitalPsi]lm, 
			HMhphc
		}, newHeads}
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


ChangeContext[HM\[ScriptCapitalA]Ins]


ChangeContext/@{
			HMZeroFunction, 
			HM\[ScriptCapitalA]Ins, HM\[ScriptCapitalA]Int, HM\[ScriptCapitalA]MR, HM\[ScriptCapitalA]IMR,
			HM\[CapitalPhi]Ins, HM\[CapitalPhi]Int, HM\[CapitalPhi]MR, HM\[CapitalPhi]IMR, HM\[CapitalPsi]lm, 
			HMhphc
}


Module[
	{SymDALIDir = NotebookDirectory[]//ParentDirectory[#, 2]&, fileName, file},
	
	
	fileName = FileNameJoin[{
		SymDALIDir, 
		"LibraryResources",
		$SystemID,
		"DerivativeRules/IMRPhenomHM/MMA/Defs.mx"
	}];
	
	file = CreateFile[fileName];
	
	DumpSave[
		file,
		{
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HMZeroFunction,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[ScriptCapitalA]Ins,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[ScriptCapitalA]Int,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[ScriptCapitalA]MR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[ScriptCapitalA]IMR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[CapitalPhi]Ins,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[CapitalPhi]Int,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[CapitalPhi]MR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[CapitalPhi]IMR,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HM\[CapitalPsi]lm,
			FelipeBarbosa`SymDALI`DALICoefficients`Private`HMhphc
		}
	]
	
]


DownValues@FelipeBarbosa`SymDALI`DALICoefficients`Private`HMZeroFunction


FelipeBarbosa`SymDALI`DALICoefficients`Private`HMZeroFunction//ClearAll


Get["/Users/felipe/Documents/GitHub/SymDALI/LibraryResources/MacOSX-ARM64/DerivativeRules/IMRPhenomHM/MMA/Defs.mx"]


DownValues@FelipeBarbosa`SymDALI`DALICoefficients`Private`HMZeroFunction
