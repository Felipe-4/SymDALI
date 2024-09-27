(* ::Package:: *)

(*SetOptions[EvaluationNotebook[], DefaultNewCellStyle->"Code"];
SetOptions[EvaluationNotebook[], WindowElements->{"MemoryMonitor","VerticalScrollBar","MenuBar", "HorizontalScrollBar"}]
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


(* ::Section:: *)
(*Main*)


BeginPackage["FelipeBarbosa`SymDALI`IMRPhenomD`"]


PhenomCoeff
Waveform


Begin["`Private`"]


G = (Quantity[("GravitationalConstant")/("SpeedOfLight")^3]//UnitConvert[#, ("Seconds")/("SolarMass")]&)[[1]];
c = (Quantity["SpeedOfLight"]//UnitConvert[#,("Megaparsecs")/("Seconds")]&)[[1]]//N;


(* ::Section::Closed:: *)
(*Helper Functions*)


PhenomCoeff[\[Eta]_, \[Chi]Pn_, \[Lambda]_?VectorQ] := With[
    {vec = {1, \[Eta], \[Eta]^2}, diff = \[Chi]Pn - 1},
    
    \[Lambda][[1;;2]] . vec[[1;;2]] + diff (\[Lambda][[3;;5]] . vec) + diff^2 (\[Lambda][[6;;8]] . vec) + diff^3 (\[Lambda][[9;;11]] . vec)
];

PheomCoeff[args___] := Throw[$Failed, failTag[PhenomCoeff]]


Erad[\[Eta]_ ,S_] := With[{
	EradNS= 0.0559745 \[Eta] + 0.580951 \[Eta]^2-0.960673 \[Eta]^3+3.35241 \[Eta]^4},
	
	(EradNS (1+ (-0.00303023` -2.00661` \[Eta]+7.70506` \[Eta]^2) S))/(1+(-0.67144` -1.47569` \[Eta]+7.30468` \[Eta]^2) S)];

Erad[args___] := Throw[$Failed, failTag[PhenomCoeff]]

aeff[\[Eta]_, S_]= S + 2 Sqrt[3.] \[Eta]+(-0.085` S+0.102` S^2-1.355` S^3-0.868` S^4) \[Eta]-4.399` \[Eta]^2+(-5.837` S-2.097` S^2+4.109` S^3+2.064` S^4) \[Eta]^2+9.397` \[Eta]^3-13.181` \[Eta]^4;


FDAMP[x_] = 0.014055787148647485` (1 - E^(8.150574465175337` (-1.0479012727408963`+x)));
FRD[x_] = 0.059184339277624086` +E^(4.494993811335756` (-1.625570103110631`+x))+0.014137205950116518` x;


SpecialFrequencies[\[Eta]_,  S_, \[Gamma]2_, \[Gamma]3_]  := Module[
    {fRD, fDamp, \[Lambda]vecs, fPeak},
    
    fRD = FRD[aeff[\[Eta], S]]/(1 - Erad[\[Eta], S/(1 - 2 \[Eta]) ]);
    
    fDamp = FDAMP[aeff[\[Eta], S]]/(1- Erad[\[Eta], S/(1 - 2 \[Eta])]);
   
   fPeak = Abs[fRD + (fDamp \[Gamma]3 (Sqrt[1- \[Gamma]2^2]-1))/\[Gamma]2]; (*Ripple does not use this for \[Gamma]2>1, why?*)
   
   
   {fPeak, fRD, fDamp} (*Dimensionless frequencies*)
]


SetAttributes[unitStep,Listable]

unitStep[x_/;x<0] :=0
unitStep[x_/;x==0] := 1/2
unitStep[x_/;x>0] :=1


(* ::Section::Closed:: *)
(*Phase*)


UncorrectedMRPhase[
	f_/;(VectorQ@f || NumberQ@f),
	fRD_?NumericQ,
	fdamp_?NumericQ, 
	\[Eta]_?NumericQ, 
	\[Alpha]_/;(VectorQ@\[Alpha]&&Length@\[Alpha]==5)] := \[Alpha][[1;;4]] . {f,-(1/f),(4 f^(3/4))/3, ArcTan[(f- \[Alpha][[-1]] fRD)/fdamp]}/\[Eta];


UncorrectedIntPhase[
	f_/;(VectorQ@f || NumberQ@f),
	\[Eta]_?NumberQ,
	\[Beta]_/;(VectorQ@\[Beta] && Length@\[Beta] == 3)] := \[Beta] . {f,Log[f],-(1/(3 f^3))}/\[Eta];


(*Function of the adimensional frequency*)

PhasePNCoefficients[\[Eta]_?NumericQ, \[Chi]s_?NumericQ, \[Chi]a_?NumericQ] := With[ (*Check elements by varying i here: Extract[DownValues[PNCoefficients], {1,2,2,i}, Hold]*)
    {\[Delta] = Sqrt[1 - 4 \[Eta]]},
    (*DEFINED WITHOUT FREQUENCY CONTRIBUTIONS BECAUSE OF VECTOR CONTRACTION*)
    List[
        1,
        0,
        3715/756 + (55 \[Eta])/9,
        -16 \[Pi] + (113 \[Delta] \[Chi]a)/3 + (113/3 - (76 \[Eta])/3) \[Chi]s,
        15293365/508032 + (27145 \[Eta])/504 + (3085 \[Eta]^2)/72 + (-405/8 + 200 \[Eta]) \[Chi]a^2 - (405 \[Delta] \[Chi]a \[Chi]s)/4 + (-405/8 + (5 \[Eta])/2) \[Chi]s^2,
        (38645 \[Pi]/756 - (65 \[Pi] \[Eta])/9 + \[Delta] (-732985/2268 - (140 \[Eta])/9) \[Chi]a + (-732985/2268 + (24260 \[Eta])/81 + (340 \[Eta]^2)/9) \[Chi]s),
        11583231236531/4694215680 - (6848 EulerGamma)/21 - (640 \[Pi]^2)/3 + (-15737765635/3048192 + (2255 \[Pi]^2)/12) \[Eta] + (76055 \[Eta]^2)/1728 - (127825 \[Eta]^3)/1296 + (2270 \[Pi] \[Delta] \[Chi]a)/3 + (2270 \[Pi]/3 - 520 \[Pi] \[Eta]) \[Chi]s,
        77096675 \[Pi]/254016 + (378515 \[Pi] \[Eta])/1512 - (74045 \[Pi] \[Eta]^2)/756 + \[Delta] (-25150083775/3048192 + (26804935 \[Eta])/6048 - (1985 \[Eta]^2)/48) \[Chi]a + (-25150083775/3048192 + (10566655595 \[Eta])/762048 - (1042165 \[Eta]^2)/3024 + (5345 \[Eta]^3)/36) \[Chi]s
    ]
];

InspiralPhase[f_?VectorQ, \[Eta]_?NumberQ, \[Chi]s_?NumberQ, \[Chi]a_?NumberQ, \[Sigma]_/;(VectorQ@\[Sigma] && Length@\[Sigma]==4)] := Module[
    {
        TaylorF2 = 3/(128 \[Eta]) PhasePNCoefficients[\[Eta], \[Chi]s, \[Chi]a] . {1/(f^(5/3) \[Pi]^(5/3)),1/(f^(4/3) \[Pi]^(4/3)),1/(f \[Pi]),1/(f^(2/3) \[Pi]^(2/3)),1/(f^(1/3) \[Pi]^(1/3)),ConstantArray[1, Length@f],f^(1/3) \[Pi]^(1/3),f^(2/3) \[Pi]^(2/3)}
    },
    TaylorF2 += 3/(128 \[Eta]) (PhasePNCoefficients[\[Eta], \[Chi]s, \[Chi]a][[6]]*Log[\[Pi] f] - (6848 Log[64 \[Pi] f])/63  (\[Pi] f)^(1/3) );  
     
     -\[Pi]/4  + TaylorF2 + \[Sigma] . {f, 3/4 f^(4/3), 3/5 f^(5/3), 1/2 f^2}/\[Eta]
];


(*This is intend to be used by \[Delta]\[Beta]0 ONLY*)

auxiliarInspiralPhase[f_?NumberQ, \[Eta]_?NumberQ, \[Chi]s_?NumberQ, \[Chi]a_?NumberQ,  \[Sigma]_/;(VectorQ@\[Sigma] && Length@\[Sigma]==4)] := Module[
    {
        TaylorF2 = 3/(128 \[Eta]) PhasePNCoefficients[\[Eta], \[Chi]s, \[Chi]a] . {1/(f^(5/3) \[Pi]^(5/3)),1/(f^(4/3) \[Pi]^(4/3)),1/(f \[Pi]),1/(f^(2/3) \[Pi]^(2/3)),1/(f^(1/3) \[Pi]^(1/3)), 1, f^(1/3) \[Pi]^(1/3),f^(2/3) \[Pi]^(2/3)}
    },
    TaylorF2 += 3/(128 \[Eta]) (PhasePNCoefficients[\[Eta], \[Chi]s, \[Chi]a][[6]]*Log[\[Pi] f] - (6848 Log[64 \[Pi] f])/63  (\[Pi] f)^(1/3));  
     
     -\[Pi]/4  + TaylorF2 + \[Sigma] . {f, 3/4 f^(4/3), 3/5 f^(5/3), 1/2 f^2}/\[Eta]
];


PNLOG5[\[Eta]_?NumericQ,\[Chi]s_?NumericQ, \[Chi]a_?NumericQ] = (38645 \[Pi]/756 - (65 \[Pi] \[Eta])/9 + Sqrt[1-4 \[Eta]] (-732985/2268 - (140 \[Eta])/9) \[Chi]a + (-732985/2268 + (24260 \[Eta])/81 + (340 \[Eta]^2)/9) \[Chi]s); (*DEFINED APART OF THE LOG[\[Pi] f] term*)
PNLOG6 = (6848)/63; (*DEFINED APART OF THE LOG[64 \[Pi] f] term*)


\[Delta]\[Beta]1[\[Eta]_?NumericQ, \[Chi]s_?NumericQ, \[Chi]a_?NumericQ, \[Beta]_/;(VectorQ@\[Beta]&&Length@\[Beta] == 3), \[Sigma]_/;(VectorQ@\[Sigma] && Length@\[Sigma]==4)] := Module[
	{dInspiralTaylor, dIntermediate, dInspiralPhenom},
	
	dInspiralTaylor = 3/(128 \[Eta]) (
		PhasePNCoefficients[\[Eta], \[Chi]s, \[Chi]a] . {-(5/(3 f^(8/3) \[Pi]^(5/3))),-(4/(3 f^(7/3) \[Pi]^(4/3))),-(1/(f^2 \[Pi])),-(2/(3 f^(5/3) \[Pi]^(2/3))),-(1/(3 f^(4/3) \[Pi]^(1/3))), 0, \[Pi]^(1/3)/(3 f^(2/3)),(2 \[Pi]^(2/3))/(3 f^(1/3))} - 
		 PNLOG6 Log[64 \[Pi] 18/1000]  \[Pi]^(1/3)/(3 f^(2/3)) + 
		 PNLOG5[\[Eta], \[Chi]s, \[Chi]a] (18/1000)^-1 - PNLOG6 (18/1000)^-1 (\[Pi] 18/1000)^(1/3)
	)/.f->18/1000; 
		
	dInspiralPhenom = \[Eta]^-1 \[Sigma] . {1, f^(1/3), f^(2/3), f}/.f->18/1000;
	
	dIntermediate= (\[Eta]^-1 \[Beta] . {1, 1/f, f^-4})/.f->18/1000;
	
	dInspiralTaylor  +  dInspiralPhenom  - dIntermediate 
];

\[Delta]\[Beta]0[\[Eta]_?NumericQ, \[Chi]s_?NumericQ, \[Chi]a_?NumericQ, \[Beta]_/;(VectorQ@\[Beta]&&Length@\[Beta] == 3), \[Sigma]_/;(VectorQ@\[Sigma] && Length@\[Sigma]==4)] := Module[
	{Inspiral, Intermediate, \[Delta]\[Beta]1correction},
	Inspiral = auxiliarInspiralPhase[18/1000, \[Eta], \[Chi]s, \[Chi]a, \[Sigma]];
	Intermediate = UncorrectedIntPhase[18/1000, \[Eta], \[Beta]];
	\[Delta]\[Beta]1correction = \[Delta]\[Beta]1[\[Eta], \[Chi]s, \[Chi]a, \[Beta], \[Sigma]] 18/1000;
	
	Inspiral - Intermediate - \[Delta]\[Beta]1correction
]

\[Delta]\[Alpha]1[fRd_, fdamp_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Alpha]_, \[Beta]_, \[Sigma]_]/;(
And@@(NumberQ/@{fRd, fdamp, \[Eta], \[Chi]s, \[Chi]a}) && And@@(VectorQ/@{\[Alpha], \[Beta], \[Sigma]}) && Length/@{\[Alpha], \[Beta], \[Sigma]} == {5,3,4}
):= Module[
	{MR, Intermediate, \[Beta]1correction},
	
	MR = (\[Eta]^-1 \[Alpha][[1;;4]] . {1, f^-2, f^(-1/4), fdamp/(fdamp^2+ (f - \[Alpha][[-1]] fRd)^2)})/.f->(fRd/2);
	
	Intermediate = (\[Eta]^-1 \[Beta] . {1, 1/f, f^-4})/.f->fRd/2;
	
	\[Beta]1correction = \[Delta]\[Beta]1[\[Eta], \[Chi]s, \[Chi]a, \[Beta], \[Sigma]];
	
	-MR + Intermediate + \[Beta]1correction
]


\[Delta]\[Alpha]0[fRd_, fdamp_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Alpha]_, \[Beta]_, \[Sigma]_]/;(
And@@(NumberQ/@{fRd, fdamp, \[Eta], \[Chi]s, \[Chi]a}) && And@@(VectorQ/@{\[Alpha], \[Beta], \[Sigma]}) && Length/@{\[Alpha], \[Beta], \[Sigma]} == {5,3,4}
) := Module[
	{MR, Intermediate, \[Beta]1correction, \[Beta]0correction, \[Alpha]1correction},
	
	MR = UncorrectedMRPhase[fRd/2, fRd, fdamp, \[Eta], \[Alpha]];
	Intermediate = UncorrectedIntPhase[fRd/2, \[Eta], \[Beta]];
	
	\[Beta]1correction = \[Delta]\[Beta]1[\[Eta], \[Chi]s, \[Chi]a, \[Beta], \[Sigma]] fRd/2;
	\[Beta]0correction = \[Delta]\[Beta]0[\[Eta], \[Chi]s, \[Chi]a, \[Beta], \[Sigma]];
	\[Alpha]1correction = \[Delta]\[Alpha]1[fRd, fdamp, \[Eta], \[Chi]s, \[Chi]a, \[Alpha], \[Beta], \[Sigma]] fRd/2;
	
	-MR - \[Alpha]1correction + Intermediate + \[Beta]0correction + \[Beta]1correction
]


(*
This does not have dimension of time bcs f is dimensionless frequency. The right answer is given by multipliyng the
function by \[ScriptCapitalM] G/c^3 
*)
t0[f_,fRd_, fdamp_, \[Eta]_, \[Alpha]_/;(VectorQ@\[Alpha]&&Length@\[Alpha]==5)]/;(
	And@@(NumberQ/@{f, fRd, fdamp, \[Eta]}) 
) := (\[Alpha][[1;;4]] . {1,f^-2, f^(-1/4), fdamp/(fdamp^2 + (f-\[Alpha][[-1]] fRd)^2)}) \[Eta]^-1;


TotalPhase[f_, frd_, fdamp_, \[Eta]_, \[Chi]s_, \[Chi]a_, \[Alpha]_, \[Beta]_, \[Sigma]_] := (

	InspiralPhase[f, \[Eta], \[Chi]s, \[Chi]a, \[Sigma]] unitStep[0.018 - f] + 
	(
	   UncorrectedIntPhase[f, \[Eta], \[Beta]] + \[Delta]\[Beta]0[\[Eta],\[Chi]s,\[Chi]a,\[Beta],\[Sigma]] + \[Delta]\[Beta]1[\[Eta], \[Chi]s, \[Chi]a, \[Beta], \[Sigma]] f
	) unitStep[(f - 0.018) (0.5 frd - f)] +
	(
		UncorrectedMRPhase[f,frd,fdamp,\[Eta], \[Alpha]] + 
		\[Delta]\[Alpha]0[frd,fdamp,\[Eta],\[Chi]s,\[Chi]a,\[Alpha], \[Beta],\[Sigma]] +
		\[Delta]\[Alpha]1[frd,fdamp,\[Eta], \[Chi]s, \[Chi]a, \[Alpha],\[Beta],\[Sigma]] f
	) unitStep[f-0.5frd]
	
);


(*JUST TO BE USED BY \[Phi]ref*)

auxiliarTotalPhase[f_, frd_?NumericQ, fdamp_?NumericQ, \[Eta]_?NumericQ, \[Chi]s_?NumericQ, \[Chi]a_?NumericQ, \[Alpha]_List, \[Beta]_List, \[Sigma]_List] := (

	auxiliarInspiralPhase[f, \[Eta], \[Chi]s, \[Chi]a, \[Sigma]] unitStep[0.018 - f] + 
	(
	   UncorrectedIntPhase[f, \[Eta], \[Beta]] + \[Delta]\[Beta]0[\[Eta],\[Chi]s,\[Chi]a,\[Beta],\[Sigma]] + \[Delta]\[Beta]1[\[Eta], \[Chi]s, \[Chi]a, \[Beta], \[Sigma]] f
	) unitStep[(f - 0.018) (0.5 frd - f)] +
	(
		UncorrectedMRPhase[f,frd,fdamp,\[Eta], \[Alpha]] + 
		\[Delta]\[Alpha]0[frd,fdamp,\[Eta],\[Chi]s,\[Chi]a,\[Alpha], \[Beta],\[Sigma]] +
		\[Delta]\[Alpha]1[frd,fdamp,\[Eta], \[Chi]s, \[Chi]a, \[Alpha],\[Beta],\[Sigma]] f
	) unitStep[f-0.5frd]
	
)//N;



(*Phase[f, variables ];//AbsoluteTiming == {0.002362`,Null}*)


(* ::Input:: *)
(**)


(* ::Section::Closed:: *)
(*Amplitude*)


MRAmplitude[f_, fRd_?NumericQ,fdamp_,\[Eta]_, \[Gamma]1_,\[Gamma]2_,\[Gamma]3_] = (\[Gamma]1 (\[Gamma]3 fdamp) Exp[-((\[Gamma]2 (f-fRd))/(\[Gamma]3 fdamp))])/((f-fRd)^2+\[Gamma]3^2 fdamp^2);


AmplitudePNCoeffs[\[Eta]_?NumericQ,\[Chi]s_,\[Chi]a_]= With[
	{\[Delta]=Sqrt[1-4 \[Eta]]},
	{
		1,
		0,
		-(323/224)+(451 \[Eta])/168,
		(27 \[Delta] \[Chi]a)/8+(27/8-(11 \[Eta])/6) \[Chi]s,
		-(27312085/8128512)-(1975055 \[Eta])/338688+(105271 \[Eta]^2)/24192+(-(81/32)+8 \[Eta]) \[Chi]a^2-(81 \[Delta] \[Chi]a \[Chi]s)/16+(-(81/32)+(17 \[Eta])/8) \[Chi]s^2,
		-(1/64) (85 \[Pi])+(85 \[Pi] \[Eta])/16+\[Delta] (285197/16128-(1579 \[Eta])/4032) \[Chi]a+(285197/16128-(15317 \[Eta])/672-(2227 \[Eta]^2)/1008) \[Chi]s,
		-(177520268561/8583708672)+(545384828789/5007163392-(205 \[Pi]^2)/48) \[Eta]-(3248849057 \[Eta]^2)/178827264+(34473079 \[Eta]^3)/6386688+(1614569/64512-(1873643 \[Eta])/16128+(2167 \[Eta]^2)/42) \[Chi]a^2+((31 \[Pi])/12-(7 \[Pi] \[Eta])/3) \[Chi]s+(1614569/64512-(61391 \[Eta])/1344+(57451 \[Eta]^2)/4032) \[Chi]s^2+\[Delta] \[Chi]a ((31 \[Pi])/12+(1614569/32256-(165961 \[Eta])/2688) \[Chi]s)
	}
];

InspiralAmplitude[f_,\[Eta]_?NumericQ,\[Chi]s_,\[Chi]a_, \[Rho]_List] := With[
	{ TF2f = {ConstantArray[1, Length@f],f^(1/3) \[Pi]^(1/3),f^(2/3) \[Pi]^(2/3),f \[Pi],f^(4/3) \[Pi]^(4/3),f^(5/3) \[Pi]^(5/3),f^2 \[Pi]^2} },
	
	AmplitudePNCoeffs[\[Eta], \[Chi]s, \[Chi]a] . TF2f  +  \[Rho] . {f^(7/3), f^(8/3), f^3}
];


auxiliarInspiralAmplitude[f_,\[Eta]_?NumericQ,\[Chi]s_,\[Chi]a_, \[Rho]_List/;Length@\[Rho]==3] := With[
	{ TF2f = {1,f^(1/3) \[Pi]^(1/3),f^(2/3) \[Pi]^(2/3),f \[Pi],f^(4/3) \[Pi]^(4/3),f^(5/3) \[Pi]^(5/3),f^2 \[Pi]^2} },
	
	AmplitudePNCoeffs[\[Eta], \[Chi]s, \[Chi]a] . TF2f  +  \[Rho] . {f^(7/3), f^(8/3), f^3}
];


(*Solutions to (22)-(26) system. Note that v_i \[And] d_i are the values and derivatives of amplitudes Amr/A0 \[And] Ainsp/A0 *)
m  = List[
	 (f^#&/@Range[0,4])/.f->f1,
	 (f^#&/@Range[0,4])/.f->f2,
	 (f^#&/@Range[0,4])/.f->f3,
	 D[f^(-7/6) (f^#&/@Range[0,4]), f]/.f->f1,
	 D[f^(-7/6) (f^#&/@Range[0,4]), f]/.f->f3
];

List[
	\[Delta]0[f1_,f2_,f3_,v1_,v2_,v3_,d1_,d3_],
    \[Delta]1[f1_,f2_,f3_,v1_,v2_,v3_,d1_,d3_],
    \[Delta]2[f1_,f2_,f3_,v1_,v2_,v3_,d1_,d3_],
    \[Delta]3[f1_,f2_,f3_,v1_,v2_,v3_,d1_,d3_],
    \[Delta]4[f1_,f2_,f3_,v1_,v2_,v3_,d1_,d3_]
    
] = LinearSolve[m, {v1,v2, v3, -7/6 f1^(-13/6) v1 + f1^(-7/6) d1,- 7/6 f3^(-13/6) v3 + f3^(-7/6) d3}]//Simplify;

Clear[m]


d1[\[Eta]_?NumericQ, \[Chi]s_?NumericQ, \[Chi]a_?NumericQ, \[Rho]_List] := Module[
	{f,v},
	v = {0, \[Pi]^(1/3)/(3 f^(2/3)),(2 \[Pi]^(2/3))/(3 f^(1/3)),\[Pi],4/3 f^(1/3) \[Pi]^(4/3),5/3 f^(2/3) \[Pi]^(5/3),2 f \[Pi]^2};
	
	(AmplitudePNCoeffs[\[Eta], \[Chi]s, \[Chi]a] . v/.f->14/1000) +  (\[Rho] . {(7 f^(4/3))/3,(8 f^(5/3))/3,3 f^2}/.f->14/1000)
];


d3[fpeak_?NumericQ,fRd_,fdamp_,\[Eta]_,\[Gamma]1_,\[Gamma]2_,\[Gamma]3_]:=-((2 E^(-(((fpeak-fRd) \[Gamma]2)/(fdamp \[Gamma]3))) fdamp (fpeak-fRd) \[Gamma]1 \[Gamma]3)/((fpeak-fRd)^2+fdamp^2 \[Gamma]3^2)^2)-(E^(-(((fpeak-fRd) \[Gamma]2)/(fdamp \[Gamma]3))) \[Gamma]1 \[Gamma]2)/((fpeak-fRd)^2+fdamp^2 \[Gamma]3^2);


IntAmp[f_,  F3_, V1_, V2_, V3_, D1_, D3_] := Module[
	{\[CapitalDelta]0,\[CapitalDelta]1,\[CapitalDelta]2,\[CapitalDelta]3,\[CapitalDelta]4},
    
    {\[CapitalDelta]0, \[CapitalDelta]1, \[CapitalDelta]2, \[CapitalDelta]3, \[CapitalDelta]4} = List[
	    \[Delta]0[14/1000, (14/1000 + F3)/2 ,F3, V1,V2,V3, D1,D3],
        \[Delta]1[14/1000,(14/1000 + F3)/2 , F3, V1,V2,V3,D1,D3],
        \[Delta]2[14/1000,(14/1000 + F3)/2 , F3, V1,V2,V3,D1,D3],
        \[Delta]3[14/1000,(14/1000 + F3)/2 , F3, V1,V2,V3,D1,D3],
        \[Delta]4[14/1000,(14/1000 + F3)/2 , F3, V1,V2,V3,D1,D3]
    ];
    
    {\[CapitalDelta]0, \[CapitalDelta]1, \[CapitalDelta]2, \[CapitalDelta]3, \[CapitalDelta]4} . (f^#&/@Range[0,4])
	];


TotalAmplitude[f_, frd_?NumericQ, fdamp_,  \[Eta]_, \[Chi]s_, \[Chi]a_, \[Gamma]1_, \[Gamma]2_, \[Gamma]3_, \[Rho]_, F3_,V1_,V2_,V3_,D1_,D3_] := (
	unitStep[0.014 - f] InspiralAmplitude[f, \[Eta], \[Chi]s, \[Chi]a, \[Rho]] + 
	unitStep[(F3 - f) (f - 0.014)] IntAmp[f, F3,V1,V2,V3,D1,D3] + 
	unitStep[(f-F3) (0.2 - f)] MRAmplitude[f, frd, fdamp, \[Eta], \[Gamma]1, \[Gamma]2, \[Gamma]3] 
)


(*TotalAmplitude[flist, variables];//AbsoluteTiming//ScientificForm == {1.519\[Times]10^(-3),Null}*)


(* ::Section::Closed:: *)
(*FullWaveform*)


Waveform[freq_, \[Eta]_, \[ScriptCapitalM]_, \[Chi]s_, \[Chi]a_, dL_,cos\[Iota]_, \[Phi]c_, tc_, freqref_]/; \[ScriptCapitalM] G freq[[1]] >= 0.2 := ConstantArray[0, Length@freq]


Waveform[freq_, \[Eta]_, \[ScriptCapitalM]_, \[Chi]s_, \[Chi]a_, dL_,cos\[Iota]_, \[Phi]c_, tc_, freqref_] := Module[

	{\[Rho], V2, \[Gamma]1, \[Gamma]2, \[Gamma]3, \[Sigma], \[Beta], \[Alpha], fpeak, frd, fdamp, V1, V3, D1, D3, \[ScriptCapitalA], \[Phi], f = \[ScriptCapitalM] G freq, \[Phi]ref, \[ScriptT]0, fref = G \[ScriptCapitalM] freqref, scalar\[Delta]\[Phi], 
	\[ScriptCapitalA]0RippleWrong, \[ScriptH],
	S = 1/4 (1+ Sqrt[1-4 \[Eta]])^2 (\[Chi]s+\[Chi]a) +  1/4 (1-Sqrt[1-4 \[Eta]])^2 (\[Chi]s-\[Chi]a),
	\[Chi]Pn = (1-76 \[Eta]/113) \[Chi]s + Sqrt[1-4 \[Eta]] \[Chi]a
	},
	
	\[Rho] = PhenomCoeff[\[Eta], \[Chi]Pn, #]&/@(PhenomDTableV[[1;;3]]);
	V2 = PhenomCoeff[\[Eta], \[Chi]Pn, #]&@PhenomDTableV[[4]];
	{\[Gamma]1, \[Gamma]2, \[Gamma]3} = PhenomCoeff[\[Eta], \[Chi]Pn, #]&/@(PhenomDTableV[[5;;7]]);
	\[Sigma] = PhenomCoeff[\[Eta], \[Chi]Pn, #]&/@(PhenomDTableV[[8;;11]]);
	\[Beta] = PhenomCoeff[\[Eta], \[Chi]Pn, #]&/@(PhenomDTableV[[12;;14]]);
	\[Alpha] = PhenomCoeff[\[Eta], \[Chi]Pn, #]&/@(PhenomDTableV[[15;;19]]);
	
	
	{fpeak, frd, fdamp} = SpecialFrequencies[\[Eta], S, \[Gamma]2, \[Gamma]3];
	
	{V1, V3, D1, D3} = {
		auxiliarInspiralAmplitude[14/1000, \[Eta], \[Chi]s, \[Chi]a, \[Rho]] ,
		MRAmplitude[fpeak, frd, fdamp, \[Eta], \[Gamma]1, \[Gamma]2, \[Gamma]3],
		d1[\[Eta], \[Chi]s, \[Chi]a, \[Rho]],
		d3[fpeak, frd, fdamp, \[Eta], \[Gamma]1, \[Gamma]2, \[Gamma]3]
		};
	
	\[ScriptCapitalA] = TotalAmplitude[f, frd, fdamp, \[Eta], \[Chi]s, \[Chi]a, \[Gamma]1, \[Gamma]2, \[Gamma]3, \[Rho], fpeak, V1, V2, V3, D1, D3];
	
	\[Phi] = TotalPhase[f, frd, fdamp, \[Eta], \[Chi]s, \[Chi]a, \[Alpha], \[Beta], \[Sigma]];
	
	\[Phi]ref = TotalPhase[{fref,fref}, frd, fdamp, \[Eta], \[Chi]s, \[Chi]a, \[Alpha], \[Beta], \[Sigma]]//Last;
	
	\[ScriptT]0 = t0[fpeak, frd, fdamp, \[Eta], \[Alpha]];
	
	scalar\[Delta]\[Phi] = - \[Phi]ref - 2 \[Phi]c + \[ScriptT]0 fref;
	
	\[Phi] += scalar\[Delta]\[Phi] - (\[ScriptT]0) f + 2 \[Pi] tc freq; 
	
	\[ScriptCapitalA]0RippleWrong = Sqrt[5/6 \[Eta]] (c f^(-7/6))/(2 \[Pi]^(2/3)) (\[ScriptCapitalM] G)^2/dL; (* G = G/c^3 DO NOT FORGET!!*)
	
	\[ScriptH] = \[ScriptCapitalA]0RippleWrong \[ScriptCapitalA] Exp[-I \[Phi]];
	
	{
		\[ScriptH] 1/2 (1 + cos\[Iota]^2),
		- I \[ScriptH] cos\[Iota]
	}
	
	
]


(* (c)/(1) (\[ScriptCapitalM] G)^2/dL/.{G->Quantity["GravitationalConstant"],
\[ScriptCapitalM]-> Quantity["SolarMass"], dL-> Quantity["SolarRadius"], c-> Quantity["SpeedOfLight"]}//UnitDimensions*)


(* ::Section::Closed:: *)
(*DATA*)


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


End[];
EndPackage[];
