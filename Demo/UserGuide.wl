(* ::Package:: *)

Quit


$HistoryLength=1;


PacletDirectoryLoad[NotebookDirectory[]//ParentDirectory[#,2]&];


(* ::Text:: *)
(*I DIDN' T SET UP PHENOMPV2 YET, SO DON' T TRY TO RUN THOSE TESTS*)


<<FelipeBarbosa`SymDALI`


(* ::Section::Closed:: *)
(*Util functions*)


(* ::Subsection::Closed:: *)
(*Comparison functions*)


(*function to calculate the relative differences between numbers*)
Clear@RelativeDiff

Attributes[RelativeDiff] = {Listable};
RelativeDiff[x_,y_]/;x==0 &&y==0 := 0
RelativeDiff[x_, y_]/; x==0 && y!=0 := 1
RelativeDiff[x_, y_]/; y==0 && x!=0 := 1

RelativeDiff[x_, y_]/; x!=0 &&y!=0 := With[
	{diff = x-y},
	Max[Abs@{diff/x, diff/y}]
]



(*function to compare complex vectors and generate plots*)
Clear@CompareComplexVectors

CompareComplexVectors[{v1_, v2_}, {name1_, name2_}, dataRange_]/;(
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
		PlotRange->{MinMax[dataRange], All},
		DataRange->MinMax[dataRange], 
		Frame->True,
		Background->White
	];
	imPlot = ListLinePlot[
		{Imv1, Imv2},
		PlotLegends->Placed[{name1, name2}, {Right, Top}],
		ImageSize->Medium,
		PlotLabel->Style["Im", Black],
		PlotRange->{MinMax[dataRange], All},
		DataRange->MinMax[dataRange],
		Frame->True,
		Background->White
		];
		
	reDiffPlot = ListLinePlot[
		Rediff,
		ImageSize->Medium,
		PlotLabel->Style["RelativeDiff Re part", Black],
		PlotRange->All,
		PlotRange->{MinMax[dataRange], All},
		DataRange->MinMax[dataRange],
		Frame->True,
		ScalingFunctions->"Log10",
		Background->White
	];
		
	imDiffPlot = ListLinePlot[
		Imdiff,
		ImageSize->Medium,
		PlotLabel->Style["RelativeDiff Im part", Black],
		PlotRange->{MinMax[dataRange], All},
		DataRange->MinMax[dataRange],
		GridLinesStyle->Directive[Red, 13, Dashed],
		Frame->True,
		ScalingFunctions->"Log10",
		Background->White
	];
	
	Grid[{{rePlot, imPlot}, {reDiffPlot, imDiffPlot}}]


]


(* ::Subsection::Closed:: *)
(*NGrad \[And] NFisher*)


NGrad//Clear

NGrad[f_, vars_, ni_, nf_] := Module[
	{h = 1. 10^-7, dummy, Point1, Point2, denominator},
	
	Point1 = Table[
		dummy = vars;
		dummy[[i]] = If[vars[[i]]==0, h, vars[[i]] + h vars[[i]]];
		dummy,
		{i, ni, nf}
	];
	
	
	Point2 = ConstantArray[vars, (nf-ni+1)];
	
	denominator = Table[If[vars[[i]]==0, h, h vars[[i]]], {i, ni, nf}];

	(f@@@Point1 - f@@@Point2)/denominator
]


NFisher[f_, vars_, ni_, nf_, Psd_List] := Module[
	{igrad},
	
	igrad = NGrad[f, vars, ni, nf]//Transpose;
	
	4 0.125 Sum[
		Conjugate[igrad[[i]]]\[TensorProduct]igrad[[i]]/Psd[[i]], 
		{i, Length@Psd}
	]//Re
]


(* ::Section::Closed:: *)
(*Test Pattern Functions against lal: *)


(* ::Text:: *)
(*Define the lal fpfc for H1:*)


DeleteObject/@ExternalSessions[];
Clear@python

python = StartExternalSession["Python"];

ExternalEvaluate[python, "
from lal import antenna
import numpy as np
"]

lalFpFc = ExternalFunction[python, "

def fpfc(ra, dec, psi):
	res = antenna.AntennaResponse('H1', ra, dec, psi=psi, times=630696086.1999)
	return np.array([res.plus, res.cross])
"]


(* ::Text:: *)
(*We work with GPS time  630696086.1999 bcs the corresponding GMST is small:*)


ExternalEvaluate["Python","
import lal

# GPS time as float (e.g., time of GW150914)
gps_time_float = 630696086.1999

# Step 1: Convert to LIGOTimeGPS object
gps_time = lal.LIGOTimeGPS(gps_time_float)

# Step 2: Compute GMST in radians
gmst_rad = lal.GreenwichMeanSiderealTime(gps_time)

print(\"GMST (radians):\", gmst_rad)

"]


Clear@test

test := Module[
	{r1, r2, \[Theta], \[Phi],\[Alpha], \[Delta], \[Psi], MMA},
	\[Delta] = RandomReal[{-\[Pi]/2,\[Pi]/2}];
	\[Psi] = RandomReal[{0, \[Pi]}];
	\[Alpha] = RandomReal[{0, 2 \[Pi]}];
	
	\[Theta] = \[Pi]/2-\[Delta]; 
	\[Phi] = \[Alpha] - (-2.821265599576199 10^-6);
	
	MMA = PatternFunctions[
		{10},
		\[Delta], \[Alpha], \[Psi],
		2.821265599576199 10^-6,
		{0,0,0}, (*The lal definition does not include the time delay so we set the position to 0.*)
		DetectorTensor["H1"]
	];
	
	
	
	r1 = RelativeDiff[
		MMA[[1]]//Last,
		(lalFpFc[\[Alpha], \[Delta], \[Psi]]//Normal)[[1]]
	];
	
	r2 = RelativeDiff[
		MMA[[2]]//Last,
		(lalFpFc[\[Alpha], \[Delta], \[Psi]]//Normal)[[2]]
	];
	
	{r1,r2}
]


test//ScientificForm


(* ::Text:: *)
(*If you get some error from the Python function just run again. I don't know why, but sometimes python external functions generate errors inside ```Do``` or ```Table```*)


list = Table[test, 10^4]; 


(* ::Text:: *)
(*The points where the transition from negative to positive happens have the highest relative differences, bcs one implementation may go to zero slightly faster than the other*)


Count[list[[All,1]], x_/;x<=10^-3]//Echo[#, "# of points with F+ agreement better than 0.1%:"]&;
Count[list[[All,2]], x_/;x<=10^-3]//Echo[#, "# of points with Fx agreement better than 0.1%:"]&;
list[[All,1]]//Sort//ListPlot[#, PlotRange->All, ScalingFunctions->"Log10"]&
list[[All,2]]//Sort//ListPlot[#, PlotRange->All, ScalingFunctions->"Log10"]&


Clear@list


(* ::Section::Closed:: *)
(*Test IMRPhenomD against lal:*)


(* ::Subsection::Closed:: *)
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


(* ::Subsection:: *)
(*Test hp and hc*)


(* ::Text:: *)
(*It is easier to just run the section "functions to generate plots" and go to "check plots" to see the plots*)


(* ::Subsubsection:: *)
(*functions to generate plots*)


Clear@testhp

testhp := Module[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA, lalD, lalR, MMAR, lalIm, MMAIm, fmax,\[Eta],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], 
		 Rediff, Imdiff, f1, sp, irrelevant, rePlot, imPlot, reDiffPlot, imDiffPlot, mc, \[Delta], \[Chi]s, \[Chi]a
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	\[Eta] = (m1 m2)/(m1+m2)^2;
	\[Delta] = Sqrt[1 - 4 \[Eta]];
	mc = (m1+m2) \[Eta]^(3/5);
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	\[Chi]s = (s1z+s2z)/2;  \[Chi]a = (s1z-s2z)/2;
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	fmax =  0.2/(G (m1+m2))//Round;
	f1 = 0.014/(G (m1+m2));
	
	
	MMA =  hphcIMRPhenomD[f, mc, \[Delta], \[Chi]s, \[Chi]a, \[Iota], 0, \[Phi]Ref, 10^3][[1]];
	lalD = lal[m1, m2, s1z,s2z, 1,\[Iota], \[Phi]Ref, 1., 10., fmax, 10., "IMRPhenomD"][[1]];
	
	(*
		lal starts at f=0 and goes to the nearest fp = 2^k, where fp>= fmax.
		This line insures the same f points for lal and MMA
	*)
	lalD = lalD[[11;;Length[MMA]+10]]; 
	
	
	CompareComplexVectors[{lalD, MMA}, {"lal", "MMA"}, f]
	
	
]


Clear@testhc

testhc := Module[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA, lalD, lalR, MMAR, lalIm, MMAIm, fmax,\[Eta],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], 
		 Rediff, Imdiff, f1, sp, irrelevant, rePlot, imPlot, reDiffPlot, imDiffPlot, mc, \[Delta], \[Chi]s, \[Chi]a
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	\[Eta] = (m1 m2)/(m1+m2)^2; 
	\[Delta] = Sqrt[1-4 \[Eta]];
	mc = (m1+m2) \[Eta]^(3/5);
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	
	\[Chi]s = (s1z+s2z)/2;  \[Chi]a = (s1z-s2z)/2;
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	fmax =  0.2/(G (m1+m2))//Round;
	f1 = 0.014/(G (m1+m2));
	
	
	MMA =  hphcIMRPhenomD[f, mc, \[Delta], \[Chi]s, \[Chi]a, \[Iota], 0, \[Phi]Ref, 10^3][[2]];
	lalD = lal[m1, m2, s1z,s2z, 1,\[Iota],\[Phi]Ref, 1., 10., fmax, 10., "IMRPhenomD"][[2]];
	
	(*
		lal starts at f=0 and goes to the nearest fp = 2^k, where fp>= fmax.
		This line insures the same f points for lal and MMA
	*)
	lalD = lalD[[11;;Length[MMA]+10]]; 
	
	
	CompareComplexVectors[{lalD, MMA}, {"lal", "MMA"}, f]
	
	
]




(* ::Subsubsection::Closed:: *)
(*check plots: *)


(* ::Text:: *)
(*The functions generates random values of {m1, m2, \[Eta], s1z, s2z, \[Iota], \[Phi]Ref} to do the analysis*)


(* ::Text:: *)
(*You can see that the highest relative differences are on points where the functions cross zero*)


testhp


(* ::Text:: *)
(*Curiously lal uses "-i Cos[\[Iota]]" instead of "i Cos[\[Iota]]" for hc like GWFAST.*)


testhc


Clear[testhc, testhp]


(* ::Section::Closed:: *)
(*Test IMRPhenomPv2 against lal:*)


(* ::Subsection::Closed:: *)
(*lal hphc function*)


Clear[python]
DeleteObject/@ExternalSessions[];

python = StartExternalSession["Python"]
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
def lal_hp(m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, dL, iota, phiRef, deltaF, f_min, f_max, f_ref, appr):

    approximant = lalsim.SimInspiralGetApproximantFromString(appr)
    
    dL_m = dL * 3.08568 * (10**22)  # convert megaparsecs to meters

    # solar mass to Kg
    m1_kg = m1 * 1.9884 * (10**30)
    m2_kg = m2 * 1.9884 * (10**30)
    
    hp, hc = lalsim.SimInspiralChooseFDWaveform(
        m1_kg,
        m2_kg,
        s1x,#s1x,
        s1y,#s1y,
        s1z,
        s2x,#s2x,
        s2y,#s2y,
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


(* ::Subsection::Closed:: *)
(*Test hp and hc*)


(* ::Text:: *)
(*It is easier to just run the section "functions to generate plots" and go to "check plots" to see the plots*)


(* ::Subsubsection::Closed:: *)
(*functions to generate plots*)


Clear@RandomSpin
RandomSpin[] := Module[
	{abs, spin}, 
	abs=2; 
	While[
		abs>1, 
		spin = RandomReal[{-1,1},3];
		abs = Norm[spin];
	];
	spin
];


Clear@testhp

testhp := Module[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA, lalD, lalR, MMAR, lalIm, MMAIm, fmax,\[Eta],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], 
		 Rediff, Imdiff, f1, sp, irrelevant, rePlot, imPlot, reDiffPlot, imDiffPlot
	},
	
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	\[Eta] = (m1 m2)/(m1+m2)^2;
	{s1x, s1y, s1z} = RandomSpin[];
	{s2x, s2y, s2z} = RandomSpin[];
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	fmax =  0.2/(G (m1+m2))//Round;
	f1 = 0.014/(G (m1+m2));
	
	
	MMA =  hphcIMRPhenomPv2[f, m1+m2, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], 0, \[Phi]Ref, 10^-3][[1]];
	lalD = lal[m1, m2, s1x, s1y, s1z,s2x, s2y, s2z, 1,\[Iota],\[Phi]Ref, 1., 10., fmax, 10., "IMRPhenomPv2"][[1]];
	
	(*
		lal starts at f=0 and goes to the nearest fp = 2^k, where fp>= fmax.
		This line insures the same f points for lal and MMA
	*)
	lalD = lalD[[11;;Length[MMA]+10]]; 
	
	
	CompareComplexVectors[{lalD, MMA}, {"lal", "MMA"}, f]
	
	
]


Clear@testhc

testhc := Module[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA, lalD, lalR, MMAR, lalIm, MMAIm, fmax,\[Eta],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], 
		 Rediff, Imdiff, f1, sp, irrelevant, rePlot, imPlot, reDiffPlot, imDiffPlot
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	\[Eta] = (m1 m2)/(m1+m2)^2;
	{s1x, s1y, s1z} = RandomSpin[];
	{s2x, s2y, s2z} = RandomSpin[];
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	fmax =  0.2/(G (m1+m2))//Round;
	f1 = 0.014/(G (m1+m2));
	
	
	MMA =  hphcIMRPhenomPv2[f, m1+m2, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], 0, \[Phi]Ref, 10^-3][[2]];
	lalD = lal[m1, m2, s1x, s1y, s1z,s2x, s2y, s2z, 1,\[Iota],\[Phi]Ref, 1., 10., fmax, 10., "IMRPhenomPv2"][[2]];
	
	(*
		lal starts at f=0 and goes to the nearest fp = 2^k, where fp>= fmax.
		This line insures the same f points for lal and MMA
	*)
	lalD = lalD[[11;;Length[MMA]+10]]; 
	
	CompareComplexVectors[{lalD, MMA}, {"lal", "MMA"}, f]
	
	
	
	
]


(* ::Subsubsection::Closed:: *)
(*check plots: *)


(* ::Text:: *)
(*The functions generates random values of {m1, m2, \[Eta], s1z, s2z, \[Iota], \[Phi]Ref} to do the analysis*)


(* ::Text:: *)
(*You can see that the highest relative differences are on points where the functions cross zero*)


testhp


testhc


Clear[testhc, testhp]


DeleteObject/@ExternalSessions[]
Clear@python


(* ::Section::Closed:: *)
(*Test numerical x symbolic Fisher matrices: *)


Clear@Test


(* ::Subsection::Closed:: *)
(*IMRPhenomD test*)


ClearAll@strain;
strain[
	\[Theta]_, \[Phi]_, \[Psi]_, Mc_, \[Delta]_, \[Chi]s_, \[Chi]a_, \[Iota]_, invdL_,tc_, \[Phi]ref_
] := Module[
	{FpFc, fvec = Range[20, 1024, 0.125], hphc},
	
	FpFc = PatternFunctions[fvec, \[Pi]/2-\[Theta], \[Phi], \[Psi], 0, DetectorVertex["H1"], DetectorTensor["H1"]];
	
	hphc = hphcIMRPhenomD[
		fvec, Mc, \[Delta], \[Chi]s, \[Chi]a, \[Iota], tc, \[Phi]ref, invdL
	];
	
	FpFc[[1]] hphc[[1]] + FpFc[[2]] hphc[[2]]
]


psd = (ASD["ET-D"]@Range[20, 1024, 0.125])^2;


ivars ={"\[ScriptCapitalM]c", "\[Delta]", "\[Chi]s", "\[Chi]a",
	"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "1/dL", 
	"tc", "\[Phi]ref"};


H1 = <|
	"ASD" -> ASD["ET-D"],
	"Position" -> DetectorVertex["H1"],
	"DetectorTensor" -> DetectorTensor["H1"]
|>;


Test["IMRPhenomD"] := Module[
	{M, \[Eta], s1z, s2z, \[Iota],Mc, \[Delta],\[Chi]s, \[Chi]a, tc, \[Phi]ref, dL, \[Theta], \[Phi], \[Psi], numerical, fp, Symbolic},
	
	M = RandomReal[{20, 120}];
	\[Eta] = RandomReal[{0.001, 0.2499}]; (*exactly at 0.25 is problematic for numeric derivatives*)
	Mc = M \[Eta]^(3/5); \[Delta] = Sqrt[1-4 \[Eta]];
	
	{s1z, s2z, tc} = RandomReal[{-1,1},3];
	\[Chi]s = (s1z+s2z)/2;  \[Chi]a = (s1z-s2z)/2;
	
	{\[Theta], \[Iota], \[Psi]} = RandomReal[{0, \[Pi]}, 3];
	{\[Phi]ref, \[Phi]} = RandomReal[{0, 2 \[Pi]},2];
	dL = RandomReal[{2, 10}];
	
	
	
	fp = Association@(Thread@Rule[ivars, {Mc,\[Delta],\[Chi]s, \[Chi]a,\[Iota],\[Theta],\[Phi],\[Psi],1/dL,tc,\[Phi]ref}]);
	
	numerical = NFisher[strain, {\[Theta],\[Phi],\[Psi],Mc,\[Delta],\[Chi]s, \[Chi]a,\[Iota],1/dL,tc,\[Phi]ref}, 1, 11, psd];
	
	Symbolic = DALITensors["IMRPhenomD", fp, {H1}, 1, "fmin"->20, "fmax"->1024, "res"->8033]//QuietEcho;
	Symbolic = ArrayReshape[Symbolic, {11,11}];
	
	
	{numerical[[9,10]] ,numerical[[9,11]], numerical[[10, 9]], numerical[[11, 9]] }//Echo;
	{Symbolic[[9,10]], Symbolic[[9,11]],Symbolic[[10, 9]], Symbolic[[11, 9]]}//Echo;
	(*
	these should be exactly zero since they are purelly imaginary contributions in the tensor product,
	we set them to zero, so they don't contaminate the comparison. The Symbolic implementation is usually closer
	to zero by some orders of magnitude
	*)
	numerical[[9,10]] = numerical[[9,11]] = numerical[[10, 9]] = numerical[[11,9]] = 0;
	Symbolic[[9,10]] = Symbolic[[9,11]] =Symbolic[[10,9]] = Symbolic[[11,9]] = 0;
	
	RelativeDiff[numerical, Symbolic]

]


a = Test["IMRPhenomD"]//UpperTriangularize;

a//MatrixForm

ListPlot[Flatten[a]//Sort, ScalingFunctions->"Log10", PlotRange->All]
(*the high differences seem to be always on the dL column*)


(* ::Subsection::Closed:: *)
(*IMRPhenomPv2 test*)


ClearAll@strain;
strain[\[Theta]_, \[Phi]_, \[Psi]_, m1_, m2_, s1x_,s1y_,s1z_,s2x_,s2y_,s2z_, \[Iota]_, dL_, tc_,\[Phi]ref_] := Module[
	{FpFc, fvec = Range[20, 1024, 0.125], hphc},
	
	FpFc = PatternFunctions[fvec, \[Pi]/2-\[Theta], \[Phi], \[Psi], 0, DetectorVertex["H1"], DetectorTensor["H1"]];
	
	hphc = hphcIMRPhenomPv2[fvec, m1+m2, (m1 m2)/(m1+m2)^2, s1x,s1y,s1z,s2x,s2y,s2z, \[Iota], tc, \[Phi]ref, dL];
	
	FpFc[[1]] hphc[[1]] + FpFc[[2]] hphc[[2]]
]


ivars = {"m1", "m2", "s1x","s1y","s1z","s2x", "s2y", "s2z",
	"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "dL", 
	"tc", "\[Phi]ref"};


psd = (ASD["ET-D"]@Range[20, 1024, 0.125])^2;


H1 = <|
	"ASD" -> ASD["ET-D"],
	"Position" -> DetectorVertex["H1"],
	"DetectorTensor" -> DetectorTensor["H1"]
|>;


RandomSpin[] := Module[
	{norm, s},
	norm=2;
	While[norm>1,s =RandomReal[{-1,1},3 ]; norm=Norm[s]];
	s
]


Test["IMRPhenomPv2"] := Module[
	{m1, m2, s1x,s1y,s1z,s2x,s2y,s2z, \[Iota], tc, \[Phi]ref, dL, \[Theta], \[Phi], \[Psi], numerical, fp, Symbolic},
	
	{m1, m2} = RandomReal[{20, 80}, 2]//ReverseSort;
	
	{s1x,s1y,s1z} = RandomSpin[];
	{s2x,s2y,s2z} = RandomSpin[];
	
	{tc} = RandomReal[{-1,1},1];
	
	{\[Theta], \[Iota], \[Psi]} = RandomReal[{0, \[Pi]}, 3];
	{\[Phi]ref, \[Phi]} = RandomReal[{0, 2 \[Pi]},2];
	dL = RandomReal[{2, 10}];
	
	fp = Association@(Thread@Rule[ivars, {m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Iota],\[Theta],\[Phi],\[Psi],dL,tc,\[Phi]ref}]);
	
	numerical = NFisher[strain, {\[Theta],\[Phi],\[Psi],m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Iota],dL,tc,\[Phi]ref}, 1, 15, psd];
	
	
	Symbolic = FisherMatrix["IMRPhenomPv2", fp, {H1}, "fmin"->20, "fmax"->1024, "res"->8032]//QuietEcho;
	
	{numerical[[-3,-2]] , numerical[[-2, -3]] }//Echo;
	{Symbolic[[-3,-2]], Symbolic[[-2, -3]]}//Echo;
	(*
	these should be exactly zero since they are purelly imaginary contributions in the tensor product,
	we set them to zero, so they don't contaminate the comparison. The Symbolic implementation is usually closer
	to zero by some orders of magnitude
	*)
	numerical[[-3,-2]]  = numerical[[-2, -3]] = 0;
	Symbolic[[-3,-2]]  =Symbolic[[-2, -3]] = 0;
	
	RelativeDiff[numerical, Symbolic]

]


a = Test["IMRPhenomPv2"]//UpperTriangularize;

a//MatrixForm

ListPlot[Flatten[a]//Sort, ScalingFunctions->"Log10", PlotRange->All]
(*the high differences seem to be always on the dL column*)




(* ::Subsection::Closed:: *)
(*IMRPhenomPv2 test with all \[Delta]pi*)


ClearAll@strain;
strain[
	\[Theta]_, \[Phi]_, \[Psi]_, 
	m1_, m2_, s1x_,s1y_,s1z_,s2x_,s2y_,s2z_, 
	\[Iota]_, dL_, tc_,\[Phi]ref_,
	\[Delta]\[CurlyPhi]minus2_, \[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_
] := Module[
	{FpFc, fvec = Range[20, 1024, 0.125], hphc, opts},
	
	FpFc = PatternFunctions[fvec, \[Pi]/2-\[Theta], \[Phi], \[Psi], 0, DetectorVertex["H1"], DetectorTensor["H1"]];
	
	opts = Thread@Rule[
		{"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"},
		{\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4}
	];
	
	
	hphc = hphcIMRPhenomPv2[
		fvec, m1+m2, (m1 m2)/(m1+m2)^2, s1x,s1y,s1z,s2x,s2y,s2z, \[Iota], tc, \[Phi]ref, dL,
		Sequence@@opts
	];
	
	FpFc[[1]] hphc[[1]] + FpFc[[2]] hphc[[2]]
]


ivars = {
	"m1", "m2", "s1x","s1y","s1z","s2x", "s2y", "s2z",
	"\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "dL", 
	"tc", "\[Phi]ref",
	"\[Delta]\[CurlyPhi]-2","\[Delta]\[CurlyPhi]0","\[Delta]\[CurlyPhi]1","\[Delta]\[CurlyPhi]2","\[Delta]\[CurlyPhi]3","\[Delta]\[CurlyPhi]4","\[Delta]\[CurlyPhi]5l","\[Delta]\[CurlyPhi]6","\[Delta]\[CurlyPhi]6l","\[Delta]\[CurlyPhi]7","\[Delta]\[Beta]2","\[Delta]\[Beta]3","\[Delta]\[Alpha]2","\[Delta]\[Alpha]3","\[Delta]\[Alpha]4"
};


psd = (ASD["ET-D"]@Range[20, 1024, 0.125])^2;


H1 = <|
	"ASD" -> ASD["ET-D"],
	"Position" -> DetectorVertex["H1"],
	"DetectorTensor" -> DetectorTensor["H1"]
|>;


RandomSpin[] := Module[
	{norm, s},
	norm=2;
	While[norm>1,s =RandomReal[{-1,1},3 ]; norm=Norm[s]];
	s
]


Test["IMRPhenomPv2"] := Module[
	{
		m1, m2, s1x,s1y,s1z,s2x,s2y,s2z, \[Iota], tc, \[Phi]ref, dL, \[Theta], \[Phi], \[Psi], numerical, fp, Symbolic, 
		\[Delta]ps
	},
	
	{m1, m2} = RandomReal[{20, 80}, 2]//ReverseSort;
	
	{s1x,s1y,s1z} = RandomSpin[];
	{s2x,s2y,s2z} = RandomSpin[];
	
	{tc} = RandomReal[{-1,1},1];
	
	{\[Theta], \[Iota], \[Psi]} = RandomReal[{0, \[Pi]}, 3];
	{\[Phi]ref, \[Phi]} = RandomReal[{0, 2 \[Pi]},2];
	dL = RandomReal[{2, 10}];
	\[Delta]ps = RandomReal[{-10,10}, 15];
	
	fp = Association@(Thread@Rule[
		ivars,
		{m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Iota],\[Theta],\[Phi],\[Psi],dL,tc,\[Phi]ref,Sequence@@\[Delta]ps}
	]);
	
	
	numerical = NFisher[strain, {\[Theta],\[Phi],\[Psi],m1,m2,s1x,s1y,s1z,s2x,s2y,s2z,\[Iota],dL,tc,\[Phi]ref, Sequence@@\[Delta]ps}, 1, 30, psd];
	(*Echo[numerical];*)
	Symbolic = FisherMatrix["IMRPhenomPv2", fp, {H1}, "fmin"->20, "fmax"->1024, "res"->8033]//QuietEcho;
	
	numerical[[13, 14;;30]]//Echo;
	Symbolic[[13, 14;;30]]//Echo;
	(*
	these should be exactly zero since they are purelly imaginary contributions in the tensor product,
	we set them to zero, so they don't contaminate the comparison. The Symbolic implementation is usually closer
	to zero by some orders of magnitude
	*)
	Do[
		numerical[[13,i]] = 0; Symbolic[[13,i]] = 0; 
		numerical[[i,13]] =0;  Symbolic[[i,13]] =0,
		{i, 14, 30}
	];
	
	RelativeDiff[numerical, Symbolic]

]


a = Test["IMRPhenomPv2"]//UpperTriangularize;

Echo[{Position[Max[a]]@a, Max@a}];

a//MatrixForm

ListPlot[Flatten[a]//Sort, ScalingFunctions->"Log10", PlotRange->All]
(*the high differences seem to be always on the dL column*)


(* ::Section:: *)
(*Testing GWFAST vs Symbolic Fishers:*)


(* ::Text:: *)
(*OBSERVE THAT DELTA_T ENTERS WITH OPPOSITE SIGN IN THE DEFINITIONS NOW, (BCS bilby USES THE OPPPOSITE SIGN OF GWFAST. THIS WILL CAUSE DISAGREEMENTS ON THE SKYPOSITION-ELEMENTS)*)


(* ::Subsection:: *)
(*Fisher defs: *)


(* ::Subsubsection::Closed:: *)
(*GWFAST*)


DeleteObject/@ExternalSessions[]
Clear@python
python = StartExternalSession[{
	"Python", (*you should change the evaluator to your installation:*)
	(*"Evaluator"-> "/home/cosmo-ufes/anaconda3/envs/gwfast_env/bin/python"*)
	"Evaluator" -> "/Users/felipe/anaconda3/envs/GWFAST/bin/python"
}];


(*Standard stuff:*)
ExternalEvaluate[python,
"
import os
import sys

import copy
import numpy as onp
from astropy.cosmology import Planck18

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd())))
sys.path.append(SCRIPT_DIR)
import gwfast.gwfastGlobals as glob



"
]


(* ::Text:: *)
(*I will be using L1 and H1 with Aplus Design:*)


ExternalEvaluate[python,
(*Replace the path by your own: *)
"
alldetectors = copy.deepcopy(glob.detectors)
LVdetectors = {det:alldetectors[det] for det in ['L1', 'H1', 'Virgo']}
#LVdetectors['L1']['psd_path'] = '/home/cosmo-ufes/Documentos/GitHub/GWFORECAST/Data/ASD/AplusDesign.txt'
LVdetectors['L1']['psd_path'] = '/Users/felipe/Documents/GitHub/GWFORECAST/Data/ASD/AplusDesign.txt'
LVdetectors['H1']['psd_path'] = LVdetectors['L1']['psd_path']
"]


ExternalEvaluate[python,
"
from gwfast.signal import GWSignal
from gwfast.network import DetNet
from gwfast.waveforms import IMRPhenomD
from fisherTools import CovMatr, compute_localization_region, check_covariance, fixParams
"
]


(* ::Text:: *)
(*Function to define parameters : *)


ExternalEvaluate[python, "
myLVSignals = {}

for d in ['L1','H1']:

    myLVSignals[d] = GWSignal(IMRPhenomD(),
                psd_path=LVdetectors[d]['psd_path'],
                detector_shape = LVdetectors[d]['shape'],
                det_lat= LVdetectors[d]['lat'],
                det_long=LVdetectors[d]['long'],
                det_xax=LVdetectors[d]['xax'],
                verbose=False,
                useEarthMotion = False,
                fmin= 20.,
                fmax = 1024,
                IntTablePath=None,
                jitCompileDerivs=True
    )


myLVNet = DetNet(myLVSignals)
myLVNet.verbose=False
"]


PythonFisher = ExternalFunction[python, 
"
def Fisher(vec):
    Mc, dL, theta, phi, iota, psi, eta, phic, chi1z, chi2z = vec
    res = {
        'Mc': onp.array([Mc]),
        'dL': onp.array([dL]),
        'theta': onp.array([theta]),
        'phi': onp.array([phi]),
        'iota': onp.array([iota]) + onp.pi, #Accounts for the fact that GWFAST uses ''+I*Cos[iota]''
        'psi': onp.array([psi]),
        'tGPS': onp.array([61094.012]),
        'eta': onp.array([eta]),
        'Phicoal': onp.array([phic]),
        'chi1z': onp.array([chi1z]),
        'chi2z': onp.array([chi2z]),
    }
    #sys.stdout = open(os.devnull, 'w')
    fm = myLVNet.FisherMatr(res, res=1000)
    #sys.stdout = sys.__stdout__
    return fm
"]


(* ::Text:: *)
(*Order of elements in the Fisher matrix of IMRPhenomD: *)


ExternalEvaluate[python, "IMRPhenomD().ParNums"]


(* ::Subsubsection:: *)
(*MMA: *)


H1 = <|
	"Position"-> DetectorVertex["H1"],
	"DetectorTensor"-> DetectorTensor["H1"],
	"ASD" -> ASD["L1H1-O5"]
|>;

L1 = <|
	"Position"-> DetectorVertex["L1"],
	"DetectorTensor"-> DetectorTensor["L1"],
	"ASD" -> ASD["L1H1-O5"]
|>;


EmptyAssociation["Aligned"]


MMAFisherMatrix[fp_Association] := Module[
	{M, chirp, eta,delta, J, fm, value, gwfast, invdL, dL},
	
	{M, delta, invdL} = fp/@{"M", "\[Delta]", "1/dL"};
	eta = (1-delta^2)/4;
	dL = 1/invdL;
	chirp = M eta^(3/5);
	
	J = DiagonalMatrix[ConstantArray[1, 11]];
	
	(*
		to change \[Delta]->\[Eta], (D[\[Delta][\[Eta]], \[Eta]])
	*)
	J[[5,5]] = -(2/Sqrt[1-4 eta]);
	
	(*(\[Chi]s, \[Chi]a)->(s1z, s2z)*)
	{J[[6,6]], J[[6,7]], J[[7,6]], J[[7,7]]} = {1/2,1/2,1/2,-(1/2)};
	
	(*
		to change tc -> -tc
		2 \[Phi]ref -> - \[Phi]c
	*)
	J[[-2,-2]] = -1;
	J[[-1,-1]] = -1/2;
	
	(*change 1/dL->dL*)
	J[[9,9]] = -dL^-2;
	
	
	
	fm = DALITensors["IMRPhenomD", fp, {L1, H1}, 1, "fmin"->20, "fmax"->1024, "res"->8033]//QuietEcho;
	
	fm = ArrayReshape[fm, {11,11}];
	
	fm = J\[Transpose] . fm . J;
	
	Clear[value];
	
	MapThread[
		(value[#1] = #2)&,
	
		{
			Flatten[{\[Theta],\[Phi],\[Psi],Mc,\[Eta],s1z,s2z,\[Iota], dL, tc, \[Phi]ref}\[TensorProduct]{\[Theta],\[Phi],\[Psi],Mc,\[Eta],s1z,s2z,\[Iota], dL, tc, \[Phi]ref}],
			Flatten[fm]
		}
	];
	
	gwfast = Flatten[{Mc, \[Eta], dL, \[Theta], \[Phi], \[Iota], \[Psi], tc, \[Phi]ref, s1z, s2z}\[TensorProduct]{Mc, \[Eta], dL, \[Theta], \[Phi], \[Iota], \[Psi], tc, \[Phi]ref, s1z, s2z}];
	
	ArrayReshape[value/@gwfast, {11,11}]
	
]


(* ::Subsection::Closed:: *)
(*Test*)


(* ::Text:: *)
(*I believe my implementation of FpFc to be slightly different from GWFast`s one. So I should report the relative differences on the fisher matrix projected onto {M, \[Eta], \[Chi]1, \[Chi]2, tc, \[Phi]c, dl}*)


Clear@test


test := Module[
	{M, Mc, dL, theta, phi, iota, psi, tc, eta, delta,chis,chia, phic, chi1z, chi2z, pyvec, fp, python, MMA},
	
	M = RandomReal[{10, 120}];
	eta = RandomReal[{0.01, 0.2499}];
	delta = Sqrt[1 - 4 eta];
	Mc = M eta^(3/5);
	
	dL = RandomReal[{2, 10}];
	{theta, psi} = RandomReal[{0, \[Pi]},2];
	iota = RandomReal[{0, \[Pi]}];
	{phi, phic} = RandomReal[{0, 2 \[Pi]},2];
	{chi1z, chi2z} = RandomReal[{-1,1},2];
	{chis, chia} = {(chi1z+chi2z)/2, (chi1z-chi2z)/2};
	
	tc = 1126259462.;
	
	pyvec = {Mc, dL, theta, phi, iota, psi, eta, phic, chi1z, chi2z};
	
	fp =  (Thread@Rule[
		{"\[ScriptCapitalM]c","\[Delta]","\[Chi]s","\[Chi]a","\[Iota]","\[Theta]","\[Phi]","\[Psi]","1/dL","tc","\[Phi]ref"},
		{Mc, delta, chis, chia, iota, theta, phi, psi, 1/dL, tc, phic}
	])//Association;
	
	MMA = MMAFisherMatrix[fp];
	python = Normal[PythonFisher[pyvec][[All, All, 1]]];
	
	
	Echo[{MMA[[-3, 3]],MMA[[3,8]]}, "Should be 0 (MMA): "];
	
	Echo[{python[[-3, 3]],python[[3,8]]}, "Should be 0 (python): "];
	(*Enforce 0 on components that should be zero ([tc, dL] and [phiref, dL]), so they don't contaminate the comparison:*)
	MMA[[3, -3]] = MMA[[-3, 3]] = MMA[[3,8]] = MMA[[8,3]] = 0;
	python[[3, -3]] = python[[-3, 3]] = python[[3,8]] = python[[8,3]] = 0;
	
	
	RelativeDiff[python, MMA]
	
]


a = UpperTriangularize[test];
Echo[{Position[Max@a]@a, Max@a}];
a//MatrixForm

rd = Extract[a, SymmetrizedIndependentComponents[{11,11}, Symmetric[All]]];

ListPlot[Sort[rd], PlotRange->All, ScalingFunctions->"Log10"]


(* ::Section:: *)
(*Testing SNRs against GWFAST:*)


(* ::Subsection:: *)
(*SNR defs: *)


(* ::Subsubsection:: *)
(*GWFAST*)


DeleteObject/@ExternalSessions[]
Clear@python
python = StartExternalSession[{
	"Python", (*you should change the evaluator to your installation:*)
	"Evaluator"-> "/Users/felipe/anaconda3/envs/GWFAST/bin/python"
	(*"Evaluator"-> "/home/cosmo-ufes/anaconda3/envs/gwfast_env/bin/python"*)
}];


(*Standard stuff:*)
ExternalEvaluate[python,
"
import os
import sys

import copy
import numpy as onp
from astropy.cosmology import Planck18

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd())))
sys.path.append(SCRIPT_DIR)
import gwfast.gwfastGlobals as glob



"
]


(* ::Text:: *)
(*I will be using L1 and H1 with Aplus Design:*)


ExternalEvaluate[python,
(*Replace the path by your own: *)
"
alldetectors = copy.deepcopy(glob.detectors)
LVdetectors = {det:alldetectors[det] for det in ['L1', 'H1', 'ETSL', 'ETMRL45d', 'ETMRLpar']}

#LVdetectors['L1']['psd_path'] = '/home/cosmo-ufes/Documentos/GitHub/GWFORECAST/Data/ASD/AplusDesign.txt'
LVdetectors['L1']['psd_path'] = '/Users/felipe/Documents/GitHub/GWFORECAST/Data/ASD/AplusDesign.txt'

LVdetectors['H1']['psd_path'] = LVdetectors['L1']['psd_path']
LVdetectors['ETSL']['psd_path'] = LVdetectors['L1']['psd_path']
LVdetectors['ETMRL45d']['psd_path'] = LVdetectors['L1']['psd_path']
LVdetectors['ETMRLpar']['psd_path'] = LVdetectors['L1']['psd_path']
"]


ExternalEvaluate[python,
"
from gwfast.signal import GWSignal
from gwfast.network import DetNet
from gwfast.waveforms import IMRPhenomD
from fisherTools import CovMatr, compute_localization_region, check_covariance, fixParams
"
]


(* ::Text:: *)
(*Function to define parameters : *)


ExternalEvaluate[python, "
myLVSignals = {}

for d in ['L1','H1', 'ETSL', 'ETMRL45d', 'ETMRLpar']:

    myLVSignals[d] = GWSignal(IMRPhenomD(),
                psd_path=LVdetectors[d]['psd_path'],
                detector_shape = LVdetectors[d]['shape'],
                det_lat= LVdetectors[d]['lat'],
                det_long=LVdetectors[d]['long'],
                det_xax=LVdetectors[d]['xax'],
                verbose=False,
                useEarthMotion = False,
                fmin= 20.,
                fmax = 1024,
                IntTablePath=None,
                jitCompileDerivs=False
    )


myLVNet = DetNet(myLVSignals)
myLVNet.verbose=False
"]


PythonSNR = ExternalFunction[python, 
"
def SNR(vec):
    Mc, dL, theta, phi, iota, psi, eta, phic, chi1z, chi2z = vec
    res = {
        'Mc': onp.array([Mc]),
        'dL': onp.array([dL]),
        'theta': onp.array([theta]),
        'phi': onp.array([phi]),
        'iota': onp.array([iota]),
        'psi': onp.array([psi]),
        'tGPS': onp.array([61094.012]),
        'eta': onp.array([eta]),
        'Phicoal': onp.array([phic]),
        'chi1z': onp.array([chi1z]),
        'chi2z': onp.array([chi2z]),
    }
    sys.stdout = open(os.devnull, 'w')

    fm = myLVNet.SNR(res, res=1000, return_all=True)

    sys.stdout = sys.__stdout__

    return [
        fm['L1'].item(),
        fm['H1'].item(),
        fm['ETSL'].item(),
        fm['ETMRL45d'].item(),
		fm['ETMRLpar'].item()
    ]

"]


(* ::Text:: *)
(*Order of elements in the Fisher matrix of IMRPhenomD: *)


(* ::Subsubsection:: *)
(*MMA: *)


H1 = <|
	"Position"-> DetectorVertex["H1"],
	"DetectorTensor"-> DetectorTensor["H1"],
	"ASD" -> ASD["L1H1-O5"]
|>;

L1 = <|
	"Position"-> DetectorVertex["L1"],
	"DetectorTensor"-> DetectorTensor["L1"],
	"ASD" -> ASD["L1H1-O5"]
|>;


ETS = <|
	"Position"->DetectorVertex["ET-S"], 
	"DetectorTensor"->DetectorTensor["ET-S-L"],
	"ASD"->ASD["L1H1-O5"]
|>;

ETMR = <|
	"Position"->DetectorVertex["ET-MR"], 
	"DetectorTensor"->DetectorTensor["ET-MR-L-45"],
	"ASD"->ASD["L1H1-O5"]
|>;

ETMR2 = <|
	"Position"->DetectorVertex["ET-MR"], 
	"DetectorTensor"->DetectorTensor["ET-MR-L-0"],
	"ASD"->ASD["L1H1-O5"]
|>;


MMASNR[fp_Association] := Module[
	{M, chirp, eta, J, fm, value, gwfast},
	
	SNR["IMRPhenomD", fp, {L1, H1, ETS, ETMR, ETMR2}, "fmin"->20, "fmax"->1024, "res"->8033, "AllSNRs"->True]
]


fp  = Thread@Rule[
	{"\[ScriptCapitalM]c","\[Delta]","\[Chi]s","\[Chi]a", "\[Iota]", "\[Theta]","\[Phi]","\[Psi]", "1/dL", "tc","\[Phi]ref"},
	{20, 0.23, 0.1, -0.2, 2, 3, 1.4, 1.2, 2, 0.2, 0}
]//Association


(* ::Subsection:: *)
(*Test*)


(* ::Text:: *)
(**)


Clear@test


test := Module[
	{M, Mc, dL, theta, phi, iota,delta, psi,chis, chia, tc, eta, phic, chi1z, chi2z, pyvec, fp, python, MMA},
	
	M = RandomReal[{10, 120}];
	eta = RandomReal[{0.01, 0.2499}];
	delta = Sqrt[1-4 eta];
	Mc = M eta^(3/5);
	
	dL = RandomReal[{0.05, 10}];
	{theta, psi} = RandomReal[{0, \[Pi]},2];
	iota = RandomReal[{0, \[Pi]}];
	{phi, phic} = RandomReal[{0, 2 \[Pi]},2];
	{chi1z, chi2z} = RandomReal[{-1,1},2];
	{chis, chia} = {(chi1z+chi2z)/2, (chi1z-chi2z)/2};
	
	tc = 61094.012;
	
	(*GWFAST uses the opposite sign for Cos(iota) in hc, that is why "iota+\[Pi]" below*)
	pyvec = {Mc, dL, theta, phi, iota+\[Pi], psi, eta, phic, chi1z, chi2z};
	
	fp =  (Thread@Rule[
		{"\[ScriptCapitalM]c","\[Delta]","\[Chi]s","\[Chi]a","\[Iota]","\[Theta]","\[Phi]","\[Psi]","1/dL","tc","\[Phi]ref"},
		{Mc, delta, chis, chia, iota, theta, phi, psi, 1/dL, tc, phic}
	])//Association;
	
	MMA = MMASNR[fp];
	
	python = PythonSNR[pyvec];
	RelativeDiff[python, MMA]
	
]


(* ::Text:: *)
(*Agreement overall better than 5%*)


test


Table[test, 1000]//MinMax


(* ::Section::Closed:: *)
(*Testing Population functionality*)


Names["FelipeBarbosa`SymDALI`Population`*"]


(* ::Subsection::Closed:: *)
(*ComovingDistance*)


l1 = ExternalEvaluate["Python", "
import numpy as np

from astropy.cosmology import Planck18 as cosmo

z = np.linspace(0, 20, 10000) 
d_c = cosmo.comoving_distance(z)
d_c.value
"]//Normal;


Clear@dc
dc = ComovingDistance["\[CapitalOmega]m" -> 0.30966] (*astropy uses this for \[CapitalOmega]m*)


l2 = dc@Range[0, 20, 20/9999.];


With[
	{l1 = RelativeDiff[l1, l2]},
	(*The first element is 1 bcs astropy evaluates to 0 at z=0 and MMA to 10^-21, but otherwhise the agreement is nice*)
	ListPlot[l1, ScalingFunctions->"Log10", Background->White, PlotRange->All, DataRange->{0, 20}]
]


(* ::Text:: *)
(*We can try with different cosmological parameters: *)


l1 = ExternalEvaluate["Python", "
import numpy as np

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.2, Tcmb0=2.725)


z = np.linspace(0, 20, 10000) 
d_c = cosmo.comoving_distance(z)
d_c.value
"]//Normal;


Clear@dc
dc = ComovingDistance["\[CapitalOmega]m" -> 0.2, "H0" -> 70]


l2 = dc@Range[0, 20, 20/9999];


(* ::Text:: *)
(*Same situation:*)


With[
	{l1 = RelativeDiff[l1, l2]},
	
	ListPlot[l1, ScalingFunctions->"Log10", Background->White, PlotRange->All, DataRange->{0, 20}]
]


Clear[l1, l2, dc]


(* ::Subsection::Closed:: *)
(*DifferentialComovingVolume*)


l1 = ExternalEvaluate["Python", "
import numpy as np

from astropy.cosmology import Planck18 as cosmo

z = np.linspace(0, 20, 10000) 
d_c = cosmo.differential_comoving_volume(z)
d_c.value
"]//Normal;


Clear@dVdz
dVdz = DifferentialComovingVolume["\[CapitalOmega]m"-> 0.30966] (*Again, adjusting to the value used by astropy:*)


l2 = dVdz@Range[0, 20, 20/9999];

l2 = l2/(4 \[Pi]); (*astropy returns Mpc/steradian*) 


With[
	{l1 = RelativeDiff[l1, l2]},
	(*The first element is 1 bcs astropy evaluates to 0 at z=0 and MMA to 10^-26, but otherwhise the agreement is nice*)
	ListPlot[l1, ScalingFunctions->"Log10", Background->White, PlotRange->All, DataRange->{0, 20}]
]


(* ::Text:: *)
(*Again, we can test with a different cosmology: *)


l1 = ExternalEvaluate["Python", "
import numpy as np

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.2, Tcmb0=2.725)

z = np.linspace(0, 20, 10000) 
d_c = cosmo.differential_comoving_volume(z)
d_c.value
"]//Normal;


Clear@dVdz

dVdz = DifferentialComovingVolume["\[CapitalOmega]m"->0.2, "H0"->70]


l2 = dVdz@Range[0, 20, 20/9999];
l2 = l2/(4 \[Pi]);


With[
	{l1 = RelativeDiff[l1, l2]},
	(*The first element is 1 bcs astropy evaluates to 0 at z=0 and MMA to 10^-26, but otherwhise the agreement is nice*)
	ListPlot[l1, ScalingFunctions->"Log10", Background->White, PlotRange->All, DataRange->{0, 20}]
]


Clear[l1, l2, dVdz]


(* ::Subsection:: *)
(*MadauDickinsonProfile*)


l1 = ExternalEvaluate["Python", "
import numpy as np
import gwpopulation.models.redshift as gwpop_redshift
import matplotlib.pyplot as plt

z = np.linspace(0, 20,10000)

# Madau-Dickinson star formation rate density from GWPopulation
sfr = gwpop_redshift.MadauDickinsonRedshift.psi_of_z(
	self = gwpop_redshift.MadauDickinsonRedshift, 
	redshift = z, 
	gamma=2.7, 
	kappa=5.7, 
	z_peak=2
)
sfr
"]//Normal;


l2 = MadauDickinsonProfile[Range[0,20, 20/9999], "R0"->1];


With[ (*overall good agreement:*)
	{l = RelativeDiff[l1, l2]//Sort},
	ListPlot[l//Sort, Background->White, ScalingFunctions->"Log10", DataRange->{0, 20}, PlotRange->All]
]


(* ::Text:: *)
(*We can check other parameters: *)


l1 = ExternalEvaluate["Python", "
import numpy as np
import gwpopulation.models.redshift as gwpop_redshift
import matplotlib.pyplot as plt

z = np.linspace(0, 20,10000)

# Madau-Dickinson star formation rate density from GWPopulation
sfr = gwpop_redshift.MadauDickinsonRedshift.psi_of_z(
	self = gwpop_redshift.MadauDickinsonRedshift, 
	redshift = z, 
	gamma=3.7, 
	kappa=5.7, 
	z_peak=3
)
sfr
"]//Normal;


l2 = MadauDickinsonProfile[Range[0,20, 20/9999], "R0"->1, "\[Alpha]z"-> 3.7, "\[Beta]z"->2, "zp"->3];


With[ (*overall good agreement:*)
	{l = RelativeDiff[l1, l2]//Sort},
	ListPlot[l//Sort, Background->White, ScalingFunctions->"Log10", DataRange->{0, 20}, PlotRange->All]
]


Clear[l1, l2]


?MadauDickinsonProfile


NIntegrate[]


(* ::Subsection:: *)
(*PowerLaw + Peak*)


<<FelipeBarbosa`SymDALI`


Options[PowerLawPlusPeak]


PowerLawPlusPeak[m1, q]


P[m1_] := NIntegrate[PowerLawPlusPeak[m1, q], {q, 0.1,1}]


DP[m1_] := (P[m1 + m1 10^-6]  - P[m1])/(m1 10^-6);


NIntegrate[MadauDickinsonProfile[z], {z,0,2}]


Plot[P[m1], {m1, 5.1, 87}, ScalingFunctions->"Log10"]//Quiet





(* ::Section::Closed:: *)
(*Testing divergence of the Derivative*)


FelipeBarbosa`GWFORECAST`Fisher`Private`NRules["IMRPhenomD"][FelipeBarbosa`GWFORECAST`Fisher`Private`\[CapitalPhi]IMR][[3,3]]


d\[CapitalPhi][f_, M_, \[Eta]_, \[Chi]1_, \[Chi]2_] := FelipeBarbosa`GWFORECAST`Fisher`Private`NRules["IMRPhenomD"][FelipeBarbosa`GWFORECAST`Fisher`Private`\[CapitalPhi]IMR][[3,3]][
	{f M 4.93 10^-6}, 4.93 10^-6 M 10, \[Eta], \[Chi]1, \[Chi]2,
	Sequence@@ConstantArray[0, 15]
]


d\[ScriptCapitalA][f_, M_, \[Eta]_, \[Chi]1_, \[Chi]2_] := FelipeBarbosa`GWFORECAST`Fisher`Private`NRules["IMRPhenomD"][FelipeBarbosa`GWFORECAST`Fisher`Private`\[ScriptA]IMR][[3, 3]][
 {f}, M, \[Eta], \[Chi]1, \[Chi]2, 0
]


FelipeBarbosa`GWFORECAST`Fisher`Private`NRules["IMRPhenomD"][FelipeBarbosa`GWFORECAST`Fisher`Private`\[ScriptA]IMR][[1,2,0]]


\[ScriptCapitalA][f_, M_, \[Eta]_, \[Chi]1_, \[Chi]2_] := FelipeBarbosa`GWFORECAST`Fisher`Private`NRules["IMRPhenomD"][FelipeBarbosa`GWFORECAST`Fisher`Private`\[ScriptA]IMR][[1,2,0]][
 {f}, M, \[Eta], \[Chi]1, \[Chi]2, 0
]


Manipulate[
Plot[
	d\[CapitalPhi][f, M, \[Eta], \[Chi]1, \[Chi]2], 
	{\[Eta], 0.24999, 0.25},
	PlotRange->All
], {f, 11, 1024}, {M, 10, 100}, {\[Chi]1, -1, 1}, {\[Chi]2, -1, 1}]


Manipulate[
Plot[
	d\[ScriptCapitalA][f, M, \[Eta], \[Chi]1, \[Chi]2][[1]]//Re, 
	{\[Eta], 0.24999, 0.25},
	PlotRange->All
], {f, 11, 1024}, {M, 10, 100}, {\[Chi]1, -1, 1}, {\[Chi]2, -1, 1}]


Manipulate[
Plot[
	{Re[d\[ScriptCapitalA][f, M, \[Eta], \[Chi]1, \[Chi]2][[1]]], Re[\[ScriptCapitalA][f, M, \[Eta], \[Chi]1, \[Chi]2][[1]]] d\[CapitalPhi][f, M, \[Eta], \[Chi]1, \[Chi]2]}, 
	{\[Eta], 0.24, 0.25},
	
	PlotLegends->{"da", "dphi"}
], {f, 11, 1024}, {M, 10, 100}, {\[Chi]1, -1, 1}, {\[Chi]2, -1, 1}]


(* ::Section::Closed:: *)
(*Testing ProcessDALITensors:*)


vars = {"\[ScriptCapitalM]c", "\[Delta]", "\[Chi]s", "\[Chi]a","\[Iota]", "\[Theta]", "\[Phi]", "\[Psi]", "1/dL", "tc", "\[Phi]ref"};

fp = MapThread[
	Rule,
	{
		vars,
		{2.3, 0.2, 0.9, -0.9, 2.3, 3.1, 2.1, 2.9, 1/3., 0.1, 6.}
	}
]//Association;


H1 = <|
"Position"->DetectorVertex["H1"],
"DetectorTensor"->DetectorTensor["H1"],
"ASD"-> ASD["L1H1-O5"]
|>;
L1 = <|
"Position"->DetectorVertex["L1"],
"DetectorTensor"->DetectorTensor["L1"],
"ASD"-> ASD["L1H1-O5"]
|>;
V1 = <|
"Position"->DetectorVertex["V1"],
"DetectorTensor"->DetectorTensor["V1"],
"ASD"-> ASD["V1-O5"]
|>;


(* ::Input:: *)
(*dali = DALITensors["IMRPhenomD", fp, {H1, L1, V1}, 3, "res"->8000, "fmin"->10, "fmax"->1024];*)


(* ::Subsection:: *)
(*Test Function*)


LITensorProduct[vec1_?VectorQ, vec2_?VectorQ] := Module[
	{l1 = Length@vec1, l2 = Length@vec2},
	If[
		l1===l2, 
		Table[vec1[[i]]*vec2[[j]], {i,l1}, {j,i,l1}]//Flatten,
		vec1\[TensorProduct]vec2
	]
]


MAKEDALI\[CapitalDelta]p[vec_] := Module[
	{dim = Length@vec, \[CapitalDelta]p, DALIorder},
	
	\[CapitalDelta]p[1] = vec;
	\[CapitalDelta]p[2] = Extract[vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim}, Symmetric[All]]];
	\[CapitalDelta]p[3] = Extract[vec\[TensorProduct]vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim, dim}, Symmetric[All]]];
	
	DALIorder = Do[
		Sow[#, j]&@(LITensorProduct[\[CapitalDelta]p[j], \[CapitalDelta]p[i]]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	DALIorder//Flatten
]


Clear@Test


Test[dim_] := Module[
	{LIComponents, LIValues, SymmetricGradients, rules, \[CapitalDelta]p, DALIscheme1, DALIscheme2, scalar1, dali\[CapitalDelta]p, scalar2},
	
	(LIComponents[#] = SymmetrizedIndependentComponents[ConstantArray[dim, #],Symmetric[All]])&/@Range[3];
	
	(LIValues[#] = RandomReal[{1,10}, Length[LIComponents[#]]])&/@Range[3];
	
	(rules[#] = MapThread[Rule, {LIComponents[#], LIValues[#]}])&/@Range[3];
	
	(SymmetricGradients[#] = SymmetrizedArray[rules[#], ConstantArray[dim, #], Symmetric[All]])&/@Range[3];
	
	DALIscheme1 = Do[
		Sow[#, j]&@(Flatten[LIValues[j]\[TensorProduct]LIValues[i]]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	
	DALIscheme2 = Do[
		Sow[#, j]&@(SymmetricGradients[j]\[TensorProduct]SymmetricGradients[i]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	scalar1 = Table[
		Dot[c[i,j]*DALIscheme2[[i,j]], Sequence@@ConstantArray[\[CapitalDelta]p, i+j]],
		{i,3},
		{j,1,i}
	]//Flatten[#, 1]&;
	
	\[CapitalDelta]p = RandomReal[{1,2}, dim];
	
	dali\[CapitalDelta]p = MAKEDALI\[CapitalDelta]p[\[CapitalDelta]p];
	
	scalar2 = Flatten[ProcessDALITensors[DALIscheme1]] . dali\[CapitalDelta]p;
	
	((Plus@@scalar1) -  scalar2)/Min[{scalar1, scalar2}]//Abs
	
]


(* ::Subsection::Closed:: *)
(*Calculate*)


MaxMemoryUsed[Test[9]]


(*hard to go beyond dim=9 in my computer, too much ram to run this.*)


(* ::Text:: *)
(*I am finding agreement to numerical precision pretty much*)


Table[Test[9], 5]


<<FelipeBarbosa`SymDALI`


(* ::Section::Closed:: *)
(*Testing that only totally symmetric part of a tensor survives the contraction*)


LITensorProduct[vec1_?VectorQ, vec2_?VectorQ] := Module[
	{l1 = Length@vec1, l2 = Length@vec2},
	If[
		l1===l2, 
		Table[vec1[[i]]*vec2[[j]], {i,l1}, {j,i,l1}]//Flatten,
		vec1\[TensorProduct]vec2
	]
]


MAKEDALI\[CapitalDelta]p//ClearAll

MAKEDALI\[CapitalDelta]p[vec_] := Module[
	{dim = Length@vec, \[CapitalDelta]p, DALIorder},
	
	\[CapitalDelta]p[1] = vec;
	\[CapitalDelta]p[2] = Extract[vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim}, Symmetric[All]]];
	\[CapitalDelta]p[3] = Extract[vec\[TensorProduct]vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim, dim}, Symmetric[All]]];
	
	DALIorder = Do[
		Sow[#, j]&@(LITensorProduct[\[CapitalDelta]p[j], \[CapitalDelta]p[i]]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	DALIorder//Flatten
]


MAKEDALI\[CapitalDelta]p//ClearAll

MAKEDALI\[CapitalDelta]p[vec_] := Module[
	{dim = Length@vec, \[CapitalDelta]p, DALIorder},
	
	\[CapitalDelta]p[1] = vec;
	\[CapitalDelta]p[2] = Extract[vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim}, Symmetric[All]]];
	\[CapitalDelta]p[3] = Extract[vec\[TensorProduct]vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim, dim}, Symmetric[All]]];
	
	DALIorder = Do[
		Sow[#, j]&@(LITensorProduct[\[CapitalDelta]p[j], \[CapitalDelta]p[i]]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	DALIorder//Flatten
]


Clear@Test
Test[dim_] := Module[
	{LIComponents, LIValues, SymmetricGradients, rules, \[CapitalDelta]p, DALIscheme1, DALIscheme2, scalar1, dali\[CapitalDelta]p, scalar2},
	
	(LIComponents[#] = SymmetrizedIndependentComponents[ConstantArray[dim, #],Symmetric[All]])&/@Range[3];
	
	(LIValues[#] = RandomReal[{1,10}, Length[LIComponents[#]]])&/@Range[3];
	
	(rules[#] = MapThread[Rule, {LIComponents[#], LIValues[#]}])&/@Range[3];
	
	(SymmetricGradients[#] = SymmetrizedArray[rules[#], ConstantArray[dim, #], Symmetric[All]])&/@Range[3];
	
	(*Way I am calculating the DALI Likelihood*)
	DALIscheme1 = Do[
		Sow[#, j]&@(Flatten[LIValues[j]\[TensorProduct]LIValues[i]]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	(*Divide by c[i,j] bcs ProcessDALITensors multiplies by it and you want to include only the multiplicities:*)
	DALIscheme1 = Table[
		DALIscheme1[[i,j]]/c[i,j],
		{i, 1, 3},
		{j,1,i}
	];
		
	
	(*The second scheme is to take only the totally symmetric part of the tensor products and contract 
	with the totally symmetric part of the Deltaps*)
	DALIscheme2 = Do[
		Sow[#, j]&@(SymmetricGradients[j]\[TensorProduct]SymmetricGradients[i]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	scalar1 = Table[
		Dot[(*Symmetrize will take the Totally symmetric part of the tensor*)
			Symmetrize[DALIscheme2[[i,j]], Symmetric[All]], 
			Sequence@@ConstantArray[\[CapitalDelta]p, i+j]
		],
		{i,3},
		{j,1,i}
	]//Flatten[#, 1]&;
	
	\[CapitalDelta]p = RandomReal[{1,2}, dim];
	
	dali\[CapitalDelta]p = MAKEDALI\[CapitalDelta]p[\[CapitalDelta]p];
	
	scalar2 = Flatten[ProcessDALITensors[DALIscheme1]] . dali\[CapitalDelta]p;
	
	((Plus@@scalar1) -  scalar2)/Min[{scalar1, scalar2}]//Abs
	
]


Table[Test[8], 10]


(* ::Section::Closed:: *)
(*Testing "TaylorForm"*)


(* ::Subsection::Closed:: *)
(*Test Function*)


LITensorProduct[vec1_?VectorQ, vec2_?VectorQ] := Module[
	{l1 = Length@vec1, l2 = Length@vec2},
	If[
		l1===l2, 
		Table[vec1[[i]]*vec2[[j]], {i,l1}, {j,i,l1}]//Flatten,
		vec1\[TensorProduct]vec2
	]
]


MAKEDALI\[CapitalDelta]p[vec_] := Module[
	{dim = Length@vec, \[CapitalDelta]p, DALIorder},
	
	\[CapitalDelta]p[1] = vec;
	\[CapitalDelta]p[2] = Extract[vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim}, Symmetric[All]]];
	\[CapitalDelta]p[3] = Extract[vec\[TensorProduct]vec\[TensorProduct]vec, SymmetrizedIndependentComponents[{dim, dim, dim}, Symmetric[All]]];
	
	DALIorder = Do[
		Sow[#, j]&@(LITensorProduct[\[CapitalDelta]p[j], \[CapitalDelta]p[i]]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	DALIorder//Flatten
]


MAKELI\[CapitalDelta]p[vec_] := Module[
	{\[CapitalDelta]p2, \[CapitalDelta]p3, \[CapitalDelta]p4, \[CapitalDelta]p5, \[CapitalDelta]p6, dim  = Length@vec},
	
	\[CapitalDelta]p2 = Table[
		vec[[i]]*vec[[j]],
		{i,dim},
		{j, i, dim}
	]//Flatten;
	
	\[CapitalDelta]p3 = Table[
		vec[[i]]*vec[[j]]*vec[[k]],
		{i,dim},
		{j, i, dim},
		{k, j, dim}
	]//Flatten;
	
	\[CapitalDelta]p4 = Table[
		vec[[i]]*vec[[j]]*vec[[k]]*vec[[l]],
		{i,dim},
		{j, i, dim},
		{k, j, dim},
		{l, k, dim}
	]//Flatten;
	
	\[CapitalDelta]p5 = Table[
		vec[[i]]*vec[[j]]*vec[[k]]*vec[[l]]*vec[[s]],
		{i,dim},
		{j, i, dim},
		{k, j, dim},
		{l, k, dim},
		{s, l, dim}
	]//Flatten;
	
	\[CapitalDelta]p6 = Table[
		vec[[i]]*vec[[j]]*vec[[k]]*vec[[l]]*vec[[s]]*vec[[q]],
		{i,dim},
		{j, i, dim},
		{k, j, dim},
		{l, k, dim},
		{s, l, dim},
		{q, s, dim}
	]//Flatten;
	
	Join[\[CapitalDelta]p2, \[CapitalDelta]p3, \[CapitalDelta]p4, \[CapitalDelta]p5, \[CapitalDelta]p6]
]


Clear@Test


Test[dim_] := Module[
	{LIComponents, LIValues, \[CapitalDelta]p, DALIscheme1,taylor,  scalar1, dali\[CapitalDelta]p, scalar2},
	
	(LIComponents[#] = SymmetrizedIndependentComponents[ConstantArray[dim, #],Symmetric[All]])&/@Range[3];
	
	(LIValues[#] = RandomReal[{1,10}, Length[LIComponents[#]]])&/@Range[3];
	
	DALIscheme1 = Do[
		Sow[#, j]&@(Flatten[LIValues[j]\[TensorProduct]LIValues[i]]),
		{i, 3},
		{j, i, 3}
	]//Reap//Last;
	
	taylor = TaylorForm[DALIscheme1]//EchoTiming;
	
	\[CapitalDelta]p = RandomReal[{1,2}, dim];
	
	dali\[CapitalDelta]p = MAKEDALI\[CapitalDelta]p[\[CapitalDelta]p];
	
	scalar2 = Flatten[ProcessDALITensors[DALIscheme1]] . dali\[CapitalDelta]p;
	
	scalar1 = Flatten[taylor] . MAKELI\[CapitalDelta]p[\[CapitalDelta]p];
	
	(scalar1 -  scalar2)/Min[{scalar1, scalar2}]//Abs
	
]


(* ::Subsection::Closed:: *)
(*Calculate*)


(* ::Text:: *)
(*I believe it takes so long bcs of ```SymmetrizedArray``` and ```Symmetrize``` functionality in Mathematica. They seem to be intrinsically slow. Not much I can I guess.*)


MaxMemoryUsed[Test[12]]


(* ::Text:: *)
(*I am finding agreement to numerical precision pretty much*)


Test[9]//EchoTiming


Table[Test[12], 5]


(* ::Section::Closed:: *)
(*Test \[Delta]\[CurlyPhi] against modified RippleGW: *)


(* ::Subsection::Closed:: *)
(*Ripple Phase: *)


DeleteObject/@ExternalSessions[];
Clear@python
python = StartExternalSession[{
	"Python",
	"Evaluator"-> "/Users/felipe/anaconda3/envs/modified_ripplegw/bin/python"
}];

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


(* ::Text:: *)
(* Phase (f : Array, theta : Array, coeffs : Array, transition_freqs : Array) -> Array :*)
(*  *)


(* ::Subsection::Closed:: *)
(*MMA def*)


Private`D\[CapitalPhi]IMR


(* ::Subsection::Closed:: *)
(*Test*)


Test := Module[
	{f, M, \[Chi]1, \[Chi]2,pos, m1, m2,\[Eta],tc, \[Phi]c, Ripple, MMA, \[Omega], \[Theta]in, \[Theta]ex, coeffs,  \[Delta], \[Chi]s, \[Chi]a,
	G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]],diff, transition,
	ringdown,\[Omega]ref,\[Delta]\[CurlyPhi]s},
	
	\[Delta]\[CurlyPhi]s = RandomReal[{-1,1}, 10];
	
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

	Ripple = Argument[f, Join[\[Theta]in, \[Delta]\[CurlyPhi]s], \[Theta]ex, coeffs, 20.];
	
	
	MMA = Private`D\[CapitalPhi]IMR[\[Omega],\[Omega]ref, \[Delta], \[Chi]s, \[Chi]a, Sequence@@\[Delta]\[CurlyPhi]s, 0,0,0,0,0];
	
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
			PlotRange->All, ImageSize->Medium, Background->White,
			ScalingFunctions->"Log10"
		],
		ListLinePlot[
			{Take[Ripple, pos], Take[MMA, pos]}, 
			GridLines->{{{0.018,Red}, {ringdown/2, Red}}, None},
			PlotRange->All,
			PlotLegends->{"Python", "MMA"}, ImageSize->Medium, Background->White
		]
	}
	

]


Test
