(* ::Package:: *)

$HistoryLength=0;


PacletDirectoryLoad[NotebookDirectory[]//ParentDirectory[#,2]&];


<<FelipeBarbosa`SymDALI`


v = Range[20., 1024., 0.125] UnitConvert["GravitationalConstant"/("SpeedOfLight")^3, "Seconds"/"SolarMass"][[1]] 30;


vref = UnitConvert["GravitationalConstant"/("SpeedOfLight")^3, "Seconds"/"SolarMass"][[1]] 30 20;


NRules["IMRPhenomD"]["\[CapitalPhi]IMR"][[20,3]][
	v, vref, 0.2, 0.1, -0.9,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
];//RepeatedTiming//ScientificForm


0.00105732421875 (*Highest function evaluation time*)


(* ::Text:: *)
(*Considering 11 heavy functions (4 derivatives of \[ScriptA]IMR , 5 of \[CapitalPhi]IMR, plus both functions without derivatives) I should aim for derivative times on the core approximant smaller than 12 ms, with 8 000 frequency points. So, 15 ms Fisher matrix calculation for 5/6 detectors should be feasible.*)


RepeatedTiming[v*v;]


(* ::Section::Closed:: *)
(*Util functions*)


(*function to calculate the relative differences between numbers*)
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


(* ::Subsection::Closed:: *)
(*Test hp and hc*)


(* ::Text:: *)
(*It is easier to just run the section "functions to generate plots" and go to "check plots" to see the plots*)


(* ::Subsubsection::Closed:: *)
(*functions to generate plots*)


Clear@testhp

testhp := Module[
	{
		m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, \[Phi]Ref, \[Iota], f, MMA, lalD, lalR, MMAR, lalIm, MMAIm, fmax,\[Eta],
		G = UnitConvert[("GravitationalConstant")/("SpeedOfLight")^3, ("Seconds")/("SolarMass")][[1]], 
		 Rediff, Imdiff, f1, sp, irrelevant, rePlot, imPlot, reDiffPlot, imDiffPlot
	},
	
	{m1, m2} = ReverseSort[RandomReal[{20, 100},2]];
	\[Eta] = (m1 m2)/(m1+m2)^2;
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	fmax =  0.2/(G (m1+m2))//Round;
	f1 = 0.014/(G (m1+m2));
	
	
	MMA =  hphcIMRPhenomD[f, m1+m2, \[Eta], s1z, s2z, \[Iota], 0, \[Phi]Ref, 1][[1]];
	lalD = lal[m1, m2, s1z,s2z, 1,\[Iota],\[Phi]Ref, 1., 10., fmax, 10., "IMRPhenomD"][[1]];
	
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
	{s1z, s2z} = RandomReal[{-1,1}, 2];
	
	\[Phi]Ref = RandomReal[{0, 2 \[Pi]}];
	\[Iota] = RandomReal[{0, \[Pi]}];
	
	f = Range[10., 0.2/(G (m1+m2)), 1.];
	fmax =  0.2/(G (m1+m2))//Round;
	f1 = 0.014/(G (m1+m2));
	
	
	MMA =  hphcIMRPhenomD[f, m1+m2, \[Eta], s1z, s2z, \[Iota], 0, \[Phi]Ref, 1][[2]];
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
	
	
	MMA =  hphcIMRPhenomPv2[f, m1+m2, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], 0, \[Phi]Ref, 1][[1]];
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
	
	
	MMA =  hphcIMRPhenomPv2[f, m1+m2, \[Eta], s1x, s1y, s1z, s2x, s2y, s2z, \[Iota], 0, \[Phi]Ref, 1][[2]];
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
