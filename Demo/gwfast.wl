(* ::Package:: *)

(* ::Input:: *)
(*Quit*)


(* ::Input:: *)
(*$HistoryLength=1;*)


(* ::Section:: *)
(*Python FM*)


python = StartExternalSession["Python"];


ExternalEvaluate[python, "
import jax
import os
import sys

import copy
import numpy as onp
from astropy.cosmology import Planck18

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd())))
sys.path.append(SCRIPT_DIR)
import gwfast.gwfastGlobals as glob
"]


ExternalEvaluate[python, "
alldetectors = copy.deepcopy(glob.detectors)
print('All available detectors are: '+str(list(alldetectors.keys())))

# select only LIGO and Virgo
LVdetectors = {det:alldetectors[det] for det in ['L1', 'H1', 'Virgo']}
print('Using detectors '+str(list(LVdetectors.keys())))
"]


ExternalEvaluate[python,"
# We use the O2 psds
LVdetectors['L1']['psd_path'] = os.path.join(glob.detPath, 'LVC_O1O2O3', '2017-08-06_DCH_C02_L1_O2_Sensitivity_strain_asd.txt')
LVdetectors['H1']['psd_path'] = os.path.join(glob.detPath, 'LVC_O1O2O3', '2017-06-10_DCH_C02_H1_O2_Sensitivity_strain_asd.txt')
"]


ExternalEvaluate[python, "
from gwfast.waveforms import IMRPhenomD
from gwfast.signal import GWSignal
from gwfast.network import DetNet
from fisherTools import CovMatr, compute_localization_region, check_covariance, fixParams
"]


ExternalEvaluate[python, "
myLVSignals = {}

for d in ['L1', 'H1']:

    myLVSignals[d] = GWSignal(IMRPhenomD(),
                psd_path=LVdetectors[d]['psd_path'],
                detector_shape = LVdetectors[d]['shape'],
                det_lat= LVdetectors[d]['lat'],
                det_long=LVdetectors[d]['long'],
                det_xax=LVdetectors[d]['xax'],
                verbose=True,
                useEarthMotion = False,
                fmin= 20.,
                fmax = 1024.,
                IntTablePath=None,
                jitCompileDerivs=True)

myLVNet = DetNet(myLVSignals)
"]


ExternalEvaluate[python, "
from gwfastUtils import GPSt_to_LMST

# Median values of the posterior samples for all the parameters,
# except psi and the coalescence phase that are set to 0

z = onp.array([0.0980])
tGPS = onp.array([61094.012])

test = {'Mc':onp.array([1.1859])*(1.+z),
            'dL':onp.array([1]), #Planck18.luminosity_distance(z).value/1000,
            'theta':onp.array([0.4080839999999999]),
            'phi':onp.array([3.4461599999999994]),
            'iota':onp.array([2.545065595974997]),
            'psi':onp.array([2.3]),
            'tcoal':onp.array([0.]), # GMST is LMST computed at long = 0\[Degree]
            'eta':onp.array([0.24786618323504223]),
            'Phicoal':onp.array([0]),
            'chi1z':onp.array([0.005136138323169717]),
            'chi2z':onp.array([0.003235146993487445])
        }

print('Parameters for test are:')
test
"]


ExternalEvaluate[python,"
import time
time0 = time.time()
totF = myLVNet.FisherMatr(test,  use_chi1chi2=True, computeAnalyticalDeriv=True, res = 1000)
time1 = time.time()
print(time1-time0)
print('The computed Fisher matrix has shape %s'%str(totF.shape))
"]


py = ExternalEvaluate[python, "totF"]//Normal;
py = ArrayReshape[py, {11,11}]


(* ::Section:: *)
(*Comparison*)


NotebookDirectory[]//ParentDirectory


PacletDirectoryLoad[NotebookDirectory[]//ParentDirectory[#,2]&];
<<FelipeBarbosa`SymDALI`


{DSymRules, DNRules} = DerivativeRulesLoad["Detectors"];
{SymRules,NRules} = DerivativeRulesLoad["IMRPhenomD"];
NRules = Join[NRules, DNRules];
SymRules = Join[SymRules, DSymRules];


Block[
{name =NotebookDirectory[]//ParentDirectory, rule},
name = FileNameJoin[{name, "/LibraryResources", $SystemID, "DerivativeRules/IMRPhenomD/Defs.wdx"}];
rule =Import[name];
h[Mc_,\[Eta]_,\[Chi]1_,\[Chi]2_, \[Delta]\[CurlyPhi]minus2_,\[Delta]\[CurlyPhi]0_,\[Delta]\[CurlyPhi]1_,\[Delta]\[CurlyPhi]2_,\[Delta]\[CurlyPhi]3_,\[Delta]\[CurlyPhi]4_,\[Delta]\[CurlyPhi]5_,\[Delta]\[CurlyPhi]5l_,\[Delta]\[CurlyPhi]6_,\[Delta]\[CurlyPhi]6l_,\[Delta]\[CurlyPhi]7_,\[Delta]\[Beta]2_,\[Delta]\[Beta]3_,\[Delta]\[Alpha]2_,\[Delta]\[Alpha]3_,\[Delta]\[Alpha]4_, fref_, f_] = rule[[2]]//.({ M->Mc \[Eta]^(-3/5)});
]


auxh[tc_, \[Phi]c_, dL_, gpsTime_, f_]  = Exp[I (2 \[Pi] f (tc + gpsTime) - \[Phi]c)]/dL;


Clear[s1];

s1[\[Theta]_, \[Phi]_, \[Psi]_, \[Iota]_, p1_,p2_,p3_,D11_,D12_,D13_,D22_,D23_,D33_, f_] =  S[
	f,\[Theta],\[Phi],  \[Psi], -Cos[\[Iota]],p1,p2,p3,D11,D12,D13,D22,D23,D33 (*
	The signal adopted for hc in GWFAST seems wrong to me, but I can't check for sure bcs 
	I can't find stuff there. My convention is easy to check and it seems right
*)
]


desktopPrefix = "/home/cosmo-ufes/anaconda3/lib/python3.11";
laptopPrefix = "/Users/felipe/anaconda3/envs/GWFAST/lib/python3.10";




Block[
	{
		asd = Import[FileNameJoin[{desktopPrefix, "/site-packages/psds/LVC_O1O2O3/2017-08-06_DCH_C02_L1_O2_Sensitivity_strain_asd.txt"}], "Data"],
		pos1, pos2
	},
	pos1 = FirstPosition[asd[[All,1]], x_/;x>=20.]//Last;
	pos2 = FirstPosition[asd[[All,1]], x_/;x>=1024.]//Last;
	PSD["L1"] = (asd[[pos1;;pos2,2 ]])^2;
]



Block[
	{
		asd = Import[FileNameJoin[{desktopPrefix, "/site-packages/psds/LVC_O1O2O3/2017-06-10_DCH_C02_H1_O2_Sensitivity_strain_asd.txt"}], "Data"], 
		pos1, pos2
	},
	pos1 = FirstPosition[asd[[All,1]], x_/;x>=20.]//Last;
	pos2 = FirstPosition[asd[[All,1]], x_/;x>=1024.]//Last;
	PSD["H1"] = (asd[[pos1;;pos2,2 ]])^2;
]


DetectorTensor["H1"] = With[
    {nx= {-0.2239, 0.7998, 0.5569}, ny = {-0.9140, 0.0261, -0.4049}},
    (nx\[TensorProduct]nx - ny\[TensorProduct]ny)/2
];

DetectorTensor["L1"] = With[
   {nx = {\[Minus]0.9546,\[Minus]0.1416,\[Minus]0.2622}, ny = {+0.2977,\[Minus]0.4879,\[Minus]0.8205} },
   0.5 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];

Vertex["H1"] = {-2.16141492636 10^6,  -3.83469517889 10^6   , 4.60035022664 10^6};
Vertex["L1"] = {-7.42760447238 10^4, -5.49628371971 10^6 ,  3.22425701744 10^6 };


(*Data from https://iopscience.iop.org/article/10.3847/1538-4357/ac4164*)
Module[
	{l},
	l = SymmetrizedIndependentComponents[{3,3}, Symmetric[All]];
	H1data = Join[Vertex["H1"], Extract[DetectorTensor["H1"], l]]; 
	L1data = Join[Vertex["L1"], Extract[DetectorTensor["L1"], l]];
]


detecVars = {\[Theta],\[Phi],\[Psi],\[Iota], p1,p2,p3,D11,D12,D13,D22,D23,D33,f, 4};
hvars = {Mc,\[Eta],\[Chi]1,\[Chi]2,\[Delta]\[CurlyPhi]minus2,\[Delta]\[CurlyPhi]0,\[Delta]\[CurlyPhi]1,\[Delta]\[CurlyPhi]2,\[Delta]\[CurlyPhi]3,\[Delta]\[CurlyPhi]4,\[Delta]\[CurlyPhi]5,\[Delta]\[CurlyPhi]5l,\[Delta]\[CurlyPhi]6,\[Delta]\[CurlyPhi]6l,\[Delta]\[CurlyPhi]7,\[Delta]\[Beta]2,\[Delta]\[Beta]3,\[Delta]\[Alpha]2,\[Delta]\[Alpha]3,\[Delta]\[Alpha]4, fref, f, 4};
auxhvars = {tc, \[Phi]c, dL, gpsTime, f, 3};


ExternalEvaluate["Python","
from astropy.time import Time
# GPS time for GW150914
gps_time = 61094.012# 1187008882.4

# Convert GPS time to UTC and then to GMST in radians
time = Time(gps_time, format='gps', scale='utc')
gmst_radians = time.sidereal_time('mean', 'greenwich').radian

print(\"GMST in radians:\", gmst_radians)
"]


Module[
	{\[Theta] = 0.4080839999999999, \[Phi] = 3.4461599999999994, gmst =6.675374052715275 10^-7, \[Psi]=2.3, \[Iota] = 2.545065595974997},
	
	fpL1 = {\[Theta], \[Phi], \[Psi], \[Iota], Sequence@@L1data};
	fpH1 =   {\[Theta], \[Phi], \[Psi], \[Iota], Sequence@@H1data};
]

fph = {1.3021182,0.24786618323504223, 0.005136138323169717,  0.003235146993487445, Sequence@@ConstantArray[0,16], 20};
fpauxh = {0, 0, 1, 61094.012};


expr = Hold[
	{s1, h, auxh}, 
	2,
	{
		{detecVars, hvars, auxhvars},
		{fpL1,fpH1, fph, fpauxh},
		3
	},
	{20.,1024.,0.125},
	{PSD["L1"], PSD["H1"]},
	SymRules,
	NRules
];


<<FelipeBarbosa`SymDALI`


MemoryConstrained[
EchoTiming[res = iGWDALICoefficients@@expr];,
12 10^9
]


ArrayReshape[res[[3,3]], {286, 286}]//SymmetricMatrixQ


contexts = ToString/@$ContextPath;
contexts = (#<>"*")&/@contexts;
AppendTo[contexts, "FelipeBarbosa`SymDALI`DALICoefficients`Private`*"];


Dynamic[Refresh[
	Length[Names[#]]&/@{
		"FelipeBarbosa`SymDALI`DALICoefficients`Private`*",
		"FelipeBarbosa`SymDALI`DALICoefficients`*",
		"FelipeBarbosa`SymDALI`*",
		"FelipeBarbosa`SymDALI`Private`*",
		"Global`*"
	},
	UpdateInterval->1
]]


Do[
	iGWDALICoefficients@@expr//QuietEcho,
	{i, 20}
]


after = ToExpression/@Names["FelipeBarbosa`SymDALI`DALICoefficients`Private`*"];


Complement[after, before]


Clear[value];

MapThread[
	(value[#1] = #2)&,
	
	{
		Flatten[{\[Theta],\[Phi],\[Psi],\[Iota],Mc,\[Eta],\[Chi]1,\[Chi]2, tc,\[Phi]c,dL}\[TensorProduct]{\[Theta],\[Phi],\[Psi],\[Iota],Mc,\[Eta],\[Chi]1,\[Chi]2, tc,\[Phi]c,dL}],
		res[[1,1]]
	}
];


gwfast = Flatten[{Mc, \[Eta], dL, \[Theta], \[Phi], \[Iota], \[Psi], tc, \[Phi]c, \[Chi]1, \[Chi]2}\[TensorProduct]{Mc,\[Eta], dL, \[Theta], \[Phi], \[Iota], \[Psi], tc, \[Phi]c, \[Chi]1, \[Chi]2}];


fisher = value/@gwfast;


fish  = ArrayReshape[fisher, {11,11}];


ratio = fish/py;


ratio//Round//MatrixForm


(* ::Text:: *)
(*-> the matrix elements (dL, tc) and (dl,\[Phi]c) should be identically zero, in linux neither gwfast nor SymDALI get zero for them, but SymDALI gets smaller numbers that is why these elements are the ones farther from 1. in the ratio.*)


ratio//Round//MatrixForm


(* ::Section::Closed:: *)
(*Testing some Likelihood calculation:*)


<<FelipeBarbosa`SymDALI`


(* ::Subsection:: *)
(*DALI*)


Module[
	{theta = 0.4080839999999999, phi = 3.4461599999999994,psi=2.3, iota = 2.545065595974997},
	fiducial = {
		{\[Theta], theta},
		{\[Phi], phi},
		{\[Psi], psi}, 
		{\[Iota], iota}, 
		{Mc,1.3021182}, 
		{\[Eta], 0.24786618323504223}, 
		{\[Chi]1,0.005136138323169717}, 
		{\[Chi]2,0.003235146993487445},
		{tc,0},
		{\[Phi]c, 0},
		{dL, 1}
	};
	
]


MemoryConstrained[
	poly = EchoTiming@TaylorForm[res, fiducial],
	12 10^9
]


Cpoly =  List[
	fiducial[[All,1]],
	poly,
	CompilationTarget->"C",
	RuntimeOptions->"Speed",
	RuntimeAttributes->{Listable},
	Parallelization->True
];

Cpoly = Compile@@Cpoly//EchoTiming;


<<CompiledFunctionTools`


Cpoly@@fiducial[[All, 2]]//AbsoluteTiming//ScientificForm


fiducial[[All,1]]


Pattern[#, Blank[]]&/@{\[Theta],\[Phi],\[Psi],\[Iota],Mc,\[Eta],\[Chi]1,\[Chi]2,tc,\[Phi]c,dL}


t[\[Theta]_?NumericQ,\[Phi]_,\[Psi]_,\[Iota]_,Mc_,\[Eta]_,\[Chi]1_,\[Chi]2_,tc_,\[Phi]c_,dL_] := Cpoly[\[Theta],\[Phi],\[Psi],\[Iota],Mc,\[Eta],\[Chi]1,\[Chi]2,tc,\[Phi]c,dL]


NMaximize[
	{
		t[\[Theta],\[Phi],\[Psi],\[Iota],Mc,\[Eta],\[Chi]1,\[Chi]2,tc,\[Phi]c,dL],
		0<=\[Theta]<=\[Pi],
		0<=\[Phi]<= 2 \[Pi],
		0<=\[Psi]<=\[Pi],
		0<=\[Iota]<=\[Pi],
		0.1<=Mc<=10,
		0.01<=\[Eta] <=0.25,
		-1 <=\[Chi]1<=1,
		-1<=\[Chi]2<=1,
		-0.1 <=tc<=0.1,
		0<=\[Phi]c<=2 \[Pi], 
		0.05<=dL<= 2
		},
	{\[Theta],\[Phi],\[Psi],\[Iota],Mc,\[Eta],\[Chi]1,\[Chi]2,tc,\[Phi]c,dL},
	Method->"NelderMead"
	
]


list = {
		0<=\[Theta]<=\[Pi],
		0<=\[Phi]<= 2 \[Pi],
		0<=\[Psi]<=\[Pi],
		0<=\[Iota]<=\[Pi],
		0.1<=Mc<=10,
		0.01<=\[Eta] <=0.25,
		-1 <=\[Chi]1<=1,
		-1<=\[Chi]2<=1,
		-0.1 <=tc<=0.1,
		0<=\[Phi]c<=2 \[Pi], 
		0.05<=dL<= 2
}


0<=\[Theta]<=\[Pi]//FullForm


prior = UniformDistribution@(Riffle[list[[All, 1]], list[[All, -1]]]//Partition[#, 2]&);


ResourceFunction["MonteCarloSample"]


chain = ResourceFunction["MonteCarloSample"][
	Cpoly,
	{prior, 100},
	100 10^5, 
	Method -> "EMCEE",
	"LogProbability"->True,
	"ListableFunction"->True	
]//EchoTiming;


python  = StartExternalSession["Python"]


chainR = ArrayReshape[chain, {100, 10^5,11}];
chainR = Transpose[chainR];
Dimensions@chainR


ExternalEvaluate[python, "
import numpy as np
import emcee
chain = np.array(<*chainR*>)
emcee.autocorr.integrated_time(chain)
"]


(* ::Input:: *)
(*RawArray["Real64",{704.3030259615037, 692.8524762880326, 709.3715457798795, 725.7979823311317, 26.485292820552143`, 1267.8179815791796`, 818.6093566364847, 800.6315454805862, 630.6211969533861, 715.4810554775202, 747.0122477537724}]//Normal*)


(* ::Input:: *)
(*{704.3030259615037`,692.8524762880326`,709.3715457798795`,725.7979823311317`,26.485292820552143`,1267.8179815791796`,818.6093566364847`,800.6315454805862`,630.6211969533861`,715.4810554775202`,747.0122477537724`}//Max*)


ExternalEvaluate[python, "
import numpy as np
import emcee
chain = np.array(<*chainR*>)
emcee.autocorr.integrated_time(chain)

"]


ExternalEvaluate[python,"
import numpy as np
def chain_reduce(chain, burn, thin):
	a = chain[burn:, :, :]
	a = a[::thin, :, :]
	nw, ns, d = a.shape
	a = a.reshape(nw*ns, d)
	return a
"]


1267.8179815791796`*2
1267.8179815791796`/2


ExternalEvaluate[python, "
import matplotlib.pyplot as plt
import corner
chain = np.array(<*chainR*>)

flat =  chain_reduce(chain, 4000, 1000)
fig = corner.corner(
    flat,
    
     
    levels = (0.68,0.954),
    
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={\"fontsize\": 12},
);

plt.savefig('test.pdf')
"]


(* ::Subsection::Closed:: *)
(*Full Likelihood*)


phase[Mc_,\[Eta]_,\[Chi]1_,\[Chi]2_, f_] = NRules["\[CapitalPhi]IMR"][[1,2,0]]@@{
	(4.9254664969309`3.6105383994801805*^-6 f Mc)/\[Eta]^(3/5),
	\[Eta],
	\[Chi]1, 
	\[Chi]2, 
	(4.9254664969309`3.6105383994801805*^-6 20 Mc)/\[Eta]^(3/5),
	Sequence@@ConstantArray[0,16]
}


amp[Mc_, \[Eta]_, \[Chi]1_, \[Chi]2_,f_] = NRules["\[ScriptA]IMR"][[1,2,0]]@@{f, Mc/\[Eta]^(3/5), \[Eta], \[Chi]1, \[Chi]2}


Detector[f_,\[Theta]_,\[Phi]_,\[Psi]_,cos\[Iota]_,p1_,p2_,p3_,D11_,D12_,D13_,D22_,D23_,D33_] = NRules["S"][[1,2]]


With[{
	f = Range[20., 1024., 0.125], 
	Mc = 1.30212, \[Eta] = 0.247866, \[Chi]1 = 0.005136138323169717`, \[Chi]2=0.003235146993487445`
	}, 
	
	signal["L1"] = Detector[
		f,
		0.4080839999999999`,
		3.4461599999999994`,
		2.3,
		Cos[2.545065595974997`],
		Sequence@@L1data
	]*auxh[0,0,1,61094.012, f]*Exp[-I phase[Mc,\[Eta], \[Chi]1, \[Chi]2,f]]*amp[Mc,\[Eta], \[Chi]1, \[Chi]2,f];
	
	signal["H1"] = Detector[
		f,
		0.4080839999999999`,
		3.4461599999999994`,2.3,
		Cos[2.545065595974997`],
		Sequence@@H1data
	]*auxh[0,0, 1, 61094.012, f]*Exp[-I phase[Mc,\[Eta], \[Chi]1, \[Chi]2,f]]*amp[Mc,\[Eta], \[Chi]1, \[Chi]2,f];
]


ClearAll[\[ScriptCapitalL]]
Attributes[\[ScriptCapitalL]] = {Listable};

\[ScriptCapitalL][\[Theta]_,\[Phi]_,\[Psi]_,\[Iota]_,Mc_,\[Eta]_,\[Chi]1_,\[Chi]2_,tc_,\[Phi]c_,dL_] := Module[
	{f = Range[20., 1024., 0.125], gpsTime = 61094.012, det1, det2, corewv, wv1, wv2, \[ScriptCapitalL]1, \[ScriptCapitalL]2},
	
	det1 = Detector[f,\[Theta],\[Phi],\[Psi],Cos[\[Iota]], Sequence@@L1data];
	det2 = Detector[f,\[Theta],\[Phi],\[Psi],Cos[\[Iota]], Sequence@@H1data];
	corewv = Exp[-I phase[Mc,\[Eta],\[Chi]1,\[Chi]2,f]]*amp[Mc,\[Eta],\[Chi]1,\[Chi]2,f]*auxh[tc,\[Phi]c,dL,gpsTime,f];
	
	wv1 = det1 corewv;
	wv2 = det2 corewv;
	
	\[ScriptCapitalL]1 = (Conjugate[signal["L1"] - wv1]*(signal["L1"] - wv1))/PSD["L1"]//Total;
	\[ScriptCapitalL]2 = (Conjugate[signal["H1"] - wv2]*(signal["H1"] - wv2))/PSD["H1"]//Total;
	
	-4 0.125 Re[\[ScriptCapitalL]1 + \[ScriptCapitalL]2]
]


prior//RandomVariate


(* ::Input:: *)
(*\[ScriptCapitalL][2.0032079985179942`,4.792234631557569`,0.9239304230491383`,2.5410568923347423`,7.7392922758502`,0.05666512378571305`,-0.8034343344162505`,0.2500306471029532`,-0.09683769310566084`,3.3242712284477736`,1.2739310035255853`]//AbsoluteTiming*)


SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber"->1}]


CloseKernels[]
LaunchKernels[2]


chainFull = ResourceFunction["MonteCarloSample"][
	\[ScriptCapitalL],
	{prior, 100},
	100 10^5, 
	Method -> "EMCEE", 
	"LogProbability"->True,
	"ListableFunction"->True,
	Parallelization->True
]//EchoTiming;


Export["matrix.mat", chainFull]


chainFull = Import["matrix.mat"];


chainFull = chainFull[[1]];
chainFull//Dimensions


chainR = ArrayReshape[chainFull, {100, 10^5,11}];
chainR = Transpose[chainR];
Dimensions@chainR


python = StartExternalSession["Python"]


ExternalEvaluate[python, "
import numpy as np
import emcee
chain = np.array(<*chainR*>)
emcee.autocorr.integrated_time(chain)
"]
