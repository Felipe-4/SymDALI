(* ::Package:: *)

BeginPackage["FelipeBarbosa`SymDALI`Detectors`"]


Fp::usage = "Pattern Function \!\(\*SubscriptBox[\(F\), \(+\)]\)";
Fx::usage= "Pattern Function \!\(\*SubscriptBox[\(F\), \(x\)]\)";


Begin["`Private`"];


LSO::usage="LSO[M] gives the GW frequency at the last stable circular orbit for
a test mass in a Schwarzschild spacetime of mass M
M[SolarMass]";


LSO[M_] = With[{G = (Quantity[("GravitationalConstant")/("SpeedOfLight")^3]//UnitConvert[#, ("Seconds")/("SolarMass")]&)[[1]]}, (6^(3/2) M G \[Pi])^-1];


DetectorTensor::usage = "Normal vectors taken from https://arxiv.org/pdf/gr-qc/0008066 which seem to agree with the ones in
https://www.ligo.org/scientists/GW100916/detectors.txt";


DetectorTensor["H1"] = With[
    {nx= {-0.2239, 0.7998, 0.5569}, ny = {-0.9140, 0.0261, -0.4049}},
    (nx\[TensorProduct]nx - ny\[TensorProduct]ny)/2
];

DetectorTensor["L1"] = With[
   {nx = {\[Minus]0.9546,\[Minus]0.1416,\[Minus]0.2622}, ny = {+0.2977,\[Minus]0.4879,\[Minus]0.8205} },
   0.5 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];

DetectorTensor["V1"] = With[
    {nx = {\[Minus]0.7005,+0.2085,+0.6826}, ny = {\[Minus]0.0538,\[Minus]0.9691,+0.2408} },
	0.5 (nx\[TensorProduct]nx - ny\[TensorProduct]ny)
];


(* ::Section:: *)
(*Antenna Pattern Functions*)


f1[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_, D33_] = Module[
	{expr},
	expr = (
		(D22 - D11 Cos[\[Theta]]^2) Cos[\[Phi]]^2
		- D33 Sin[\[Theta]]^2
		+ Sin[2 \[Theta]] (D13 Cos[\[Phi]] + D23 Sin[\[Phi]])
		+ (D11 - D22 Cos[\[Theta]]^2) Sin[\[Phi]]^2
		- D12 Sin[2 \[Phi]] (1+Cos[\[Theta]]^2)
	);
	
	expr = (expr)//.{Cos[d_]^2 :> (1+ Cos[2 d])/2, Sin[d_]^2 :> (1- Cos[2 d])/2 };
	
	
	expr//Simplify
] 


f2[\[Theta]_, \[Phi]_, D11_, D12_, D13_, D22_, D23_] = (
	2 Sin[\[Theta]] (D13 Sin[\[Phi]] - D23 Cos[\[Phi]]) +
	Cos[\[Theta]] (2 D12 Cos[2 \[Phi]] + (D22-D11) Sin[2 \[Phi]])
)//Simplify


Fp[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	Cos[2 \[Psi]] f1[\[Theta], \[Phi], D11, D12, D13, D22, D23, D33] + 
	Sin[2 \[Psi]] f2[\[Theta], \[Phi], D11, D12, D13, D22, D23]
);

Fx[\[Theta]_, \[Phi]_,\[Psi]_, D11_, D12_, D13_, D22_, D23_, D33_] = (
	- Sin[2 \[Psi]] f1[\[Theta], \[Phi], D11, D12, D13, D22, D23, D33] +
	Cos[2 \[Psi]] f2[\[Theta], \[Phi], D11, D12, D13, D22, D23]
);


End[];


EndPackage[]
