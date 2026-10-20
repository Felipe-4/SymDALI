(* ::Package:: *)

PacletDirectoryLoad["/Users/felipe/Documents/GitHub"];


<<FelipeBarbosa`SymDALI`


fp = <|
    "dec" -> -0.883787, "\[Phi]" -> 6.18412, "\[Psi]" -> 0.471846,
    "\[ScriptCapitalM]c" -> 190.258, "q" -> 0.935815,
    "s1z" -> 0.107562, "s2z" -> -0.0254365, "\[Iota]" -> 2.60021,
    "1/dL" -> 0.0284797, "tc" -> 61093.87906420612, "\[Phi]ref" -> 5.93578
|>;

ET1 = <|
        "Position"-> DetectorVertex["ET-S"],
        "DetectorTensor"->DetectorTensor["ET-S-L"],
        "ASD"->ASD["ET-15-lfhf"]
|>;

ET2 = <|
        "Position"->DetectorVertex["ET-MR"],
        "DetectorTensor"->DetectorTensor["ET-MR-L-45"],
        "ASD"->ASD["ET-15-lfhf"]
|>;

CE = <|
        "Position"-> DetectorVertex["CE-NM"],
        "DetectorTensor"->DetectorTensor["CE-NM"],
        "ASD"->ASD["CE-40"]
|>;

(*Given the total mass of the injection after 93 Hz the waveform is 0*)
fcut = 93.;

(*We take Delta_f to be 1/Delta_t == 1/2 for this signal*)
res = Length@Range[5., 93., 1/2.] + 1;

tensors = DALITensors["IMRPhenomHM", fp, {ET1, ET2, CE}, 3, "fmin"->5, "fmax"->93, "res"->res];


taylor = TaylorForm[tensors];


SetDirectory[DirectoryName[$InputFileName]]


Export["taylor.txt", Flatten[taylor], "Table"];
