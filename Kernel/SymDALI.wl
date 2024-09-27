(* ::Package:: *)

BeginPackage["FelipeBarbosa`SymDALI`"];
ClearAll[us\[Theta]];
us\[Theta];
Begin["`Private`"]
Attributes@us\[Theta] = {Listable};
us\[Theta][x_]/; x<0 = 0; us\[Theta][x_]/; x==0 = 1/2; us\[Theta][x_]/; x>0 = 1;
Derivative[n_][us\[Theta]][x_]/;n>=1 := 0
End[];


EndPackage[];
<<FelipeBarbosa`SymDALI`DALICoefficients`;
<<FelipeBarbosa`SymDALI`DerivativeTools`;
<<FelipeBarbosa`SymDALI`Detectors`;
<<FelipeBarbosa`SymDALI`IMRPhenomD`;
<<FelipeBarbosa`SymDALI`PNCoefficients`;
<<FelipeBarbosa`SymDALI`DALIPolynomial`;
