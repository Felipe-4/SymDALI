(* ::Package:: *)

System`$D::usage="$D[{n__}, symbol][y__] is equivalent to Derivative[n__][symbol][y__].";
System`TagRule::usage="TagRule[f, g[f[x_]], 3] is equivalent in spirit to f/: g[f[x_]] = 3. At some point I will fix the formating so that you 
can type f/: g[f[x_]] -> 3. For now the syntax is TagRule[tag, lhs, rhs]";

BeginPackage["FelipeBarbosa`SymDALI`"];
ClearAll[us\[Theta]];
us\[Theta];
Begin["`Private`"]
Attributes@us\[Theta] = {Listable};
us\[Theta][x_]/; x<0 = 0; us\[Theta][x_]/; x==0 = 1/2; us\[Theta][x_]/; x>0 = 1;
Unprotect[Derivative, $D]; Clear[Derivative, $D];

Derivative[n_][us\[Theta]][x_]/;n>=1 := 0; $D[n_, us\[Theta]][x_] := 0

Protect[Derivative, $D];
End[];


EndPackage[];


<<FelipeBarbosa`SymDALI`DALICoefficients`;
<<FelipeBarbosa`SymDALI`DerivativeTools`;
<<FelipeBarbosa`SymDALI`Detectors`;
<<FelipeBarbosa`SymDALI`Population`;
<<FelipeBarbosa`SymDALI`DALIPolynomial`;
<<FelipeBarbosa`SymDALI`Utils`;
