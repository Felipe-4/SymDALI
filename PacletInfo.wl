(* ::Package:: *)

PacletObject[
  <|
    "Name" -> "FelipeBarbosa/SymDALI",
    "Description" -> "Implementation of the DALI algorithm (Derivative Approximation for LIkelihoods)",
	"Creator" -> "Felipe Barbosa",
    "Version" -> "1.0.0",
    "WolframVersion" -> "14.1+",
    "PublisherID" -> "FelipeBarbosa",
    "License" -> "MIT",
    "PrimaryContext" -> "FelipeBarbosa`SymDALI`",
    "DocumentationURL" -> "https://resources.wolframcloud.com/PacletRepository/resources",
    "Extensions" -> {
      {
        "Kernel",
        "Root" -> "Kernel",
        "Context" -> {
            {"FelipeBarbosa`SymDALI`","SymDALI.wl"},
            {"FelipeBarbosa`SymDALI`DALICoefficients`", "DALICoefficients.wl"}, 
            {"FelipeBarbosa`SymDALI`Detectors`", "Detectors.wl"}, 
            {"FelipeBarbosa`SymDALI`IMRPhenomD`", "IMRPhenomD.wl"}, 
            {"PNCoefficients`", "PNCoefficients.wl"}, 
            {"FelipeBarbosa`SymDALI`DerivativeTools`", "DerivativeTools.wl"}, 
            {"FelipeBarbosa`SymDALI`DALIPolynomial`", "DALIPolynomial.wl"},
            {"FelipeBarbosa`SymDALI`PNCoefficients`","PNCoefficients.wl" }
        }
      },
      {
        "Documentation",
        "Root" -> "Documentation",
        "Language" -> "English"
      },
      {"LibraryLink"}
    }
  |>
]
