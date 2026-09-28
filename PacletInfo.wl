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
            {"FelipeBarbosa`SymDALI`DerivativeTools`", "DerivativeTools.wl"}, 
            {"FelipeBarbosa`SymDALI`DALIPolynomial`", "DALIPolynomial.wl"},
            {"FelipeBarbosa`SymDALI`Population`", "Population.wl"},
            {"FelipeBarbosa`SymDALI`Utils`", "Utils.wl"}
        }
      },
      {
          "Asset",
          "Root"->"Assets", 
          "Assets"->{
              {"FpFc", "Derivatives/Detectors/Defs.mx"},
              {"PhenomD", "Derivatives/IMRPhenomD/Defs.mx"},
              {"PhenomHM", "Derivatives/IMRPhenomHM/Defs.mx"},
              {"ET_D",  "ASD/ET-0000A-18_ETDSensitivityCurveTxtFile.txt"},
              {"CE-20", "ASD/cosmic_explorer_20km_strain.txt"},
              {"CE-40-lf", "ASD/cosmic_explorer_40km_lf_strain.txt"},
              {"CE-20-pm", "ASD/cosmic_explorer_20km_pm_strain.txt"},
              {"CE-40", "ASD/cosmic_explorer_strain.txt"},
              {"ET-10", "ASD/18213_ET10kmcolumns.txt"},
              {"ET-15", "ASD/18213_ET15kmcolumns.txt"},
              {"ET-20", "ASD/18213_ET20kmcolumns.txt"},
              {"Aplus", "ASD/AplusDesign.txt"},
              {"kagra", "ASD/kagra_80Mpc.txt"},
              {"V_o5", "ASD/avirgo_O5low_NEW.txt"}
          }
      }
    }
  |>
]
