# Contributing to SymDALI

This document is intended for collaborators who want to extend or maintain the codebase. It assumes familiarity with Wolfram Language and a basic understanding of GW parameter estimation.

## Repository Layout

```
SymDALI/
├── Kernel/           — Wolfram Language source packages
├── Rules/            — Pre-computed symbolic and numerical derivative rules
├── LibraryResources/ — Platform-specific compiled C libraries
├── Data/ASD/         — Noise curve data files
├── Documentation/    — Reference pages (Mathematica notebook format)
└── Demo/             — Usage examples
```

## Core Architecture

Understanding the data flow is essential before contributing.

### Gradient computation pipeline

The central computation is the DALI tensor network sum over frequency bins:

$$\mathcal{F}^{(i,j)}_{\mu_1 \cdots \mu_i \nu_1 \cdots \nu_j} = 4 \Delta f \sum_k \frac{\text{Re}\left[\partial_{\mu_1 \cdots \mu_i} \tilde{s}^*(f_k)\, \partial_{\nu_1 \cdots \nu_j} \tilde{s}(f_k)\right]}{S_n(f_k)}$$

where $\tilde{s}(f) = F_+ h_+(f) + F_\times h_\times(f)$ is the detector signal. The computation is factored:

1. **Derivative rules loading** (`DerivativeTools.wl`): Pre-computed rules for amplitude $\mathcal{A}$, phase $\Phi$, and pattern functions $F_{+,\times}$ are loaded from disk as compiled C `LibraryFunction` objects.

2. **Gradient computation** (`DALICoefficients.wl`): `iGenGrads` iteratively applies `TakeGrad` to build gradients up to the requested order, exploiting Mathematica's `SymmetrizedIndependentComponents` to work with only the linearly independent tensor components.

3. **Tensor contraction** (`DALICoefficients.wl`): `GenDaliTerm` and `GenDaliList` perform the frequency-domain integral via the compiled `CiGenDaliTerm` function.

4. **Post-processing** (`DALIPolynomial.wl`): `ProcessDALITensors` and `TaylorForm` convert the raw output into forms convenient for sampling.

### Derivative rules format

The `Rules/` directory stores pre-computed derivatives in two forms:

- **`SymRules`**: Associations mapping `$D[{n,...}, f] -> number` for derivatives that evaluate to constants (e.g., phase derivatives that do not depend on frequency). Stored as `.wdx` or `.mx` files.
- **`NRules`**: Associations mapping `$D[{n,...}, f][args_] -> CompiledFunction[...]` for derivatives that are functions of the arguments. The `CompiledFunction` wraps a `LibraryFunction` loaded from a `.dylib` file.

The helper `$D[{n1,n2,...}, f][x,y,...]` denotes the mixed partial derivative $\partial^{n_1+n_2+\cdots} f / \partial x_1^{n_1} \partial x_2^{n_2} \cdots$ evaluated at $(x,y,\ldots)$.

## How to Add a New Detector

### 1. Add detector tensor and vertex

In `Kernel/Detectors.wl`, add entries to the `DetectorTensor` and `DetectorVertex` associations. Use `ArmDirection` and `FromGeocentricCoordinates` to compute the tensors from geographic coordinates:

```mathematica
Module[
  {lat = ... Degree, lon = ... Degree, n1, n2},
  {n1, n2} = ArmDirection[{lat, lon}, {{omega1, psi1}, {omega2, psi2}}];
  DetectorTensor["MyDet"] = (n1⊗n1 - n2⊗n2)/2;
]

DetectorVertex["MyDet"] = FromGeocentricCoordinates[{h, lat, lon}];
```

The arm direction parameters `{omega_i, psi_i}` follow the WGS-84 convention from `LALDetectors.h`: `omega` is the altitude and `psi` is the azimuth of the arm.

### 2. Add a noise curve

Place a two-column text file (frequency [Hz], ASD [1/√Hz]) in `Data/ASD/`, then add an interpolation in `Kernel/Detectors.wl`:

```mathematica
data = Import[FileNameJoin[{ASDdir, "my_detector_asd.txt"}], "Data"] // N;
ASD["MyDet"] = Interpolation[data];
```

No changes to any other file are needed.

## How to Add a New Waveform Model

Adding a new waveform is the most involved type of contribution. The following steps are required.

### 1. Implement the waveform function

The waveform must be expressed as a product of an amplitude factor and a phase factor, separated so that the amplitude can be computed independently. Study `ihphcIMRPhenomD` in `Kernel/DALICoefficients.wl` for the expected structure.

The function should take all physical parameters as separate scalar arguments (not a list or association), followed by a reference frequency and the evaluation frequency as the last two arguments.

### 2. Generate derivative rules

The most expensive step. Use `DerivativeRules` from `Kernel/DerivativeTools.wl` to compute all symbolic derivatives and compile them to C:

```mathematica
(* Example: generate derivatives of the amplitude function *)
ampRules = DerivativeRules[
  {myAmp, {param1, param2, ...}, 3},  (* name, vars, max order *)
  $Block[
    {{internalHead1, internalHead2}},
    internalHead1[x_] = ...; (* internal function definitions *)
    myAmpExpression[param1, param2, ...]
  ],
  "KeepDefs" -> {localVar -> expression, ...}
]
```

For waveforms with many parameters, this can take hours. The demo notebooks in `Demo/` and `Documentation/English/Tutorials/` contain worked examples for IMRPhenomD.

### 3. Save the derivative rules

Use `SaveSymRules` and `SaveNRules` (or the serialization utilities in `DerivativeTools.wl`) to persist the rules to `Rules/YourModel/` and to `LibraryResources/<platform>/DerivativeRules/YourModel/`.

### 4. Register the model in `DALICoefficients.wl`

Add entries in the `Module` block that loads derivative rules at package initialization, and implement the corresponding `ihphcYourModel` function. Register the model name in `headhphc`, `iihphcvars`, `isAligned`, and `iconvert` lookup tables. Add validation in `RetrieveFiducial` and `iDALITensors`.

### 5. Add a public waveform evaluation function in `Utils.wl`

Following the pattern of `hphcIMRPhenomD`, implement a user-facing function `hphcYourModel` with input validation and `Options` for any deviation parameters.

## How to Generate New Derivative Rules for an Existing Model

If you want to extend an existing model to a higher DALI order or to additional parameters, you need to generate the additional derivative rules and add them to the existing rule sets.

The notebook `Demo/DerivativeRules_test.nb` and `Demo/manual_C_functions.nb` show examples of how this was done for IMRPhenomD. The key function is:

```mathematica
DerivativeRules[{name, vars, n}, expr, opts]
```

Important options:
- `"IncludeZeroDerivative" -> True` — include the zeroth derivative (function value itself)
- `"KeepDefs" -> {var -> expr, ...}` — local variable substitutions to apply after differentiation
- `"Parallel" -> {k}` — use `k` parallel kernels

## Platform Support

All compiled C libraries are currently provided only for `MacOSX-ARM64`. The `DerivativeRulesLoad` function in `DerivativeTools.wl` uses `SystemInformation["Kernel", "SystemID"]` to locate the appropriate library directory. To add support for a new platform:

1. Recompile the `.dylib` files on the target platform using the C source in `LibraryResources/MacOSX-ARM64/Polynomial_Functions/DALI.c` and by regenerating the waveform derivative rules.
2. Place the resulting libraries under `LibraryResources/<SystemID>/`.

The `SystemID` string for common platforms is:
- `"MacOSX-ARM64"` — macOS Apple Silicon
- `"MacOSX-x86-64"` — macOS Intel
- `"Linux-x86-64"` — Linux 64-bit

## Code Conventions

- **Wolfram Language style**: Follow the existing patterns. Internal (non-exported) functions are prefixed with a lowercase `i` (e.g., `iGenGrads`, `iDALITensors`).
- **Error handling**: Use `Throw[$Failed, failTag[functionName]]` for internal errors and `Throw["message"]` for user-facing errors. The top-level `DALITensors` wraps everything in `Catch`.
- **Protection**: All public symbols are `Protect`ed after definition. Use `Unprotect` before modifying, then `Protect` again.
- **Package contexts**: Each submodule lives in its own context under `FelipeBarbosa`SymDALI``. Keep package dependencies explicit via the `BeginPackage` second argument.
- **No comments on the obvious**: Comment only non-obvious derivations, physical assumptions, or known bugs/workarounds.

## Testing

There is no automated test suite. Validation is done by cross-checking DALI tensor components against independent numerical differentiation and against results from [gwfast](https://arxiv.org/abs/2203.02670) for the Fisher matrix limit. Before submitting changes:

1. Run the `Demo/UserGuide.wl` examples end to end.
2. For Fisher matrix changes, verify that `DALITensors[..., 1][[1,1]]` agrees with an independent Fisher matrix calculation (e.g., numerical finite differences of the log-likelihood).
3. For waveform changes, compare `hphcYourModel` output against a reference implementation over a broad parameter range.

## Opening Issues and Pull Requests

This repository is currently in a pre-publication state. Please open an issue to discuss any proposed changes before starting significant development work, to avoid duplicating ongoing efforts.
