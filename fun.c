#include "math.h"

#include "WolframRTL.h"

static WolframCompileLibrary_Functions funStructCompile;

static const mint UnitIncrements[3] = {1, 1, 1};

static UnaryMathFunctionPointer FP0;

static BinaryMathFunctionPointer FP1;

static UnaryMathFunctionPointer FP2;

static LibraryFunctionPointer FP3;

static MArgument FPA[4];


static mint I0_0;

static mint I0_1;

static mbool initialize = 1;

#include "fun.h"

DLLEXPORT int Initialize_fun(WolframLibraryData libData)
{
if( initialize)
{
funStructCompile = libData->compileLibraryFunctions;
I0_0 = (mint) 2;
I0_1 = (mint) 6;
FP0 = funStructCompile->getUnaryMathFunction(1, 3);/*  Sin  */
if( FP0 == 0)
{
return LIBRARY_FUNCTION_ERROR;
}
FP1 = funStructCompile->getBinaryMathFunction(257, 3, 3);/*  Plus  */
if( FP1 == 0)
{
return LIBRARY_FUNCTION_ERROR;
}
FP2 = funStructCompile->getUnaryMathFunction(2, 3);/*  Cos  */
if( FP2 == 0)
{
return LIBRARY_FUNCTION_ERROR;
}
FP3 = funStructCompile->getFunctionCallPointer("DotVV");
if( FP3 == 0)
{
return LIBRARY_FUNCTION_ERROR;
}
initialize = 0;
}
return 0;
}

DLLEXPORT void Uninitialize_fun(WolframLibraryData libData)
{
if( !initialize)
{
initialize = 1;
}
}

DLLEXPORT int fun(WolframLibraryData libData, mreal A1, MTensor A2, mreal *Res)
{
mreal R0_0;
mreal R0_1;
mreal R0_2;
MTensor* T0_0;
MTensor* T0_1;
MTensorInitializationData Tinit;
mreal *P1;
MArgument FPA[4];
int err = 0;
Tinit = funStructCompile->GetInitializedMTensors(libData, 1);
T0_1 = MTensorInitializationData_getTensor(Tinit, 0);
R0_0 = A1;
T0_0 = &A2;
{
mint S0 = FP0((void*) (&R0_1), (void*) (&R0_0), 1, UnitIncrements, 6);/*  Sin  */
err = S0 == 0 ? 0 : LIBRARY_NUMERICAL_ERROR;
if( err)
{
goto error_label;
}
}
R0_2 = (mreal) I0_0;
{
mint S0 = FP1((void*) (&R0_1), (void*) (&R0_1), (void*) (&R0_2), 1, UnitIncrements, 6);/*  Plus  */
err = S0 == 0 ? 0 : LIBRARY_NUMERICAL_ERROR;
if( err)
{
goto error_label;
}
}
{
mint S0 = FP2((void*) (&R0_2), (void*) (&R0_0), 1, UnitIncrements, 6);/*  Cos  */
err = S0 == 0 ? 0 : LIBRARY_NUMERICAL_ERROR;
if( err)
{
goto error_label;
}
}
{
mint S0[1] = {2};
err = funStructCompile->MTensor_allocate(T0_1, 3, 1, S0);
if( err)
{
goto error_label;
}
P1 = MTensor_getRealDataMacro(*T0_1);
P1[0] = R0_2;
P1[1] = R0_1;
}
MArgument_getMTensorAddress(FPA[0]) = T0_1;
MArgument_getMTensorAddress(FPA[1]) = T0_0;
MArgument_getIntegerAddress(FPA[2]) = &I0_1;
MArgument_getRealAddress(FPA[3]) = &R0_2;
err = FP3(libData, 3, FPA, FPA[3]);/*  DotVV  */
if( err)
{
goto error_label;
}
*Res = R0_2;
error_label:
funStructCompile->ReleaseInitializedMTensors(Tinit);
funStructCompile->WolframLibraryData_cleanUp(libData, 1);
return err;
}

