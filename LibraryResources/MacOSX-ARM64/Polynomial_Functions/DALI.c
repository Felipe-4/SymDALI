/* Include required header */
#include "WolframLibrary.h"

typedef struct {
    mreal *data;
    mint   length;
} Piece;

/* Return the version of Library Link */
DLLEXPORT mint WolframLibrary_getVersion() {
    return WolframLibraryVersion;
}

/* Initialize Library */
DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData) {
    return LIBRARY_NO_ERROR;
}

/* Uninitialize Library */
DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData) {
    return;
}

int TensorProduct(WolframLibraryData libData, MTensor V1, MTensor V2, MTensor *result) {
    mint l1,l2;
    mreal *V1Pointer, *V2Pointer, *resultPointer;

    V1Pointer = libData->MTensor_getRealData(V1);
    V2Pointer = libData->MTensor_getRealData(V2);

    l1 = libData->MTensor_getFlattenedLength(V1);
    l2 = libData->MTensor_getFlattenedLength(V2);
    
    if (l1==l2) {
        mint dims[1] = {l1*(l1+1)/2};
        int err = libData->MTensor_new(MType_Real,1, dims, result);
        if (err) return err;

        resultPointer = libData->MTensor_getRealData(*result);
        
        mint i,j;
        mint idx=0;
        for(i=0;i<l1;i++) {
            #pragma omp simd
            for (j=i;j<l1;j++) {
                resultPointer[idx++] = V1Pointer[i]*V2Pointer[j];
            }
        }


    } else {
         mint dims[1] = {l1*l2};
        int err = libData->MTensor_new(MType_Real,1,dims, result);
        if (err) return err;

        resultPointer = libData->MTensor_getRealData(*result);

        mint i, j;
        mint idx=0;
        for (i = 0; i < l1; i++) {
            #pragma omp simd
            for (j = 0; j < l2; j++) {
                resultPointer[idx++] = V1Pointer[i] * V2Pointer[j];
            }
        }

    }
    
    return LIBRARY_NO_ERROR;
}

int Deltap2(WolframLibraryData libData, MTensor tensor, MTensor *result) {
    mint len;
    mreal *tensorPointer, *resultPointer;

    len = libData->MTensor_getFlattenedLength(tensor);
    mint dims[1] = {len*(len+1)/2};
    
    //create result
    int err = libData->MTensor_new(MType_Real, 1, dims, result);
    if (err) return err;

    tensorPointer = libData->MTensor_getRealData(tensor);
    resultPointer = libData->MTensor_getRealData(*result);
    
    mint i,j;
    mint idx=0;
    for (i=0; i<len; i++) {
        #pragma omp simd
        for (j=i; j<len; j++) {
           resultPointer[idx++] =  tensorPointer[i]*tensorPointer[j];
        }
    }
    return LIBRARY_NO_ERROR;
}


int Deltap3(WolframLibraryData libData, MTensor tensor, MTensor *result) {
    mint len;
    mreal *tensorPointer, *resultPointer;

    len = libData->MTensor_getFlattenedLength(tensor);
    mint dims[1] = {len*(len+1)*(len+2)/6};
    
    //create result
    int err = libData->MTensor_new(MType_Real, 1, dims, result);
    if (err) return err;

    tensorPointer = libData->MTensor_getRealData(tensor);
    resultPointer = libData->MTensor_getRealData(*result);
    
    mint i,j,k;
    mint idx=0;
    for (i=0; i<len; i++) {
        for (j=i; j<len; j++) {
            #pragma omp simd
            for (k=j;k<len;k++) {
                resultPointer[idx++] = tensorPointer[i]*tensorPointer[j]*tensorPointer[k];
            }
        }
    }
    return LIBRARY_NO_ERROR;
}

int Deltap4(WolframLibraryData libData, MTensor tensor, MTensor *result) {
    mint len;
    mreal *tensorPointer, *resultPointer;

    len = libData->MTensor_getFlattenedLength(tensor);
    mint dims[1] = {len*(len+1)*(len+2)*(len+3)/24};
    
    //create result
    int err = libData->MTensor_new(MType_Real, 1, dims, result);
    if (err) return err;

    tensorPointer = libData->MTensor_getRealData(tensor);
    resultPointer = libData->MTensor_getRealData(*result);
    
    mint i,j,k,p;
    mint idx=0;
    for (i=0; i<len; i++) {
        for (j=i; j<len; j++) {
            for (k=j;k<len;k++) {
                #pragma omp simd
                for (p=k;p<len;p++) {
                    resultPointer[idx++] = tensorPointer[i]*tensorPointer[j]*tensorPointer[k]*tensorPointer[p];
                }
                
            }
        }
    }
    return LIBRARY_NO_ERROR;
}


int Deltap5(WolframLibraryData libData, MTensor tensor, MTensor *result) {
    mint len;
    mreal *tensorPointer, *resultPointer;

    len = libData->MTensor_getFlattenedLength(tensor);
    mint dims[1] = {len*(len+1)*(len+2)*(len+3)*(len+4)/120};
    
    //create result
    int err = libData->MTensor_new(MType_Real, 1, dims, result);
    if (err) return err;

    tensorPointer = libData->MTensor_getRealData(tensor);
    resultPointer = libData->MTensor_getRealData(*result);
    
    mint i,j,k,p,q;
    mint idx=0;
    for (i=0; i<len; i++) {
        for (j=i; j<len; j++) {
            for (k=j;k<len;k++) {
                for (p=k;p<len;p++) {
                    #pragma omp simd
                    for (q=p;q<len;q++) {
                        resultPointer[idx++] = tensorPointer[i]*tensorPointer[j]*tensorPointer[k]*tensorPointer[p]*tensorPointer[q];
                    }
                }
                
            }
        }
    }
    return LIBRARY_NO_ERROR;
}



DLLEXPORT int DALI_Delta_Ps(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) {
    MTensor DeltaP;
    mint Order;

    DeltaP = MArgument_getMTensor(Args[0]);
    Order = MArgument_getInteger(Args[1]);
    int err;
    
    if (Order==1) {
        // declare t11
        MTensor t11;
        err = Deltap2(libData, DeltaP, &t11);
        if (err) return err;
        MArgument_setMTensor(Res, t11);

    } else if (Order==2) {
        // declare t12 and t22 and t11
        MTensor t11, t12, t22;
        //calculate t11
        err = Deltap2(libData, DeltaP, &t11); if (err) return err;
        //calculate t12
        err = TensorProduct(libData, DeltaP, t11, &t12); if (err) return err;
        //calculate t22
        err = TensorProduct(libData, t11, t11, &t22); if (err) return err;
        
        // 2) Build an array of (pointer, length) pairs
        Piece pieces[3];
        pieces[0].data   = libData->MTensor_getRealData(t11);
        pieces[0].length = libData->MTensor_getFlattenedLength(t11);
        pieces[1].data   = libData->MTensor_getRealData(t12);
        pieces[1].length = libData->MTensor_getFlattenedLength(t12);
        pieces[2].data   = libData->MTensor_getRealData(t22);
        pieces[2].length = libData->MTensor_getFlattenedLength(t22);

        // 3) Allocate the final result
        mint totalLen = pieces[0].length + pieces[1].length + pieces[2].length;
        mint dims[1] = { totalLen };
        MTensor result;
        err = libData->MTensor_new(MType_Real, 1, dims, &result);
        if (err) {
            libData->MTensor_free(t11);
            libData->MTensor_free(t12);
            libData->MTensor_free(t22);
            return err;
        }

        mreal *resPtr = libData->MTensor_getRealData(result);

        // 4) One generic copy loop
        mint offset = 0;
        for (int k = 0; k < 3; ++k) {
            mreal *src = pieces[k].data;
            mint   len = pieces[k].length;
            #pragma omp simd
            for (mint i = 0; i < len; ++i) {
                resPtr[offset + i] = src[i];
            }
            offset += len;
        }

        // 5) Clean up temporaries
        libData->MTensor_free(t11);
        libData->MTensor_free(t12);
        libData->MTensor_free(t22);

        MArgument_setMTensor(Res, result);
        return LIBRARY_NO_ERROR;
    } else if (Order==3) {
         // declare t11, t12, t22, t13, t23, t33
        MTensor t11, t12, t22, t13, t23, t33, P3;

        //calculate p3
        err = Deltap3(libData, DeltaP, &P3); if (err) return err;
        
        //calculate t11
        err = Deltap2(libData, DeltaP, &t11); if (err) return err;
        
        //calculate t12
        err = TensorProduct(libData, DeltaP, t11, &t12); if (err) return err;
        
        //calculate t22
        err = TensorProduct(libData, t11, t11, &t22); if (err) return err;

        //calculate t13
        err = TensorProduct(libData, DeltaP, P3, &t13); if (err) return err;

        //calculate t23
        err = TensorProduct(libData, t11, P3, &t23); if (err) return err;

        //calculate t33
        err = TensorProduct(libData, P3, P3, &t33); if (err) return err;
        
        // 2) Build an array of (pointer, length) pairs
        Piece pieces[6];
        pieces[0].data   = libData->MTensor_getRealData(t11);
        pieces[0].length = libData->MTensor_getFlattenedLength(t11);
        
        pieces[1].data   = libData->MTensor_getRealData(t12);
        pieces[1].length = libData->MTensor_getFlattenedLength(t12);
        
        pieces[2].data   = libData->MTensor_getRealData(t22);
        pieces[2].length = libData->MTensor_getFlattenedLength(t22);
        
        pieces[3].data   = libData->MTensor_getRealData(t13);
        pieces[3].length = libData->MTensor_getFlattenedLength(t13);

        pieces[4].data   = libData->MTensor_getRealData(t23);
        pieces[4].length = libData->MTensor_getFlattenedLength(t23);

        pieces[5].data   = libData->MTensor_getRealData(t33);
        pieces[5].length = libData->MTensor_getFlattenedLength(t33);

        // 3) Allocate the final result
        mint totalLen = pieces[0].length + pieces[1].length + pieces[2].length + pieces[3].length +pieces[4].length +pieces[5].length;
        mint dims[1] = { totalLen };
        MTensor result;
        err = libData->MTensor_new(MType_Real, 1, dims, &result);
        
        if (err) {
            libData->MTensor_free(t11);
            libData->MTensor_free(t12);
            libData->MTensor_free(t22);
            libData->MTensor_free(t13);
            libData->MTensor_free(t23);
            libData->MTensor_free(t33);
            return err;
        }

        mreal *resPtr = libData->MTensor_getRealData(result);

        // 4) One generic copy loop
        mint offset = 0;
        for (int k = 0; k < 6; ++k) {
            mreal *src = pieces[k].data;
            mint   len = pieces[k].length;
            #pragma omp simd
            for (mint i = 0; i < len; ++i) {
                resPtr[offset + i] = src[i];
            }
            offset += len;
        }

        // 5) Clean up temporaries
        libData->MTensor_free(P3);
        libData->MTensor_free(t11);
        libData->MTensor_free(t12);
        libData->MTensor_free(t22);
        libData->MTensor_free(t13);
        libData->MTensor_free(t23);
        libData->MTensor_free(t33);

        MArgument_setMTensor(Res, result);
        return LIBRARY_NO_ERROR;

    }  else if (Order==4) {
         // declare t11, t12, ..., t44
        MTensor t11, t12, t22, t13, t23, t33,t14, t24, t34, t44, P3, P4;

        //calculate p3
        err = Deltap3(libData, DeltaP, &P3); if (err) return err;
        
        //calculate p4
        err = Deltap4(libData, DeltaP, &P4); if (err) return err;
        
        //calculate t11
        err = Deltap2(libData, DeltaP, &t11); if (err) return err;
        
        //calculate t12
        err = TensorProduct(libData, DeltaP, t11, &t12); if (err) return err;
        
        //calculate t22
        err = TensorProduct(libData, t11, t11, &t22); if (err) return err;

        //calculate t13
        err = TensorProduct(libData, DeltaP, P3, &t13); if (err) return err;

        //calculate t23
        err = TensorProduct(libData, t11, P3, &t23); if (err) return err;

        //calculate t33
        err = TensorProduct(libData, P3, P3, &t33); if (err) return err;

        //calculate t14
        err = TensorProduct(libData, DeltaP, P4, &t14); if (err) return err;
        
        //calculate t24
        err = TensorProduct(libData, t11, P4, &t24); if (err) return err;

        //calculate t34
        err = TensorProduct(libData, P3, P4, &t34); if (err) return err;

        //calculate t44
        err = TensorProduct(libData, P4, P4, &t44); if (err) return err;
        
        // 2) Build an array of (pointer, length) pairs
        Piece pieces[10];
        pieces[0].data   = libData->MTensor_getRealData(t11);
        pieces[0].length = libData->MTensor_getFlattenedLength(t11);
        
        pieces[1].data   = libData->MTensor_getRealData(t12);
        pieces[1].length = libData->MTensor_getFlattenedLength(t12);
        
        pieces[2].data   = libData->MTensor_getRealData(t22);
        pieces[2].length = libData->MTensor_getFlattenedLength(t22);
        
        pieces[3].data   = libData->MTensor_getRealData(t13);
        pieces[3].length = libData->MTensor_getFlattenedLength(t13);

        pieces[4].data   = libData->MTensor_getRealData(t23);
        pieces[4].length = libData->MTensor_getFlattenedLength(t23);

        pieces[5].data   = libData->MTensor_getRealData(t33);
        pieces[5].length = libData->MTensor_getFlattenedLength(t33);

        pieces[6].data   = libData->MTensor_getRealData(t14);
        pieces[6].length = libData->MTensor_getFlattenedLength(t14);

        pieces[7].data   = libData->MTensor_getRealData(t24);
        pieces[7].length = libData->MTensor_getFlattenedLength(t24);

        pieces[8].data   = libData->MTensor_getRealData(t34);
        pieces[8].length = libData->MTensor_getFlattenedLength(t34);

        pieces[9].data   = libData->MTensor_getRealData(t44);
        pieces[9].length = libData->MTensor_getFlattenedLength(t44);

        // 3) Allocate the final result
        mint totalLen = pieces[0].length + pieces[1].length + pieces[2].length + pieces[3].length +pieces[4].length +pieces[5].length + pieces[6].length+pieces[7].length+pieces[8].length + pieces[9].length;
        mint dims[1] = { totalLen };
        MTensor result;
        err = libData->MTensor_new(MType_Real, 1, dims, &result);
        
        if (err) {
            libData->MTensor_free(t11);
            libData->MTensor_free(t12);
            libData->MTensor_free(t22);
            libData->MTensor_free(t13);
            libData->MTensor_free(t23);
            libData->MTensor_free(t33);
            libData->MTensor_free(t14);
            libData->MTensor_free(t24);
            libData->MTensor_free(t34);
            libData->MTensor_free(t44);
            return err;
        }

        mreal *resPtr = libData->MTensor_getRealData(result);

        // 4) One generic copy loop
        mint offset = 0;
        for (int k = 0; k < 10; ++k) {
            mreal *src = pieces[k].data;
            mint   len = pieces[k].length;
            #pragma omp simd
            for (mint i = 0; i < len; ++i) {
                resPtr[offset + i] = src[i];
            }
            offset += len;
        }

        // 5) Clean up temporaries
        libData->MTensor_free(P3);
        libData->MTensor_free(P4);
        libData->MTensor_free(t11);
        libData->MTensor_free(t12);
        libData->MTensor_free(t22);
        libData->MTensor_free(t13);
        libData->MTensor_free(t23);
        libData->MTensor_free(t33);
        libData->MTensor_free(t14);
        libData->MTensor_free(t24);
        libData->MTensor_free(t34);
        libData->MTensor_free(t44);
        

        MArgument_setMTensor(Res, result);
        return LIBRARY_NO_ERROR;

    }else if (Order==5) {
         // declare t11, t12, ..., t55
        MTensor t11, t12, t22, t13, t23, t33,t14, t24, t34, t44, t15, t25, t35, t45, t55, P3, P4, P5;

        //calculate p3
        err = Deltap3(libData, DeltaP, &P3); if (err) return err;
        
        //calculate p4
        err = Deltap4(libData, DeltaP, &P4); if (err) return err;
        
        //calculate p5
        err = Deltap5(libData, DeltaP, &P5); if (err) return err;
        
        //calculate t11
        err = Deltap2(libData, DeltaP, &t11); if (err) return err;
        
        //calculate t12
        err = TensorProduct(libData, DeltaP, t11, &t12); if (err) return err;
        
        //calculate t22
        err = TensorProduct(libData, t11, t11, &t22); if (err) return err;

        //calculate t13
        err = TensorProduct(libData, DeltaP, P3, &t13); if (err) return err;

        //calculate t23
        err = TensorProduct(libData, t11, P3, &t23); if (err) return err;

        //calculate t33
        err = TensorProduct(libData, P3, P3, &t33); if (err) return err;

        //calculate t14
        err = TensorProduct(libData, DeltaP, P4, &t14); if (err) return err;
        
        //calculate t24
        err = TensorProduct(libData, t11, P4, &t24); if (err) return err;

        //calculate t34
        err = TensorProduct(libData, P3, P4, &t34); if (err) return err;

        //calculate t44
        err = TensorProduct(libData, P4, P4, &t44); if (err) return err;

        // calculate t15 
        err = TensorProduct(libData, DeltaP, P5, &t15); if (err) return err;

        // calculate t25 
        err = TensorProduct(libData, t11, P5, &t25); if (err) return err;

        // calculate t35 
        err = TensorProduct(libData, P3, P5, &t35); if (err) return err;

        // calculate t45 
        err = TensorProduct(libData, P4, P5, &t45); if (err) return err;

        // calculate t55
        err = TensorProduct(libData, P5, P5, &t55); if (err) return err;
        
        // 2) Build an array of (pointer, length) pairs
        Piece pieces[15];
        pieces[0].data   = libData->MTensor_getRealData(t11);
        pieces[0].length = libData->MTensor_getFlattenedLength(t11);
        
        pieces[1].data   = libData->MTensor_getRealData(t12);
        pieces[1].length = libData->MTensor_getFlattenedLength(t12);
        
        pieces[2].data   = libData->MTensor_getRealData(t22);
        pieces[2].length = libData->MTensor_getFlattenedLength(t22);
        
        pieces[3].data   = libData->MTensor_getRealData(t13);
        pieces[3].length = libData->MTensor_getFlattenedLength(t13);

        pieces[4].data   = libData->MTensor_getRealData(t23);
        pieces[4].length = libData->MTensor_getFlattenedLength(t23);

        pieces[5].data   = libData->MTensor_getRealData(t33);
        pieces[5].length = libData->MTensor_getFlattenedLength(t33);

        pieces[6].data   = libData->MTensor_getRealData(t14);
        pieces[6].length = libData->MTensor_getFlattenedLength(t14);

        pieces[7].data   = libData->MTensor_getRealData(t24);
        pieces[7].length = libData->MTensor_getFlattenedLength(t24);

        pieces[8].data   = libData->MTensor_getRealData(t34);
        pieces[8].length = libData->MTensor_getFlattenedLength(t34);

        pieces[9].data   = libData->MTensor_getRealData(t44);
        pieces[9].length = libData->MTensor_getFlattenedLength(t44);

        pieces[10].data   = libData->MTensor_getRealData(t15);
        pieces[10].length = libData->MTensor_getFlattenedLength(t15);

        pieces[11].data   = libData->MTensor_getRealData(t25);
        pieces[11].length = libData->MTensor_getFlattenedLength(t25);

        pieces[12].data   = libData->MTensor_getRealData(t35);
        pieces[12].length = libData->MTensor_getFlattenedLength(t35);

        pieces[13].data   = libData->MTensor_getRealData(t45);
        pieces[13].length = libData->MTensor_getFlattenedLength(t45);

        pieces[14].data   = libData->MTensor_getRealData(t55);
        pieces[14].length = libData->MTensor_getFlattenedLength(t55);


        // 3) Allocate the final result
        mint totalLen=0;
        for (mint dummy=0; dummy<15; dummy++) {totalLen +=pieces[dummy].length;}

        mint dims[1] = { totalLen };
        MTensor result;
        err = libData->MTensor_new(MType_Real, 1, dims, &result);
        
        if (err) {
            libData->MTensor_free(t11);
            libData->MTensor_free(t12);
            libData->MTensor_free(t22);
            libData->MTensor_free(t13);
            libData->MTensor_free(t23);
            libData->MTensor_free(t33);
            libData->MTensor_free(t14);
            libData->MTensor_free(t24);
            libData->MTensor_free(t34);
            libData->MTensor_free(t44);
            libData->MTensor_free(t15);
            libData->MTensor_free(t25);
            libData->MTensor_free(t35);
            libData->MTensor_free(t45);
            libData->MTensor_free(t55);
            return err;
        }

        mreal *resPtr = libData->MTensor_getRealData(result);

        // 4) One generic copy loop
        mint offset = 0;
        for (int k = 0; k < 15; ++k) {
            mreal *src = pieces[k].data;
            mint   len = pieces[k].length;
            #pragma omp simd
            for (mint i = 0; i < len; ++i) {
                resPtr[offset + i] = src[i];
            }
            offset += len;
        }

        // 5) Clean up temporaries
        libData->MTensor_free(P3);
        libData->MTensor_free(P4);
        libData->MTensor_free(P5);
        libData->MTensor_free(t11);
        libData->MTensor_free(t12);
        libData->MTensor_free(t22);
        libData->MTensor_free(t13);
        libData->MTensor_free(t23);
        libData->MTensor_free(t33);
        libData->MTensor_free(t14);
        libData->MTensor_free(t24);
        libData->MTensor_free(t34);
        libData->MTensor_free(t44);
        libData->MTensor_free(t15);
        libData->MTensor_free(t25);
        libData->MTensor_free(t35);
        libData->MTensor_free(t45);
        libData->MTensor_free(t55);
        

        MArgument_setMTensor(Res, result);
        return LIBRARY_NO_ERROR;

    } else {
        return LIBRARY_FUNCTION_ERROR; 
    }

    return LIBRARY_NO_ERROR;

}
