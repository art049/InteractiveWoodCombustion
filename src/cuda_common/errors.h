/* cf. Jason Sanders, Edward Kandrot. CUDA by Example: An Introduction to General-Purpose GPU Programming */
// I found this copyright notice off of jiekebo's repo:
/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation.
 * Any use, reproduction, disclosure, or distribution of this software
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA)
 * associated with this source code for terms and conditions that govern
 * your use of this NVIDIA software.
 *
 */
#ifndef __ERRORS_H__
#define __ERRORS_H__
#include <stdio.h>
#include <cuda_runtime.h>
#include <cstdlib>

inline void HandleError(cudaError_t err,
                        const char *file,
                        int line) {
    if (err != cudaSuccess) {
        printf( "Cuda error : %s in %s at line %d\n", cudaGetErrorString(err), 
                file, line );
        exit( EXIT_FAILURE );
    }
}

#define HANDLE_ERROR( err) (HandleError( err, __FILE__, __LINE__ ))



#endif // __ERRORS_H__
