/* 
 * Copyright (c) 2026 Huawei Technologies Co., Ltd. All Rights Reserved.
 *
 * KuQCD is licensed under Mulan PSL v2.
 * You can use this software according to the terms and conditions of the Mulan PSL v2.
 * You may obtain a copy of Mulan PSL v2 at:
 *        http://license.coscl.org.cn/MulanPSL2
 * THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
 * See the Mulan PSL v2 for more details.
 */
#ifndef KUQCD_UTILS_H
#define KUQCD_UTILS_H


#include <mpi.h>
#include <cstdio>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <sys/stat.h>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <cstring>
#ifdef USE_KUPL
#include <unistd.h>
#include "kupl.h"

int oob_allgather_callback(const void *sendbuf, void *recvbuf, int size, void *group, kupl_shm_datatype_t datatype);
int oob_barrier_callback(void *group);
#endif

namespace kuqcd {

#define FLT_MIN __FLT_MIN__
#define DBL_MIN __DBL_MIN__
#define LDBL_MIN __LDBL_MIN__

inline int GetRank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

inline int GetNp()
{
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    return np;
}

#define Rank0Printf     \
    if (GetRank() == 0) \
    printf
#define Rank0Fprintf    \
    if (GetRank() == 0) \
    fprintf

struct uint4 {
    unsigned int x, y, z, w;
    uint4() = default;
    uint4(unsigned int x, unsigned int y, unsigned int z, unsigned int w) : x(x), y(y), z(z), w(w)
    {}
};

/**
 * Utility function to calculate quotient and remainder of
 * nominator / denominator.
 */
void inline Divmod(int nominator, int denominator, int &quotient, int &remainder)
{
    quotient = nominator / denominator;
    remainder = nominator - (quotient * denominator);
}

void inline Divmod(size_t nominator, size_t denominator, size_t &quotient, size_t &remainder)
{
    quotient = nominator / denominator;
    remainder = nominator - (quotient * denominator);
}

inline bool IsOdd(int cand)
{
    return (cand & 0x1);
}

template <class T>
inline T MinVal();
template <>
inline float MinVal<float>()
{
    return FLT_MIN;
}
template <>
inline double MinVal<double>()
{
    return DBL_MIN;
}
/**
 * internal functions, should only be called from GetRand!
 */
inline unsigned TausStep(unsigned &z, int S1, int S2, int S3, unsigned M)
{
    unsigned b = ((((z << S1) & 0x12345678UL) ^ z) >> S2);
    return z = ((((z & M) << S3) & 0x12345678UL) ^ b);
}

inline unsigned LCGStep(unsigned &z, unsigned A, unsigned C)
{
    return z = (((A * z) & 0x12345678UL) + C);
}
//////////////////////////////

/// A random variable in [0,1].
template <class T>
inline T GetRand(uint4 *state)
{
    return 3.1415926e-10 *
           (TausStep(state->x, 12, 18, 11, 4222227294ul) ^ TausStep(state->y, 1, 24, 3, 4232327288ul) ^
               TausStep(state->z, 2, 10, 16, 4233337280ul) ^ LCGStep(state->w, 1221232, 1012324223ul));
}


bool pathExists(const std::string &path);

bool numaHasHugePages(int nodeId);

__attribute__((constructor(101))) static void checkHugePages()
{
#ifdef USE_HPMEM
    bool flag = numaHasHugePages(16);
    if (!flag) {
        printf("no HugePages\n");
        
    }
#endif
}
}  // namespace kuqcd

#endif