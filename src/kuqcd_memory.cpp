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
#include "kuqcd_memory.h"
#include <omp.h>
#include <cstdlib>
#include <cstring>
#include "kuqcd_utils.h"

#ifdef USE_KUPL
#include "kupl.h"
#endif

#ifdef USE_HPMEM
#include <sys/mman.h>
#endif

namespace kuqcd {

KMemory::KMemory(const KMemory &mem)
{
    InitMem(mem.size);
    Copy(mem);
}

KMemory::KMemory(size_t sz)
{
    InitMem(sz);
}

KMemory &KMemory::operator=(const KMemory &mem)
{
    Copy(mem);
    return *this;
}

void KMemory::InitMem(size_t sz)
{
    size = sz;
    if (sz > 0) {
#ifdef USE_HPMEM
        rawPtr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#else
        rawPtr = std::aligned_alloc(64, sz);
#endif
    if (rawPtr == nullptr) {
        throw(std::bad_alloc());
    }
    } else {
        rawPtr = nullptr;
    }
}

bool KMemory::Alloc(size_t sz)
{
    bool resize = size < sz;
    if (resize) {
        if (size > 0) {
            Free();
        }
        InitMem(sz);
        return true;
    } else {
        return false;
    }
}

void KMemory::Copy(const kuqcd::KMemory &mem)
{
    if (size != mem.size) {
        Rank0Fprintf(stderr, "memory copy size error! %s %d\n", __FILE__, __LINE__);
    }
    if (size > 0) {
        char *ptr = static_cast<char *>(rawPtr);
        char *ptr1 = static_cast<char *>(mem.rawPtr);
#pragma omp parallel for
        for (size_t i = 0; i < size; ++i) {
            ptr[i] = ptr1[i];
        }
    }
}

void KMemory::Clear()
{
    if (size > 0) {
        char *ptr = static_cast<char *>(rawPtr);
#pragma omp parallel for
        for (size_t i = 0; i < size; ++i) {
            ptr[i] = 0;
        }
    }
}

void KMemory::Free()
{
    if (rawPtr != nullptr) {
#ifdef USE_HPMEM
        munmap(rawPtr, size);
#else
        free(rawPtr);
#endif
        rawPtr = nullptr;
    }
    size = 0;
}

HaloMemory::HaloMemory(size_t sz)
{
    InitHaloMem(sz);
}

HaloMemory::HaloMemory(const kuqcd::HaloMemory &hmem)
{
    InitHaloMem(hmem.size);
    if (size > 0) {
        char *ptr = static_cast<char *>(this->rawPtr);
        char *ptr1 = static_cast<char *>(hmem.rawPtr);
#pragma omp parallel for
        for (size_t i = 0; i < size; ++i) {
            ptr[i] = ptr1[i];
        }
    }
}

void HaloMemory::InitHaloMem(size_t sz)
{
    this->size = sz;
    if (sz > 0) {
#ifdef USE_KUPL_COMM
        int np;
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        oob_cbs_h = &oob_cbs;
        oob_cbs_h->oob_allgather = oob_allgather_callback;
        oob_cbs_h->oob_barrier = oob_barrier_callback;
        kupl_shm_comm_create(np, rank, getpid(), oob_cbs_h, (void *)comm, &kupl_comm);
        kupl_shm_win_alloc(sz, kupl_comm, &(this->rawPtr), &win);
#elif USE_HPMEM
        rawPtr = mmap(NULL, sz, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#else
        this->rawPtr = std::aligned_alloc(64, sz);
#endif
    if (this->rawPtr == nullptr) {
        throw(std::bad_alloc());
    }
    } else {
        this->rawPtr = nullptr;
    }
}

bool HaloMemory::Alloc(size_t sz)
{
    bool resize = this->size < sz;
    if (resize) {
        if (this->size > 0) {
            Free();
        }
        InitHaloMem(sz);
        return true;
    } else {
        return false;
    }
}

void HaloMemory::Free()
{
    if (this->rawPtr != nullptr) {
#ifdef USE_KUPL_COMM
        kupl_shm_win_free(win);
        kupl_shm_comm_destroy(kupl_comm);
#elif USE_HPMEM
        munmap(this->rawPtr, size);
#else
        free(this->rawPtr);
#endif
        this->rawPtr = nullptr;
    }
    this->size = 0;
}

Index::Index(const Index &other)
{
    Init(other.length);
#pragma omp parallel for
    for (size_t i = 0; i < length; ++i) {
        idx[i] = other.idx[i];
    }
}

Index::Index(size_t len)
{
    Init(len);
}

Index &Index::operator=(const Index &other)
{
    if (length != other.length) {
        Rank0Fprintf(stderr, "memory copy size error! %s %d\n", __FILE__, __LINE__);
    }
    if (length > 0) {
#pragma omp parallel for
        for (size_t i = 0; i < length; ++i) {
            idx[i] = other.idx[i];
        }
    }
    return *this;
}

void Index::Init(size_t len)
{
    length = len;
    if (len > 0) {
        idx = static_cast<long *>(std::aligned_alloc(64, len * sizeof(long)));
    if (idx == nullptr) {
        throw(std::bad_alloc());
    }
    } else {
        idx = nullptr;
    }
}

bool Index::Alloc(size_t len)
{
    bool resize = length < len;
    if (resize) {
        if (len > 0) {
            Free();
        }
        Init(len);
        return true;
    } else {
        return false;
    }
}

void Index::Free()
{
    if (idx != nullptr) {
        free(idx);
        idx = nullptr;
    }
}

}  // namespace kuqcd