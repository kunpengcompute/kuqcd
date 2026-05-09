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
#ifndef KUQCD_MEMORY_H
#define KUQCD_MEMORY_H

#include <cstdlib>
#include "kuqcd_utils.h"

#ifdef USE_KUPL
#include "kupl.h"
#include <mpi.h>
#endif

namespace kuqcd {

class KMemory {
public:
    size_t size = 0;
    void *rawPtr = nullptr;
    KMemory() = default;
    KMemory(const KMemory &mem);
    KMemory(size_t sz);
    virtual KMemory &operator=(const KMemory &mem);

    void InitMem(size_t sz);
    virtual bool Alloc(size_t sz);
    void Copy(const KMemory &mem);
    void Clear();
    virtual void Free();

    virtual ~KMemory()
    {
        Free();
    }
};

class HaloMemory : public KMemory {
public:
#ifdef USE_KUPL_COMM
    kupl_shm_oob_cb_t oob_cbs;
    kupl_shm_oob_cb_h oob_cbs_h;
    kupl_shm_comm_h kupl_comm;
    kupl_shm_win_h win;
    MPI_Comm comm = MPI_COMM_WORLD;
#endif
    HaloMemory() = default;
    HaloMemory(size_t sz);
    HaloMemory(const HaloMemory &hmem);

    void InitHaloMem(size_t sz);
    virtual bool Alloc(size_t sz);
    virtual void Free();
    virtual ~HaloMemory()
    {
        Free();
    }
};

class Index {
public:
    size_t length = 0;
    long *idx = nullptr;
    Index() = default;
    Index(const Index &other);
    Index(size_t len);
    Index &operator=(const Index &other);
    void Init(size_t len);
    bool Alloc(size_t len);
    void Free();
    ~Index()
    {
        Free();
    }
};

} // namespace kuqcd
#endif
