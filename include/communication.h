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
#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <cstdlib>
#include <cstdint>
#include <string>
#include <mpi.h>
#include "space_time.h"
#include "spinor.h"
#include "su3.h"
#include "kuqcd_memory.h"
#include "latt_base.h"

namespace kuqcd {

void Root2All(int &value);
void Root2All(int64_t &value);
void Root2All(size_t &value);
void Root2All(bool &value);
void Root2All(float &value);
void Root2All(double &value);

int Reduce(int in);
uint32_t Reduce(uint32_t in);
float Reduce(float in);
double Reduce(double in);

template <typename FieldT>
class Communicator {
public:
    int myRank = 0;
    int nbRank = 0;
    int tagSend;
    int tagRecv;
    size_t length = 0;
    MPI_Request stateSend, stateRecv;
    HaloMemory memSend;
    FieldT *dataSend = nullptr;
#ifndef USE_KUPL_COMM
    HaloMemory memRecv;
    FieldT *dataRecv = nullptr;
#endif
    Communicator() = default;
    Communicator(const Communicator &comm);
    void Init(int nbr, size_t len, size_t size, int ts, int tr);
#ifndef USE_KUPL_COMM
    void Send();
    void Recv();
    void Finish(bool recv);
#endif
    void Clear();
    void Print();
    ~Communicator()
    {
        dataSend = nullptr;
#ifndef USE_KUPL_COMM
        dataRecv = nullptr;
#endif
    }
};

#define DECLARE_EXTERN_COMMUNICATOR_TEMPLATE(T, Nbatch)                       \
    extern template class Communicator<ColorSpinor<T, Nbatch>>;               \
    extern template class Communicator<SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO>>; \
    extern template class Communicator<SU3<T, Nbatch, KUQCD_GAUGE_RECON_12>>; \
    extern template class Communicator<SU3<T, Nbatch, KUQCD_GAUGE_RECON_U14>>

DECLARE_EXTERN_COMMUNICATOR_TEMPLATE(float, 1);
DECLARE_EXTERN_COMMUNICATOR_TEMPLATE(double, 1);
DECLARE_EXTERN_COMMUNICATOR_TEMPLATE(float, 16);
DECLARE_EXTERN_COMMUNICATOR_TEMPLATE(double, 8);

/* only for gaugefield */
class CommBase {
public:
    LattBase &latBase;

    /* neigbhor rank and buffer length for all geometry */
    /* +X, -X, +Y, -Y, +Z, -Z, +T, -T */
    std::array<std::array<int, 2>, 4> rankHyperPlane;
    std::array<size_t, 4> lenHyperPlane;
    /* +X+Y, -X+Y, +X-Y, -X-Y, +X+Z, -X+Z, +X-Z, -X-Z,
     * +X+T, -X+T, +X-T, -X-T, +Y+Z, -Y+Z, +Y-Z, -Y-Z,
     * +Y+T, -Y+T, +Y-T, -Y-T, +Z+T, -Z+T, +Z-T, -Z-T */
    std::array<std::array<int, 4>, 6> rankPlane;
    std::array<size_t, 6> lenPlane;
    /* +X+Y+Z, -X+Y+Z, +X-Y+Z, -X-Y+Z, +X+Y-Z, -X+Y-Z, +X-Y-Z, -X-Y-Z,
     * +X+Y+T, -X+Y+T, +X-Y+T, -X-Y+T, +X+Y-T, -X+Y-T, +X-Y-T, -X-Y-T,
     * +X+Z+T, -X+Z+T, +X-Z+T, -X-Z+T, +X+Z-T, -X+Z-T, +X-Z-T, -X-Z-T,
     * +Y+Z+T, -Y+Z+T, +Y-Z+T, -Y-Z+T, +Y+Z-T, -Y+Z-T, +Y-Z-T, -Y-Z-T */
    std::array<std::array<int, 8>, 4> rankStripe;
    std::array<size_t, 4> lenStripe;
    /* +X+Y+Z+T, -X+Y+Z+T, +X-Y+Z+T, -X-Y+Z+T,
     * +X+Y-Z+T, -X+Y-Z+T, +X-Y-Z+T, -X-Y-Z+T,
     * +X+Y+Z-T, -X+Y+Z-T, +X-Y+Z-T, -X-Y+Z-T,
     * +X+Y-Z-T, -X+Y-Z-T, +X-Y-Z-T, -X-Y-Z-T */
    std::array<int, 16> rankCorner;
    size_t lenCorner;

    Index nbIdx;

    std::array<std::array<Index, 2>, 4> haloIdxSend1;
    std::array<std::array<Index, 2>, 4> haloIdxRecv1;
    std::array<std::array<Index, 4>, 6> haloIdxSend2;
    std::array<std::array<Index, 4>, 6> haloIdxRecv2;
    std::array<std::array<Index, 8>, 4> haloIdxSend3;
    std::array<std::array<Index, 8>, 4> haloIdxRecv3;
    std::array<Index, 16> haloIdxSend4;
    std::array<Index, 16> haloIdxRecv4;

    CommBase() = delete;
    CommBase(const CommBase &comm);
    CommBase(LattBase &latBase);
    CommBase &operator=(const CommBase &comm) = delete;

    void InitNbRank();
    void InitNbRank1(const Position &pos, const Lattice &grid, int mu);
    void InitNbRank2(const Position &pos, const Lattice &grid, int idx, int mu, int nu);
    void InitNbRank3(const Position &pos, const Lattice &grid, int idx, int mu, int nu, int rho);
    void InitNbRank4(const Position &pos, const Lattice &grid);
    void InitNbIdx();
    void InitNbIdxParity(KuQCDSiteSubset site_subset);
    void InitHaloIdx();
    void InitHaloIdx1(int mu);
    void InitHaloIdx2(int idx, int mu, int nu);
    void InitHaloIdx3(int idx, int mu, int nu, int rho);
    void InitHaloIdx4();
};

}  // namespace kuqcd

#endif
