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
#include <cstring>

#include <mpi.h>

#include "communication.h"
#include "logger.h"
#include "kuqcd_utils.h"
#include "index_utils.h"

namespace kuqcd {

void Root2All(int &value)
{
    MPI_Bcast(&value, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
void Root2All(int64_t &value)
{
    MPI_Bcast(&value, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
}
void Root2All(size_t &value)
{
    MPI_Bcast(&value, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
}
void Root2All(bool &value)
{
    MPI_Bcast(&value, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
}
void Root2All(float &value)
{
    MPI_Bcast(&value, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
}
void Root2All(double &value)
{
    MPI_Bcast(&value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

int Reduce(int in)
{
    int recv;
    MPI_Allreduce(&in, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return recv;
}
uint32_t Reduce(uint32_t in)
{
    uint32_t recv;
    MPI_Allreduce(&in, &recv, 1, MPI_UINT32_T, MPI_SUM, MPI_COMM_WORLD);
    return recv;
}
float Reduce(float in)
{
    float recv;
    MPI_Allreduce(&in, &recv, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    return recv;
}
double Reduce(double in)
{
    double recv;
    MPI_Allreduce(&in, &recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return recv;
}

template <typename FieldT>
Communicator<FieldT>::Communicator(const Communicator<FieldT> &comm)
    : myRank(comm.myRank), nbRank(comm.nbRank), tagSend(comm.tagSend), tagRecv(comm.tagRecv), length(comm.length),
      stateSend(comm.stateSend), stateRecv(comm.stateRecv), memSend(comm.memSend)
{
    dataSend = static_cast<FieldT *>(memSend.rawPtr);
#ifndef USE_KUPL_COMM
    memRecv.InitHaloMem(comm.memRecv.size);
    memRecv = comm.memRecv;
    dataRecv = static_cast<FieldT *>(memRecv.rawPtr);
#endif
}

template <typename FieldT>
void Communicator<FieldT>::Init(int nbr, size_t len, size_t size, int ts, int tr)
{
    myRank = GetRank();
    nbRank = nbr;
    length = len;
    tagSend = ts;
    tagRecv = tr;
    memSend.InitHaloMem(size);
    dataSend = static_cast<FieldT *>(memSend.rawPtr);
#ifndef USE_KUPL_COMM
    memRecv.InitHaloMem(size);
    dataRecv = static_cast<FieldT *>(memRecv.rawPtr);
#endif
}

#ifndef USE_KUPL_COMM

template <typename FieldT>
void Communicator<FieldT>::Send()
{
    MPI_Issend(memSend.rawPtr, memSend.size, MPI_CHAR, nbRank, tagSend, MPI_COMM_WORLD, &stateSend);
}

template <typename FieldT>
void Communicator<FieldT>::Recv()
{
    MPI_Irecv(memRecv.rawPtr, memSend.size, MPI_CHAR, nbRank, tagRecv, MPI_COMM_WORLD, &stateRecv);
}

template <typename FieldT>
void Communicator<FieldT>::Finish(bool recv)
{
    if (recv) {
        MPI_Wait(&stateRecv, MPI_STATUS_IGNORE);
        if (stateRecv != MPI_REQUEST_NULL) {
            Rank0Fprintf(stderr, "receive state error!\n");
        }
    } else {
        MPI_Wait(&stateSend, MPI_STATUS_IGNORE);
        if (stateSend != MPI_REQUEST_NULL) {
            Rank0Fprintf(stderr, "send state error!\n");
        }
    }
}

#endif

template <typename FieldT>
void Communicator<FieldT>::Clear()
{
    memSend.Clear();
#ifndef USE_KUPL_COMM
    memRecv.Clear();
#endif
}

template <typename FieldT>
void Communicator<FieldT>::Print()
{
    Rank0Printf("%d %d %d %d %zu\n", myRank, nbRank, tagSend, tagRecv, length);
}

#define INSTANTIATE_COMMUNICATOR_TEMPLATE(T, Nbatch)                   \
    template class Communicator<ColorSpinor<T, Nbatch>>;               \
    template class Communicator<SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO>>; \
    template class Communicator<SU3<T, Nbatch, KUQCD_GAUGE_RECON_12>>; \
    template class Communicator<SU3<T, Nbatch, KUQCD_GAUGE_RECON_U14>>
INSTANTIATE_COMMUNICATOR_TEMPLATE(float, 1);
INSTANTIATE_COMMUNICATOR_TEMPLATE(double, 1);
INSTANTIATE_COMMUNICATOR_TEMPLATE(float, 16);
INSTANTIATE_COMMUNICATOR_TEMPLATE(double, 8);

CommBase::CommBase(const CommBase &comm)
    : latBase(comm.latBase), rankHyperPlane(comm.rankHyperPlane), lenHyperPlane(comm.lenHyperPlane),
      rankPlane(comm.rankPlane), lenPlane(comm.lenPlane), rankStripe(comm.rankStripe), lenStripe(comm.lenStripe),
      rankCorner(comm.rankCorner), lenCorner(comm.lenCorner), nbIdx(comm.nbIdx), haloIdxSend1(comm.haloIdxSend1),
      haloIdxRecv1(comm.haloIdxRecv1), haloIdxSend2(comm.haloIdxSend2), haloIdxRecv2(comm.haloIdxRecv2),
      haloIdxSend3(comm.haloIdxSend3), haloIdxRecv3(comm.haloIdxRecv3), haloIdxSend4(comm.haloIdxSend4),
      haloIdxRecv4(comm.haloIdxRecv4)
{}

CommBase::CommBase(LattBase &latBase_) : latBase(latBase_)
{
    InitNbRank();
    InitNbIdx();
    InitHaloIdx();
}

void CommBase::InitNbRank()
{
    Lattice grid(latBase.gridDim[0], latBase.gridDim[1], latBase.gridDim[2], latBase.gridDim[3]);
    Position pos;
    if (latBase.rankOrder == KUQCD_XYZT_ORDER) {
        pos = GetPosition(latBase.myRank, grid);
    } else {
        pos = GetAntiPosition(latBase.myRank, grid);
    }

    /* hyper plane */
    for (int mu = 0; mu < 4; ++mu) {
        InitNbRank1(pos, grid, mu);
    }

    /* plane */
    int idx = 0;
    for (int mu = 0; mu < 3; ++mu) {
        for (int nu = mu + 1; nu < 4; ++nu) {
            InitNbRank2(pos, grid, idx, mu, nu);
            idx++;
        }
    }

    /* stripe */
    idx = 0;
    for (int mu = 0; mu < 2; ++mu) {
        for (int nu = mu + 1; nu < 3; ++nu) {
            for (int rho = nu + 1; rho < 4; ++rho) {
                InitNbRank3(pos, grid, idx, mu, nu, rho);
                ++idx;
            }
        }
    }

    /* corner */
    InitNbRank4(pos, grid);
}

void CommBase::InitNbRank1(const Position &pos, const Lattice &grid, int mu)
{
    if (latBase.gridDim[mu] > 1) {
        Position myPos;
        Position dPos(0, 0, 0, 0);
        size_t surf = 1;
        for (int dir = 0; dir < 4; ++dir) {
            if (dir == mu) {
                surf *= latBase.haloDepth;
            } else {
                surf *= latBase.latSubBulkFL[dir];
            }
        }
        lenHyperPlane[mu] = surf;

        int i = 0;
        for (int dx = 1; dx >= -1; dx -= 2) {
            dPos[mu] = dx;
            myPos = pos + dPos;
            myPos.Normalize(mu, grid);
            if (latBase.rankOrder == KUQCD_XYZT_ORDER) {
                rankHyperPlane[mu][i] = myPos.Index(grid);
            } else {
                rankHyperPlane[mu][i] = myPos.IndexAnti(grid);
            }
            ++i;
        }
    } else {
        for (int i = 0; i < 2; ++i) {
            rankHyperPlane[mu][i] = latBase.myRank;
        }
        lenHyperPlane[mu] = 0;
    }
}

void CommBase::InitNbRank2(const Position &pos, const Lattice &grid, int idx, int mu, int nu)
{
    if ((latBase.gridDim[mu] > 1) && (latBase.gridDim[nu] > 1)) {
        Position myPos;
        Position dPos(0, 0, 0, 0);
        size_t surf = 1;
        for (int dir = 0; dir < 4; ++dir) {
            if ((dir == mu) || (dir == nu)) {
                surf *= latBase.haloDepth;
            } else {
                surf *= latBase.latSubBulkFL[dir];
            }
        }
        lenPlane[idx] = surf;

        int i = 0;
        for (int dy = 1; dy >= -1; dy -= 2) {
            for (int dx = 1; dx >= -1; dx -= 2) {
                dPos[mu] = dx;
                dPos[nu] = dy;
                myPos = pos + dPos;
                myPos.Normalize(mu, grid);
                myPos.Normalize(nu, grid);
                if (latBase.rankOrder == KUQCD_XYZT_ORDER) {
                    rankPlane[idx][i] = myPos.Index(grid);
                } else {
                    rankPlane[idx][i] = myPos.IndexAnti(grid);
                }
                ++i;
            }
        }
    } else {
        for (int i = 0; i < 4; ++i) {
            rankPlane[idx][i] = latBase.myRank;
        }
        lenPlane[idx] = 0;
    }
}

void CommBase::InitNbRank3(const Position &pos, const Lattice &grid, int idx, int mu, int nu, int rho)
{
    if ((latBase.gridDim[mu] > 1) && (latBase.gridDim[nu] > 1) && (latBase.gridDim[rho] > 1)) {
        Position myPos;
        Position dPos(0, 0, 0, 0);
        size_t surf = 1;
        for (int dir = 0; dir < 4; ++dir) {
            if ((dir == mu) || (dir == nu) || (dir == rho)) {
                surf *= latBase.haloDepth;
            } else {
                surf *= latBase.latSubBulkFL[dir];
            }
        }
        lenStripe[idx] = surf;

        int i = 0;
        for (int dz = 1; dz >= -1; dz -= 2) {
            for (int dy = 1; dy >= -1; dy -= 2) {
                for (int dx = 1; dx >= -1; dx -= 2) {
                    dPos[mu] = dx;
                    dPos[nu] = dy;
                    dPos[rho] = dz;
                    myPos = pos + dPos;
                    myPos.Normalize(mu, grid);
                    myPos.Normalize(nu, grid);
                    myPos.Normalize(rho, grid);
                    if (latBase.rankOrder == KUQCD_XYZT_ORDER) {
                        rankStripe[idx][i] = myPos.Index(grid);
                    } else {
                        rankStripe[idx][i] = myPos.IndexAnti(grid);
                    }
                    ++i;
                }
            }
        }
    } else {
        for (int i = 0; i < 8; ++i) {
            rankStripe[idx][i] = latBase.myRank;
        }
        lenStripe[idx] = 0;
    }
}

void CommBase::InitNbRank4(const Position &pos, const Lattice &grid)
{
    if ((latBase.gridDim[0] > 1) && (latBase.gridDim[1] > 1) && (latBase.gridDim[2] > 1) && (latBase.gridDim[3] > 1)) {
        Position myPos;
        Position dPos(0, 0, 0, 0);
        lenCorner = latBase.haloDepth * latBase.haloDepth * latBase.haloDepth * latBase.haloDepth;

        for (int dt = 1; dt >= -1; dt -= 2) {
            for (int dz = 1; dz >= -1; dz -= 2) {
                for (int dy = 1; dy >= -1; dy -= 2) {
                    for (int dx = 1; dx >= -1; dx -= 2) {
                        dPos[3] = dt;
                        dPos[2] = dz;
                        dPos[1] = dy;
                        dPos[0] = dx;
                        myPos = pos + dPos;
                        myPos.Normalize(0, grid);
                        myPos.Normalize(1, grid);
                        myPos.Normalize(2, grid);
                        myPos.Normalize(3, grid);
                    }
                }
            }
        }
    } else {
        for (int i = 0; i < 16; ++i) {
            rankCorner[i] = latBase.myRank;
        }
        lenCorner = 0;
    }
}

void CommBase::InitNbIdx()
{
    nbIdx.Init(8 * latBase.lenBulk);
    InitNbIdxParity(KUQCD_EVEN);
    InitNbIdxParity(KUQCD_ODD);
}

void CommBase::InitNbIdxParity(KuQCDSiteSubset site_subset)
{
    size_t offsetnb;
    size_t offseti;
    offseti = (site_subset == KUQCD_EVEN) ? 0 : latBase.latSubBulkCB.V;
    offsetnb = (site_subset == KUQCD_EVEN) ? latBase.latSubFullCB.V : 0;
#pragma omp parallel for
    for (size_t i = 0; i < latBase.latSubBulkCB.V; ++i) {
        size_t j = i + offseti;
        // sub bulk CB position
        Position pos = GetPosition(i, latBase.latSubBulkCB);
        // sub bulk FL position
        pos = GetPosFromCBPos(pos, site_subset);
        // sub full FL position
        pos = GetFullPosFromBulkPos(pos, latBase.haloDepth, latBase.gridDim);
        for (int dir = 0; dir < 4; ++dir) {
            Position dx(0, 0, 0, 0);
            // forward neigbhor
            dx[dir] = latBase.haloDepth;
            Position myPos = pos + dx;
            // sub full FL position
            myPos.Normalize(dir, latBase.latSubFullFL);
            // sub full CB position
            myPos = GetCBPos(myPos);
            nbIdx.idx[dir + 8 * j] = myPos.Index(latBase.latSubFullCB) + offsetnb;

            // backward neigbhor
            myPos = pos - dx;
            // sub full FL position
            myPos.Normalize(dir, latBase.latSubFullFL);
            // sub full CB position
            myPos = GetCBPos(myPos);
            nbIdx.idx[dir + 8 * j] = myPos.Index(latBase.latSubFullCB) + offsetnb;
        }
    }
}

void CommBase::InitHaloIdx()
{
    for (int mu = 0; mu < 4; ++mu) {
        InitHaloIdx1(mu);
    }
}

void CommBase::InitHaloIdx1(int mu)
{
    int split = 0;
    for (auto i : latBase.gridDim) {
        if (i > 1)
            split++;
    }
    KuQCDSiteSubset parityshift = ((split & 0x1) && (latBase.haloDepth & 0x1)) ? KUQCD_ODD : KUQCD_EVEN;
    size_t haloLen = lenHyperPlane[mu];
    Lattice latFullCB = latBase.latSubFullCB;
    Lattice latFullFL = latBase.latSubFullFL;
    Lattice latBulkFL = latBase.latSubBulkFL;
    if (latBase.gridDim[mu] > 1) {
        for (int i = 0; i < 2; ++i) {
            haloIdxSend1[mu][i].Alloc(haloLen);
            haloIdxRecv1[mu][i].Alloc(haloLen);
        }
        Lattice latHalo(latBase.latSubBulkFL);
        latHalo[mu] = latBase.haloDepth;
        latHalo.ComputeVolume();
        if (latHalo.V != haloLen) {
            Rank0Fprintf(stderr, "halo length error %zu %zu! %s %d\n", latHalo.V, haloLen, __FILE__, __LINE__);
        }
#pragma omp parallel for
        for (size_t i = 0; i < haloLen; ++i) {
            Position posSend = GetPosition(i, latHalo);
            posSend = GetFullPosFromBulkPos(posSend, latBase.haloDepth, latBase.gridDim);
            Position posRecv = posSend;
            posRecv[mu] -= latBase.haloDepth;
            Position posCB;
            // dir -mu send index
            size_t offset = (posSend.Parity() == parityshift) ? 0 : latBase.latSubFullCB.V;
            posCB = GetCBPos(posSend);
            haloIdxSend1[mu][1].idx[i] = posCB.Index(latFullCB) + offset;
            // dir +mu send index
            posSend[mu] += latBulkFL[mu] - latBase.haloDepth;
            offset = (posSend.Parity() == parityshift) ? 0 : latBase.latSubFullCB.V;
            posCB = GetCBPos(posSend);
            haloIdxSend1[mu][0].idx[i] = posCB.Index(latFullCB) + offset;
            // dir -mu recv index
            offset = (posRecv.Parity() == parityshift) ? 0 : latBase.latSubFullCB.V;
            posCB = GetCBPos(posRecv);
            haloIdxRecv1[mu][1].idx[i] = posCB.Index(latFullCB) + offset;
            // dir +mu recv index
            posRecv[mu] += latFullFL[mu] - latBase.haloDepth;
            offset = (posRecv.Parity() == parityshift) ? 0 : latBase.latSubFullCB.V;
            posCB = GetCBPos(posRecv);
            haloIdxRecv1[mu][0].idx[i] = posCB.Index(latFullCB) + offset;
        }
    }
}

}  // namespace kuqcd