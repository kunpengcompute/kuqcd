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
#include <cassert>

#include "dirac_hisq.h"

#include "index_utils.h"
#include "kuqcd_math.h"
#include "logger.h"

#ifdef USE_KUPL
#include "kupl.h"
#endif

#ifdef SVEBITS
#include "arm_sve.h"
#endif

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
DiracHISQ<T, Nbatch, ReconS, ReconN>::DiracHISQ(KuQCDDslashType diracType, const LinkST &gaugeSmeared_,
    const LinkNT &gaugeNaik_, const double mass_, VecT *tmp_, T naikEps)
    : Dirac<VecT>(diracType), gaugeSmeared(gaugeSmeared_), gaugeNaik(gaugeNaik_), tmp(tmp_),
      latBase(gaugeSmeared_.latBase), mass(mass_), mass2(mass_ * mass_),
      c3000((-1. / 48.0) * (1.0 + (double)naikEps)), gridDim(gaugeSmeared.latBase.gridDim)
{
    InitComm();
    InitNbIdx();
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::InitComm()
{
    np = GetNp();
    myRank = GetRank();

    Lattice grid(gridDim[0], gridDim[1], gridDim[2], gridDim[3]);
    Position pos;
    if (latBase.rankOrder == KUQCD_XYZT_ORDER) {
        pos = GetPosition(myRank, grid);
    } else {
        pos = GetAntiPosition(myRank, grid);
    }
    int locSize = sizeof(SpinorT);

    for (int dir = 0; dir < 4; ++dir) {
        if (gridDim[dir] > 1) {
            size_t srf = 1;
            for (int mu = 0; mu < 4; ++mu) {
                if (mu == dir)
                    continue;
                srf *= latBase.latSubBulkFL[mu];
            }
            if (this->diracType == KUQCD_HISQ_DIRAC) {
                surf[dir] = srf;
            } else if (this->diracType == KUQCD_HISQPC_DIRAC) {
                surf[dir] = srf / 2;
            } else {
                rootLogger.Error("dirac type should be hisq or hisq pc!");
            }

            Position myPos;
            Position dPos(0, 0, 0, 0);
            dPos[dir] = 1;
            myPos = pos + dPos;
            myPos.Normalize(dir, grid);
            if (latBase.rankOrder == KUQCD_XYZT_ORDER) {
                nbRank[dir] = myPos.Index(grid);
            } else {
                nbRank[dir] = myPos.IndexAnti(grid);
            }

            myPos = pos - dPos;
            myPos.Normalize(dir, grid);
            if (latBase.rankOrder == KUQCD_XYZT_ORDER) {
                nbRank[dir + 4] = myPos.Index(grid);
            } else {
                nbRank[dir + 4] = myPos.IndexAnti(grid);
            }
        } else {
            surf[dir] = 0;
            nbRank[dir] = myRank;
            nbRank[dir + 4] = myRank;
        }
        comm[dir].Init(nbRank[dir], surf[dir], surf[dir] * locSize, dir, dir + 4);
        comm[dir + 4].Init(nbRank[dir + 4], surf[dir], surf[dir] * locSize, dir + 4, dir);
        comm3[dir].Init(nbRank[dir], surf[dir] * 3, 3 * surf[dir] * locSize, dir, dir + 4);
        comm3[dir + 4].Init(nbRank[dir + 4], surf[dir] * 3, 3 * surf[dir] * locSize, dir + 4, dir);
    }
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::InitNbIdx()
{
    nb1Idx.Init(8 * latBase.lenBulk);
    nb3Idx.Init(8 * latBase.lenBulk);
    for (int dir = 0; dir < 4; ++dir) {
        edge1Idx[dir].Init(surf[dir]);
        edge1Idx[dir + 4].Init(surf[dir]);
        edge3Idx[dir].Init(3 * surf[dir]);
        edge3Idx[dir + 4].Init(3 * surf[dir]);
    }
    InitNbIdxParity(KUQCD_EVEN);
    InitNbIdxParity(KUQCD_ODD);
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::InitNbIdxParity(KuQCDSiteSubset site_subset)
{
    int offsetnb;
    int offseti;
    offseti = (site_subset == KUQCD_EVEN) ? 0 : latBase.latSubBulkCB.V;
    if (this->diracType == KUQCD_HISQ_DIRAC) {
        offsetnb = (site_subset == KUQCD_EVEN) ? latBase.latSubBulkCB.V : 0;
    } else if (this->diracType == KUQCD_HISQPC_DIRAC) {
        offsetnb = 0;
        // rootLogger.Warn("hisq pc dirac is under development!");
    } else {
        rootLogger.Error("dirac type should be hisq or hisq pc!");
    }

#pragma omp parallel for
    for (size_t i = 0; i < latBase.latSubBulkCB.V; ++i) {
        size_t j = i + offseti;
        // sub bulk CB position
        Position pos = GetPosition(i, latBase.latSubBulkCB);
        //        Rank0Printf("%5zu sub cb pos=%2d %2d %2d %2d ", i, pos.x, pos.y, pos.z, pos.t);
        // sup bulk CB position
        pos = GetSupPosFromSubPos(pos, latBase.latSubBulkCB, gridDim, myRank);
        // sup bulk FL position
        pos = GetPosFromCBPos(pos, site_subset);
        for (int dir = 0; dir < 4; ++dir) {
            Position dx(0, 0, 0, 0);
            int offsethalo;
            if (this->diracType == KUQCD_HISQ_DIRAC) {
                offsethalo = (site_subset == KUQCD_ODD) ? 0 : surf[dir] / 2;
            } else {
                offsethalo = 0;
            }
            // forward neigbhor
            dx[dir] = 1;
            Position myPos = pos + dx;
            // sup bulk FL position
            myPos.Normalize(dir, latBase.latSupBulkFL);
            //            Rank0Printf("dir=%1d sup fl pos=%2d %2d %2d %2d ", dir, myPos.x, myPos.y, myPos.z, myPos.t);
            // sup bulk CB position
            myPos = GetCBPos(myPos);
            int rank = GetRankFromSupPos(myPos, latBase.latSubBulkCB, gridDim);
            // sub bulk CB position
            myPos = GetSubPosFromSupPos(myPos, latBase.latSubBulkCB);
            int subIdx;
            if (rank == myRank) {
                subIdx = myPos.Index(latBase.latSubBulkCB);
                nb1Idx.idx[dir + (j << 3)] = subIdx + offsetnb;
                //                Rank0Printf("idx=%5d offsetnb=%5d %zu, %ld\n", subIdx, offsetnb, dir + (i<<3),
                //                nb1Idx.idx[dir + (i<<3)]);
            } else {
                subIdx = GetBufIdx1(myPos, dir, latBase.latSubBulkCB);
                nb1Idx.idx[dir + (j << 3)] = -subIdx - offsethalo - 1;
                edge1Idx[dir].idx[subIdx + offsethalo] = i;
            }
            // backward neigbhor
            myPos = pos - dx;
            // sup bulk FL position
            myPos.Normalize(dir, latBase.latSupBulkFL);
            //  sup bulk CB position
            myPos = GetCBPos(myPos);
            rank = GetRankFromSupPos(myPos, latBase.latSubBulkCB, gridDim);
            // sub bulk CB position
            myPos = GetSubPosFromSupPos(myPos, latBase.latSubBulkCB);
            if (rank == myRank) {
                subIdx = myPos.Index(latBase.latSubBulkCB);
                nb1Idx.idx[dir + 4 + (j << 3)] = subIdx + offsetnb;
            } else {
                subIdx = GetBufIdx1(myPos, dir, latBase.latSubBulkCB);
                nb1Idx.idx[dir + 4 + (j << 3)] = -subIdx - offsethalo - 1;
                edge1Idx[dir + 4].idx[subIdx + offsethalo] = i;
            }

            // forward next next nearest neigbhor
            offsethalo *= 3;
            dx[dir] = 3;
            myPos = pos + dx;
            // sup bulk FL position
            myPos.Normalize(dir, latBase.latSupBulkFL);
            // sup bulk CB position
            myPos = GetCBPos(myPos);
            rank = GetRankFromSupPos(myPos, latBase.latSubBulkCB, gridDim);
            // sub bulk CB position
            myPos = GetSubPosFromSupPos(myPos, latBase.latSubBulkCB);
            if (rank == myRank) {
                subIdx = myPos.Index(latBase.latSubBulkCB);
                nb3Idx.idx[dir + (j << 3)] = subIdx + offsetnb;
            } else {
                subIdx = GetBufIdx3(myPos, dir, latBase.latSubBulkCB);
                nb3Idx.idx[dir + (j << 3)] = -subIdx - offsethalo - 1;
                edge3Idx[dir].idx[subIdx + offsethalo] = i;
            }
            // backward next next nearest neigbhor
            myPos = pos - dx;
            // sup bulk FL position
            myPos.Normalize(dir, latBase.latSupBulkFL);
            // sup bulk CB position
            myPos = GetCBPos(myPos);
            rank = GetRankFromSupPos(myPos, latBase.latSubBulkCB, gridDim);
            // sub bulk CB position
            myPos = GetSubPosFromSupPos(myPos, latBase.latSubBulkCB);
            if (rank == myRank) {
                subIdx = myPos.Index(latBase.latSubBulkCB);
                nb3Idx.idx[dir + 4 + (j << 3)] = subIdx + offsetnb;
            } else {
                subIdx = GetBufIdx3(myPos, dir + 4, latBase.latSubBulkCB);
                nb3Idx.idx[dir + 4 + (j << 3)] = -subIdx - offsethalo - 1;
                edge3Idx[dir + 4].idx[subIdx + offsethalo] = i;
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::PackHalo(const VecT &vec, int dir)
{
    assert(comm[dir].dataSend != nullptr);
    int d = (dir < 4) ? dir : (dir - 4);
    size_t offseti, offsetg;
    if (this->diracType == KUQCD_HISQ_DIRAC) {
        offseti = latBase.latSubBulkCB.V;
        offsetg = 0;
    } else if (this->diracType == KUQCD_HISQPC_DIRAC) {
        offseti = 0;
        if (vec.latBase.site_subset == KUQCD_EVEN) {
            offsetg = 0;
        } else if (vec.latBase.site_subset == KUQCD_ODD) {
            offsetg = latBase.latSubBulkCB.V;
        } else {
            rootLogger.Error("site_subset of spinor for hisq pc dirac should be even or odd!");
        }
    } else {
        rootLogger.Error("dirac type should be hisq or hisq pc!");
    }
#pragma omp parallel for
    for (size_t i = 0; i < surf[d]; ++i) {
        size_t offset = (i < surf[d] / 2) ? 0 : offseti;
        size_t ni = edge1Idx[dir].idx[i] + offset;
        if (dir < 4) {
            SpinorT &v = vec.data[vec.latBase.myIdx.idx[ni]];
            GaugeST &U = gaugeSmeared.data[d + 4 * gaugeSmeared.latBase.myIdx.idx[ni + offsetg]];
            SpinorT &res = comm[dir].dataSend[i];
            UDaggMultV(U, v, res);
        } else {
            comm[dir].dataSend[i] = vec.data[vec.latBase.myIdx.idx[ni]];
        }

        for (int depth = 0; depth < 3; ++depth) {
            size_t l = depth + 3 * i;
            ni = edge3Idx[dir].idx[l] + offset;
            if (dir < 4) {
                SpinorT &v = vec.data[vec.latBase.myIdx.idx[ni]];
                GaugeNT &U = gaugeNaik.data[d + 4 * gaugeNaik.latBase.myIdx.idx[ni + offsetg]];
                SpinorT &res = comm3[dir].dataSend[l];
                UDaggMultV(U, v, res);
            } else {
                comm3[dir].dataSend[l] = vec.data[vec.latBase.myIdx.idx[ni]];
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::UnpackHalo(const VecT &vec, int dir)
{
    assert(comm[dir].dataSend != nullptr);
    int d = (dir < 4) ? dir : (dir - 4);
    size_t offseti, offsetg;
    if (this->diracType == KUQCD_HISQ_DIRAC) {
        offseti = latBase.latSubBulkCB.V;
        offsetg = 0;
    } else if (this->diracType == KUQCD_HISQPC_DIRAC) {
        offseti = 0;
        if (vec.latBase.site_subset == KUQCD_EVEN) {
            offsetg = 0;
        } else if (vec.latBase.site_subset == KUQCD_ODD) {
            offsetg = latBase.latSubBulkCB.V;
        } else {
            rootLogger.Error("site_subset of spinor for hisq pc dirac should be even or odd!");
        }
    } else {
        rootLogger.Error("dirac type should be hisq or hisq pc!");
    }

#ifdef USE_KUPL_COMM
    SpinorT *ghost1;
    SpinorT *ghost3;
    void *remote1;
    void *remote3;
    int remoteCommId = (dir < 4) ? (dir + 4) : (dir - 4);
    kupl_shm_win_query(this->comm[remoteCommId].memSend.win, this->comm[dir].nbRank, &remote1);
    kupl_shm_win_query(this->comm3[remoteCommId].memSend.win, this->comm[dir].nbRank, &remote3);
    ghost1 = static_cast<SpinorT *>(remote1);
    ghost3 = static_cast<SpinorT *>(remote3);
#endif

#pragma omp parallel for
    for (size_t i = 0; i < surf[d]; ++i) {
        size_t offset = (i < surf[d] / 2) ? offseti : 0;
        size_t ni = edge1Idx[dir].idx[i] + offset;
#ifdef USE_KUPL_COMM
        SpinorT *ghost = ghost1;
#else
        SpinorT *ghost = comm[dir].dataRecv;
#endif
        if (dir < 4) {
            SpinorT &v = ghost[i];
            GaugeST &U = gaugeSmeared.data[d + 4 * gaugeSmeared.latBase.myIdx.idx[ni + offsetg]];
            SpinorT &res = vec.data[vec.latBase.myIdx.idx[ni]];
            UMultV(U, v, res, static_cast<T>(0.5));
        } else {
            vec.data[vec.latBase.myIdx.idx[ni]] -= ghost[i] * static_cast<T>(0.5);
        }

        for (int depth = 0; depth < 3; ++depth) {
            size_t i3 = 3 * i;
            size_t l = depth + i3;
            size_t k = 2 - depth + i3;
            ni = edge3Idx[dir].idx[k] + offset;
#ifdef USE_KUPL_COMM
            SpinorT *ghost = ghost3;
#else
            ghost = comm3[dir].dataRecv;
#endif
            if (dir < 4) {
                SpinorT &v = ghost[l];
                GaugeNT &U = gaugeNaik.data[d + 4 * gaugeNaik.latBase.myIdx.idx[ni + offsetg]];
                SpinorT &res = vec.data[vec.latBase.myIdx.idx[ni]];
                UMultV(U, v, res, c3000);
            } else {
                vec.data[vec.latBase.myIdx.idx[ni]] -= ghost[l] * static_cast<T>(c3000);
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::StartComm(const VecT &vec, bool update)
{
    for (int dir = 7; dir >= 0; --dir) {
        if (gridDim[(dir < 4) ? dir : (dir - 4)] > 1) {
#ifndef USE_KUPL_COMM
            if (update) {
                comm[dir].Recv();
                comm3[dir].Recv();
            }
#endif
            PackHalo(vec, dir);
#ifndef USE_KUPL_COMM
            if (update) {
                comm[dir].Send();
                comm3[dir].Send();
            }
#endif
        }
    }
#ifdef USE_KUPL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::EndComm(const VecT &vec, bool update)
{
    for (int dir = 0; dir < 8; ++dir) {
        if (gridDim[(dir < 4) ? dir : (dir - 4)] > 1) {
#ifndef USE_KUPL_COMM
            if (update) {
                comm[dir].Finish(true);
                comm3[dir].Finish(true);
            }
#endif
            UnpackHalo(vec, dir);
#ifndef USE_KUPL_COMM
            if (update) {
                comm[dir].Finish(false);
                comm3[dir].Finish(false);
            }
#endif
        }
    }
#ifdef USE_KUPL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::Dslash(const VecT &src, VecT &dst, bool update)
{
    if (np > 1) {
        StartComm(src, update);
    }
    bool p = dst.latBase.site_subset == KUQCD_EVEN;
    size_t offseti;
    size_t offsetg;
    if (this->diracType == KUQCD_HISQPC_DIRAC) {
        offseti = p ? 0 : latBase.latSubBulkCB.V;
        offsetg = p ? latBase.latSubBulkCB.V : 0;
    } else {
        offseti = 0;
        offsetg = 0;
    }

#pragma omp parallel for
    for (size_t i = 0; i < dst.latBase.lenBulk; ++i) {
        size_t j = i + offseti;

        bool isdsteven =
            (src.latBase.site_subset == KUQCD_ODD) || ((src.latBase.site_subset == KUQCD_ALL) && (j < latBase.latSubBulkCB.V));
        size_t xyzt = (j < latBase.latSubBulkCB.V) ? j : (j - latBase.latSubBulkCB.V);
        size_t yzt = xyzt / latBase.latSubBulkCB.nx;
        size_t zt = yzt / latBase.latSubBulkCB.ny;
        size_t y = yzt - zt * latBase.latSubBulkCB.ny;
        size_t t = zt / latBase.latSubBulkCB.nz;
        size_t z = zt - t * latBase.latSubBulkCB.nz;
        bool ismatch = isdsteven == (!((y + z + t) & 0x1));

        SpinorT &res = dst.data[dst.latBase.myIdx.idx[i]];
        res = SpinorT();
        FwdDslash<0>(j, src, res, offsetg, ismatch);
        FwdDslash<1>(j, src, res, offsetg);
        FwdDslash<2>(j, src, res, offsetg);
        FwdDslash<3>(j, src, res, offsetg);
        BwdDslash<0>(j, src, res, offsetg, ismatch);
        BwdDslash<1>(j, src, res, offsetg);
        BwdDslash<2>(j, src, res, offsetg);
        BwdDslash<3>(j, src, res, offsetg);

    }

    if (np > 1) {
        EndComm(dst, update);
    }
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::MdaggM(const VecT &src, VecT &dst, bool update)
{
    Dslash(src, *tmp, update);
    Dslash(*tmp, dst, update);
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
template <int dir>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::FwdDslash(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    long ni = nb1Idx.idx[dir + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[i];
        res += UMultV(gaugeSmeared.data[dir + 4 * igg], src.data[isp] * static_cast<T>(0.5));
    }

    ni = nb3Idx.idx[dir + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[i];
        res += UMultV(gaugeNaik.data[dir + 4 * igg], src.data[isp] * c3000);
    }
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
template <int dir>
void DiracHISQ<T, Nbatch, ReconS, ReconN>::BwdDslash(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    long ni = nb1Idx.idx[dir + 4 + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[ni + offsetg];
        res -= UDaggMultV(gaugeSmeared.data[dir + 4 * igg], src.data[isp] * static_cast<T>(0.5));
    }

    ni = nb3Idx.idx[dir + 4 + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        res -= UDaggMultV(gaugeNaik.data[dir + 4 * igg], src.data[isp] * c3000);
    }
}

#ifndef SVEBITS

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::FwdDslash<0>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int Nbatch = 16;

    size_t igg = gaugeSmeared.latBase.myIdx.idx[i];
    if (ismatch) {
        long ni = nb1Idx.idx[8 * i];
        size_t isp = src.latBase.myIdx.idx[ni];
        SpinorT &v1 = src.data[isp];
        res += UMultV(gaugeSmeared.data[4 * igg], v1 * static_cast<float>(0.5));

        ni = nb3Idx.idx[8 * i];
        isp = src.latBase.myIdx.idx[ni];
        SpinorT &v3 = src.data[isp];
        SpinorT v;
        for (int nc = 0; nc < 3; ++nc) {
            for (int j = 0; j < Nbatch - 1; ++j) {
                v.c[nc].z[j] = v1.c[nc].z[j + 1];
                v.c[nc].z[j + Nbatch] = v1.c[nc].z[j + 1 + Nbatch];
            }
            v.c[nc].z[Nbatch - 1] = v3.c[nc].z[0];
            v.c[nc].z[Nbatch - 1 + Nbatch] = v3.c[nc].z[Nbatch];
        }
        res += UMultV(gaugeNaik.data[4 * igg], v * c3000);
    } else {
        long fni = nb1Idx.idx[8 * i];
        long bni = nb1Idx.idx[4 + 8 * i];
        SpinorT &fv = src.data[src.latBase.myIdx.idx[fni]];
        SpinorT &bv = src.data[src.latBase.myIdx.idx[bni]];
        SpinorT v;
        for (int nc = 0; nc < 3; ++nc) {
            for (int j = 0; j < Nbatch - 1; ++j) {
                v.c[nc].z[j] = bv.c[nc].z[j + 1];
                v.c[nc].z[j + Nbatch] = bv.c[nc].z[j + 1 + Nbatch];
            }
            v.c[nc].z[Nbatch - 1] = fv.c[nc].z[0];
            v.c[nc].z[Nbatch - 1 + Nbatch] = fv.c[nc].z[Nbatch];
        }
        res += UMultV(gaugeSmeared.data[4 * igg], v * static_cast<float>(0.5));

        for (int nc = 0; nc < 3; ++nc) {
            for (int j = 0; j < Nbatch - 2; ++j) {
                v.c[nc].z[j] = bv.c[nc].z[j + 2];
                v.c[nc].z[j + Nbatch] = bv.c[nc].z[j + 2 + Nbatch];
            }
            v.c[nc].z[Nbatch - 2] = fv.c[nc].z[0];
            v.c[nc].z[Nbatch - 1] = fv.c[nc].z[1];
            v.c[nc].z[Nbatch - 2 + Nbatch] = fv.c[nc].z[Nbatch];
            v.c[nc].z[Nbatch - 1 + Nbatch] = fv.c[nc].z[Nbatch + 1];
        }
        res += UMultV(gaugeNaik.data[4 * igg], v * c3000);
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::BwdDslash<0>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int Nbatch = 16;

    if (ismatch) {
        long fni = nb1Idx.idx[8 * i];
        long bni = nb1Idx.idx[4 + 8 * i];
        size_t fisp = src.latBase.myIdx.idx[fni];
        size_t bisp = src.latBase.myIdx.idx[bni];
        size_t figg = gaugeSmeared.latBase.myIdx.idx[fni + offsetg];
        size_t bigg = gaugeSmeared.latBase.myIdx.idx[bni + offsetg];
        SpinorT &fv = src.data[fisp];
        SpinorT &bv = src.data[bisp];
        GaugeST &fus = gaugeSmeared.data[4 * figg];
        GaugeST &bus = gaugeSmeared.data[4 * bigg];
        SpinorT v;
        GaugeST u;
        for (int nc = 0; nc < 3; ++nc) {
            v.c[nc].z[0] = bv.c[nc].z[Nbatch - 1];
            v.c[nc].z[Nbatch] = bv.c[nc].z[Nbatch - 1 + Nbatch];
            for (int j = 0; j < Nbatch - 1; ++j) {
                v.c[nc].z[j + 1] = fv.c[nc].z[j];
                v.c[nc].z[j + 1 + Nbatch] = fv.c[nc].z[j + Nbatch];
            }
        }
        for (int n = 0; n < 9; ++n) {
            u.m[n].z[0] = bus.m[n].z[Nbatch - 1];
            u.m[n].z[Nbatch] = bus.m[n].z[Nbatch - 1 + Nbatch];
            for (int j = 0; j < Nbatch - 1; ++j) {
                u.m[n].z[j + 1] = fus.m[n].z[j];
                u.m[n].z[j + 1 + Nbatch] = fus.m[n].z[j + Nbatch];
            }
        }
        res -= UDaggMultV(u, v * static_cast<float>(0.5));

        figg = gaugeNaik.latBase.myIdx.idx[fni + offsetg];
        bigg = gaugeNaik.latBase.myIdx.idx[bni + offsetg];
        GaugeST &fun = gaugeNaik.data[4 * figg];
        GaugeST &bun = gaugeNaik.data[4 * bigg];
        for (int nc = 0; nc < 3; ++nc) {
            v.c[nc].z[0] = bv.c[nc].z[Nbatch - 2];
            v.c[nc].z[1] = bv.c[nc].z[Nbatch - 1];
            v.c[nc].z[Nbatch] = bv.c[nc].z[Nbatch - 2 + Nbatch];
            v.c[nc].z[Nbatch + 1] = bv.c[nc].z[Nbatch - 1 + Nbatch];
            for (int j = 0; j < Nbatch - 2; ++j) {
                v.c[nc].z[j + 2] = fv.c[nc].z[j];
                v.c[nc].z[j + 2 + Nbatch] = fv.c[nc].z[j + Nbatch];
            }
        }
        for (int n = 0; n < 9; ++n) {
            u.m[n].z[0] = bun.m[n].z[Nbatch - 2];
            u.m[n].z[1] = bun.m[n].z[Nbatch - 1];
            u.m[n].z[Nbatch] = bun.m[n].z[Nbatch - 2 + Nbatch];
            u.m[n].z[Nbatch + 1] = bun.m[n].z[Nbatch - 1 + Nbatch];
            for (int j = 0; j < Nbatch - 2; ++j) {
                u.m[n].z[j + 2] = fun.m[n].z[j];
                u.m[n].z[j + 2 + Nbatch] = fun.m[n].z[j + Nbatch];
            }
        }
        res -= UDaggMultV(u, v * c3000);
    } else {
        long ni = nb1Idx.idx[4 + 8 * i];
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v1 = src.data[isp];
        GaugeST &u1s = gaugeSmeared.data[4 * igg];
        res -= UDaggMultV(u1s, v1 * static_cast<float>(0.5));

        igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        GaugeST &u1n = gaugeNaik.data[4 * igg];
        ni = nb3Idx.idx[4 + 8 * i];
        isp = src.latBase.myIdx.idx[ni];
        igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v3 = src.data[isp];
        GaugeST &u3 = gaugeNaik.data[4 * igg];
        SpinorT v;
        GaugeST u;
        for (int nc = 0; nc < 3; ++nc) {
            v.c[nc].z[0] = v3.c[nc].z[Nbatch - 1];
            v.c[nc].z[Nbatch] = v3.c[nc].z[Nbatch - 1 + Nbatch];
            for (int j = 0; j < Nbatch - 1; ++j) {
                v.c[nc].z[j + 1] = v1.c[nc].z[j];
                v.c[nc].z[j + 1 + Nbatch] = v1.c[nc].z[j + Nbatch];
            }
        }
        for (int n = 0; n < 9; ++n) {
            u.m[n].z[0] = u3.m[n].z[Nbatch - 1];
            u.m[n].z[Nbatch] = u3.m[n].z[Nbatch - 1 + Nbatch];
            for (int j = 0; j < Nbatch - 1; ++j) {
                u.m[n].z[j + 1] = u1n.m[n].z[j];
                u.m[n].z[j + 1 + Nbatch] = u1n.m[n].z[j + Nbatch];
            }
        }
        res -= UDaggMultV(u, v * c3000);
    }
}

#else

#if SVEBITS == 512

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::FwdDslash<0>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int Nbatch = 16;

    svbool_t pg = svptrue_b32();
    svbool_t pgSp;
    svfloat32_t src1;
    svfloat32_t src2;
    svfloat32_t tmpsrc;
    svfloat32_t tmp1;
    svfloat32_t tmp2;
    svfloat32_t res1;
    svfloat32_t res2;
    svfloat32_t u1;
    svfloat32_t u2;
    svfloat32_t tmpc;

    size_t igg = gaugeSmeared.latBase.myIdx.idx[i];
    if (ismatch) {
        long ni = nb1Idx.idx[8 * i];
        size_t isp = src.latBase.myIdx.idx[ni];
        SpinorT &v1 = src.data[isp];
        GaugeST &Us = gaugeSmeared.data[4 * igg];
        tmpc = svdup_n_f32(0.5);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &v1.c[n].z[0]);
                src2 = svld1_f32(pg, &v1.c[n].z[Nbatch]);
                u1 = svld1_f32(pg, &Us.m[n + 3 * m].z[0]);
                u2 = svld1_f32(pg, &Us.m[n + 3 * m].z[Nbatch]);
                tmp1 = svmla_f32_m(pg, svmls_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmla_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmla_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmla_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }

        ni = nb3Idx.idx[8 * i];
        isp = src.latBase.myIdx.idx[ni];
        SpinorT &v3 = src.data[isp];
        GaugeNT &Un = gaugeNaik.data[4 * igg];
        pgSp = svwhilegt_b32_s32(16, 1);
        tmpc = svdup_n_f32(c3000);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &v1.c[n].z[0]);
                tmpsrc = svld1_f32(pg, &v3.c[n].z[0]);
                src1 = svsplice_f32(pgSp, src1, tmpsrc);
                src2 = svld1_f32(pg, &v1.c[n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &v3.c[n].z[Nbatch]);
                src2 = svsplice_f32(pgSp, src2, tmpsrc);
                u1 = svld1_f32(pg, &Un.m[n + 3 * m].z[0]);
                u2 = svld1_f32(pg, &Un.m[n + 3 * m].z[Nbatch]);
                tmp1 = svmla_f32_m(pg, svmls_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmla_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmla_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmla_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }
    } else {
        long fni = nb1Idx.idx[8 * i];
        long bni = nb1Idx.idx[4 + 8 * i];
        SpinorT &fv = src.data[src.latBase.myIdx.idx[fni]];
        SpinorT &bv = src.data[src.latBase.myIdx.idx[bni]];
        GaugeST &Us = gaugeSmeared.data[4 * igg];
        pgSp = svwhilegt_b32_s32(16, 1);
        tmpc = svdup_n_f32(0.5);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &bv.c[n].z[0]);
                tmpsrc = svld1_f32(pg, &fv.c[n].z[0]);
                src1 = svsplice_f32(pgSp, src1, tmpsrc);
                src2 = svld1_f32(pg, &bv.c[n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &fv.c[n].z[Nbatch]);
                src2 = svsplice_f32(pgSp, src2, tmpsrc);
                u1 = svld1_f32(pg, &Us.m[n + 3 * m].z[0]);
                u2 = svld1_f32(pg, &Us.m[n + 3 * m].z[Nbatch]);
                tmp1 = svmla_f32_m(pg, svmls_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmla_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmla_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmla_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }

        GaugeNT &Un = gaugeNaik.data[4 * igg];
        pgSp = svwhilegt_b32_s32(16, 2);
        tmpc = svdup_n_f32(c3000);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &bv.c[n].z[0]);
                tmpsrc = svld1_f32(pg, &fv.c[n].z[0]);
                src1 = svsplice_f32(pgSp, src1, tmpsrc);
                src2 = svld1_f32(pg, &bv.c[n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &fv.c[n].z[Nbatch]);
                src2 = svsplice_f32(pgSp, src2, tmpsrc);
                u1 = svld1_f32(pg, &Un.m[n + 3 * m].z[0]);
                u2 = svld1_f32(pg, &Un.m[n + 3 * m].z[Nbatch]);
                tmp1 = svmla_f32_m(pg, svmls_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmla_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmla_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmla_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::BwdDslash<0>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int Nbatch = 16;

    svbool_t pg = svptrue_b32();
    svbool_t pgSp;
    svfloat32_t src1;
    svfloat32_t src2;
    svfloat32_t tmpsrc;
    svfloat32_t tmp1;
    svfloat32_t tmp2;
    svfloat32_t res1;
    svfloat32_t res2;
    svfloat32_t u1;
    svfloat32_t u2;
    svfloat32_t tmpc;

    if (ismatch) {
        long fni = nb1Idx.idx[8 * i];
        long bni = nb1Idx.idx[4 + 8 * i];
        size_t fisp = src.latBase.myIdx.idx[fni];
        size_t bisp = src.latBase.myIdx.idx[bni];
        size_t figg = gaugeSmeared.latBase.myIdx.idx[fni + offsetg];
        size_t bigg = gaugeSmeared.latBase.myIdx.idx[bni + offsetg];
        SpinorT &fv = src.data[fisp];
        SpinorT &bv = src.data[bisp];
        GaugeST &fus = gaugeSmeared.data[4 * figg];
        GaugeST &bus = gaugeSmeared.data[4 * bigg];

        pgSp = svwhilegt_b32_s32(16, 15);
        tmpc = svdup_n_f32(0.5);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &fv.c[n].z[0]);
                tmpsrc = svld1_f32(pg, &bv.c[n].z[0]);
                src1 = svsplice_f32(pgSp, tmpsrc, src1);
                src2 = svld1_f32(pg, &fv.c[n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &bv.c[n].z[Nbatch]);
                src2 = svsplice_f32(pgSp, tmpsrc, src2);
                u1 = svld1_f32(pg, &fus.m[m + 3 * n].z[0]);
                tmpsrc = svld1_f32(pg, &bus.m[m + 3 * n].z[0]);
                u1 = svsplice_f32(pgSp, tmpsrc, u1);
                u2 = svld1_f32(pg, &fus.m[m + 3 * n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &bus.m[m + 3 * n].z[Nbatch]);
                u2 = svsplice_f32(pgSp, tmpsrc, u2);
                tmp1 = svmla_f32_m(pg, svmla_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmls_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmls_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmls_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }

        figg = gaugeNaik.latBase.myIdx.idx[fni + offsetg];
        bigg = gaugeNaik.latBase.myIdx.idx[bni + offsetg];
        GaugeST &fun = gaugeNaik.data[4 * figg];
        GaugeST &bun = gaugeNaik.data[4 * bigg];

        pgSp = svwhilegt_b32_s32(16, 14);
        tmpc = svdup_n_f32(c3000);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &fv.c[n].z[0]);
                tmpsrc = svld1_f32(pg, &bv.c[n].z[0]);
                src1 = svsplice_f32(pgSp, tmpsrc, src1);
                src2 = svld1_f32(pg, &fv.c[n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &bv.c[n].z[Nbatch]);
                src2 = svsplice_f32(pgSp, tmpsrc, src2);
                u1 = svld1_f32(pg, &fun.m[m + 3 * n].z[0]);
                tmpsrc = svld1_f32(pg, &bun.m[m + 3 * n].z[0]);
                u1 = svsplice_f32(pgSp, tmpsrc, u1);
                u2 = svld1_f32(pg, &fun.m[m + 3 * n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &bun.m[m + 3 * n].z[Nbatch]);
                u2 = svsplice_f32(pgSp, tmpsrc, u2);
                tmp1 = svmla_f32_m(pg, svmla_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmls_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmls_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmls_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }
    } else {
        long ni = nb1Idx.idx[4 + 8 * i];
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v1 = src.data[isp];
        GaugeST &u1s = gaugeSmeared.data[4 * igg];
        tmpc = svdup_n_f32(0.5);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &v1.c[n].z[0]);
                src2 = svld1_f32(pg, &v1.c[n].z[Nbatch]);
                u1 = svld1_f32(pg, &u1s.m[m + 3 * n].z[0]);
                u2 = svld1_f32(pg, &u1s.m[m + 3 * n].z[Nbatch]);
                tmp1 = svmla_f32_m(pg, svmla_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmls_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmls_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmls_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }

        igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        GaugeST &u1n = gaugeNaik.data[4 * igg];
        ni = nb3Idx.idx[4 + 8 * i];
        isp = src.latBase.myIdx.idx[ni];
        igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v3 = src.data[isp];
        GaugeST &u3 = gaugeNaik.data[4 * igg];
        pgSp = svwhilegt_b32_s32(16, 15);
        tmpc = svdup_n_f32(c3000);
#pragma unroll
        for (int m = 0; m < 3; ++m) {
            tmp1 = svdup_n_f32(0);
            tmp2 = svdup_n_f32(0);
#pragma unroll
            for (int n = 0; n < 3; ++n) {
                src1 = svld1_f32(pg, &v1.c[n].z[0]);
                tmpsrc = svld1_f32(pg, &v3.c[n].z[0]);
                src1 = svsplice_f32(pgSp, tmpsrc, src1);
                src2 = svld1_f32(pg, &v1.c[n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &v3.c[n].z[Nbatch]);
                src2 = svsplice_f32(pgSp, tmpsrc, src2);
                u1 = svld1_f32(pg, &u1n.m[m + 3 * n].z[0]);
                tmpsrc = svld1_f32(pg, &u3.m[m + 3 * n].z[0]);
                u1 = svsplice_f32(pgSp, tmpsrc, u1);
                u2 = svld1_f32(pg, &u1n.m[m + 3 * n].z[Nbatch]);
                tmpsrc = svld1_f32(pg, &u3.m[m + 3 * n].z[Nbatch]);
                u2 = svsplice_f32(pgSp, tmpsrc, u2);
                tmp1 = svmla_f32_m(pg, svmla_f32_m(pg, tmp1, u2, src2), u1, src1);
                tmp2 = svmla_f32_m(pg, svmls_f32_m(pg, tmp2, u2, src1), u1, src2);
            }
            res1 = svld1_f32(pg, &res.c[m].z[0]);
            res2 = svld1_f32(pg, &res.c[m].z[Nbatch]);
            res1 = svmls_f32_m(pg, res1, tmp1, tmpc);
            res2 = svmls_f32_m(pg, res2, tmp2, tmpc);
            svst1_f32(pg, &res.c[m].z[0], res1);
            svst1_f32(pg, &res.c[m].z[Nbatch], res2);
        }
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::FwdDslash<1>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int dir = 1;
    long ni = nb1Idx.idx[dir + 8 * i];

    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[i];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeSmeared.data[dir + 4 * igg];
        UMultV(U, v, res, 0.5f);
    }

    ni = nb3Idx.idx[dir + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[i];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeNaik.data[dir + 4 * igg];
        UMultV(U, v, res, c3000);
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::BwdDslash<1>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int dir = 1;

    long ni = nb1Idx.idx[dir + 4 + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeSmeared.data[dir + 4 * igg];
        UDaggMultV(U, v, res, 0.5f);
    }

    ni = nb3Idx.idx[dir + 4 + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeNaik.data[dir + 4 * igg];
        UDaggMultV(U, v, res, c3000);
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::FwdDslash<2>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int dir = 2;
    long ni = nb1Idx.idx[dir + 8 * i];

    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[i];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeSmeared.data[dir + 4 * igg];
        UMultV(U, v, res, 0.5f);
    }

    ni = nb3Idx.idx[dir + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[i];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeNaik.data[dir + 4 * igg];
        UMultV(U, v, res, c3000);
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::BwdDslash<2>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int dir = 2;
    long ni = nb1Idx.idx[dir + 4 + 8 * i];

    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeSmeared.data[dir + 4 * igg];
        UDaggMultV(U, v, res, 0.5f);
    }

    ni = nb3Idx.idx[dir + 4 + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeNaik.data[dir + 4 * igg];
        UDaggMultV(U, v, res, c3000);
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::FwdDslash<3>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int dir = 3;
    long ni = nb1Idx.idx[dir + 8 * i];

    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[i];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeSmeared.data[dir + 4 * igg];
        UMultV(U, v, res, 0.5f);

    }

    ni = nb3Idx.idx[dir + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[i];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeNaik.data[dir + 4 * igg];
        UMultV(U, v, res, c3000);
    }
}

template <>
template <>
void DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>::BwdDslash<3>(
    size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch)
{
    int dir = 3;
    long ni = nb1Idx.idx[dir + 4 + 8 * i];


    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeSmeared.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeSmeared.data[dir + 4 * igg];
        UDaggMultV(U, v, res, 0.5f);
    }

    ni = nb3Idx.idx[dir + 4 + 8 * i];
    if (ni >= 0) {
        size_t isp = src.latBase.myIdx.idx[ni];
        size_t igg = gaugeNaik.latBase.myIdx.idx[ni + offsetg];
        SpinorT &v = src.data[isp];
        GaugeST &U = gaugeNaik.data[dir + 4 * igg];
        UDaggMultV(U, v, res, c3000);
    }
}

#endif

#endif

template class DiracHISQ<float, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;
template class DiracHISQ<double, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;
template class DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;
template class DiracHISQ<double, 8, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;

template class DiracHISQ<float, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;
template class DiracHISQ<double, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;
template class DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;
template class DiracHISQ<double, 8, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;

}  // namespace kuqcd
