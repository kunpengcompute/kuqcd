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
#include "latt_base.h"
#include "omp.h"
#include "index_utils.h"

namespace kuqcd {

LattBase::LattBase(const LattBase &latBase)
{
    latSupBulkFL = latBase.latSupBulkFL;
    latSupBulkCB = latBase.latSupBulkCB;
    latSubFullFL = latBase.latSubFullFL;
    latSubFullCB = latBase.latSubFullCB;
    latSubBulkFL = latBase.latSubBulkFL;
    latSubBulkCB = latBase.latSubBulkCB;
    lenFull = latBase.lenFull;
    lenBulk = latBase.lenBulk;
    gridDim = latBase.gridDim;
    haloDepth = latBase.haloDepth;
    Nbatch = latBase.Nbatch;
    myRank = latBase.myRank;
    np = latBase.np;
    prec = latBase.prec;
    gauge_order = latBase.gauge_order;
    rankOrder = latBase.rankOrder;
    site_subset = latBase.site_subset;
    myIdx = latBase.myIdx;

    if (site_subset == KUQCD_ALL) {
        latFull = &latSubFullFL;
        latBulk = &latSubBulkFL;
    } else {
        latFull = &latSubFullCB;
        latBulk = &latSubBulkCB;
    }
}

LattBase::LattBase(const Lattice &lat, const std::array<unsigned short, 4> &gridDim_, unsigned short depth,
    KuQCDGaugeOrder order_, KuQCDRankOrder rankOrder_, KuQCDPrecision prec_, KuQCDSiteSubset parity_)
{
    gridDim = gridDim_;
    haloDepth = depth;
    gauge_order = order_;
    rankOrder = rankOrder_;
    prec = prec_;
    site_subset = parity_;
    myRank = GetRank();
    np = GetNp();
    bool pc = (site_subset == KUQCD_EVEN) || (site_subset == KUQCD_ODD);

    unsigned short n[4];
    unsigned short nFull[4];

    int tmp = 1;
    for (int i = 0; i < 4; ++i) {
        tmp *= gridDim[i];
    }
    if (tmp != np) {
        Rank0Fprintf(stderr, "%d != %d * %d * %d * %d!\n", np, gridDim[0], gridDim[1], gridDim[2], gridDim[3]);
    }

    if ((gauge_order == KUQCD_BATCH_ORDER) && (gridDim[0] > 1)) {
        Rank0Fprintf(
            stderr, "if using batch gauge_order the x direction should not communicate! %s %d\n", __FILE__, __LINE__);
    }

    if (gauge_order == KUQCD_BATCH_ORDER) {
        if (prec == KUQCD_DOUBLE) {
            Nbatch = 8;
        } else if (prec == KUQCD_SINGLE) {
            Nbatch = 16;
        } else {
            Rank0Fprintf(stderr, "error precision! %s %d\n", __FILE__, __LINE__);
        }
        if ((lat[0] % Nbatch) != 0) {
            Rank0Fprintf(
                stderr, "x direction lattice %d cannot divided by %d! %s %d\n", lat[0], Nbatch, __FILE__, __LINE__);
        }
    } else if (gauge_order == KUQCD_NORMAL_ORDER) {
        Nbatch = 1;
    } else {
        Rank0Fprintf(stderr, "error gauge_order! %s %d\n", __FILE__, __LINE__);
    }

    for (int i = 0; i < 4; i++) {
        if ((lat[i] % gridDim[i]) != 0) {
            Rank0Fprintf(stderr,
                "dir %d lattice %d cannot be splitted into %d part! %s %d\n",
                i,
                lat[i],
                gridDim[i],
                __FILE__,
                __LINE__);
        }
        n[i] = lat[i] / gridDim[i];
        n[i] = (i == 0) ? (n[i] / Nbatch) : n[i];

        if (gridDim[i] > 1) {
            if (n[i] < depth) {
                Rank0Fprintf(stderr,
                    "dir %d lattice %d is smaller than halo depth %d! %s %d\n",
                    i,
                    n[i],
                    depth,
                    __FILE__,
                    __LINE__);
            }
            nFull[i] = n[i] + 2 * depth;
        } else {
            nFull[i] = n[i];
        }
    }

    if (pc && (n[0] & 1)) {
        Rank0Fprintf(stderr, "x direction lattice %d should be even! %s %d\n", n[0], __FILE__, __LINE__);
    }

    latSupBulkFL = Lattice(lat.nx / Nbatch, lat.ny, lat.nz, lat.nt);
    latSupBulkCB = Lattice(lat.nx / Nbatch / 2, lat.ny, lat.nz, lat.nt);
    latSubFullFL = Lattice(nFull[0], nFull[1], nFull[2], nFull[3]);
    latSubFullCB = Lattice(nFull[0] / 2, nFull[1], nFull[2], nFull[3]);
    latSubBulkFL = Lattice(n[0], n[1], n[2], n[3]);
    latSubBulkCB = Lattice(n[0] / 2, n[1], n[2], n[3]);

    if (pc) {
        lenFull = latSubFullCB.V;
        lenBulk = latSubBulkCB.V;
        latFull = &latSubFullCB;
        latBulk = &latSubBulkCB;
    } else {
        lenFull = latSubFullFL.V;
        lenBulk = latSubBulkFL.V;
        latFull = &latSubFullFL;
        latBulk = &latSubBulkFL;
    }

    InitIdx();
}

LattBase &LattBase::operator=(const LattBase &latBase)
{
    latSupBulkFL = latBase.latSupBulkFL;
    latSupBulkCB = latBase.latSupBulkCB;
    latSubFullFL = latBase.latSubFullFL;
    latSubFullCB = latBase.latSubFullCB;
    latSubBulkFL = latBase.latSubBulkFL;
    latSubBulkCB = latBase.latSubBulkCB;
    lenFull = latBase.lenFull;
    lenBulk = latBase.lenBulk;
    gridDim = latBase.gridDim;
    haloDepth = latBase.haloDepth;
    Nbatch = latBase.Nbatch;
    prec = latBase.prec;
    gauge_order = latBase.gauge_order;
    rankOrder = latBase.rankOrder;
    site_subset = latBase.site_subset;
    myIdx = latBase.myIdx;

    if (site_subset == KUQCD_ALL) {
        latFull = &latSubFullFL;
        latBulk = &latSubBulkFL;
    } else {
        latFull = &latSubFullCB;
        latBulk = &latSubBulkCB;
    }

    return *this;
}

void LattBase::InitIdx()
{
    myIdx.Alloc(lenBulk);

#pragma omp parallel for
    for (size_t i = 0; i < latSubBulkCB.V; ++i) {
        // sub bulk CB position
        Position pos = GetPosition(i, latSubBulkCB);
        //        Rank0Printf("i=%5zu sub bulk cb pos=%2d %2d %2d %2d\n", i, pos.x, pos.y, pos.z, pos.t);
        if (site_subset == KUQCD_ALL) {
            /* even site_subset */
            // sub bulk FL position
            Position posP = GetPosFromCBPos(pos, KUQCD_EVEN);
            //            Rank0Printf("even sub bulk FL pos=%2d %2d %2d %2d ", posP.x, posP.y, posP.z, posP.t);
            // sub full FL position
            posP = GetFullPosFromBulkPos(posP, haloDepth, gridDim);
            //            Rank0Printf("sub full FL pos=%2d %2d %2d %2d ", posP.x, posP.y, posP.z, posP.t);
            // sub full FL index
            size_t idx = posP.Index(latSubFullFL);
            // sub full CB index
            idx = idx / 2;
            myIdx.idx[i] = idx;

            /* odd site_subset */
            // sub bulk FL position
            posP = GetPosFromCBPos(pos, KUQCD_ODD);
            //            Rank0Printf("odd sub bulk FL pos=%2d %2d %2d %2d ", posP.x, posP.y, posP.z, posP.t);
            // sub full FL position
            posP = GetFullPosFromBulkPos(posP, haloDepth, gridDim);
            //            Rank0Printf("sub full FL pos=%2d %2d %2d %2d ", posP.x, posP.y, posP.z, posP.t);
            // sub full FL index
            idx = posP.Index(latSubFullFL);
            // sub full CB index
            idx = idx / 2 + latSubFullFL.V / 2;
            myIdx.idx[i + latSubBulkCB.V] = idx;
        } else {
            // sub full CB position
            pos = GetFullPosFromBulkPos(pos, haloDepth, gridDim);
            size_t idx = pos.Index(latSubFullCB);
            myIdx.idx[i] = idx;
        }
    }
}

}  // namespace kuqcd