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
#ifndef LATT_BASE_H
#define LATT_BASE_H

#include "space_time.h"
#include "kuqcd_memory.h"

namespace kuqcd {

class LattBase {
public:
    Lattice latSupBulkFL;  // supper bulk flat lattice
    Lattice latSupBulkCB;  // supper bulk checkerboard lattice
    Lattice latSubFullFL;  // sub full flat lattice
    Lattice latSubFullCB;  // sub full checkerboard lattice
    Lattice latSubBulkFL;  // sub bulk flat lattice
    Lattice latSubBulkCB;  // sub bulk checkerboard lattice
    Lattice *latFull;
    Lattice *latBulk;
    size_t lenFull;
    size_t lenBulk;
    std::array<unsigned short, 4> gridDim{1, 1, 1, 1};
    unsigned short haloDepth;
    unsigned int Nbatch;
    int myRank;
    int np;
    KuQCDPrecision prec;
    KuQCDGaugeOrder gauge_order;
    KuQCDRankOrder rankOrder;
    KuQCDSiteSubset site_subset;
    Index myIdx;

    LattBase() = default;
    LattBase(const LattBase &latBase);
    LattBase(const Lattice &lat, const std::array<unsigned short, 4> &gridDim_, unsigned short depth, KuQCDGaugeOrder order_,
        KuQCDRankOrder rankOrder_, KuQCDPrecision prec_, KuQCDSiteSubset parity_ = KUQCD_ALL);
    LattBase &operator=(const LattBase &latBase);
    void InitIdx();
};

}  // namespace kuqcd

#endif