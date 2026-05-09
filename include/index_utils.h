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
#ifndef INDEX_UTILS_H
#define INDEX_UTILS_H

#include "space_time.h"

/*************************************
 * four level lattice:
 * supper lattice
 * checkerboard supper lattice
 * sub lattice
 * checkerboard sub lattice
 *************************************/

namespace kuqcd {

// On lattice get position from index
Position GetPosition(const unsigned int i, const Lattice &lat);

// On lattice get position from index, t-z-y-x order
Position GetAntiPosition(const unsigned int i, const Lattice &lat);

// On lattice get index from position
inline unsigned int GetIndex(const Position &pos, const Lattice &lat)
{
    return pos.Index(lat);
}

// get (CB) supper lattice position from (CB) sub lattice position
Position GetSupPosFromSubPos(const Position &pos, const Lattice &subLat, const std::array<unsigned short, 4> &gridDim,
                             const int myRank);

// get (CB) sub lattice position from (CB) supper lattice position
Position GetSubPosFromSupPos(const Position &pos, const Lattice &subLat);

// get supper (sub) lattice position from CB supper (sub) lattice position
Position GetPosFromCBPos(const Position &pos, const KuQCDSiteSubset site_subset);

// get CB supper (sub) lattice position from supper (sub) lattice position
Position GetCBPos(const Position &pos);

// get (CB) full lattice position from (CB) bulk lattice position
Position GetFullPosFromBulkPos(const Position &pos, int haloDepth, const std::array<unsigned short, 4> &gridDim);

// get supper lattice index from sub lattice index, do not consider site_subset
unsigned int GetSupIdxFromSubIdx(const unsigned int i, const Lattice &subLat, const Lattice &supLat,
                                 const std::array<unsigned short, 4> &gridDim, const int myRank);

// get sub lattice index from supper lattice index, do not consider site_subset
unsigned int GetSubIdxFromSupIdx(const unsigned int i, const Lattice &subLat, const Lattice &supLat);

// get supper (sub) lattice index from checkerboard supper (sub) index
unsigned int GetSupIdx(const unsigned int i, const KuQCDSiteSubset site_subset, const Lattice &latCB);

// get checkerboard supper (sub) index from supper (sub) lattice
inline unsigned int GetCBIdx(const unsigned int i)
{
    return i >> 0x1;
};

// get rank from supper lattice index
unsigned int GetRankFromSupIdx(const unsigned int i, const Lattice &subLat, const Lattice &supLat);

// get rank from supper lattice position
unsigned int GetRankFromSupPos(const Position &pos, const Lattice &subLat,
                               const std::array<unsigned short, 4> &gridDim);

// get buffer index from sub lattice position, nearest neighbor
unsigned int GetBufIdx1(const Position &pos, const int dir, const Lattice &lat);

// get buffer index from sub lattice position, next next nearest neighbor
unsigned int GetBufIdx3(const Position &pos, const int dir, const Lattice &lat);

} // namespace kuqcd

#endif
