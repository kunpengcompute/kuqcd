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
#include <cstdio>

#include "index_utils.h"
namespace kuqcd {

Position GetPosition(const unsigned int i, const Lattice &lat)
{
    unsigned int x;
    unsigned int y;
    unsigned int z;
    unsigned int t;
    unsigned int offx;
    unsigned int offy;
    offx = i / lat.nx;
    offy = offx / lat.ny;
    x = i - offx * lat.nx;
    y = offx - offy * lat.ny;
    t = offy / lat.nz;
    z = offy - t * lat.nz;
    Position p(x, y, z, t);
    return p;
}

Position GetAntiPosition(const unsigned int i, const Lattice &lat)
{
    unsigned int x;
    unsigned int y;
    unsigned int z;
    unsigned int t;
    unsigned int offt;
    unsigned int offz;
    offt = i / lat.nt;
    offz = offt / lat.nz;
    t = i - offt * lat.nt;
    z = offt - offz * lat.nz;
    x = offz / lat.ny;
    y = offz - x * lat.ny;
    Position p(x, y, z, t);
    return p;
}

Position GetSupPosFromSubPos(
    const Position &pos, const Lattice &subLat, const std::array<unsigned short, 4> &gridDim, const int myRank)
{
    unsigned int gridx;
    unsigned int gridy;
    unsigned int gridz;
    unsigned int gridt;
    unsigned int offx;
    unsigned int offy;
    Position p;

    offx = myRank / gridDim[0];
    offy = offx / gridDim[1];
    gridx = myRank - offx * gridDim[0];
    gridy = offx - offy * gridDim[1];
    gridt = offy / gridDim[2];
    gridz = offy - gridt * gridDim[2];

    p[0] = pos[0] + gridx * subLat[0];
    p[1] = pos[1] + gridy * subLat[1];
    p[2] = pos[2] + gridz * subLat[2];
    p[3] = pos[3] + gridt * subLat[3];

    return p;
}

Position GetSubPosFromSupPos(const Position &pos, const Lattice &subLat)
{
    unsigned int grid;
    Position p;

    for (int i = 0; i < 4; ++i) {
        grid = pos[i] / subLat[i];
        p[i] = pos[i] - grid * subLat[i];
    }
    return p;
}

Position GetPosFromCBPos(const Position &pos, const KuQCDSiteSubset site_subset)
{
    Position p;
    p.x = (pos.x << 1) + ((pos.y + pos.z + pos.t + (int)site_subset) & 1);
    p.y = pos.y;
    p.z = pos.z;
    p.t = pos.t;

    return p;
}

Position GetCBPos(const Position &pos)
{
    Position p;
    p.x = pos.x >> 1;
    p.y = pos.y;
    p.z = pos.z;
    p.t = pos.t;

    return p;
}

Position GetFullPosFromBulkPos(const Position &pos, int haloDepth, const std::array<unsigned short, 4> &gridDim)
{
    Position p(pos);
    for (int dir = 0; dir < 4; ++dir) {
        if (gridDim[dir] > 1) {
            p[dir] += haloDepth;
        }
    }
    return p;
}

unsigned int GetSupIdxFromSubIdx(const unsigned int i, const Lattice &subLat, const Lattice &supLat,
    const std::array<unsigned short, 4> &gridDim, const int myRank)
{
    unsigned int x;
    unsigned int y;
    unsigned int z;
    unsigned int t;
    unsigned int gridx;
    unsigned int gridy;
    unsigned int gridz;
    unsigned int gridt;
    unsigned int offx;
    unsigned int offy;
    offx = i / subLat.nx;
    offy = offx / subLat.ny;
    x = i - offx * subLat.nx;
    y = offx - offy * subLat.ny;
    t = offy / subLat.nz;
    z = offy - t * subLat.nz;

    offx = myRank / gridDim[0];
    offy = offx / gridDim[1];
    gridx = myRank - offx * gridDim[0];
    gridy = offx - offy * gridDim[1];
    gridt = offy / gridDim[2];
    gridz = offy - gridt * gridDim[2];

    x += gridx * subLat.nx;
    y += gridy * subLat.ny;
    z += gridz * subLat.nz;
    t += gridt * subLat.nt;

    return x + supLat.nx * (y + supLat.ny * (z + supLat.nz * t));
}

unsigned int GetSubIdxFromSupIdx(const unsigned int i, const Lattice &subLat, const Lattice &supLat)
{
    unsigned int x;
    unsigned int y;
    unsigned int z;
    unsigned int t;
    unsigned int offx;
    unsigned int offy;
    offx = i / supLat.nx;
    offy = offx / supLat.ny;
    x = i - offx * supLat.nx;
    y = offx - offy * supLat.ny;
    t = offy / supLat.nz;
    z = offy - t * supLat.nz;

    offx = x / subLat.nx;
    x -= offx * subLat.nx;
    offx = y / subLat.ny;
    y -= offx * subLat.ny;
    offx = z / subLat.nz;
    z -= offx * subLat.nz;
    offx = t / subLat.nt;
    t -= offx * subLat.nt;

    return x + subLat.nx * (y + subLat.ny * (z + subLat.nz * t));
}

unsigned int GetSupIdx(const unsigned int i, const KuQCDSiteSubset site_subset, const Lattice &latCB)
{
    unsigned int y;
    unsigned int z;
    unsigned int t;
    unsigned int offx;
    unsigned int offy;
    offx = i / latCB.nx;
    offy = offx / latCB.ny;
    y = offx - offy * latCB.ny;
    t = offy / latCB.nz;
    z = offy - t * latCB.nz;

    return 2 * i + ((y + z + t + site_subset) & 1);
}

unsigned int GetRankFromSupIdx(const unsigned int i, const Lattice &subLat, const Lattice &supLat)
{
    unsigned int x;
    unsigned int y;
    unsigned int z;
    unsigned int t;
    unsigned int offx;
    unsigned int offy;
    unsigned int gridx;
    unsigned int gridy;
    unsigned int gridz;
    unsigned int gridt;
    offx = i / supLat.nx;
    offy = offx / supLat.ny;
    x = i - offx * supLat.nx;
    y = offx - offy * supLat.ny;
    t = offy / supLat.nz;
    z = offy - t * supLat.nz;

    gridx = x / subLat.nx;
    gridy = y / subLat.ny;
    gridz = z / subLat.nz;
    gridt = t / subLat.nt;

    return gridx + subLat.nx * (gridy + subLat.ny * (gridz + subLat.nz * gridt));
}

unsigned int GetRankFromSupPos(const Position &pos, const Lattice &subLat, const std::array<unsigned short, 4> &gridDim)
{
    unsigned int gridx;
    unsigned int gridy;
    unsigned int gridz;
    unsigned int gridt;
    gridx = pos[0] / subLat.nx;
    gridy = pos[1] / subLat.ny;
    gridz = pos[2] / subLat.nz;
    gridt = pos[3] / subLat.nt;

    return gridx + gridDim[0] * (gridy + gridDim[1] * (gridz + gridDim[2] * gridt));
}

unsigned int GetBufIdx1(const Position &pos, const int dir, const Lattice &lat)
{
    unsigned int idx = 0;
    switch (dir) {
        case 0:
            idx = pos.y + lat.ny * (pos.z + lat.nz * pos.t);
            break;
        case 1:
            idx = pos.x + lat.nx * (pos.z + lat.nz * pos.t);
            break;
        case 2:
            idx = pos.x + lat.nx * (pos.y + lat.ny * pos.t);
            break;
        case 3:
            idx = pos.x + lat.nx * (pos.y + lat.ny * pos.z);
            break;
        default:
            Rank0Fprintf(stderr, "error direction! %d %d %s\n", dir, __LINE__, __FILE__);
    }
    return idx;
}

unsigned int GetBufIdx3(const Position &pos, const int dir, const Lattice &lat)
{
    unsigned int idx = 0;
    short x = pos.x;
    short y = pos.y;
    short z = pos.z;
    short t = pos.t;
    short nx = lat.nx;
    short ny = lat.ny;
    short nz = lat.nz;
    short nt = lat.nt;
    short depth = -1;
    int d;
    if (dir < 4) {
        d = dir;
        if (pos[d] == 0) {
            depth = 0;
        } else if (pos[d] == 1) {
            depth = 1;
        } else if (pos[d] == 2) {
            depth = 2;
        } else {
            Rank0Fprintf(stderr, "error buffer depth %d! %s %d\n", depth, __FILE__, __LINE__);
            Rank0Fprintf(stderr, "pos = %d %d %d %d dir = %d\n", x, y, z, t, dir);
        }
    } else {
        d = dir - 4;
        if (pos[d] == lat[d] - 1) {
            depth = 0;
        } else if (pos[d] == lat[d] - 2) {
            depth = 1;
        } else if (pos[d] == lat[d] - 3) {
            depth = 2;
        } else {
            Rank0Fprintf(stderr, "error buffer depth %d! %s %d\n", depth, __FILE__, __LINE__);
            Rank0Fprintf(stderr, "pos = %d %d %d %d dir = %d\n", x, y, z, t, dir);
        }
    }

    switch (d) {
        case 0:
            x = depth;
            nx = 3;
            idx = x + static_cast<unsigned int>(nx) *
                          (y + static_cast<unsigned int>(ny) * (z + static_cast<unsigned int>(nz) * t));
            break;
        case 1:
            y = depth;
            ny = 3;
            idx = y + static_cast<unsigned int>(ny) *
                          (x + static_cast<unsigned int>(nx) * (z + static_cast<unsigned int>(nz) * t));
            break;
        case 2:
            z = depth;
            nz = 3;
            idx = z + static_cast<unsigned int>(nz) *
                          (x + static_cast<unsigned int>(nx) * (y + static_cast<unsigned int>(ny) * t));
            break;
        case 3:
            t = depth;
            nt = 3;
            idx = t + static_cast<unsigned int>(nt) *
                          (x + static_cast<unsigned int>(nx) * (y + static_cast<unsigned int>(ny) * z));
            break;
        default:
            Rank0Fprintf(stderr, "error buffer depth %d! %s %d\n", depth, __FILE__, __LINE__);
            Rank0Fprintf(stderr, "pos = %d %d %d %d dir = %d\n", x, y, z, t, dir);
            break;
    }
    return idx;
}

}  // namespace kuqcd