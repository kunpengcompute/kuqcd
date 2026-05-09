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
#include <cstdio>
#include <mpi.h>
#include <iomanip>

#include "space_time.h"
#include "kuqcd_utils.h"

namespace kuqcd {

Lattice &Lattice::operator=(const Lattice &lat)
{
    nx = lat.nx;
    ny = lat.ny;
    nz = lat.nz;
    nt = lat.nt;
    V = lat.V;
    return *this;
}

unsigned short &Lattice::operator[](int dir)
{
    switch (dir) {
        case 0:
            return nx;
            break;
        case 1:
            return ny;
            break;
        case 2:
            return nz;
            break;
        case 3:
            return nt;
            break;
        default:
            perror("error direction!");
            return nx;
    }
}

const unsigned short &Lattice::operator[](int dir) const
{
    switch (dir) {
        case 0:
            return nx;
            break;
        case 1:
            return ny;
            break;
        case 2:
            return nz;
            break;
        case 3:
            return nt;
            break;
        default:
            perror("error direction!");
            return nx;
    }
}

size_t Lattice::ComputeVolume()
{
    V = nx * ny * nz * nt;
    return V;
}

void Lattice::Print() const
{
    Rank0Printf("%d %d %d %d %zu\n", nx, ny, nz, nt, V);
}

std::ostream &operator<<(std::ostream &os, const Lattice &lat)
{
    os << std::setw(3) << lat.nx;
    os << std::setw(3) << lat.ny;
    os << std::setw(3) << lat.nz;
    os << std::setw(3) << lat.nt;
    os << std::setw(9) << lat.V;
    return os;
}

Position &Position::operator=(const Position &pos)
{
    x = pos.x;
    y = pos.y;
    z = pos.z;
    t = pos.t;
    return *this;
}

short &Position::operator[](int dir)
{
    switch (dir) {
        case 0:
            return x;
            break;
        case 1:
            return y;
            break;
        case 2:
            return z;
            break;
        case 3:
            return t;
        default:
            perror("error direction!");
            return x;
    }
}

const short &Position::operator[](int dir) const
{
    switch (dir) {
        case 0:
            return x;
            break;
        case 1:
            return y;
            break;
        case 2:
            return z;
            break;
        case 3:
            return t;
        default:
            perror("error direction!");
            return x;
    }
}

bool Position::operator==(const Position &pos) const
{
    return ((x == pos.x) && (y == pos.y) && (z == pos.z) && (t == pos.t));
}

bool Position::operator!=(const Position &pos) const
{
    return ((x != pos.x) || (y != pos.y) || (z != pos.z) || (t != pos.t));
}

Position Position::operator+(const Position &pos) const
{
    Position res = *this;
    for (int i = 0; i < 4; ++i) {
        res[i] += pos[i];
    }
    return res;
}

Position Position::operator-(const Position &pos) const
{
    Position res = *this;
    for (int i = 0; i < 4; ++i) {
        res[i] -= pos[i];
    }
    return res;
}

Position &Position::operator+=(const Position &pos)
{
    for (int i = 0; i < 4; ++i) {
        (*this)[i] += pos[i];
    }
    return *this;
}
Position &Position::operator-=(const Position &pos)
{
    for (int i = 0; i < 4; ++i) {
        (*this)[i] -= pos[i];
    }
    return *this;
}

void Position::Normalize(int dir, const Lattice &lat)
{
    assert((dir >= 0) && (dir < 4));
    (*this)[dir] = (*this)[dir] % lat[dir];
    while ((*this)[dir] >= lat[dir]) {
        (*this)[dir] -= lat[dir];
    }
    while ((*this)[dir] < 0) {
        (*this)[dir] += lat[dir];
    }
}

KuQCDSiteSubset Position::Parity() const
{
    return (KuQCDSiteSubset)((x + y + z + t) & 1);
}

unsigned int Position::Index(const Lattice &lat) const
{
    assert((x < lat.nx) && (x >= 0));
    assert((y < lat.ny) && (y >= 0));
    assert((z < lat.nz) && (z >= 0));
    assert((t < lat.nt) && (t >= 0));
    return x + lat.nx * (y + lat.ny * (z + lat.nz * t));
}

unsigned int Position::IndexAnti(const Lattice &lat) const
{
    assert((x < lat.nx) && (x >= 0));
    assert((y < lat.ny) && (y >= 0));
    assert((z < lat.nz) && (z >= 0));
    assert((t < lat.nt) && (t >= 0));
    return t + lat.nt * (z + lat.nz * (y + lat.ny * x));
}

void Position::Print() const
{
    Rank0Printf("%d %d %d %d\n", x, y, z, t);
}

std::ostream &operator<<(std::ostream &os, const Position &pos)
{
    os << std::setw(3) << pos.x;
    os << std::setw(3) << pos.y;
    os << std::setw(3) << pos.z;
    os << std::setw(3) << pos.t;
    return os;
}

} // namespace kuqcd
