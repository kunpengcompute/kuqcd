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
#ifndef SPACE_TIME_H
#define SPACE_TIME_H


#include <array>
#include <iostream>
#include "kuqcd_enum.h"
#include "kuqcd_utils.h"

namespace kuqcd {

class Lattice {
public:
    unsigned short nx = 0;
    unsigned short ny = 0;
    unsigned short nz = 0;
    unsigned short nt = 0;
    size_t V = 0;  // nx*ny*nz*nt

    Lattice() = default;
    Lattice(const Lattice &lat) : nx(lat.nx), ny(lat.ny), nz(lat.nz), nt(lat.nt), V(lat.V)
    {}
    Lattice(const unsigned short x, const unsigned short y, const unsigned short z, const unsigned short t)
        : nx(x), ny(y), nz(z), nt(t)
    {
        V = (size_t)nx * ny * nz * nt;
    }

    Lattice &operator=(const Lattice &lat);

    unsigned short &operator[](int dir);
    const unsigned short &operator[](int dir) const;
    size_t ComputeVolume();
    void Print() const;
};

std::ostream &operator<<(std::ostream &, const Lattice &);

class Position {
public:
    short x = 0;
    short y = 0;
    short z = 0;
    short t = 0;

    Position() = default;
    Position(const Position &pos) : x(pos.x), y(pos.y), z(pos.z), t(pos.t)
    {}
    Position(const short x_, const short y_, const short z_, const short t_) : x(x_), y(y_), z(z_), t(t_)
    {}
    Position &operator=(const Position &pos);

    short &operator[](int dir);
    const short &operator[](int dir) const;

    bool operator==(const Position &pos) const;
    bool operator!=(const Position &pos) const;

    Position operator+(const Position &pos) const;
    Position operator-(const Position &pos) const;
    Position &operator+=(const Position &pos);
    Position &operator-=(const Position &pos);

    void Normalize(int dir, const Lattice &lat);
    KuQCDSiteSubset Parity() const;
    unsigned int Index(const Lattice &lat) const;
    unsigned int IndexAnti(const Lattice &lat) const;
    void Print() const;
};

std::ostream &operator<<(std::ostream &, const Position &);

}  // namespace kuqcd

#endif