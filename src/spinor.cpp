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
#include <cmath>
#include <iomanip>
#include "spinor.h"

namespace kuqcd {

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch>::ColorSpinor(const ColorSpinor<T, Nbatch> &other)
{
    for (int i = 0; i < 3; ++i) {
        c[i] = other.c[i];
    }
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch>::ColorSpinor(const Complex<T, Nbatch> &other)
{
    for (int i = 0; i < 3; ++i) {
        c[i] = other;
    }
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch>::ColorSpinor(const Complex<T, Nbatch> &red, const Complex<T, Nbatch> &green,
                                    const Complex<T, Nbatch> &blue)
{
    c[0] = red;
    c[1] = green;
    c[2] = blue;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch>::ColorSpinor(const Complex<T, Nbatch> ComponentsArray[3])
{
    for (int i = 0; i < 3; ++i) {
        c[i] = ComponentsArray[i];
    }
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch>::ColorSpinor(T real, T imag)
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < Nbatch; ++j) {
            c[i].z[j] = real;
            c[i].z[j + Nbatch] = imag;
        }
    }
}

// Element access - safe version with bounds checking
template <typename T, int Nbatch>
Complex<T, Nbatch> &ColorSpinor<T, Nbatch>::ColorComponent(int i)
{
    if ((i >= 3) || (i < 0)) {
        throw std::out_of_range("ColorSpinor index out of range. Must be 0-2.");
    }
    return c[i];
}

template <typename T, int Nbatch>
const Complex<T, Nbatch> &ColorSpinor<T, Nbatch>::ColorComponent(int i) const
{
    if ((i >= 3) || (i < 0)) {
        throw std::out_of_range("ColorSpinor index out of range. Must be 0-2.");
    }
    return c[i];
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> &ColorSpinor<T, Nbatch>::operator=(const ColorSpinor<T, Nbatch> &other)
{
    for (int i = 0; i < 3; ++i) {
        c[i] = other.c[i];
    }
    return *this;
}

// Assignment operators
template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> &ColorSpinor<T, Nbatch>::operator+=(const ColorSpinor<T, Nbatch> &other)
{
    c[0] += other.c[0];
    c[1] += other.c[1];
    c[2] += other.c[2];
    return *this;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> &ColorSpinor<T, Nbatch>::operator-=(const ColorSpinor<T, Nbatch> &other)
{
    c[0] -= other.c[0];
    c[1] -= other.c[1];
    c[2] -= other.c[2];
    return *this;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> &ColorSpinor<T, Nbatch>::operator*=(const ColorSpinor<T, Nbatch> &other)
{
    // Component-wise multiplication
    c[0] *= other.c[0];
    c[1] *= other.c[1];
    c[2] *= other.c[2];
    return *this;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> &ColorSpinor<T, Nbatch>::operator*=(T scalar)
{
    c[0] *= scalar;
    c[1] *= scalar;
    c[2] *= scalar;
    return *this;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> &ColorSpinor<T, Nbatch>::operator/=(T scalar)
{
    c[0] /= scalar;
    c[1] /= scalar;
    c[2] /= scalar;
    return *this;
}

// Unary operators
template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> ColorSpinor<T, Nbatch>::operator-() const
{
    ColorSpinor<T, Nbatch> a(*this);
    a.c[0] = -a.c[0];
    a.c[1] = -a.c[1];
    a.c[2] = -a.c[2];
    return a;
}

template <typename T, int Nbatch>
void ColorSpinor<T, Nbatch>::Conj()
{
    c[0].Conj();
    c[1].Conj();
    c[2].Conj();
}

template <typename T, int Nbatch>
ColorSpinor<T, 1> ColorSpinor<T, Nbatch>::GetBatch(int i) const
{
    ColorSpinor<T, 1> a;
    a.c[0] = this->c[0].GetBatch(i);
    a.c[1] = this->c[1].GetBatch(i);
    a.c[2] = this->c[2].GetBatch(i);
    return a;
}

template <typename T, int Nbatch>
void ColorSpinor<T, Nbatch>::RandGauss(uint4 *state)
{
    T radius0;
    T radius1;
    T radius2;
    T phi0;
    T phi1;
    T phi2;

    phi0 = 2.0 * M_PI * GetRand<T>(state);
    phi1 = 2.0 * M_PI * GetRand<T>(state);
    phi2 = 2.0 * M_PI * GetRand<T>(state);

    radius0 = GetRand<T>(state);
    radius0 = radius0 + (1.0 - radius0) * MinVal<T>(); // exclude 0 from random numbers!
    radius0 = sqrt(-1.0 * log(radius0));
    radius1 = GetRand<T>(state);
    radius1 = radius1 + (1.0 - radius1) * MinVal<T>(); // exclude 0 from random numbers!
    radius1 = sqrt(-1.0 * log(radius1));
    radius2 = GetRand<T>(state);
    radius2 = radius2 + (1.0 - radius2) * MinVal<T>(); // exclude 0 from random numbers!
    radius2 = sqrt(-1.0 * log(radius2));

    c[0] = Complex<T, Nbatch>(radius0 * cos(phi0), radius0 * sin(phi0));
    c[1] = Complex<T, Nbatch>(radius1 * cos(phi1), radius1 * sin(phi1));
    c[2] = Complex<T, Nbatch>(radius2 * cos(phi2), radius2 * sin(phi2));
}

template <typename T, int Nbatch>
void ColorSpinor<T, Nbatch>::Clear()
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < (Nbatch << 1); ++j) {
            c[i].z[j] = 0;
        }
    }
}

template <typename T, int Nbatch>
void ColorSpinor<T, Nbatch>::PrintBatch(int i) const
{
    for (int nc = 0; nc < 3; ++nc) {
        Rank0Printf("%13.6e %13.6e ", c[nc].z[i], c[nc].z[i + Nbatch]);
    }
    Rank0Printf("\n");
}

template <typename T, int Nbatch>
T Norm2(const ColorSpinor<T, Nbatch> &a)
{
    T nrm = 0;
    for (int nc = 0; nc < 3; ++nc) {
        for (int i = 0; i < 2 * Nbatch; ++i) {
            nrm += a.c[nc].z[i] * a.c[nc].z[i];
        }
    }
    return nrm;
}

template <typename T, int Nbatch>
T Norm(const ColorSpinor<T, Nbatch> &a)
{
    return std::sqrt(Norm2(a));
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> Conj(const ColorSpinor<T, Nbatch> &a)
{
    ColorSpinor<T, Nbatch> b(a);
    b.c[0].Conj();
    b.c[1].Conj();
    b.c[2].Conj();
    return b;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> Dot(const ColorSpinor<T, Nbatch> &bra, const ColorSpinor<T, Nbatch> &ket)
{
    Complex<T, Nbatch> res = Conj(bra.c[0]) * ket.c[0];
    res += Conj(bra.c[1]) * ket.c[1];
    res += Conj(bra.c[2]) * ket.c[2];
    return res;
}

template <typename T, int Nbatch>
T RealDot(const ColorSpinor<T, Nbatch> &bra, const ColorSpinor<T, Nbatch> &ket)
{
    return RealSum(Dot(bra, ket));
}

// Non-member binary operators
template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator+(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator-(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs)
{
    lhs *= rhs;
    return lhs;
}

// Scalar multiplication/division
template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator*(T scalar, ColorSpinor<T, Nbatch> rhs)
{
    rhs *= scalar;
    return rhs;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch> lhs, T scalar)
{
    lhs *= scalar;
    return lhs;
}

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator/(ColorSpinor<T, Nbatch> lhs, T scalar)
{
    lhs /= scalar;
    return lhs;
}

// Comparison operators
template <typename T, int Nbatch>
bool operator==(const ColorSpinor<T, Nbatch> &lhs, const ColorSpinor<T, Nbatch> &rhs)
{
    return lhs.c[0] == rhs.c[0] && lhs.c[1] == rhs.c[1] && lhs.c[2] == rhs.c[2];
}

template <typename T, int Nbatch>
bool operator!=(const ColorSpinor<T, Nbatch> &lhs, const ColorSpinor<T, Nbatch> &rhs)
{
    return !(lhs == rhs);
}

// Stream output
template <typename T, int Nbatch>
std::ostream &operator<<(std::ostream &os, const ColorSpinor<T, Nbatch> &cs)
{
    os << "ColorSpinor(\n";
    os << "  red:   " << cs.c[0] << "\n";
    os << "  green: " << cs.c[1] << "\n";
    os << "  blue:  " << cs.c[2] << "\n";
    os << ")";
    return os;
}

// Explicit instantiations for specific template parameters
#define INSTANTIATE_COLOR_SPINOR_TEMPLATE(T, Nbatch)                                                          \
    template class ColorSpinor<T, Nbatch>;                                                                    \
    template T Norm2(const ColorSpinor<T, Nbatch> &a);                                                        \
    template T Norm(const ColorSpinor<T, Nbatch> &a);                                                         \
    template ColorSpinor<T, Nbatch> Conj(const ColorSpinor<T, Nbatch> &a);                                    \
    template Complex<T, Nbatch> Dot(const ColorSpinor<T, Nbatch> &bra, const ColorSpinor<T, Nbatch> &ket);    \
    template T RealDot(const ColorSpinor<T, Nbatch> &bra, const ColorSpinor<T, Nbatch> &ket);                 \
    template ColorSpinor<T, Nbatch> operator+(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs); \
    template ColorSpinor<T, Nbatch> operator-(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs); \
    template ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs); \
    template ColorSpinor<T, Nbatch> operator*(T scalar, ColorSpinor<T, Nbatch> rhs);                          \
    template ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch> lhs, T scalar);                          \
    template ColorSpinor<T, Nbatch> operator/(ColorSpinor<T, Nbatch> lhs, T scalar);                          \
    template bool operator==(const ColorSpinor<T, Nbatch> &lhs, const ColorSpinor<T, Nbatch> &rhs);           \
    template bool operator!=(const ColorSpinor<T, Nbatch> &lhs, const ColorSpinor<T, Nbatch> &rhs);           \
    template std::ostream &operator<<(std::ostream &os, const ColorSpinor<T, Nbatch> &cs)
// Instantiate the commonly used template parameter combinations
INSTANTIATE_COLOR_SPINOR_TEMPLATE(float, 1);
INSTANTIATE_COLOR_SPINOR_TEMPLATE(double, 1);
INSTANTIATE_COLOR_SPINOR_TEMPLATE(float, 16);
INSTANTIATE_COLOR_SPINOR_TEMPLATE(double, 8);

} // namespace kuqcd
