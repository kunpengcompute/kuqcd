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
#ifndef SPINOR_H
#define SPINOR_H

#include "kuqcd_complex.h"

namespace kuqcd {

template <typename T, int Nbatch>
class ColorSpinor {
public:
    Complex<T, Nbatch> c[3];
    ColorSpinor() = default;
    ColorSpinor(const ColorSpinor<T, Nbatch> &other);
    ColorSpinor(const Complex<T, Nbatch> &other);
    ColorSpinor(const Complex<T, Nbatch> &red, const Complex<T, Nbatch> &green, const Complex<T, Nbatch> &blue);
    explicit ColorSpinor(const Complex<T, Nbatch> ComponentsArray[3]);
    ColorSpinor(T real, T imag = 0);

    /**
     * @brief Safe element access with bounds checking
     * Use for debugging or when safety is preferred over performance
     */
    Complex<T, Nbatch> &ColorComponent(int i);
    const Complex<T, Nbatch> &ColorComponent(int i) const;

    // Assignment operators
    ColorSpinor<T, Nbatch> &operator=(const ColorSpinor<T, Nbatch> &other);
    ColorSpinor<T, Nbatch> &operator+=(const ColorSpinor<T, Nbatch> &other);
    ColorSpinor<T, Nbatch> &operator-=(const ColorSpinor<T, Nbatch> &other);
    ColorSpinor<T, Nbatch> &operator*=(const ColorSpinor<T, Nbatch> &other); // Component-wise multiplication
    ColorSpinor<T, Nbatch> &operator*=(T scalar);
    ColorSpinor<T, Nbatch> &operator/=(T scalar);

    // Unary operators
    ColorSpinor<T, Nbatch> operator-() const;

    void Conj();
    ColorSpinor<T, 1> GetBatch(int i = 0) const;

    void RandGauss(uint4 *state);

    void Clear();

    void PrintBatch(int i) const;
};

// Physics operations
template <typename T, int Nbatch>
T Norm2(const ColorSpinor<T, Nbatch> &a); // Sum of |component|^2 over all colors and batch elements

template <typename T, int Nbatch>
T Norm(const ColorSpinor<T, Nbatch> &a); // sqrt(norm_sq())

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> Conj(const ColorSpinor<T, Nbatch> &a);

// Color dot product: sum over color indices
template <typename T, int Nbatch>
Complex<T, Nbatch> Dot(const ColorSpinor<T, Nbatch> &bra, const ColorSpinor<T, Nbatch> &ket);

// Real part of color dot product
template <typename T, int Nbatch>
T RealDot(const ColorSpinor<T, Nbatch> &bra, const ColorSpinor<T, Nbatch> &ket);

// Non-member binary operators
template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator+(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs);

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator-(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs);

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch> lhs, const ColorSpinor<T, Nbatch> &rhs);

// Scalar multiplication/division
template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator*(T scalar, ColorSpinor<T, Nbatch> rhs);

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch> lhs, T scalar);

template <typename T, int Nbatch>
ColorSpinor<T, Nbatch> operator/(ColorSpinor<T, Nbatch> lhs, T scalar);

// Comparison operators
template <typename T, int Nbatch>
bool operator==(const ColorSpinor<T, Nbatch> &lhs, const ColorSpinor<T, Nbatch> &rhs);

template <typename T, int Nbatch>
bool operator!=(const ColorSpinor<T, Nbatch> &lhs, const ColorSpinor<T, Nbatch> &rhs);

// Stream output
template <typename T, int Nbatch>
std::ostream &operator<<(std::ostream &os, const ColorSpinor<T, Nbatch> &cs);

// Macro to declare extern templates
#define DECLARE_COLOR_SPINOR_EXTERN_TEMPLATE(T, Nbatch)                                                       \
    extern template class ColorSpinor<T, Nbatch>;                                                             \
    extern template T Norm2(const ColorSpinor<T, Nbatch> &);                                                  \
    extern template T Norm(const ColorSpinor<T, Nbatch> &);                                                   \
    extern template ColorSpinor<T, Nbatch> Conj(const ColorSpinor<T, Nbatch> &);                              \
    extern template Complex<T, Nbatch> Dot(const ColorSpinor<T, Nbatch> &, const ColorSpinor<T, Nbatch> &);   \
    extern template T RealDot(const ColorSpinor<T, Nbatch> &, const ColorSpinor<T, Nbatch> &);                \
    extern template ColorSpinor<T, Nbatch> operator+(ColorSpinor<T, Nbatch>, const ColorSpinor<T, Nbatch> &); \
    extern template ColorSpinor<T, Nbatch> operator-(ColorSpinor<T, Nbatch>, const ColorSpinor<T, Nbatch> &); \
    extern template ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch>, const ColorSpinor<T, Nbatch> &); \
    extern template ColorSpinor<T, Nbatch> operator*(T, ColorSpinor<T, Nbatch>);                              \
    extern template ColorSpinor<T, Nbatch> operator*(ColorSpinor<T, Nbatch>, T);                              \
    extern template ColorSpinor<T, Nbatch> operator/(ColorSpinor<T, Nbatch>, T);                              \
    extern template bool operator==(const ColorSpinor<T, Nbatch> &, const ColorSpinor<T, Nbatch> &);          \
    extern template bool operator!=(const ColorSpinor<T, Nbatch> &, const ColorSpinor<T, Nbatch> &);          \
    extern template std::ostream &operator<<(std::ostream &, const ColorSpinor<T, Nbatch> &)
// Declare extern templates for common instantiations
DECLARE_COLOR_SPINOR_EXTERN_TEMPLATE(float, 1);
DECLARE_COLOR_SPINOR_EXTERN_TEMPLATE(double, 1);
DECLARE_COLOR_SPINOR_EXTERN_TEMPLATE(float, 16);
DECLARE_COLOR_SPINOR_EXTERN_TEMPLATE(double, 8);

} // namespace kuqcd
#endif
