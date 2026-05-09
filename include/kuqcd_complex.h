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
#ifndef KUQCD_COMPLEX_H
#define KUQCD_COMPLEX_H

#include <omp.h>
#include <iostream>
#include "kuqcd_utils.h"

namespace kuqcd {

template <typename T, int Nbatch>
class Complex {
public:
    T z[Nbatch << 1] = {0};

    // Constructors
    Complex() = default;
    explicit Complex(T real, T imag = 0); // For Nbatch = 1
    explicit Complex(const T values[2 * Nbatch]);
    explicit Complex(const T realParts[Nbatch], const T imagParts[Nbatch]);

    // Copy constructor and assignment operator
    Complex(const Complex<T, Nbatch> &other);
    Complex<T, Nbatch> &operator=(T real);
    Complex<T, Nbatch> &operator=(const Complex<T, Nbatch> &other);

    /**
     * @brief Safe raw data access with bounds checking
     * Use for debugging or when safety is preferred over performance
     */
    T &DataSafe(int i);
    const T &DataSafe(int i) const;

    // Named component access - zero overhead, always inline
    constexpr T Real(int i = 0) const noexcept
    {
        return z[i];
    }
    constexpr T &Real(int i = 0) noexcept
    {
        return z[i];
    }
    constexpr T Imag(int i = 0) const noexcept
    {
        return z[Nbatch + i];
    }
    constexpr T &Imag(int i = 0) noexcept
    {
        return z[Nbatch + i];
    }

    // Member functions
    void Conj();
    Complex<T, 1> GetBatch(int i = 0) const;

    // Compound assignment operators
    Complex<T, Nbatch> &operator+=(const Complex<T, Nbatch> &other);
    Complex<T, Nbatch> &operator-=(const Complex<T, Nbatch> &other);
    Complex<T, Nbatch> &operator*=(const Complex<T, Nbatch> &other);
    Complex<T, Nbatch> &operator/=(const Complex<T, Nbatch> &other);
    Complex<T, Nbatch> &operator*=(T scalar); // Scalar multiplication
    Complex<T, Nbatch> &operator/=(T scalar); // Scalar division

    // Unary minus
    Complex<T, Nbatch> operator-() const;
};

// Non-member single operators
template <typename T, int Nbatch>
Complex<T, Nbatch> Conj(const Complex<T, Nbatch> &a);

template <typename T, int Nbatch>
T Norm2(const Complex<T, Nbatch> &a, int i); // |z_i|^2

template <typename T, int Nbatch>
T Abs(const Complex<T, Nbatch> &a, int i); // |z_i|

template <typename T, int Nbatch>
T Arg(const Complex<T, Nbatch> &a, int i); // arg(z_i)

template <typename T, int Nbatch>
T Norm2(const Complex<T, Nbatch> &a); // Sum of norms squared

template <typename T, int Nbatch>
T Abs(const Complex<T, Nbatch> &a); // Sum of norms

template <typename T, int Nbatch>
T RealSum(const Complex<T, Nbatch> &a); // Sum of real parts

// Non-member binary operators
template <typename T, int Nbatch>
Complex<T, Nbatch> operator+(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs);

template <typename T, int Nbatch>
Complex<T, Nbatch> operator-(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs);

template <typename T, int Nbatch>
Complex<T, Nbatch> operator*(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs);

template <typename T, int Nbatch>
Complex<T, Nbatch> operator/(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs);

// Scalar multiplication/division (non-member)
template <typename T, int Nbatch>
Complex<T, Nbatch> operator*(T scalar, Complex<T, Nbatch> rhs);

template <typename T, int Nbatch>
Complex<T, Nbatch> operator*(Complex<T, Nbatch> lhs, T scalar);

template <typename T, int Nbatch>
Complex<T, Nbatch> operator/(Complex<T, Nbatch> lhs, T scalar);

// Comparison operators (primarily for Nbatch = 1, or element-wise for Nbatch > 1 if needed)
template <typename T, int Nbatch>
bool operator==(const Complex<T, Nbatch> &lhs, const Complex<T, Nbatch> &rhs);

template <typename T, int Nbatch>
bool operator!=(const Complex<T, Nbatch> &lhs, const Complex<T, Nbatch> &rhs);

// Stream output
template <typename T, int Nbatch>
std::ostream &operator<<(std::ostream &os, const Complex<T, Nbatch> &a);

// Macro to declare extern templates for a given Complex type
#define DECLARE_COMPLEX_EXTERN_TEMPLATE(T, Nbatch)                                                \
    extern template class Complex<T, Nbatch>;                                                     \
    extern template Complex<T, Nbatch> Conj(const Complex<T, Nbatch> &);                          \
    extern template T Norm2(const Complex<T, Nbatch> &, int);                                     \
    extern template T Abs(const Complex<T, Nbatch> &, int);                                       \
    extern template T Arg(const Complex<T, Nbatch> &, int);                                       \
    extern template T Norm2(const Complex<T, Nbatch> &);                                          \
    extern template T Abs(const Complex<T, Nbatch> &);                                            \
    extern template T RealSum(const Complex<T, Nbatch> &);                                        \
    extern template Complex<T, Nbatch> operator+(Complex<T, Nbatch>, const Complex<T, Nbatch> &); \
    extern template Complex<T, Nbatch> operator-(Complex<T, Nbatch>, const Complex<T, Nbatch> &); \
    extern template Complex<T, Nbatch> operator*(Complex<T, Nbatch>, const Complex<T, Nbatch> &); \
    extern template Complex<T, Nbatch> operator/(Complex<T, Nbatch>, const Complex<T, Nbatch> &); \
    extern template Complex<T, Nbatch> operator*(T, Complex<T, Nbatch>);                          \
    extern template Complex<T, Nbatch> operator*(Complex<T, Nbatch>, T);                          \
    extern template Complex<T, Nbatch> operator/(Complex<T, Nbatch>, T);                          \
    extern template bool operator==(const Complex<T, Nbatch> &, const Complex<T, Nbatch> &);      \
    extern template bool operator!=(const Complex<T, Nbatch> &, const Complex<T, Nbatch> &);      \
    extern template std::ostream &operator<<(std::ostream &, const Complex<T, Nbatch> &)

// Declare extern templates for the specific instantiations
DECLARE_COMPLEX_EXTERN_TEMPLATE(float, 1);
DECLARE_COMPLEX_EXTERN_TEMPLATE(double, 1);
DECLARE_COMPLEX_EXTERN_TEMPLATE(float, 16);
DECLARE_COMPLEX_EXTERN_TEMPLATE(double, 8);

} // namespace kuqcd
#endif
