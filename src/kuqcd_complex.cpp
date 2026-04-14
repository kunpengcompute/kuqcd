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

#include "kuqcd_complex.h"


namespace kuqcd {

template <typename T, int Nbatch>
Complex<T, Nbatch>::Complex(const T real, const T imag)
{
    for (int i = 0; i < Nbatch; ++i) {
        z[i] = real;
        z[i + Nbatch] = imag;
    }
}

template <typename T, int Nbatch>
Complex<T, Nbatch>::Complex(const T values[2 * Nbatch])
{
    for (int i = 0; i < 2 * Nbatch; ++i) {
        z[i] = values[i];
    }
}

template <typename T, int Nbatch>
Complex<T, Nbatch>::Complex(const T realParts[Nbatch], const T imagParts[Nbatch])
{
    for (int i = 0; i < Nbatch; ++i) {
        z[i] = realParts[i];
        z[Nbatch + i] = imagParts[i];
    }
}

template <typename T, int Nbatch>
Complex<T, Nbatch>::Complex(const Complex &other)
{
    for (int i = 0; i < 2 * Nbatch; ++i) {
        z[i] = other.z[i];
    }
}

template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator=(T real)
{
    for (int i = 0; i < Nbatch; ++i) {
        z[i] = real;
        z[i + Nbatch] = 0;
    }
    return *this;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator=(const Complex<T, Nbatch> &other)
{
    if (this != &other) {  // Self-assignment check
        for (int i = 0; i < 2 * Nbatch; ++i) {
            z[i] = other.z[i];
        }
    }
    return *this;
}

// Safe data access with bounds checking
template <typename T, int Nbatch>
T &Complex<T, Nbatch>::DataSafe(int i)
{
    if ((i >= 2 * Nbatch) || (i < 0)) {
        throw std::out_of_range("Complex data index out of range. Must be < 2 * Nbatch and > 0.");
    }
    return z[i];
}

template <typename T, int Nbatch>
const T &Complex<T, Nbatch>::DataSafe(int i) const
{
    if ((i >= 2 * Nbatch) || (i < 0)) {
        throw std::out_of_range("Complex data index out of range. Must be < 2 * Nbatch and > 0.");
    }
    return z[i];
}

template <typename T, int Nbatch>
void Complex<T, Nbatch>::Conj()
{
    for (int i = 0; i < Nbatch; ++i) {
        z[Nbatch + i] = -z[Nbatch + i];
    }
}

template <typename T, int Nbatch>
Complex<T, 1> Complex<T, Nbatch>::GetBatch(int i) const
{
    Complex<T, 1> a(z[i], z[i + Nbatch]);
    return a;
}

// Compound assignment operators
template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator+=(const Complex<T, Nbatch> &other)
{
    for (int i = 0; i < 2 * Nbatch; ++i) {
        z[i] += other.z[i];
    }
    return *this;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator-=(const Complex<T, Nbatch> &other)
{
    for (int i = 0; i < 2 * Nbatch; ++i) {
        z[i] -= other.z[i];
    }
    return *this;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator*=(const Complex<T, Nbatch> &other)
{
    for (int i = 0; i < Nbatch; ++i) {
        T r1 = z[i];
        T i1 = z[Nbatch + i];
        T r2 = other.z[i];
        T i2 = other.z[Nbatch + i];
        z[i] = r1 * r2 - i1 * i2;
        z[Nbatch + i] = r1 * i2 + i1 * r2;
    }
    return *this;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator/=(const Complex<T, Nbatch> &other)
{
    for (int i = 0; i < Nbatch; ++i) {
        T r1 = z[i];
        T i1 = z[Nbatch + i];
        T r2 = other.z[i];
        T i2 = other.z[Nbatch + i];
        T den = r2 * r2 + i2 * i2;
        den = 1.0 / den;
        z[i] = (r1 * r2 + i1 * i2) * den;
        z[Nbatch + i] = (i1 * r2 - r1 * i2) * den;
    }
    return *this;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator*=(T scalar)
{
    for (int i = 0; i < 2 * Nbatch; ++i) {
        z[i] *= scalar;
    }
    return *this;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> &Complex<T, Nbatch>::operator/=(T scalar)
{
    scalar = 1.0 / scalar;
    for (int i = 0; i < 2 * Nbatch; ++i) {
        z[i] *= scalar;
    }
    return *this;
}

// Unary minus
template <typename T, int Nbatch>
Complex<T, Nbatch> Complex<T, Nbatch>::operator-() const
{
    Complex<T, Nbatch> res;
    for (int i = 0; i < 2 * Nbatch; ++i) {
        res.z[i] = -z[i];
    }
    return res;
}

// Non-member single operators
template <typename T, int Nbatch>
Complex<T, Nbatch> Conj(const Complex<T, Nbatch> &a)
{
    Complex<T, Nbatch> res(a);
    for (int i = 0; i < Nbatch; ++i) {
        res.z[Nbatch + i] = -res.z[Nbatch + i];
    }
    return res;
}

template <typename T, int Nbatch>
T Norm2(const Complex<T, Nbatch> &a, int i)
{
    return a.z[i] * a.z[i] + a.z[Nbatch + i] * a.z[Nbatch + i];
}

template <typename T, int Nbatch>
T Abs(const Complex<T, Nbatch> &a, int i)
{
    return std::sqrt(Norm2(a, i));
}

template <typename T, int Nbatch>
T Arg(const Complex<T, Nbatch> &a, int i)
{
    return std::atan2(a.z[Nbatch + i], a.z[i]);
}

template <typename T, int Nbatch>
T Norm2(const Complex<T, Nbatch> &a)
{
    T tmp = 0;
    for (int i = 0; i < Nbatch; ++i) {
        tmp += a.z[i] * a.z[i] + a.z[Nbatch + i] * a.z[Nbatch + i];
    }
    return tmp;
}

template <typename T, int Nbatch>
T Abs(const Complex<T, Nbatch> &a)
{
    return std::sqrt(Norm2(a));
}

template <typename T, int Nbatch>
T RealSum(const Complex<T, Nbatch> &a)
{
    T sum = T(0);
    for (int i = 0; i < Nbatch; ++i) {
        sum += a.z[i];
    }
    return sum;
}

// Non-member binary operators
template <typename T, int Nbatch>
Complex<T, Nbatch> operator+(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> operator-(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs)
{
    lhs -= rhs;
    return lhs;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> operator*(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs)
{
    lhs *= rhs;
    return lhs;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> operator/(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs)
{
    lhs /= rhs;
    return lhs;
}

// Scalar multiplication/division (non-member)
template <typename T, int Nbatch>
Complex<T, Nbatch> operator*(T scalar, Complex<T, Nbatch> rhs)
{
    rhs *= scalar;
    return rhs;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> operator*(Complex<T, Nbatch> lhs, T scalar)
{
    lhs *= scalar;
    return lhs;
}

template <typename T, int Nbatch>
Complex<T, Nbatch> operator/(Complex<T, Nbatch> lhs, T scalar)
{
    lhs /= scalar;
    return lhs;
}

// Comparison operators
template <typename T, int Nbatch>
bool operator==(const Complex<T, Nbatch> &lhs, const Complex<T, Nbatch> &rhs)
{
    T eps;
    if (sizeof(T) == sizeof(float)) {
        eps = 1.0e-5;
    } else if (sizeof(T) == sizeof(double)) {
        eps = 1.0e-10;
    }
    for (int i = 0; i < 2 * Nbatch; ++i) {
        if (std::fabs(lhs.z[i] - rhs.z[i]) / std::fabs(rhs.z[i]) > eps) {
            return false;
        }
    }
    return true;
}

template <typename T, int Nbatch>
bool operator!=(const Complex<T, Nbatch> &lhs, const Complex<T, Nbatch> &rhs)
{
    return !(lhs == rhs);
}

// Stream output
template <typename T, int Nbatch>
std::ostream &operator<<(std::ostream &os, const Complex<T, Nbatch> &a)
{
    os << "[";
    for (int i = 0; i < Nbatch; ++i) {
        os << std::scientific << std::setprecision(6) << a.Real(i);
        if (a.Imag(i) >= 0) {
            os << "+";
        }
        os << std::scientific << std::setprecision(6) << a.Imag(i) << "j";
        if (i < Nbatch - 1) {
            os << ", ";
        }
    }
    os << "]";
    os.unsetf(std::ios_base::floatfield);
    return os;
}

// Explicitly instantiate non-member functions for each required specialization
#define INSTANTIATE_COMPLEX_TEMPLATE(T, Nbatch)                                                   \
    template class Complex<T, Nbatch>;                                                            \
    template Complex<T, Nbatch> Conj(const Complex<T, Nbatch> &a);                                \
    template T Norm2(const Complex<T, Nbatch> &a, int i);                                         \
    template T Abs(const Complex<T, Nbatch> &a, int i);                                           \
    template T Arg(const Complex<T, Nbatch> &a, int i);                                           \
    template T Norm2(const Complex<T, Nbatch> &a);                                                \
    template T Abs(const Complex<T, Nbatch> &a);                                                  \
    template T RealSum(const Complex<T, Nbatch> &a);                                              \
    template Complex<T, Nbatch> operator+(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs); \
    template Complex<T, Nbatch> operator-(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs); \
    template Complex<T, Nbatch> operator*(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs); \
    template Complex<T, Nbatch> operator/(Complex<T, Nbatch> lhs, const Complex<T, Nbatch> &rhs); \
    template Complex<T, Nbatch> operator*(T scalar, Complex<T, Nbatch> rhs);                      \
    template Complex<T, Nbatch> operator*(Complex<T, Nbatch> lhs, T scalar);                      \
    template Complex<T, Nbatch> operator/(Complex<T, Nbatch> lhs, T scalar);                      \
    template bool operator==(const Complex<T, Nbatch> &lhs, const Complex<T, Nbatch> &rhs);       \
    template bool operator!=(const Complex<T, Nbatch> &lhs, const Complex<T, Nbatch> &rhs);       \
    template std::ostream &operator<<(std::ostream &os, const Complex<T, Nbatch> &a)

INSTANTIATE_COMPLEX_TEMPLATE(float, 1);
INSTANTIATE_COMPLEX_TEMPLATE(double, 1);
INSTANTIATE_COMPLEX_TEMPLATE(float, 16);
INSTANTIATE_COMPLEX_TEMPLATE(double, 8);

}  // namespace kuqcd