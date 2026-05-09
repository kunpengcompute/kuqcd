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
#ifndef SU3_H
#define SU3_H

#include "kuqcd_enum.h"
#include "kuqcd_complex.h"
#include "spinor.h"

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType Recon>
class SU3 {
public:
    Complex<T, Nbatch> m[Recon >> 1];
    SU3() = default;
    SU3(const SU3<T, Nbatch, Recon> &other);
    SU3<T, Nbatch, Recon> &operator=(const SU3<T, Nbatch, Recon> &other);

    // Matrix operations
    SU3<T, Nbatch, Recon> &operator+=(const SU3<T, Nbatch, Recon> &other);
    SU3<T, Nbatch, Recon> &operator-=(const SU3<T, Nbatch, Recon> &other);
    SU3<T, Nbatch, Recon> &operator*=(const SU3<T, Nbatch, Recon> &other);
    SU3<T, Nbatch, Recon> &operator*=(const Complex<T, Nbatch> &scalar);
    SU3<T, Nbatch, Recon> &operator*=(T scalar);
    SU3<T, Nbatch, Recon> &operator/=(T scalar);

    // Unary operators
    SU3<T, Nbatch, Recon> operator-() const;

    // SU(3) group operations
    void Conj();                      // In-place conjugate
    void Trans();                     // In-place transpose
    void Dagger();                    // In-place dagger
    Complex<T, Nbatch> Trace() const; // Matrix trace
    T RealTrace() const;              // Real part of trace

    // Special matrices
    void Identity(); // Set to identity matrix

    // SU(3) specific operations
    void Unitarize();             // Project to SU(3) group
    void Unitarize(int j);        // Project to SU(3) group for batch j
    T UnitarityDeviation() const; // Measure unitarity violation
    SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO> Reconstruct() const;

    void Print() const;
    void PrintBatch(int i) const;
    void RankPrintBatch(int rank, int i) const;
};

// Non-member binary operators
template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator+(SU3<T, Nbatch, Recon> lhs, const SU3<T, Nbatch, Recon> &rhs);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator-(SU3<T, Nbatch, Recon> lhs, const SU3<T, Nbatch, Recon> &rhs);

// Scalar multiplication
template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(T scalar, SU3<T, Nbatch, Recon> rhs);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon> lhs, T scalar);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(const Complex<T, Nbatch> &scalar, SU3<T, Nbatch, Recon> rhs);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon> lhs, const Complex<T, Nbatch> &scalar);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator/(SU3<T, Nbatch, Recon> lhs, T scalar);

// Stream output
template <typename T, int Nbatch, KuQCDReconstructType Recon>
std::ostream &operator<<(std::ostream &os, const SU3<T, Nbatch, Recon> &matrix);

// Macro to declare extern templates
#define DECLARE_SU3_EXTERN_TEMPLATE(T, Nbatch, Recon)                                                             \
    extern template class SU3<T, Nbatch, Recon>;                                                                  \
    extern template SU3<T, Nbatch, Recon> operator+(SU3<T, Nbatch, Recon>, const SU3<T, Nbatch, Recon> &);        \
    extern template SU3<T, Nbatch, Recon> operator-(SU3<T, Nbatch, Recon>, const SU3<T, Nbatch, Recon> &);        \
    extern template SU3<T, Nbatch, Recon> operator*(T, SU3<T, Nbatch, Recon>);                                    \
    extern template SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon>, T);                                    \
    extern template SU3<T, Nbatch, Recon> operator*(const Complex<T, Nbatch> &scalar, SU3<T, Nbatch, Recon> rhs); \
    extern template SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon> lhs, const Complex<T, Nbatch> &scalar); \
    extern template SU3<T, Nbatch, Recon> operator/(SU3<T, Nbatch, Recon> lhs, T scalar);                         \
    extern template std::ostream &operator<<(std::ostream &, const SU3<T, Nbatch, Recon> &)
// Declare extern templates for common instantiations
DECLARE_SU3_EXTERN_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_NO);
DECLARE_SU3_EXTERN_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_NO);
DECLARE_SU3_EXTERN_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_NO);
DECLARE_SU3_EXTERN_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_NO);

DECLARE_SU3_EXTERN_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_12);
DECLARE_SU3_EXTERN_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_12);
DECLARE_SU3_EXTERN_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_12);
DECLARE_SU3_EXTERN_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_12);

DECLARE_SU3_EXTERN_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_U14);
DECLARE_SU3_EXTERN_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_U14);
DECLARE_SU3_EXTERN_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_U14);
DECLARE_SU3_EXTERN_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_U14);

} // namespace kuqcd

#endif
