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
#ifndef KUQCD_MATH_H
#define KUQCD_MATH_H

#include "spinor.h"
#include "su3.h"

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType Recon>
ColorSpinor<T, Nbatch> UMultV(const SU3<T, Nbatch, Recon> &u, const ColorSpinor<T, Nbatch> &v);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UMultV(const SU3<T, Nbatch, Recon> &u, const ColorSpinor<T, Nbatch> &src, ColorSpinor<T, Nbatch> &res);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UMultV(const SU3<T, Nbatch, Recon> &u, const ColorSpinor<T, Nbatch> &src, ColorSpinor<T, Nbatch> &res, T c);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UMultV(const T *const u, const T *const src, T *res, T c);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
ColorSpinor<T, Nbatch> UDaggMultV(const SU3<T, Nbatch, Recon> &u, const ColorSpinor<T, Nbatch> &v);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UDaggMultV(const SU3<T, Nbatch, Recon> &u, const ColorSpinor<T, Nbatch> &src, ColorSpinor<T, Nbatch> &res);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UDaggMultV(const SU3<T, Nbatch, Recon> &u, const ColorSpinor<T, Nbatch> &src, ColorSpinor<T, Nbatch> &res, T c);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UDaggMultV(const T *const u, const T *const src, T *res, T c);

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> UMultU(const SU3<T, Nbatch, Recon> &u1, const SU3<T, Nbatch, Recon> &u2);

#ifdef SVEBITS

#if SVEBITS == 512

template <>
ColorSpinor<float, 16> UMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO> &, const ColorSpinor<float, 16> &);

template <>
void UMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO> &, const ColorSpinor<float, 16> &, ColorSpinor<float, 16> &);

template <>
void UMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO> &, const ColorSpinor<float, 16> &, ColorSpinor<float, 16> &,
            float);

template <>
void UMultV<float, 16, KUQCD_GAUGE_RECON_NO>(const float *const u, const float *const src, float *res, float c);

template <>
ColorSpinor<float, 16> UDaggMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO> &, const ColorSpinor<float, 16> &);

template <>
void UDaggMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO> &, const ColorSpinor<float, 16> &, ColorSpinor<float, 16> &);

template <>
void UDaggMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO> &, const ColorSpinor<float, 16> &, ColorSpinor<float, 16> &,
                float);

template <>
void UDaggMultV<float, 16, KUQCD_GAUGE_RECON_NO>(const float *const u, const float *const src, float *res, float c);

#endif

#endif

#define DECLARE_MATH_EXTERN_TEMPLATE(T, Nbatch, Recon)                                                                \
    extern template ColorSpinor<T, Nbatch> UMultV(const SU3<T, Nbatch, Recon> &, const ColorSpinor<T, Nbatch> &);     \
    extern template void UMultV(const SU3<T, Nbatch, Recon> &, const ColorSpinor<T, Nbatch> &,                        \
                                ColorSpinor<T, Nbatch> &);                                                            \
    extern template void UMultV(const SU3<T, Nbatch, Recon> &, const ColorSpinor<T, Nbatch> &,                        \
                                ColorSpinor<T, Nbatch> &, T);                                                         \
    extern template void UMultV<T, Nbatch, Recon>(const T *const u, const T *const src, T *res, T c);                 \
    extern template ColorSpinor<T, Nbatch> UDaggMultV(const SU3<T, Nbatch, Recon> &, const ColorSpinor<T, Nbatch> &); \
    extern template void UDaggMultV(const SU3<T, Nbatch, Recon> &, const ColorSpinor<T, Nbatch> &,                    \
                                    ColorSpinor<T, Nbatch> &);                                                        \
    extern template void UDaggMultV(const SU3<T, Nbatch, Recon> &, const ColorSpinor<T, Nbatch> &,                    \
                                    ColorSpinor<T, Nbatch> &, T);                                                     \
    extern template void UDaggMultV<T, Nbatch, Recon>(const T *const u, const T *const src, T *res, T c);             \
    extern template SU3<T, Nbatch, Recon> UMultU(const SU3<T, Nbatch, Recon> &, const SU3<T, Nbatch, Recon> &)

DECLARE_MATH_EXTERN_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_NO);
DECLARE_MATH_EXTERN_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_NO);
DECLARE_MATH_EXTERN_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_NO);
DECLARE_MATH_EXTERN_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_NO);

DECLARE_MATH_EXTERN_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_12);
DECLARE_MATH_EXTERN_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_12);
DECLARE_MATH_EXTERN_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_12);
DECLARE_MATH_EXTERN_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_12);

DECLARE_MATH_EXTERN_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_U14);
DECLARE_MATH_EXTERN_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_U14);
DECLARE_MATH_EXTERN_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_U14);
DECLARE_MATH_EXTERN_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_U14);

} // namespace kuqcd
#endif
