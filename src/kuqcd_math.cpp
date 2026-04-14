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
#include "kuqcd_math.h"

#ifdef SVEBITS
#include "arm_sve.h"
#endif

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType Recon>
ColorSpinor<T, Nbatch> UMultV(const SU3<T, Nbatch, Recon>& u, const ColorSpinor<T, Nbatch>& v)
{
    ColorSpinor<T, Nbatch> res;
    for (int raw = 0; raw < 3; ++raw) {
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                res.c[raw].z[i       ] += u.m[col + 3*raw].z[i] * v.c[col].z[i] - u.m[col + 3*raw].z[i+Nbatch] * v.c[col].z[i+Nbatch];
                res.c[raw].z[i+Nbatch] += u.m[col + 3*raw].z[i] * v.c[col].z[i+Nbatch] + u.m[col + 3*raw].z[i+Nbatch] * v.c[col].z[i];
            }
        }
    }
    return res;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UMultV(const SU3<T, Nbatch, Recon>& u, const ColorSpinor<T, Nbatch>& src, ColorSpinor<T, Nbatch>& res)
{
    res = ColorSpinor<T, Nbatch>();
    for (int raw = 0; raw < 3; ++raw) {
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                res.c[raw].z[i       ] += u.m[col + 3*raw].z[i] * src.c[col].z[i] - u.m[col + 3*raw].z[i+Nbatch] * src.c[col].z[i+Nbatch];
                res.c[raw].z[i+Nbatch] += u.m[col + 3*raw].z[i] * src.c[col].z[i+Nbatch] + u.m[col + 3*raw].z[i+Nbatch] * src.c[col].z[i];
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UMultV(const SU3<T, Nbatch, Recon>& u, const ColorSpinor<T, Nbatch>& src, ColorSpinor<T, Nbatch>& res, T c)
{
    for (int raw = 0; raw < 3; ++raw) {
        T tmp[2*Nbatch]{0.0};
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                tmp[i       ] += u.m[col + 3*raw].z[i] * src.c[col].z[i] - u.m[col + 3*raw].z[i+Nbatch] * src.c[col].z[i+Nbatch];
                tmp[i+Nbatch] += u.m[col + 3*raw].z[i] * src.c[col].z[i+Nbatch] + u.m[col + 3*raw].z[i+Nbatch] * src.c[col].z[i];
            }
        }
        for (int i = 0; i < 2*Nbatch; ++i) {
            res.c[raw].z[i] += c * tmp[i];
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UMultV(const T* const u, const T* const src, T* res, T c)
{
    for (int raw = 0; raw < 3; ++raw) {
        T tmp[2*Nbatch]{0.0};
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                tmp[i       ] += u[i+2*Nbatch*(col+3*raw)] * src[i+2*Nbatch*col] - u[i+Nbatch+2*Nbatch*(col+3*raw)] * src[i+Nbatch+2*Nbatch*col];
                tmp[i+Nbatch] += u[i+2*Nbatch*(col+3*raw)] * src[i+Nbatch+2*Nbatch*col] + u[i+Nbatch+2*Nbatch*(col+3*raw)] * src[i+2*Nbatch*col];
            }
        }
        for (int i = 0; i < 2*Nbatch; ++i) {
            res[i+2*Nbatch*raw] += c * tmp[i];
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
ColorSpinor<T, Nbatch> UDaggMultV(const SU3<T, Nbatch, Recon>& u, const ColorSpinor<T, Nbatch>& v)
{
    ColorSpinor<T, Nbatch> res;
    for (int raw = 0; raw < 3; ++raw) {
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                res.c[raw].z[i       ] += u.m[raw + 3*col].z[i] * v.c[col].z[i] + u.m[raw + 3*col].z[i+Nbatch] * v.c[col].z[i+Nbatch];
                res.c[raw].z[i+Nbatch] += u.m[raw + 3*col].z[i] * v.c[col].z[i+Nbatch] - u.m[raw + 3*col].z[i+Nbatch] * v.c[col].z[i];
            }
        }
    }
    return res;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UDaggMultV(const SU3<T, Nbatch, Recon>& u, const ColorSpinor<T, Nbatch>& src, ColorSpinor<T, Nbatch>& res)
{
    res = ColorSpinor<T, Nbatch>();
    for (int raw = 0; raw < 3; ++raw) {
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                res.c[raw].z[i       ] += u.m[raw + 3*col].z[i] * src.c[col].z[i] + u.m[raw + 3*col].z[i+Nbatch] * src.c[col].z[i+Nbatch];
                res.c[raw].z[i+Nbatch] += u.m[raw + 3*col].z[i] * src.c[col].z[i+Nbatch] - u.m[raw + 3*col].z[i+Nbatch] * src.c[col].z[i];
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UDaggMultV(const SU3<T, Nbatch, Recon>& u, const ColorSpinor<T, Nbatch>& src, ColorSpinor<T, Nbatch>& res, T c)
{
    for (int raw = 0; raw < 3; ++raw) {
        T tmp[2*Nbatch]{0.0};
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                tmp[i       ] += u.m[raw + 3*col].z[i] * src.c[col].z[i] + u.m[raw + 3*col].z[i+Nbatch] * src.c[col].z[i+Nbatch];
                tmp[i+Nbatch] += u.m[raw + 3*col].z[i] * src.c[col].z[i+Nbatch] - u.m[raw + 3*col].z[i+Nbatch] * src.c[col].z[i];
            }
        }
        for (int i = 0; i < 2*Nbatch; ++i) {
            res.c[raw].z[i] -= c * tmp[i];
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void UDaggMultV(const T* const u, const T* const src, T* res, T c)
{
    for (int raw = 0; raw < 3; ++raw) {
        T tmp[2*Nbatch]{0.0};
        for (int col = 0; col < 3; ++col) {
            for (int i = 0; i < Nbatch; ++i) {
                tmp[i       ] += u[i+2*Nbatch*(raw+3*col)] * src[i+2*Nbatch*col] + u[i+Nbatch+2*Nbatch*(raw+3*col)] * src[i+Nbatch+2*Nbatch*col];
                tmp[i+Nbatch] += u[i+2*Nbatch*(raw+3*col)] * src[i+Nbatch+2*Nbatch*col] - u[i+Nbatch+2*Nbatch*(raw+3*col)] * src[i+2*Nbatch*col];
            }
        }
        for (int i = 0; i < 2*Nbatch; ++i) {
            res[i+2*Nbatch*raw] -= c * tmp[i];
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> UMultU(const SU3<T, Nbatch, Recon>& u1, const SU3<T, Nbatch, Recon>& u2)
{
    SU3<T, Nbatch, Recon> u;
    for (int mm = 0; mm < 3; ++mm) {
        for (int kk = 0; kk < 3; ++kk) {
            for (int nn = 0; nn < 3; ++nn) {
                for (int i = 0; i < Nbatch; ++i) {
                    u.m[nn + 3*mm].z[i       ] += u1.m[kk + 3*mm].z[i] * u2.m[nn + 3*kk].z[i] - u1.m[kk + 3*mm].z[i+Nbatch] * u2.m[nn + 3*kk].z[i+Nbatch];
                    u.m[nn + 3*mm].z[i+Nbatch] += u1.m[kk + 3*mm].z[i] * u2.m[nn + 3*kk].z[i+Nbatch] + u1.m[kk + 3*mm].z[i+Nbatch] * u2.m[nn + 3*kk].z[i];
                }
            }
        }
    }
    return u;
}

#define INSTANTIATE_MATH_TEMPLATE(T, Nbatch, Recon) \
    template ColorSpinor<T, Nbatch> UMultV(const SU3<T, Nbatch, Recon>&, const ColorSpinor<T, Nbatch>&); \
    template void UMultV(const SU3<T, Nbatch, Recon>&, const ColorSpinor<T, Nbatch>&, ColorSpinor<T, Nbatch>&); \
    template void UMultV(const SU3<T, Nbatch, Recon>&, const ColorSpinor<T, Nbatch>&, ColorSpinor<T, Nbatch>&, T);\
    template void UMultV<T, Nbatch, Recon>(const T* const u, const T* const src, T* res, T c);               \
    template ColorSpinor<T, Nbatch> UDaggMultV(const SU3<T, Nbatch, Recon>&, const ColorSpinor<T, Nbatch>&);    \
    template void UDaggMultV(const SU3<T, Nbatch, Recon>&, const ColorSpinor<T, Nbatch>&, ColorSpinor<T, Nbatch>&); \
    template void UDaggMultV(const SU3<T, Nbatch, Recon>&, const ColorSpinor<T, Nbatch>&, ColorSpinor<T, Nbatch>&, T);\
    template void UDaggMultV<T, Nbatch, Recon>(const T* const u, const T* const src, T* res, T c);           \
    template SU3<T, Nbatch, Recon> UMultU(const SU3<T, Nbatch, Recon>&, const SU3<T, Nbatch, Recon>&);

INSTANTIATE_MATH_TEMPLATE(float,  1, KUQCD_GAUGE_RECON_NO)
INSTANTIATE_MATH_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_NO)
#ifndef SVEBITS
INSTANTIATE_MATH_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_NO)
#endif
INSTANTIATE_MATH_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_NO)

INSTANTIATE_MATH_TEMPLATE(float,  1, KUQCD_GAUGE_RECON_12)
INSTANTIATE_MATH_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_12)
INSTANTIATE_MATH_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_12)
INSTANTIATE_MATH_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_12)

INSTANTIATE_MATH_TEMPLATE(float,  1, KUQCD_GAUGE_RECON_U14)
INSTANTIATE_MATH_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_U14)
INSTANTIATE_MATH_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_U14)
INSTANTIATE_MATH_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_U14)

#ifdef SVEBITS

#if SVEBITS == 512

template <>
ColorSpinor<float, 16> UMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO>& u, const ColorSpinor<float, 16>& v)
{
    int Nbatch = 16;
    ColorSpinor<float, 16> res;
    svfloat32_t src1, src2, res1, res2, u1, u2;
    svbool_t pg = svptrue_b32();
    #pragma unroll
    for (int raw = 0; raw < 3; ++raw) {
        res1 = svdup_n_f32(0);
        res2 = svdup_n_f32(0);
        #pragma unroll
        for (int col = 0; col < 3; ++col) {
            src1 = svld1_f32(pg, &v.c[col].z[0]);
            src2 = svld1_f32(pg, &v.c[col].z[Nbatch]);
            u1 = svld1_f32(pg, &u.m[col + 3*raw].z[0]);
            u2 = svld1_f32(pg, &u.m[col + 3*raw].z[Nbatch]);
            res1 = svmla_f32_m(pg, svmls_f32_m(pg, res1, u2, src2), u1, src1);
            res2 = svmla_f32_m(pg, svmla_f32_m(pg, res2, u2, src1), u1, src2);
        }
        svst1_f32(pg, &res.c[raw].z[0], res1);
        svst1_f32(pg, &res.c[raw].z[Nbatch], res2);
    }
    return res;
}

template <>
void UMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO>& u, const ColorSpinor<float, 16>& src, ColorSpinor<float, 16>& res)
{
    int Nbatch = 16;
    res = ColorSpinor<float, 16>();
    svfloat32_t src1, src2, res1, res2, u1, u2;
    svbool_t pg = svptrue_b32();
    #pragma unroll
    for (int raw = 0; raw < 3; ++raw) {
        res1 = svdup_n_f32(0);
        res2 = svdup_n_f32(0);
        #pragma unroll
        for (int col = 0; col < 3; ++col) {
            src1 = svld1_f32(pg, &src.c[col].z[0]);
            src2 = svld1_f32(pg, &src.c[col].z[Nbatch]);
            u1 = svld1_f32(pg, &u.m[col + 3*raw].z[0]);
            u2 = svld1_f32(pg, &u.m[col + 3*raw].z[Nbatch]);
            res1 = svmla_f32_m(pg, svmls_f32_m(pg, res1, u2, src2), u1, src1);
            res2 = svmla_f32_m(pg, svmla_f32_m(pg, res2, u2, src1), u1, src2);
        }
        svst1_f32(pg, &res.c[raw].z[0], res1);
        svst1_f32(pg, &res.c[raw].z[Nbatch], res2);
    }
}

template <>
void UMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO>& u, const ColorSpinor<float, 16>& src, ColorSpinor<float, 16>& res, float c)
{
    int Nbatch = 16;
    svfloat32_t src1, src2, res1, res2, tmp1, tmp2, u1, u2, tmpc;
    tmpc = svdup_n_f32(c);
    svbool_t pg = svptrue_b32();
    #pragma unroll
    for (int raw = 0; raw < 3; ++raw) {
        tmp1 = svdup_n_f32(0);
        tmp2 = svdup_n_f32(0);
        #pragma unroll
        for (int col = 0; col < 3; ++col) {
            src1 = svld1_f32(pg, &src.c[col].z[0]);
            src2 = svld1_f32(pg, &src.c[col].z[Nbatch]);
            u1 = svld1_f32(pg, &u.m[col + 3*raw].z[0]);
            u2 = svld1_f32(pg, &u.m[col + 3*raw].z[Nbatch]);
            tmp1 = svmla_f32_m(pg, svmls_f32_m(pg, tmp1, u2, src2), u1, src1);
            tmp2 = svmla_f32_m(pg, svmla_f32_m(pg, tmp2, u2, src1), u1, src2);
        }
        res1 = svld1_f32(pg, &res.c[raw].z[0]);
        res2 = svld1_f32(pg, &res.c[raw].z[Nbatch]);
        res1 = svmla_f32_m(pg, res1, tmp1, tmpc);
        res2 = svmla_f32_m(pg, res2, tmp2, tmpc);
        svst1_f32(pg, &res.c[raw].z[0], res1);
        svst1_f32(pg, &res.c[raw].z[Nbatch], res2);
    }
}

template <>
void UMultV<float, 16, KUQCD_GAUGE_RECON_NO>(const float* const u, const float* const src, float* res, float c)
{
    int Nbatch = 16;
    svbool_t pg = svptrue_b32();
    svfloat32_t src1, src2, tmp1, tmp2, res1, res2, u1, u2, tmpc;
    tmpc = svdup_n_f32(c);
    #pragma unroll
    for (int m = 0; m < 3; ++m) {
        tmp1 = svdup_n_f32(0);
        tmp2 = svdup_n_f32(0);
        #pragma unroll
        for (int n = 0; n < 3; ++n) {
            src1 = svld1_f32(pg, src+2*n*Nbatch);
            src2 = svld1_f32(pg, src+(1+2*n)*Nbatch);
            u1 = svld1_f32(pg, u+2*(n+3*m)*Nbatch);
            u2 = svld1_f32(pg, u+(1+2*(n+3*m))*Nbatch);
            tmp1 = svmla_f32_m(pg, svmls_f32_m(pg, tmp1, u2, src2), u1, src1);
            tmp2 = svmla_f32_m(pg, svmla_f32_m(pg, tmp2, u2, src1), u1, src2);
        }
        res1 = svld1_f32(pg, res+2*m*Nbatch);
        res2 = svld1_f32(pg, res+(1+2*m)*Nbatch);
        res1 = svmla_f32_m(pg, res1, tmp1, tmpc);
        res2 = svmla_f32_m(pg, res2, tmp2, tmpc);
        svst1_f32(pg, res+2*m*Nbatch, res1);
        svst1_f32(pg, res+(1+2*m)*Nbatch, res2);
    }
}

template <>
ColorSpinor<float, 16> UDaggMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO>& u, const ColorSpinor<float, 16>& v)
{
    int Nbatch = 16;
    ColorSpinor<float, 16> res;
    svfloat32_t src1, src2, res1, res2, u1, u2;
    svbool_t pg = svptrue_b32();
    #pragma unroll
    for (int raw = 0; raw < 3; ++raw) {
        res1 = svdup_n_f32(0);
        res2 = svdup_n_f32(0);
        #pragma unroll
        for (int col = 0; col < 3; ++col) {
            src1 = svld1_f32(pg, &v.c[col].z[0]);
            src2 = svld1_f32(pg, &v.c[col].z[Nbatch]);
            u1 = svld1_f32(pg, &u.m[raw + 3*col].z[0]);
            u2 = svld1_f32(pg, &u.m[raw + 3*col].z[Nbatch]);
            res1 = svmla_f32_m(pg, svmla_f32_m(pg, res1, u2, src2), u1, src1);
            res2 = svmla_f32_m(pg, svmls_f32_m(pg, res2, u2, src1), u1, src2);
        }
        svst1_f32(pg, &res.c[raw].z[0], res1);
        svst1_f32(pg, &res.c[raw].z[Nbatch], res2);
    }
    return res;
}

template <>
void UDaggMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO>& u, const ColorSpinor<float, 16>& src, ColorSpinor<float, 16>& res)
{
    int Nbatch = 16;
    res = ColorSpinor<float, 16>();
    svfloat32_t src1, src2, res1, res2, u1, u2;
    svbool_t pg = svptrue_b32();
    #pragma unroll
    for (int raw = 0; raw < 3; ++raw) {
        res1 = svdup_n_f32(0);
        res2 = svdup_n_f32(0);
        #pragma unroll
        for (int col = 0; col < 3; ++col) {
            src1 = svld1_f32(pg, &src.c[col].z[0]);
            src2 = svld1_f32(pg, &src.c[col].z[Nbatch]);
            u1 = svld1_f32(pg, &u.m[raw + 3*col].z[0]);
            u2 = svld1_f32(pg, &u.m[raw + 3*col].z[Nbatch]);
            res1 = svmla_f32_m(pg, svmla_f32_m(pg, res1, u2, src2), u1, src1);
            res2 = svmla_f32_m(pg, svmls_f32_m(pg, res2, u2, src1), u1, src2);
        }
        svst1_f32(pg, &res.c[raw].z[0], res1);
        svst1_f32(pg, &res.c[raw].z[Nbatch], res2);
    }
}

template <>
void UDaggMultV(const SU3<float, 16, KUQCD_GAUGE_RECON_NO>& u, const ColorSpinor<float, 16>& src, ColorSpinor<float, 16>& res, float c)
{
    int Nbatch = 16;
    svfloat32_t src1, src2, res1, res2, tmp1, tmp2, u1, u2, tmpc;
    tmpc = svdup_n_f32(c);
    svbool_t pg = svptrue_b32();
    #pragma unroll
    for (int raw = 0; raw < 3; ++raw) {
        tmp1 = svdup_n_f32(0);
        tmp2 = svdup_n_f32(0);
        #pragma unroll
        for (int col = 0; col < 3; ++col) {
            src1 = svld1_f32(pg, &src.c[col].z[0]);
            src2 = svld1_f32(pg, &src.c[col].z[Nbatch]);
            u1 = svld1_f32(pg, &u.m[raw + 3*col].z[0]);
            u2 = svld1_f32(pg, &u.m[raw + 3*col].z[Nbatch]);
            tmp1 = svmla_f32_m(pg, svmla_f32_m(pg, tmp1, u2, src2), u1, src1);
            tmp2 = svmla_f32_m(pg, svmls_f32_m(pg, tmp2, u2, src1), u1, src2);
        }
        res1 = svld1_f32(pg, &res.c[raw].z[0]);
        res2 = svld1_f32(pg, &res.c[raw].z[Nbatch]);
        res1 = svmls_f32_m(pg, res1, tmp1, tmpc);
        res2 = svmls_f32_m(pg, res2, tmp2, tmpc);
        svst1_f32(pg, &res.c[raw].z[0], res1);
        svst1_f32(pg, &res.c[raw].z[Nbatch], res2);
    }
}

template <>
void UDaggMultV<float, 16, KUQCD_GAUGE_RECON_NO>(const float* const u, const float* const src, float* res, float c)
{
    int Nbatch = 16;
    svbool_t pg = svptrue_b32();
    svfloat32_t src1, src2, tmp1, tmp2, res1, res2, u1, u2, tmpc;
    tmpc = svdup_n_f32(c);
    #pragma unroll
    for (int m = 0; m < 3; ++m) {
        tmp1 = svdup_n_f32(0);
        tmp2 = svdup_n_f32(0);
        #pragma unroll
        for (int n = 0; n < 3; ++n) {
            src1 = svld1_f32(pg, src+2*n*Nbatch);
            src2 = svld1_f32(pg, src+(1+2*n)*Nbatch);
            u1 = svld1_f32(pg, u+2*(m+3*n)*Nbatch);
            u2 = svld1_f32(pg, u+(1+2*(m+3*n))*Nbatch);
            tmp1 = svmla_f32_m(pg, svmla_f32_m(pg, tmp1, u2, src2), u1, src1);
            tmp2 = svmla_f32_m(pg, svmls_f32_m(pg, tmp2, u2, src1), u1, src2);
        }
        res1 = svld1_f32(pg, res+2*m*Nbatch);
        res2 = svld1_f32(pg, res+(1+2*m)*Nbatch);
        res1 = svmls_f32_m(pg, res1, tmp1, tmpc);
        res2 = svmls_f32_m(pg, res2, tmp2, tmpc);
        svst1_f32(pg, res+2*m*Nbatch, res1);
        svst1_f32(pg, res+(1+2*m)*Nbatch, res2);
    }
}

#endif

#endif

}