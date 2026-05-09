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
#include <mpi.h>
#include <omp.h>
#include "spinor_field.h"

#ifdef SVEBITS
#include "arm_sve.h"
#endif

namespace kuqcd {

template <typename T, int Nbatch, int Ns>
SpinorField<T, Nbatch, Ns>::SpinorField(const SpinorField<T, Nbatch, Ns> &spinor) : LatticeField<SpinorT>(spinor)
{
}

template <typename T, int Nbatch, int Ns>
SpinorField<T, Nbatch, Ns>::SpinorField(LattBase &latBase) : LatticeField<SpinorT>(latBase, Ns * Nbatch * 6)
{
}

template <typename T, int Nbatch, int Ns>
ColorSpinor<T, Nbatch> SpinorField<T, Nbatch, Ns>::GetSpinorCopy(size_t i, int ns) const
{
    return this->data[ns + Ns * i];
}

template <typename T, int Nbatch, int Ns>
ColorSpinor<T, Nbatch> &SpinorField<T, Nbatch, Ns>::GetSpinorRef(size_t i, int ns)
{
    return this->data[ns + Ns * i];
}

template <typename T, int Nbatch, int Ns>
void SpinorField<T, Nbatch, Ns>::RandGauss(uint4 *randState)
{
#pragma omp parallel for
    for (size_t i = 0; i < this->latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int j = 0; j < Nbatch; ++j) {
                ColorSpinor<T, 1> tmp;
                tmp.RandGauss(&randState[j + Nbatch * i]);
                for (int nc = 0; nc < 3; ++nc) {
                    this->data[ns + Ns * i].c[nc].z[j] = tmp.c[nc].z[0];
                    this->data[ns + Ns * i].c[nc].z[j + Nbatch] = tmp.c[nc].z[1];
                }
            }
        }
    }
}

template <typename T, int Nbatch, int Ns>
SpinorField<T, Nbatch, Ns> operator-(const SpinorField<T, Nbatch, Ns> &vec1, const SpinorField<T, Nbatch, Ns> &vec2)
{
    SpinorField<T, Nbatch, Ns> res(vec1);
#pragma omp parallel for
    for (size_t i = 0; i < res.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < (Nbatch << 1); ++j) {
                    res.data[ns + Ns * i].c[nc].z[j] -= vec2.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }
    return res;
}

template <typename T, int Nbatch, int Ns>
SpinorField<T, Nbatch, Ns> SpinorField<T, Nbatch, Ns>::operator*(T scalar)
{
    SpinorField<T, Nbatch, Ns> res(*this);
#pragma omp parallel for
    for (size_t i = 0; i < res.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < (Nbatch << 1); ++j) {
                    res.data[ns + Ns * i].c[nc].z[j] = this->data[ns + Ns * i].c[nc].z[j] * scalar;
                }
            }
        }
    }
    return res;
}

template <typename T, int Nbatch, int Ns>
SpinorField<T, Nbatch, Ns> &SpinorField<T, Nbatch, Ns>::operator-=(const SpinorField<T, Nbatch, Ns> &vec1)
{
    for (size_t i = 0; i < vec1.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < (Nbatch << 1); ++j) {
                    this->data[ns + Ns * i].c[nc].z[j] -= vec1.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }
    return *this;
}

template <typename T, int Nbatch, int Ns>
SpinorField<T, Nbatch, Ns> SpinorField<T, Nbatch, Ns>::operator+(const SpinorField<T, Nbatch, Ns> &vec1)
{
    SpinorField<T, Nbatch, Ns> res(*this);

    for (size_t i = 0; i < vec1.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < (Nbatch << 1); ++j) {
                    res.data[ns + Ns * i].c[nc].z[j] =
                        this->data[ns + Ns * i].c[nc].z[j] + vec1.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }
    return res;
}

template <typename T, int Nbatch, int Ns>
void Axmy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a)
{
#pragma omp parallel for
    for (size_t i = 0; i < y.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < 2 * Nbatch; ++j) {
                    y.data[ns + Ns * i].c[nc].z[j] =
                        a * x.data[ns + Ns * i].c[nc].z[j] - y.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }
}

template <typename T, int Nbatch, int Ns>
void Axpy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a)
{
#pragma omp parallel for
    for (size_t i = 0; i < y.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < 2 * Nbatch; ++j) {
                    y.data[ns + Ns * i].c[nc].z[j] += a * x.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }
}

template <typename T, int Nbatch, int Ns>
void Ymax(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a)
{
#pragma omp parallel for
    for (size_t i = 0; i < y.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < 2 * Nbatch; ++j) {
                    y.data[ns + Ns * i].c[nc].z[j] -= a * x.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }
}

template <typename T, int Nbatch, int Ns>
void Axpby(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a, T b)
{
#pragma omp parallel for
    for (size_t i = 0; i < y.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < 2 * Nbatch; ++j) {
                    y.data[ns + Ns * i].c[nc].z[j] =
                        a * x.data[ns + Ns * i].c[nc].z[j] + b * y.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }
}

template <typename T, int Nbatch, int Ns>
T RealDot(const SpinorField<T, Nbatch, Ns> &bra, const SpinorField<T, Nbatch, Ns> &ket)
{
    T res = 0;
#pragma omp parallel for reduction(+ : res)
    for (size_t i = 0; i < bra.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < (Nbatch << 1); ++j) {
                    res += bra.data[ns + Ns * i].c[nc].z[j] * ket.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }

    if (GetNp() > 1) {
        T tmp = 0.0;
        if (sizeof(T) == sizeof(double)) {
            MPI_Allreduce(&res, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        } else {
            MPI_Allreduce(&res, &tmp, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        }
        res = tmp;
    }

    return res;
}

template <typename T, int Nbatch, int Ns>
T Norm2(const SpinorField<T, Nbatch, Ns> &vec)
{
    T res = 0;
#pragma omp parallel for reduction(+ : res)
    for (size_t i = 0; i < vec.latBase.lenBulk; ++i) {
        for (int ns = 0; ns < Ns; ++ns) {
            for (int nc = 0; nc < 3; ++nc) {
                for (int j = 0; j < (Nbatch << 1); ++j) {
                    res += vec.data[ns + Ns * i].c[nc].z[j] * vec.data[ns + Ns * i].c[nc].z[j];
                }
            }
        }
    }

    if (GetNp() > 1) {
        T tmp = 0.0;
        if (sizeof(T) == sizeof(double)) {
            MPI_Allreduce(&res, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        } else {
            MPI_Allreduce(&res, &tmp, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        }
        res = tmp;
    }

    return res;
}

template <typename T, int Nbatch, int Ns>
T Norm(const SpinorField<T, Nbatch, Ns> &vec)
{
    return sqrt(Norm2<T, Nbatch, Ns>(vec));
}

#ifdef SVEBITS

#if SVEBITS == 512

template <>
void Ymax<float, 16, 1>(SpinorField<float, 16, 1> &y, const SpinorField<float, 16, 1> &x, float a)
{
    int Nbatch = 16;
#pragma omp parallel for
    for (size_t i = 0; i < y.latBase.lenBulk; ++i) {
        for (int nc = 0; nc < 3; ++nc) {
            svfloat32_t tmpy, tmpx;
            svfloat32_t tmpa;
            svbool_t pg = svptrue_b32();
            tmpa = svdup_n_f32(a);
            tmpx = svld1_f32(pg, &x.data[i].c[nc].z[0]);
            tmpy = svld1_f32(pg, &y.data[i].c[nc].z[0]);
            tmpy = svmls_f32_m(pg, tmpy, tmpa, tmpx);
            svst1_f32(pg, &y.data[i].c[nc].z[0], tmpy);

            tmpx = svld1_f32(pg, &x.data[i].c[nc].z[Nbatch]);
            tmpy = svld1_f32(pg, &y.data[i].c[nc].z[Nbatch]);
            tmpy = svmls_f32_m(pg, tmpy, tmpa, tmpx);
            svst1_f32(pg, &y.data[i].c[nc].z[Nbatch], tmpy);
        }
    }
}

template <>
void Axpby<float, 16, 1>(SpinorField<float, 16, 1> &y, const SpinorField<float, 16, 1> &x, float a, float b)
{
    int Nbatch = 16;

#pragma omp parallel for
    for (size_t i = 0; i < y.latBase.lenBulk; ++i) {
        for (int nc = 0; nc < 3; ++nc) {
            svfloat32_t tmpy, tmpx, tmp;
            svfloat32_t tmpa, tmpb;
            svbool_t pg = svptrue_b32();
            tmpa = svdup_n_f32(a);
            tmpb = svdup_n_f32(b);
            tmp = svdup_n_f32(0);
            tmpx = svld1_f32(pg, &x.data[i].c[nc].z[0]);
            tmp = svmla_f32_m(pg, tmp, tmpx, tmpa);
            tmpy = svld1_f32(pg, &y.data[i].c[nc].z[0]);
            tmpy = svmla_f32_m(pg, tmp, tmpy, tmpb);
            svst1_f32(pg, &y.data[i].c[nc].z[0], tmpy);

            tmp = svdup_n_f32(0);
            tmpx = svld1_f32(pg, &x.data[i].c[nc].z[Nbatch]);
            tmp = svmla_f32_m(pg, tmp, tmpx, tmpa);
            tmpy = svld1_f32(pg, &y.data[i].c[nc].z[Nbatch]);
            tmpy = svmla_f32_m(pg, tmp, tmpy, tmpb);
            svst1_f32(pg, &y.data[i].c[nc].z[Nbatch], tmpy);
        }
    }
}
#endif

#endif

#define INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(T, Nbatch, Ns)                                          \
    template class SpinorField<T, Nbatch, Ns>;                                                         \
    template SpinorField<T, Nbatch, Ns> operator-(const SpinorField<T, Nbatch, Ns> &,                  \
                                                  const SpinorField<T, Nbatch, Ns> &);                 \
    template T Norm2(const SpinorField<T, Nbatch, Ns> &);                                              \
    template T Norm(const SpinorField<T, Nbatch, Ns> &);                                               \
    template void Axmy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);       \
    template void Axpy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);       \
    template void Ymax(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);       \
    template void Axpby(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a, T b); \
    template T RealDot(const SpinorField<T, Nbatch, Ns> &bra, const SpinorField<T, Nbatch, Ns> &ket)
INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(float, 1, 1);
INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(double, 1, 1);
INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(float, 16, 1);
INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(double, 8, 1);

INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(float, 1, 4);
INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(double, 1, 4);
INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(float, 16, 4);
INSTANTIATE_TEMPLATE_CLASS_SPINORFIELD(double, 8, 4);

} // namespace kuqcd
