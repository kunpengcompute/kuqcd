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

#include "su3.h"

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon>::SU3(const SU3<T, Nbatch, Recon> &other)
{
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i] = other.m[i];
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> &SU3<T, Nbatch, Recon>::operator=(const SU3<T, Nbatch, Recon> &other)
{
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i] = other.m[i];
    }
    return *this;
}

// Matrix operations
template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> &SU3<T, Nbatch, Recon>::operator+=(const SU3<T, Nbatch, Recon> &other)
{
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i] += other.m[i];
    }
    return *this;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> &SU3<T, Nbatch, Recon>::operator-=(const SU3<T, Nbatch, Recon> &other)
{
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i] -= other.m[i];
    }
    return *this;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> &SU3<T, Nbatch, Recon>::operator*=(const SU3<T, Nbatch, Recon> &other)
{
    // Full 3x3 matrix multiplication
    Complex<T, Nbatch> res[9];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            res[i * 3 + j] = Complex<T, Nbatch>();
            for (int k = 0; k < 3; ++k) {
                res[i * 3 + j] += m[i * 3 + k] * other.m[k * 3 + j];
            }
        }
    }
    for (int i = 0; i < 9; ++i) {
        m[i] = res[i];
    }
    return *this;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> &SU3<T, Nbatch, Recon>::operator*=(const Complex<T, Nbatch> &scalar)
{
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i] *= scalar;
    }
    return *this;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> &SU3<T, Nbatch, Recon>::operator*=(T scalar)
{
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i] *= scalar;
    }
    return *this;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> &SU3<T, Nbatch, Recon>::operator/=(T scalar)
{
    scalar = 1.0 / scalar;
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i] *= scalar;
    }
    return *this;
}

// Unary operators
template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> SU3<T, Nbatch, Recon>::operator-() const
{
    SU3<T, Nbatch, Recon> result;
    for (int i = 0; i < (Recon >> 1); ++i) {
        result.m[i] = -m[i];
    }
    return result;
}

// SU(3) group operations
template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::Conj()
{
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i].Conj();
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::Trans()
{
    for (int i = 0; i < 2; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            Complex<T, Nbatch> tmp = m[j + 3 * i];
            m[j + 3 * i] = m[i + 3 * j];
            m[i + 3 * j] = tmp;
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::Dagger()
{
    this->Trans();
    for (int i = 0; i < (Recon >> 1); ++i) {
        m[i].Conj();
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
Complex<T, Nbatch> SU3<T, Nbatch, Recon>::Trace() const
{
    Rank0Fprintf(stderr, "su3 trace is under development! %s %d\n", __FILE__, __LINE__);
    return Complex<T, Nbatch>();
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
T SU3<T, Nbatch, Recon>::RealTrace() const
{
    Rank0Fprintf(stderr, "su3 real trace is under development! %s %d\n", __FILE__, __LINE__);
    return 0.0;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::Identity()
{
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j)
                m[j + 3 * i] = Complex<T, Nbatch>(1.0);
            else
                m[j + 3 * i] = Complex<T, Nbatch>(0.0);
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::Unitarize()
{
    for (int j = 0; j < Nbatch; ++j)
        Unitarize(j);
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::Unitarize(int j)
{
    double quadnorm;
    double invnorm;
    double Cre;
    double Cim;

    quadnorm = m[0].z[j] * m[0].z[j] + m[0].z[j + Nbatch] * m[0].z[j + Nbatch] + m[1].z[j] * m[1].z[j] +
               m[1].z[j + Nbatch] * m[1].z[j + Nbatch] + m[2].z[j] * m[2].z[j] +
               m[2].z[j + Nbatch] * m[2].z[j + Nbatch];

    invnorm = 1.0 / sqrt(quadnorm);

    m[0].z[j] *= invnorm;
    m[0].z[j + Nbatch] *= invnorm;
    m[1].z[j] *= invnorm;
    m[1].z[j + Nbatch] *= invnorm;
    m[2].z[j] *= invnorm;
    m[2].z[j + Nbatch] *= invnorm;

    Cre = m[3].z[j] * m[0].z[j] + m[3].z[j + Nbatch] * m[0].z[j + Nbatch] + m[4].z[j] * m[1].z[j] +
          m[4].z[j + Nbatch] * m[1].z[j + Nbatch] + m[5].z[j] * m[2].z[j] + m[5].z[j + Nbatch] * m[2].z[j + Nbatch];

    Cim = m[3].z[j + Nbatch] * m[0].z[j] - m[3].z[j] * m[0].z[j + Nbatch] + m[4].z[j + Nbatch] * m[1].z[j] -
          m[4].z[j] * m[1].z[j + Nbatch] + m[5].z[j + Nbatch] * m[2].z[j] - m[5].z[j] * m[2].z[j + Nbatch];

    m[3].z[j] -= Cre * m[0].z[j] - Cim * m[0].z[j + Nbatch];
    m[3].z[j + Nbatch] -= Cre * m[0].z[j + Nbatch] + Cim * m[0].z[j];
    m[4].z[j] -= Cre * m[1].z[j] - Cim * m[1].z[j + Nbatch];
    m[4].z[j + Nbatch] -= Cre * m[1].z[j + Nbatch] + Cim * m[1].z[j];
    m[5].z[j] -= Cre * m[2].z[j] - Cim * m[2].z[j + Nbatch];
    m[5].z[j + Nbatch] -= Cre * m[2].z[j + Nbatch] + Cim * m[2].z[j];

    quadnorm = m[3].z[j] * m[3].z[j] + m[3].z[j + Nbatch] * m[3].z[j + Nbatch] + m[4].z[j] * m[4].z[j] +
               m[4].z[j + Nbatch] * m[4].z[j + Nbatch] + m[5].z[j] * m[5].z[j] +
               m[5].z[j + Nbatch] * m[5].z[j + Nbatch];

    invnorm = 1.0 / sqrt(quadnorm);

    m[3].z[j] *= invnorm;
    m[3].z[j + Nbatch] *= invnorm;
    m[4].z[j] *= invnorm;
    m[4].z[j + Nbatch] *= invnorm;
    m[5].z[j] *= invnorm;
    m[5].z[j + Nbatch] *= invnorm;

    if constexpr (Recon == KUQCD_GAUGE_RECON_NO) {
        m[6].z[j] = (m[1].z[j] * m[5].z[j] - m[1].z[j + Nbatch] * m[5].z[j + Nbatch]) -
                    (m[2].z[j] * m[4].z[j] - m[2].z[j + Nbatch] * m[4].z[j + Nbatch]);
        m[6].z[j + Nbatch] = -(m[1].z[j + Nbatch] * m[5].z[j] + m[1].z[j] * m[5].z[j + Nbatch]) +
                             (m[2].z[j + Nbatch] * m[4].z[j] + m[2].z[j] * m[4].z[j + Nbatch]);

        m[7].z[j] = (m[2].z[j] * m[3].z[j] - m[2].z[j + Nbatch] * m[3].z[j + Nbatch]) -
                    (m[0].z[j] * m[5].z[j] - m[0].z[j + Nbatch] * m[5].z[j + Nbatch]);
        m[7].z[j + Nbatch] = -(m[2].z[j + Nbatch] * m[3].z[j] + m[2].z[j] * m[3].z[j + Nbatch]) +
                             (m[0].z[j + Nbatch] * m[5].z[j] + m[0].z[j] * m[5].z[j + Nbatch]);

        m[8].z[j] = (m[0].z[j] * m[4].z[j] - m[0].z[j + Nbatch] * m[4].z[j + Nbatch]) -
                    (m[1].z[j] * m[3].z[j] - m[1].z[j + Nbatch] * m[3].z[j + Nbatch]);
        m[8].z[j + Nbatch] = -(m[0].z[j + Nbatch] * m[4].z[j] + m[0].z[j] * m[4].z[j + Nbatch]) +
                             (m[1].z[j + Nbatch] * m[3].z[j] + m[1].z[j] * m[3].z[j + Nbatch]);
    } else if constexpr (Recon == KUQCD_GAUGE_RECON_U14) {
        m[6].z[j] = 1.0;
        m[6].z[j + Nbatch] = 0.0;
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
T SU3<T, Nbatch, Recon>::UnitarityDeviation() const
{
    return 0.0;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO> SU3<T, Nbatch, Recon>::Reconstruct() const
{
    return SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO>();
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::Print() const
{
    for (int i = 0; i < Nbatch; ++i) {
        for (int j = 0; j < (Recon >> 1); ++j) {
            printf("%13.6e %13.6e  ", m[j].z[i], m[j].z[i + Nbatch]);
            if (j % 3 == 2) {
                printf("\n");
            }
        }
        printf("\n");
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::PrintBatch(int i) const
{
    for (int j = 0; j < (Recon >> 1); ++j) {
        Rank0Printf("%13.6e %13.6e  ", m[j].z[i], m[j].z[i + Nbatch]);
        if (j % 3 == 2) {
            Rank0Printf("\n");
        }
    }
    Rank0Printf("\n");
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void SU3<T, Nbatch, Recon>::RankPrintBatch(int rank, int i) const
{
    if (GetRank() == 0) {
        for (int j = 0; j < (Recon >> 1); ++j) {
            printf("%13.6e %13.6e  ", m[j].z[i], m[j].z[i + Nbatch]);
            if (j % 3 == 2) {
                printf("\n");
            }
        }
        printf("\n");
    }
}

// Non-member binary operators
template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator+(SU3<T, Nbatch, Recon> lhs, const SU3<T, Nbatch, Recon> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator-(SU3<T, Nbatch, Recon> lhs, const SU3<T, Nbatch, Recon> &rhs)
{
    lhs -= rhs;
    return lhs;
}

// Scalar multiplication
template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(T scalar, SU3<T, Nbatch, Recon> rhs)
{
    rhs *= scalar;
    return rhs;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon> lhs, T scalar)
{
    lhs *= scalar;
    return lhs;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(const Complex<T, Nbatch> &scalar, SU3<T, Nbatch, Recon> rhs)
{
    rhs *= scalar;
    return rhs;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon> lhs, const Complex<T, Nbatch> &scalar)
{
    lhs *= scalar;
    return lhs;
}

// Stream output
template <typename T, int Nbatch, KuQCDReconstructType Recon>
std::ostream &operator<<(std::ostream &os, const SU3<T, Nbatch, Recon> &matrix)
{
    os << "SU3<" << Recon << ">(\n";

    for (int i = 0; i < 3; ++i) {
        os << "  [";
        for (int j = 0; j < 3; ++j) {
            os << matrix.m[j + 3 * i];
            if (j < 2)
                os << ", ";
        }
        os << "]\n";
    }
    os << ")";

    return os;
}

// Explicit instantiations for specific template parameters
#define INSTANTIATE_SU3_TEMPLATE(T, Nbatch, Recon)                                                         \
    template class SU3<T, Nbatch, Recon>;                                                                  \
    template SU3<T, Nbatch, Recon> operator+(SU3<T, Nbatch, Recon> lhs, const SU3<T, Nbatch, Recon> &rhs); \
    template SU3<T, Nbatch, Recon> operator-(SU3<T, Nbatch, Recon> lhs, const SU3<T, Nbatch, Recon> &rhs); \
    template SU3<T, Nbatch, Recon> operator*(T scalar, SU3<T, Nbatch, Recon> rhs);                         \
    template SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon> lhs, T scalar);                         \
    template SU3<T, Nbatch, Recon> operator*(const Complex<T, Nbatch> &scalar, SU3<T, Nbatch, Recon> rhs); \
    template SU3<T, Nbatch, Recon> operator*(SU3<T, Nbatch, Recon> lhs, const Complex<T, Nbatch> &scalar); \
    template std::ostream &operator<<(std::ostream &os, const SU3<T, Nbatch, Recon> &matrix)
// Instantiate the commonly used template parameter combinations
INSTANTIATE_SU3_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_NO);
INSTANTIATE_SU3_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_NO);
INSTANTIATE_SU3_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_NO);
INSTANTIATE_SU3_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_NO);

INSTANTIATE_SU3_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_12);
INSTANTIATE_SU3_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_12);
INSTANTIATE_SU3_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_12);
INSTANTIATE_SU3_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_12);

INSTANTIATE_SU3_TEMPLATE(float, 1, KUQCD_GAUGE_RECON_U14);
INSTANTIATE_SU3_TEMPLATE(double, 1, KUQCD_GAUGE_RECON_U14);
INSTANTIATE_SU3_TEMPLATE(float, 16, KUQCD_GAUGE_RECON_U14);
INSTANTIATE_SU3_TEMPLATE(double, 8, KUQCD_GAUGE_RECON_U14);

} // namespace kuqcd
