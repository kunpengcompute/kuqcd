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
#ifndef SPINOR_FIELD_H
#define SPINOR_FIELD_H

#include "lattice_field.h"

namespace kuqcd {

template <typename T, int Nbatch, int Ns>
class SpinorField : public LatticeField<ColorSpinor<T, Nbatch>> {
public:
    typedef ColorSpinor<T, Nbatch> SpinorT;
    SpinorField() = delete;
    SpinorField(const SpinorField<T, Nbatch, Ns> &spinor);
    SpinorField(LattBase &latBase);
    SpinorT GetSpinorCopy(size_t i, int ns) const;
    SpinorT &GetSpinorRef(size_t, int ns);
    void RandGauss(uint4 *randState);
    inline int GetNbatch()
    {
        return Nbatch;
    }

    // SpinorField<T, Nbatch, Ns> &operator=(const SpinorField<T, Nbatch, Ns> &vec1);
    SpinorField<T, Nbatch, Ns> operator*(T scalar);
    SpinorField<T, Nbatch, Ns> &operator-=(const SpinorField<T, Nbatch, Ns> &vec1);
    SpinorField<T, Nbatch, Ns> operator+(const SpinorField<T, Nbatch, Ns> &vec1);
};

template <typename T, int Nbatch, int Ns>
SpinorField<T, Nbatch, Ns> operator-(const SpinorField<T, Nbatch, Ns> &vec1, const SpinorField<T, Nbatch, Ns> &vec2);

template <typename T, int Nbatch, int Ns>
T Norm2(const SpinorField<T, Nbatch, Ns> &vec);

template <typename T, int Nbatch, int Ns>
void Axmy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);

template <typename T, int Nbatch, int Ns>
void Axpy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);

template <typename T, int Nbatch, int Ns>
void Ymax(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);

template <typename T, int Nbatch, int Ns>
void Axpby(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a, T b);

template <typename T, int Nbatch, int Ns>
T RealDot(const SpinorField<T, Nbatch, Ns> &bra, const SpinorField<T, Nbatch, Ns> &ket);

template <typename T, int Nbatch, int Ns>
T Norm(const SpinorField<T, Nbatch, Ns> &vec);

#ifdef SVEBITS

#if SVEBITS == 512

template <>
void Ymax<float, 16, 1>(SpinorField<float, 16, 1> &y, const SpinorField<float, 16, 1> &x, float a);

template <>
void Axpby<float, 16, 1>(SpinorField<float, 16, 1> &y, const SpinorField<float, 16, 1> &x, float a, float b);
#endif

#endif

#define DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(T, Nbatch, Ns)                                              \
    extern template class SpinorField<T, Nbatch, Ns>;                                                         \
    extern template SpinorField<T, Nbatch, Ns> operator-(const SpinorField<T, Nbatch, Ns> &,                  \
                                                         const SpinorField<T, Nbatch, Ns> &);                 \
    extern template T Norm2(const SpinorField<T, Nbatch, Ns> &);                                              \
    extern template T Norm(const SpinorField<T, Nbatch, Ns> &);                                               \
    extern template void Axmy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);       \
    extern template void Axpy(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);       \
    extern template void Ymax(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a);       \
    extern template void Axpby(SpinorField<T, Nbatch, Ns> &y, const SpinorField<T, Nbatch, Ns> &x, T a, T b); \
    extern template T RealDot(const SpinorField<T, Nbatch, Ns> &bra, const SpinorField<T, Nbatch, Ns> &ket)
DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(float, 1, 1);
DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(double, 1, 1);
DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(float, 16, 1);
DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(double, 8, 1);

DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(float, 1, 4);
DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(double, 1, 4);
DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(float, 16, 4);
DECLARE_EXTERN_TEMPLATE_CLASS_SPINORFIELD(double, 8, 4);

} // namespace kuqcd
#endif
