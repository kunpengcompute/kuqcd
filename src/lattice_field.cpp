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
#include "lattice_field.h"
#include "omp.h"
#include "index_utils.h"

#ifdef USE_KUPL
#include "kupl.h"
#endif

namespace kuqcd {

template <typename FieldT>
LatticeField<FieldT>::LatticeField(const LatticeField<FieldT> &lf)
    : KMemory(lf), latBase(lf.latBase), locSize(lf.locSize), length(latBase.lenFull)
{
    data = static_cast<FieldT *>(this->rawPtr);
}

template <typename FieldT>
LatticeField<FieldT>::LatticeField(LattBase &latBase_, const int locSize_)
    : KMemory(locSize_ * latBase_.lenFull * latBase_.prec),
      latBase(latBase_),
      locSize(locSize_),
      length(latBase_.lenFull)
{
    data = static_cast<FieldT *>(this->rawPtr);
}

template <typename FieldT>
LatticeField<FieldT> &LatticeField<FieldT>::operator=(const LatticeField<FieldT> &lf)
{
    if ((locSize != lf.locSize) || (length != lf.length) || (this->size != lf.size)) {
        Rank0Fprintf(stderr, "lattice field length does not match! %s %d\n", __FILE__, __LINE__);
        Rank0Fprintf(stderr, "this: %d %zu %zu\n", locSize, length, this->size);
        Rank0Fprintf(stderr, "that: %d %zu %zu\n", lf.locSize, lf.length, lf.size);
    }
    Copy(lf);
    return *this;
}

#define INSTANTIATE_LATTICE_FIELD_TEMPLATE(T, Nbatch)                  \
    template class LatticeField<ColorSpinor<T, Nbatch>>;               \
    template class LatticeField<SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO>>; \
    template class LatticeField<SU3<T, Nbatch, KUQCD_GAUGE_RECON_12>>; \
    template class LatticeField<SU3<T, Nbatch, KUQCD_GAUGE_RECON_U14>>
INSTANTIATE_LATTICE_FIELD_TEMPLATE(float, 1);
INSTANTIATE_LATTICE_FIELD_TEMPLATE(double, 1);
INSTANTIATE_LATTICE_FIELD_TEMPLATE(float, 16);
INSTANTIATE_LATTICE_FIELD_TEMPLATE(double, 8);

template class LatticeField<uint4>;

} // namespace kuqcd
