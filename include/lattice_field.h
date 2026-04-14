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
#ifndef LATTICE_FIELD_H
#define LATTICE_FIELD_H


#include "space_time.h"
#include "kuqcd_memory.h"
#include "latt_base.h"
#include "spinor.h"
#include "su3.h"

namespace kuqcd {

template <typename FieldT>
class LatticeField : public KMemory {
public:
    LattBase &latBase;
    FieldT *data = nullptr;
    int locSize;
    size_t length;

    LatticeField() = delete;
    LatticeField(const LatticeField<FieldT> &lf);
    LatticeField(LattBase &latBase_, const int locSize_);
    LatticeField<FieldT> &operator=(const LatticeField<FieldT> &lf);
    virtual ~LatticeField()
    {
        data = nullptr;
    }
};

#define DECLARE_EXTERN_LATTICE_FIELD_TEMPLATE(T, Nbatch)                      \
    extern template class LatticeField<ColorSpinor<T, Nbatch>>;               \
    extern template class LatticeField<SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO>>; \
    extern template class LatticeField<SU3<T, Nbatch, KUQCD_GAUGE_RECON_12>>; \
    extern template class LatticeField<SU3<T, Nbatch, KUQCD_GAUGE_RECON_U14>>
DECLARE_EXTERN_LATTICE_FIELD_TEMPLATE(float, 1);
DECLARE_EXTERN_LATTICE_FIELD_TEMPLATE(double, 1);
DECLARE_EXTERN_LATTICE_FIELD_TEMPLATE(float, 16);
DECLARE_EXTERN_LATTICE_FIELD_TEMPLATE(double, 8);

extern template class LatticeField<uint4>;

}  // namespace kuqcd

#endif