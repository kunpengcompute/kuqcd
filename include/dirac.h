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
#ifndef DIRAC_H
#define DIRAC_H

#include "kuqcd_enum.h"

namespace kuqcd {

template <typename VecT>
class Dirac {
public:
    KuQCDDslashType diracType = KUQCD_INVALID_DIRAC;
    Dirac(KuQCDDslashType dt) : diracType(dt) {}
    virtual ~Dirac() {};
    virtual void Dslash(const VecT &src, VecT &dst, bool update) = 0;
    virtual void MdaggM(const VecT &src, VecT &dst, bool update) = 0;
};

} // namespace kuqcd

#endif
