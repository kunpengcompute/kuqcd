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
#ifndef CGM_H
#define CGM_H

#include <vector>
#include "inverter.h"

namespace kuqcd {

template <typename T, typename VecT>
class CGM : public Inverter<T, VecT> {
public:
    CGM(KuQCDInverterType it, Dirac<VecT> *dirac_) : Inverter<T, VecT>(it, dirac_) {}
    ~CGM() override {};
    int Invert(int noShifts, T shifts[], const VecT & src, VecT *sol[], const double tol, int maxIters) override;
};

} // namespace kuqcd

#endif