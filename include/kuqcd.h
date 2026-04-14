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

#ifndef KUQCD_H
#define KUQCD_H

#include <omp.h>
#include <mpi.h>
#include <cstdio>
#include "dirac_hisq.h"
#include "index_utils.h"
#include "inverter_cgm.h"
#include "logger.h"
#include "kuqcd_enum.h"
#include "command_line.h"

namespace kuqcd {

KuQCDStatus kuqcd_init();
KuQCDStatus kuqcd_load_gaugefield(void *h_gauge, KuQCDGaugeParam *param);
KuQCDStatus kuqcd_solve_multishift(void **_hp_x, void *_hp_b, KuQCDInvertParam *param);

void kuqcd_gauge_free();
void kuqcd_end();
}  // namespace kuqcd
#endif