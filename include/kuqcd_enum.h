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
#ifndef KUQCD_ENUM_H
#define KUQCD_ENUM_H
#include <array>
typedef enum {
    KUQCD_SUCCESS = 0,                    /**< Success */
    KUQCD_ERR_INVALID_PARAM = -1,         /**< The parameters is invalid */
    KUQCD_ERR_UNSUPPORTED = -2,           /**< The function invoked is not yet supported */
    KUQCD_ERR_NULL_POINTER = -3,          /**< NULL Pointer */
    KUQCD_ERR_NO_MEMORY = -4,             /**< Failed to allocate memory */
    KUQCD_ERR_UNKNOWN_TYPE = -5,          /**< Unknown types, such as undefined matrix store types */
    KUQCD_ERR_UNEXPECTED_CALL_ORDER = -6, /**< Unexpected call order */
    KUQCD_ERR_NO_RESOURCE = -7,           /**< Resource is not ready */
    KUQCD_ERR_UNKNOWN_ERROR = -8,         /**< Exception occurs*/
} KuQCDStatus;

typedef enum {
    KUQCD_READ,
    KUQCD_WRITE,
    KUQCD_READWRITE
} KuQCDIOMode;

typedef enum {
    KUQCD_EVEN = 0,
    KUQCD_ODD = 1,
    KUQCD_ALL = 2,
    KUQCD_INVALID_PARITY
} KuQCDSiteSubset;

typedef enum {
    KUQCD_HALF = 2,
    KUQCD_SINGLE = 4,
    KUQCD_DOUBLE = 8,
    KUQCD_INVALID_PRECISION
} KuQCDPrecision;

typedef enum {
    KUQCD_HYPERPLANE = 1,
    KUQCD_PLANE = 2,
    KUQCD_STRIPE = 4,
    KUQCD_CORNER = 8,
    KUQCD_COMMALL = 15,
    KUQCD_INVALID_HALO
} KuQCDHalo;

typedef enum {
    KUQCD_XYZT_ORDER = 0,
    KUQCD_TZYX_ORDER = 1
} KuQCDRankOrder;

typedef enum {
    KUQCD_NORMAL_ORDER = 0,
    KUQCD_BATCH_ORDER = 1,
    KUQCD_INVALID_ORDER
} KuQCDGaugeOrder;
typedef enum {
    KUQCD_EVEN_EVEN = 0,
    KUQCD_ODD_ODD = 1
} KuQCDMatPCType;
typedef enum {
    KUQCD_WILSON_LINK = 0,
    KUQCD_FAT_LINK = 1,
    KUQCD_LONG_LINK = 2
} KuQCDLinkType;

typedef enum {
    KUQCD_GAUGE_RECON_NO = 18,
    KUQCD_GAUGE_RECON_U14 = 14,
    KUQCD_GAUGE_RECON_12 = 12,
    KUQCD_GAUGE_RECON_8 = 8,
    KUQCD_INVALID_RECON = 0
} KuQCDReconstructType;

typedef enum {
    KUQCD_WILSON_DIRAC,
    KUQCD_WILSONPC_DIRAC,
    KUQCD_HISQ_DIRAC,
    KUQCD_HISQPC_DIRAC,
    KUQCD_INVALID_DIRAC
} KuQCDDslashType;

typedef enum {
    KUQCD_CG,
    KUQCD_CGM,
    KUQCD_BICGSTAB,
    KUQCD_BICGSTABM,
    KUQCD_INVALID_INVERTER
} KuQCDInverterType;

typedef struct {
    unsigned int lat[4];
    int halo_depth;
    KuQCDRankOrder od;
    KuQCDSiteSubset site_subset;
    KuQCDGaugeOrder gauge_order;
    KuQCDPrecision prec;
    KuQCDLinkType link_type;
    KuQCDReconstructType recon;
} KuQCDGaugeParam;

typedef struct {
    KuQCDDslashType dslash_type;
    KuQCDInverterType inverter_type;
    KuQCDMatPCType matpc_type;
    KuQCDPrecision prec;
    int num_shift;
    int max_iter;
    int iter;
    double tol;
    double *shifts;
} KuQCDInvertParam;

#endif
