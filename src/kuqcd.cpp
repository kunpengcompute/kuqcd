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

#include "kuqcd.h"

namespace kuqcd {

LattBase *g_latbase_even_halo0 = nullptr;
LattBase *g_latbase_odd_halo0 = nullptr;
LattBase *g_latbase_all_halo0 = nullptr;
LattBase *g_latbase_all_halo3 = nullptr;
CommBase *g_comm = nullptr;

GaugeField<float, 1, KUQCD_GAUGE_RECON_NO> *g_gauge_smeared_f1 = nullptr;
GaugeField<float, 16, KUQCD_GAUGE_RECON_NO> *g_gauge_smeared_f16 = nullptr;
GaugeField<double, 1, KUQCD_GAUGE_RECON_NO> *g_gauge_smeared_d1 = nullptr;
GaugeField<double, 8, KUQCD_GAUGE_RECON_NO> *g_gauge_smeared_d8 = nullptr;

GaugeField<float, 1, KUQCD_GAUGE_RECON_NO> *g_gauge_naik_f1 = nullptr;
GaugeField<float, 16, KUQCD_GAUGE_RECON_NO> *g_gauge_naik_f16 = nullptr;
GaugeField<double, 1, KUQCD_GAUGE_RECON_NO> *g_gauge_naik_d1 = nullptr;
GaugeField<double, 8, KUQCD_GAUGE_RECON_NO> *g_gauge_naik_d8 = nullptr;

Inverter<float, SpinorField<float, 1, 1>> *g_inverter_f1 = nullptr;
Inverter<float, SpinorField<float, 16, 1>> *g_inverter_f16 = nullptr;
Inverter<double, SpinorField<double, 1, 1>> *g_inverter_d1 = nullptr;
Inverter<double, SpinorField<double, 8, 1>> *g_inverter_d8 = nullptr;

std::array<unsigned short, 4> g_grid_dim;

// set latbase to g_latbase_all_halo3,g_latbase_all_halo0,g_latbase_even_halo0,g_latbase_odd_halo0
KuQCDStatus set_latbase(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int nt, int halo_depth,
                        KuQCDRankOrder od, KuQCDSiteSubset site_subset, KuQCDGaugeOrder gauge_order,
                        KuQCDPrecision prec, std::array<unsigned short, 4> grid_dim, LattBase **latbase)
{
    if (gauge_order < KUQCD_NORMAL_ORDER || gauge_order >= KUQCD_INVALID_ORDER) {
        rootLogger.Error("gauge_order is not valid");
        return KUQCD_ERR_INVALID_PARAM;
    }

    if (prec != KUQCD_SINGLE && prec != KUQCD_DOUBLE) {
        if (prec == KUQCD_HALF) {
            rootLogger.Error("prec is not support");
            return KUQCD_ERR_UNSUPPORTED;
        }
        rootLogger.Error("prec is not valid");
        return KUQCD_ERR_INVALID_PARAM;
    }
    if (site_subset < KUQCD_EVEN || site_subset >= KUQCD_INVALID_PARITY) {
        rootLogger.Error("site_subset is not valid");
        return KUQCD_ERR_INVALID_PARAM;
    }
    if (halo_depth == 3 && site_subset == KUQCD_ALL) {
        if (g_latbase_all_halo3 == nullptr) {
            Lattice lat_sup(nx, ny, nz, nt);
            *latbase = new LattBase(lat_sup, grid_dim, halo_depth, gauge_order, od, prec, site_subset);
            g_latbase_all_halo3 = *latbase;
        } else {
            *latbase = g_latbase_all_halo3;
        }
    }
    if (halo_depth == 0 && site_subset == KUQCD_ALL) {
        if (g_latbase_all_halo0 == nullptr) {
            Lattice lat_sup(nx, ny, nz, nt);
            *latbase = new LattBase(lat_sup, grid_dim, halo_depth, gauge_order, od, prec, site_subset);
            g_latbase_all_halo0 = *latbase;
        } else {
            *latbase = g_latbase_all_halo0;
        }
    }
    if (halo_depth == 0 && site_subset == KUQCD_EVEN) {
        if (g_latbase_even_halo0 == nullptr) {
            Lattice lat_sup(nx, ny, nz, nt);
            *latbase = new LattBase(lat_sup, grid_dim, halo_depth, gauge_order, od, prec, site_subset);
            g_latbase_even_halo0 = *latbase;
            g_latbase_odd_halo0 = new LattBase(lat_sup, grid_dim, halo_depth, gauge_order, od, prec, KUQCD_ODD);
        } else {
            *latbase = g_latbase_even_halo0;
        }
    }
    if (halo_depth == 0 && site_subset == KUQCD_ODD) {
        if (g_latbase_odd_halo0 == nullptr) {
            Lattice lat_sup(nx, ny, nz, nt);
            *latbase = new LattBase(lat_sup, grid_dim, halo_depth, gauge_order, od, prec, site_subset);
            g_latbase_odd_halo0 = *latbase;
            g_latbase_even_halo0 = new LattBase(lat_sup, grid_dim, halo_depth, gauge_order, od, prec, KUQCD_EVEN);
        } else {
            *latbase = g_latbase_odd_halo0;
        }
    }
    return KUQCD_SUCCESS;
}

template <typename T, int Nbatch, KuQCDReconstructType ReconS>
KuQCDStatus kuqcd_gaugefield_set_values(GaugeField<T, Nbatch, ReconS> *gauge, T *values)
{
    if (gauge == nullptr || values == nullptr) {
        rootLogger.Error("gauge or values is nullptr");
        return KUQCD_ERR_NULL_POINTER;
    }
    size_t vol = gauge->latBase.lenBulk;

    int gaugeS = ReconS >> 1;
    for (size_t i = 0; i < vol; ++i) {
        size_t iFull;
        iFull = gauge->latBase.myIdx.idx[i];
        for (int dir = 0; dir < 4; ++dir) {
            for (int l = 0; l < gaugeS; ++l) {
                for (int j = 0; j < Nbatch; ++j) {
                    size_t ni = j + Nbatch * i;
                    gauge->data[dir + 4 * iFull].m[l].z[j] = values[0 + 2 * (l + gaugeS * (dir + 4 * ni))];
                    gauge->data[dir + 4 * iFull].m[l].z[j + Nbatch] = values[1 + 2 * (l + gaugeS * (dir + 4 * ni))];
                }
            }
        }
    }
    return KUQCD_SUCCESS;
}

template <typename T, int Nbatch, int Ns>
KuQCDStatus kuqcd_spinorfield_set_values(SpinorField<T, Nbatch, Ns> *spinor, T *values)
{
    if (spinor == nullptr || values == nullptr) {
        rootLogger.Error("spinor or values is nullptr");
        return KUQCD_ERR_NULL_POINTER;
    }
    size_t volh = spinor->latBase.lenBulk;
    for (size_t i = 0; i < volh; ++i) {
        size_t iFull = spinor->latBase.myIdx.idx[i];
        for (int m = 0; m < 3; ++m) {
            for (int j = 0; j < Nbatch; ++j) {
                size_t ni = j + Nbatch * i;
                spinor->data[iFull].c[m].z[j] = values[0 + 2 * (m + 3 * ni)];
                spinor->data[iFull].c[m].z[j + Nbatch] = values[1 + 2 * (m + 3 * ni)];
            }
        }
    }
    return KUQCD_SUCCESS;
}

template <typename T, int Nbatch, int Ns>
KuQCDStatus kuqcd_spinorfield_get_values(SpinorField<T, Nbatch, Ns> *spinor, T *values)
{
    if (spinor == nullptr || values == nullptr) {
        rootLogger.Error("spinor or values is nullptr");
        return KUQCD_ERR_NULL_POINTER;
    }
    size_t volh = spinor->latBase.lenBulk;
    for (size_t i = 0; i < volh; ++i) {
        size_t iFull = spinor->latBase.myIdx.idx[i];
        for (int m = 0; m < 3; ++m) {
            for (int j = 0; j < Nbatch; ++j) {
                size_t ni = j + Nbatch * i;
                values[0 + 2 * (m + 3 * ni)] = spinor->data[iFull].c[m].z[j];
                values[1 + 2 * (m + 3 * ni)] = spinor->data[iFull].c[m].z[j + Nbatch];
            }
        }
    }
    return KUQCD_SUCCESS;
}

template <typename T, int Nbatch, int Ns, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
KuQCDStatus kuqcd_inverter(Inverter<T, SpinorField<T, Nbatch, Ns>> **inverter,
                           GaugeField<T, Nbatch, ReconS> *gauge_smeared, GaugeField<T, Nbatch, ReconN> *gauge_naik,
                           void **_hp_x, void *_hp_b, KuQCDInvertParam &param)
{
    if (inverter == nullptr || gauge_smeared == nullptr || gauge_naik == nullptr) {
        rootLogger.Error("inverter or gauge_smeared or gauge_naik is nullptr");
        return KUQCD_ERR_NULL_POINTER;
    }
    if (param.inverter_type != KUQCD_CGM) {
        rootLogger.Error("Inverter is not support");
        return KUQCD_ERR_UNSUPPORTED;
    }
    if (param.dslash_type == KUQCD_WILSON_DIRAC || param.dslash_type == KUQCD_WILSONPC_DIRAC) {
        rootLogger.Error("dirac type should be hisq or hisq pc!");
        return KUQCD_ERR_UNSUPPORTED;
    }
    if (param.dslash_type != KUQCD_HISQ_DIRAC && param.dslash_type != KUQCD_HISQPC_DIRAC) {
        rootLogger.Error("dirac type is error!");
        return KUQCD_ERR_INVALID_PARAM;
    }
    if (param.matpc_type != KUQCD_EVEN_EVEN && param.matpc_type != KUQCD_ODD_ODD) {
        rootLogger.Error("matpc_type is error!");
        return KUQCD_ERR_INVALID_PARAM;
    }
    if (param.num_shift <= 0) {
        rootLogger.Error("num_shift < 0!");
        return KUQCD_ERR_INVALID_PARAM;
    }
    KuQCDSiteSubset src_parity = (param.matpc_type == KUQCD_EVEN_EVEN) ? KUQCD_EVEN : KUQCD_ODD;
    T *src_value = (T *)_hp_b;
    T *dst_value[param.num_shift];
    for (int i = 0; i < param.num_shift; i++) {
        dst_value[i] = (T *)(_hp_x[i]);
    }
    int nx = gauge_smeared->latBase.latSupBulkFL.nx * Nbatch;
    int ny = gauge_smeared->latBase.latSupBulkFL.ny;
    int nz = gauge_smeared->latBase.latSupBulkFL.nz;
    int nt = gauge_smeared->latBase.latSupBulkFL.nt;

    KuQCDGaugeOrder gauge_order = gauge_smeared->latBase.gauge_order;
    KuQCDRankOrder od = gauge_smeared->latBase.rankOrder;
    LattBase *latbase = nullptr;
    KuQCDStatus rtn = set_latbase(nx, ny, nz, nt, 0, od, src_parity, gauge_order, param.prec, g_grid_dim, &latbase);
    if (rtn != KUQCD_SUCCESS) {
        return rtn;
    }
    SpinorField<T, Nbatch, Ns> *src = new SpinorField<T, Nbatch, Ns>(*latbase);

    kuqcd_spinorfield_set_values<T, Nbatch, Ns>(src, src_value);
    SpinorField<T, Nbatch, Ns> *dst[param.num_shift];
    Dirac<SpinorField<T, Nbatch, Ns>> *dirac = nullptr;
    SpinorField<T, Nbatch, Ns> *tmp = nullptr;
    T *shift = nullptr;
    for (int i = 0; i < param.num_shift; ++i) {
        dst[i] = nullptr;
    }
    for (int i = 0; i < param.num_shift; ++i) {
        dst[i] = new (std::nothrow) SpinorField<T, Nbatch, Ns>(*latbase);
        if (dst[i] == nullptr) {
            rtn = KUQCD_ERR_NO_MEMORY;
            goto mem_release;
        }
    }

    if (g_latbase_even_halo0 == nullptr && g_latbase_odd_halo0 == nullptr) {
        rootLogger.Error("no init spinorfield");
        return KUQCD_ERR_UNEXPECTED_CALL_ORDER;
    }
    if (param.matpc_type == KUQCD_ODD_ODD) {
        tmp = new (std::nothrow) SpinorField<T, Nbatch, Ns>(*g_latbase_even_halo0);
        if (tmp == nullptr) {
            rtn = KUQCD_ERR_NO_MEMORY;
            goto mem_release;
        }
    } else {
        tmp = new (std::nothrow) SpinorField<T, Nbatch, Ns>(*g_latbase_odd_halo0);
        if (tmp == nullptr) {
            rtn = KUQCD_ERR_NO_MEMORY;
            goto mem_release;
        }
    }

    dirac = new (std::nothrow)
        DiracHISQ<T, Nbatch, ReconS, ReconN>(param.dslash_type, *gauge_smeared, *gauge_naik, 0, tmp);
    if (dirac == nullptr) {
        rtn = KUQCD_ERR_NO_MEMORY;
        goto mem_release;
    }
    *inverter = new (std::nothrow) CGM<T, SpinorField<T, Nbatch, Ns>>(param.inverter_type, dirac);
    if (*inverter == nullptr) {
        rtn = KUQCD_ERR_NO_MEMORY;
        goto mem_release;
    }
    shift = (T *)malloc(param.num_shift * sizeof(T));
    if (shift == nullptr) {
        rtn = KUQCD_ERR_NO_MEMORY;
        goto mem_release;
    }
    for (int i = 0; i < param.num_shift; i++) {
        shift[i] = (T)param.shifts[i];
    }
    param.iter = (*inverter)->Invert(param.num_shift, shift, *src, dst, param.tol, param.max_iter);

    for (int i = 0; i < param.num_shift; i++) {
        kuqcd_spinorfield_get_values<T, Nbatch, Ns>(dst[i], dst_value[i]);
    }
mem_release:
    delete src;
    delete tmp;
    for (int i = 0; i < param.num_shift; ++i) {
        delete dst[i];
    }
    delete dirac;
    return rtn;
}

KuQCDStatus kuqcd_init()
{
    CmdLineParam cmdParam;
    cmdParam.ReadParam();
    g_grid_dim[0] = cmdParam.gridx;
    g_grid_dim[1] = cmdParam.gridy;
    g_grid_dim[2] = cmdParam.gridz;
    g_grid_dim[3] = cmdParam.gridt;
    if (GetNp() != g_grid_dim[0] * g_grid_dim[1] * g_grid_dim[2] * g_grid_dim[3]) {
        rootLogger.Error("g_grid_dim is not equal np");
        return KUQCD_ERR_INVALID_PARAM;
    }
    return KUQCD_SUCCESS;
}

KuQCDStatus kuqcd_load_gaugefield(void *h_gauge, KuQCDGaugeParam *param)
{
    try {
        if (h_gauge == nullptr || param == nullptr) {
            rootLogger.Error("h_gauge or param is nullptr");
            return KUQCD_ERR_NULL_POINTER;
        }

        if (param->gauge_order != KUQCD_BATCH_ORDER && param->gauge_order != KUQCD_NORMAL_ORDER) {
            rootLogger.Error("gauge_order is invalid");
            return KUQCD_ERR_INVALID_PARAM;
        }
        if (param->link_type != KUQCD_FAT_LINK && param->link_type != KUQCD_LONG_LINK) {
            rootLogger.Error("link_type is invalid");
            return KUQCD_ERR_INVALID_PARAM;
        }
        if (param->recon != KUQCD_GAUGE_RECON_NO) {
            if (param->recon == KUQCD_GAUGE_RECON_U14 || param->recon == KUQCD_GAUGE_RECON_12 ||
                param->recon == KUQCD_GAUGE_RECON_8) {
                rootLogger.Error("recon is not support");
                return KUQCD_ERR_UNSUPPORTED;
            }
            rootLogger.Error("recon is invalid");
            return KUQCD_ERR_INVALID_PARAM;
        }

        if (param->lat[0] == 0 || param->lat[1] == 0 || param->lat[2] == 0 || param->lat[3] == 0) {
            rootLogger.Error("param->lat is 0");
            return KUQCD_ERR_INVALID_PARAM;
        }

        LattBase *latbase = nullptr;
        KuQCDStatus rtn = set_latbase(param->lat[0], param->lat[1], param->lat[2], param->lat[3], param->halo_depth,
                                      param->od, param->site_subset, param->gauge_order, param->prec, g_grid_dim,
                                      &latbase);
        if (rtn != KUQCD_SUCCESS) {
            return rtn;
        }
        const KuQCDReconstructType recon = KUQCD_GAUGE_RECON_NO;
        g_comm = new CommBase(*latbase);
        if (param->gauge_order == KUQCD_BATCH_ORDER) {
            if (param->link_type == KUQCD_FAT_LINK) {
                if (param->prec == KUQCD_SINGLE) {
                    g_gauge_smeared_f16 = new GaugeField<float, 16, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<float, 16, recon>(g_gauge_smeared_f16, (float *)h_gauge);
                    return KUQCD_SUCCESS;
                }
                if (param->prec == KUQCD_DOUBLE) {
                    g_gauge_smeared_d8 = new GaugeField<double, 8, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<double, 8, recon>(g_gauge_smeared_d8, (double *)h_gauge);
                    return KUQCD_SUCCESS;
                }
            }
            if (param->link_type == KUQCD_LONG_LINK) {
                if (param->prec == KUQCD_SINGLE) {
                    g_gauge_naik_f16 = new GaugeField<float, 16, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<float, 16, recon>(g_gauge_naik_f16, (float *)h_gauge);
                    return KUQCD_SUCCESS;
                }
                if (param->prec == KUQCD_DOUBLE) {
                    g_gauge_naik_d8 = new GaugeField<double, 8, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<double, 8, recon>(g_gauge_naik_d8, (double *)h_gauge);
                    return KUQCD_SUCCESS;
                }
            }
        } else {
            if (param->link_type == KUQCD_FAT_LINK) {
                if (param->prec == KUQCD_SINGLE) {
                    g_gauge_smeared_f1 = new GaugeField<float, 1, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<float, 1, recon>(g_gauge_smeared_f1, (float *)h_gauge);
                    return KUQCD_SUCCESS;
                }
                if (param->prec == KUQCD_DOUBLE) {
                    g_gauge_smeared_d1 = new GaugeField<double, 1, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<double, 1, recon>(g_gauge_smeared_d1, (double *)h_gauge);
                    return KUQCD_SUCCESS;
                }
            }
            if (param->link_type == KUQCD_LONG_LINK) {
                if (param->prec == KUQCD_SINGLE) {
                    g_gauge_naik_f1 = new GaugeField<float, 1, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<float, 1, recon>(g_gauge_naik_f1, (float *)h_gauge);
                    return KUQCD_SUCCESS;
                }
                if (param->prec == KUQCD_DOUBLE) {
                    g_gauge_naik_d1 = new GaugeField<double, 1, recon>(*latbase, *g_comm);
                    kuqcd_gaugefield_set_values<double, 1, recon>(g_gauge_naik_d1, (double *)h_gauge);
                    return KUQCD_SUCCESS;
                }
            }
        }
    } catch (std::exception &e) {
        if (typeid(e) == typeid(std::bad_alloc)) {
            return KUQCD_ERR_NO_MEMORY;
        } else {
            return KUQCD_ERR_UNKNOWN_ERROR;
        }
    }
    return KUQCD_SUCCESS;
}

KuQCDStatus kuqcd_solve_multishift(void **_hp_x, void *_hp_b, KuQCDInvertParam *param)
{
    try {
        if (_hp_x == nullptr || _hp_b == nullptr || param == nullptr) {
            rootLogger.Error("hp_x or hp_b or param is nullptr");
            return KUQCD_ERR_NULL_POINTER;
        }
        if (param->matpc_type != KUQCD_EVEN_EVEN && param->matpc_type != KUQCD_ODD_ODD) {
            rootLogger.Error("matpc_type is error!");
            return KUQCD_ERR_INVALID_PARAM;
        }
        if (g_gauge_smeared_f16 != nullptr && g_gauge_naik_f16 != nullptr) {
            const int Nbatch = 16;
            return kuqcd_inverter<float, Nbatch, 1>(&g_inverter_f16, g_gauge_smeared_f16, g_gauge_naik_f16, _hp_x,
                                                    _hp_b, *param);
        }
        if (g_gauge_smeared_d8 != nullptr && g_gauge_naik_d8 != nullptr) {
            const int Nbatch = 8;
            return kuqcd_inverter<double, Nbatch, 1>(&g_inverter_d8, g_gauge_smeared_d8, g_gauge_naik_d8, _hp_x, _hp_b,
                                                     *param);
        }

        if (g_gauge_smeared_f1 != nullptr && g_gauge_naik_f1 != nullptr) {
            const int Nbatch = 1;
            return kuqcd_inverter<float, Nbatch, 1>(&g_inverter_f1, g_gauge_smeared_f1, g_gauge_naik_f1, _hp_x, _hp_b,
                                                    *param);
        }
        if (g_gauge_smeared_d1 != nullptr && g_gauge_naik_d1 != nullptr) {
            const int Nbatch = 1;
            return kuqcd_inverter<double, Nbatch, 1>(&g_inverter_d1, g_gauge_smeared_d1, g_gauge_naik_d1, _hp_x, _hp_b,
                                                     *param);
        }
    } catch (std::exception &e) {
        if (typeid(e) == typeid(std::bad_alloc)) {
            return KUQCD_ERR_NO_MEMORY;
        } else {
            return KUQCD_ERR_UNKNOWN_ERROR;
        }
    }
    rootLogger.Error("no init gauge_smeared and gauge_naik");
    return KUQCD_ERR_UNKNOWN_TYPE;
}

void kuqcd_gauge_free()
{
    delete g_latbase_even_halo0;
    g_latbase_even_halo0 = nullptr;
    delete g_latbase_odd_halo0;
    g_latbase_odd_halo0 = nullptr;
    delete g_latbase_all_halo0;
    g_latbase_all_halo0 = nullptr;
    delete g_latbase_all_halo3;
    g_latbase_all_halo3 = nullptr;
    delete g_comm;
    g_comm = nullptr;
    delete g_gauge_smeared_f1;
    g_gauge_smeared_f1 = nullptr;
    delete g_gauge_smeared_f16;
    g_gauge_smeared_f16 = nullptr;
    delete g_gauge_smeared_d1;
    g_gauge_smeared_d1 = nullptr;
    delete g_gauge_smeared_d8;
    g_gauge_smeared_d8 = nullptr;

    delete g_gauge_naik_f1;
    g_gauge_naik_f1 = nullptr;
    delete g_gauge_naik_f16;
    g_gauge_naik_f16 = nullptr;
    delete g_gauge_naik_d1;
    g_gauge_naik_d1 = nullptr;
    delete g_gauge_naik_d8;
    g_gauge_naik_d8 = nullptr;
}
void kuqcd_end()
{
    delete g_inverter_d1;
    g_inverter_d1 = nullptr;
    delete g_inverter_d8;
    g_inverter_d8 = nullptr;
    delete g_inverter_f1;
    g_inverter_f1 = nullptr;
    delete g_inverter_f16;
    g_inverter_f16 = nullptr;
}

} // namespace kuqcd
