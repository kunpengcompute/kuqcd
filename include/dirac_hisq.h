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
#ifndef DIRAC_HISQ_H
#define DIRAC_HISQ_H

#include <array>

#include "dirac.h"

#include "spinor_field.h"
#include "gauge_field.h"
#include "communication.h"

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType ReconS, KuQCDReconstructType ReconN>
class DiracHISQ : public Dirac<SpinorField<T, Nbatch, 1>> {
public:
    typedef ColorSpinor<T, Nbatch> SpinorT;
    typedef SpinorField<T, Nbatch, 1> VecT;
    typedef SU3<T, Nbatch, ReconS> GaugeST;
    typedef SU3<T, Nbatch, ReconN> GaugeNT;
    typedef GaugeField<T, Nbatch, ReconS> LinkST;
    typedef GaugeField<T, Nbatch, ReconN> LinkNT;
    const LinkST &gaugeSmeared;
    const LinkNT &gaugeNaik;
    VecT *tmp = nullptr;
    LattBase &latBase;

    double mass = 0;
    T mass2 = 0, c3000 = -1.0 / 48.0;

    int myRank = 0;
    int np = 0;
    std::array<unsigned short, 4> &gridDim;
    std::array<unsigned int, 4> surf{0, 0, 0, 0};
    std::array<int, 8> nbRank{0, 0, 0, 0, 0, 0, 0, 0};
    Index nb1Idx;
    Index nb3Idx;
    std::array<Index, 8> edge1Idx;
    std::array<Index, 8> edge3Idx;
    std::array<Communicator<ColorSpinor<T, Nbatch>>, 8> comm;
    std::array<Communicator<ColorSpinor<T, Nbatch>>, 8> comm3;

    DiracHISQ(KuQCDDslashType diracType, const LinkST &gaugeSmeared_, const LinkNT &gaugeNaik_, const double mass_,
        VecT *tmp_ = nullptr, T naikEps = 0);
    ~DiracHISQ() override {};
    void Dslash(const VecT &src, VecT &dst, bool update = false) override;
    void MdaggM(const VecT &src, VecT &dst, bool update = false) override;
    void InitComm();
    void InitNbIdx();
    void InitNbIdxParity(KuQCDSiteSubset site_subset);  // result site_subset
    void PackHalo(const VecT &vec, int dir);
    void UnpackHalo(const VecT &vec, int dir);
    void StartComm(const VecT &vec, bool update);
    void EndComm(const VecT &vec, bool updata);

    template <int dir>
    void FwdDslash(size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch = true);
    template <int dir>
    void BwdDslash(size_t i, const VecT &src, SpinorT &res, size_t offsetg, bool ismatch = true);
};

extern template class DiracHISQ<float, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;
extern template class DiracHISQ<double, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;
extern template class DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;
extern template class DiracHISQ<double, 8, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_NO>;

extern template class DiracHISQ<float, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;
extern template class DiracHISQ<double, 1, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;
extern template class DiracHISQ<float, 16, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;
extern template class DiracHISQ<double, 8, KUQCD_GAUGE_RECON_NO, KUQCD_GAUGE_RECON_U14>;

}  // namespace kuqcd

#endif
