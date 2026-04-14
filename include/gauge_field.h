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
#ifndef GAUGE_FIELD_H
#define GAUGE_FIELD_H

#include "lattice_field.h"
#include "communication.h"

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType Recon>
class GaugeField : public LatticeField<SU3<T, Nbatch, Recon>> {
public:
    typedef SU3<T, Nbatch, Recon> GaugeT;

    GaugeT *even;
    GaugeT *odd;
    CommBase &comm;

    std::array<std::array<Communicator<GaugeT>, 2>, 4> commHyperPlane;
    std::array<std::array<Communicator<GaugeT>, 4>, 6> commPlane;
    std::array<std::array<Communicator<GaugeT>, 8>, 4> commStripe;
    std::array<Communicator<GaugeT>, 16> commCorner;

    GaugeField() = delete;
    GaugeField(const GaugeField<T, Nbatch, Recon> &gauge);
    GaugeField(LattBase &latBase, CommBase &comm_);
    SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO> GetLinkCopy(size_t i, int dir) const;
    SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO> GetLinkDaggCopy(size_t i, int dir) const;

    void Update(KuQCDHalo halo = KUQCD_COMMALL);
    void StartComm(KuQCDHalo halo);
    void EndComm(KuQCDHalo halo);
    void PackHalo(Communicator<GaugeT> &communicator, const Index &idxSend);
    void UnpackHalo(Communicator<GaugeT> &communicator, const Index &idxRecv);
};

extern template class GaugeField<float, 1, KUQCD_GAUGE_RECON_NO>;
extern template class GaugeField<double, 1, KUQCD_GAUGE_RECON_NO>;
extern template class GaugeField<float, 16, KUQCD_GAUGE_RECON_NO>;
extern template class GaugeField<double, 8, KUQCD_GAUGE_RECON_NO>;

extern template class GaugeField<float, 1, KUQCD_GAUGE_RECON_12>;
extern template class GaugeField<double, 1, KUQCD_GAUGE_RECON_12>;
extern template class GaugeField<float, 16, KUQCD_GAUGE_RECON_12>;
extern template class GaugeField<double, 8, KUQCD_GAUGE_RECON_12>;

extern template class GaugeField<float, 1, KUQCD_GAUGE_RECON_U14>;
extern template class GaugeField<double, 1, KUQCD_GAUGE_RECON_U14>;
extern template class GaugeField<float, 16, KUQCD_GAUGE_RECON_U14>;
extern template class GaugeField<double, 8, KUQCD_GAUGE_RECON_U14>;

}  // namespace kuqcd

#endif
