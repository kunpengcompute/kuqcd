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
#include "gauge_field.h"

namespace kuqcd {

template <typename T, int Nbatch, KuQCDReconstructType Recon>
GaugeField<T, Nbatch, Recon>::GaugeField(const GaugeField<T, Nbatch, Recon> &gauge)
    : LatticeField<GaugeT>(gauge), comm(gauge.comm)
{
    even = this->data;
    odd = this->data + this->latBase.latSubFullCB.V;
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
GaugeField<T, Nbatch, Recon>::GaugeField(LattBase &latBase, CommBase &comm_)
    : LatticeField<GaugeT>(latBase, 4 * Nbatch * Recon), comm(comm_)
{
    even = this->data;
    odd = this->data + this->latBase.latSubFullCB.V;

    int offset = 0;
    int idx = 0;
    int num = 2;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            int rk = comm.rankHyperPlane[i][j];
            size_t len = comm.lenHyperPlane[i];
            int div = idx / num;
            int mod = idx - div * num;
            int ts;
            int tr;
            ts = offset + idx;
            tr = offset + (num - 1 - mod) + div * num;
            commHyperPlane[i][j].Init(rk, len, len * this->locSize * latBase.prec, ts, tr);
            idx++;
        }
    }

    offset += idx;
    idx = 0;
    num = 4;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
            int rk = comm.rankPlane[i][j];
            size_t len = comm.lenPlane[i];
            int div = idx / num;
            int mod = idx - div * num;
            int ts;
            int tr;
            ts = offset + idx;
            tr = offset + (num - 1 - mod) + div * num;
            commPlane[i][j].Init(rk, len, len * this->locSize * latBase.prec, ts, tr);
            idx++;
        }
    }

    offset += idx;
    idx = 0;
    num = 8;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 8; ++j) {
            int rk = comm.rankStripe[i][j];
            size_t len = comm.lenStripe[i];
            int div = idx / num;
            int mod = idx - div * num;
            int ts;
            int tr;
            ts = offset + idx;
            tr = offset + (num - 1 - mod) + div * num;
            commStripe[i][j].Init(rk, len, len * this->locSize * latBase.prec, ts, tr);
            idx++;
        }
    }

    offset += idx;
    idx = 0;
    num = 16;
    for (int i = 0; i < 16; ++i) {
        int rk = comm.rankCorner[i];
        size_t len = comm.lenCorner;
        int div = idx / num;
        int mod = idx - div * num;
        int ts;
        int tr;
        ts = offset + idx;
        tr = offset + (num - 1 - mod) + div * num;
        commCorner[i].Init(rk, len, len * this->locSize * latBase.prec, ts, tr);
        idx++;
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO> GaugeField<T, Nbatch, Recon>::GetLinkCopy(size_t i, int dir) const
{
    return (this->data[dir + 4 * i]).Reconstruct();
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
SU3<T, Nbatch, KUQCD_GAUGE_RECON_NO> GaugeField<T, Nbatch, Recon>::GetLinkDaggCopy(size_t i, int dir) const
{
    return (this->data[dir + 4 * i]).Reconstruct();
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void GaugeField<T, Nbatch, Recon>::Update(KuQCDHalo halo)
{
    StartComm(halo);
    EndComm(halo);
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void GaugeField<T, Nbatch, Recon>::StartComm(KuQCDHalo halo)
{
    if (halo & KUQCD_HYPERPLANE) {
        for (int i = 0; i < 4; ++i) {
            if (comm.lenHyperPlane[i] > 0) {
                for (int j = 0; j < 2; ++j) {
#ifndef USE_KUPL_COMM
                    commHyperPlane[i][j].Recv();
#endif
                    PackHalo(commHyperPlane[i][j], comm.haloIdxSend1[i][j]);
#ifndef USE_KUPL_COMM
                    commHyperPlane[i][j].Send();
#endif
                }
            }
        }
    }
    if (halo & KUQCD_PLANE) {
        for (int i = 0; i < 6; ++i) {
            if (comm.lenPlane[i] > 0) {
                for (int j = 0; j < 4; ++j) {
#ifndef USE_KUPL_COMM
                    commPlane[i][j].Recv();
#endif
                    PackHalo(commPlane[i][j], comm.haloIdxSend2[i][j]);
#ifndef USE_KUPL_COMM
                    commPlane[i][j].Send();
#endif
                }
            }
        }
    }
    if (halo & KUQCD_STRIPE) {
        for (int i = 0; i < 4; ++i) {
            if (comm.lenStripe[i] > 0) {
                for (int j = 0; j < 8; ++j) {
#ifndef USE_KUPL_COMM
                    commStripe[i][j].Recv();
#endif
                    PackHalo(commStripe[i][j], comm.haloIdxSend3[i][j]);
#ifndef USE_KUPL_COMM
                    commStripe[i][j].Send();
#endif
                }
            }
        }
    }
    if (halo & KUQCD_CORNER) {
        for (int i = 0; i < 16; ++i) {
            if (comm.lenCorner > 0) {
#ifndef USE_KUPL_COMM
                commCorner[i].Recv();
#endif
                PackHalo(commCorner[i], comm.haloIdxSend4[i]);
#ifndef USE_KUPL_COMM
                commCorner[i].Send();
#endif
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void GaugeField<T, Nbatch, Recon>::EndComm(KuQCDHalo halo)
{
    if (halo & KUQCD_HYPERPLANE) {
        for (int i = 0; i < 4; ++i) {
            if (comm.lenHyperPlane[i] > 0) {
                for (int j = 0; j < 2; ++j) {
#ifndef USE_KUPL_COMM
                    commHyperPlane[i][j].Finish(true);
#endif
                    UnpackHalo(commHyperPlane[i][j], comm.haloIdxRecv1[i][j]);
#ifndef USE_KUPL_COMM
                    commHyperPlane[i][j].Finish(false);
#endif
                }
            }
        }
    }
    if (halo & KUQCD_PLANE) {
        for (int i = 0; i < 6; ++i) {
            if (comm.lenPlane[i] > 0) {
                for (int j = 0; j < 4; ++j) {
#ifndef USE_KUPL_COMM
                    commPlane[i][j].Finish(true);
#endif
                    UnpackHalo(commPlane[i][j], comm.haloIdxRecv2[i][j]);
#ifndef USE_KUPL_COMM
                    commPlane[i][j].Finish(false);
#endif
                }
            }
        }
    }
    if (halo & KUQCD_STRIPE) {
        for (int i = 0; i < 4; ++i) {
            if (comm.lenStripe[i] > 0) {
                for (int j = 0; j < 8; ++j) {
#ifndef USE_KUPL_COMM
                    commStripe[i][j].Finish(true);
#endif
                    UnpackHalo(commStripe[i][j], comm.haloIdxRecv3[i][j]);
#ifndef USE_KUPL_COMM
                    commStripe[i][j].Finish(false);
#endif
                }
            }
        }
    }
    if (halo & KUQCD_CORNER) {
        for (int i = 0; i < 16; ++i) {
            if (comm.lenCorner > 0) {
#ifndef USE_KUPL_COMM
                commCorner[i].Finish(true);
#endif
                UnpackHalo(commCorner[i], comm.haloIdxRecv4[i]);
#ifndef USE_KUPL_COMM
                commCorner[i].Finish(false);
#endif
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void GaugeField<T, Nbatch, Recon>::PackHalo(Communicator<GaugeT> &communicator, const Index &idxSend)
{
    int ls = this->locSize / this->latBase.Nbatch;
    if ((ls != 6) && (ls != 24) && (ls != 72) && (ls != 56) && (ls != 48)) {
        Rank0Fprintf(stderr, "error local size! %d %s %d\n", ls, __FILE__, __LINE__);
    }

    if (ls == 6) {
        for (size_t i = 0; i < communicator.length; ++i) {
            communicator.dataSend[i] = this->data[idxSend.idx[i]];
        }
    } else {
        for (size_t i = 0; i < communicator.length; ++i) {
            for (int j = 0; j < 4; ++j) {
                communicator.dataSend[j + 4 * i] = this->data[j + 4 * idxSend.idx[i]];
            }
        }
    }
}

template <typename T, int Nbatch, KuQCDReconstructType Recon>
void GaugeField<T, Nbatch, Recon>::UnpackHalo(Communicator<GaugeT> &communicator, const Index &idxRecv)
{
    int ls = this->locSize / this->latBase.Nbatch;
    if ((ls != 6) && (ls != 24) && (ls != 72) && (ls != 56) && (ls != 48)) {
        Rank0Fprintf(stderr, "error local size! %d %s %d\n", ls, __FILE__, __LINE__);
    }
#ifdef USE_KUPL_COMM
    MPI_Barrier(MPI_COMM_WORLD);
    void *remote;
    GaugeT *halo;
    int ts = communicator.tagSend;
    HaloMemory *memSend;
    if (ts < 8) {
        int num = 2;
        int i = ts / num;
        int j = ts - num * i;
        j = num - 1 - j;
        memSend = &commHyperPlane[i][j].memSend;
    } else if (ts < 32) {
        ts -= 8;
        int num = 4;
        int i = ts / num;
        int j = ts - num * i;
        j = num - 1 - j;
        memSend = &commPlane[i][j].memSend;
    } else if (ts < 64) {
        ts -= 32;
        int num = 8;
        int i = ts / num;
        int j = ts - num * i;
        j = num - 1 - j;
        memSend = &commStripe[i][j].memSend;
    } else {
        ts -= 64;
        int num = 16;
        int i = num - 1 - ts;
        memSend = &commCorner[i].memSend;
    }
    kupl_shm_win_query(memSend->win, communicator.nbRank, &remote);
    halo = static_cast<GaugeT *>(remote);

    if (ls == 6) {
        for (size_t i = 0; i < communicator.length; ++i) {
            this->data[idxRecv.idx[i]] = halo[i];
        }
    } else {
        for (size_t i = 0; i < communicator.length; ++i) {
            for (int j = 0; j < 4; ++j) {
                this->data[j + 4 * idxRecv.idx[i]] = halo[j + 4 * i];
            }
        }
    }
#else
    if (ls == 6) {
        for (size_t i = 0; i < communicator.length; ++i) {
            this->data[idxRecv.idx[i]] = communicator.dataRecv[i];
        }
    } else {
        for (size_t i = 0; i < communicator.length; ++i) {
            for (int j = 0; j < 4; ++j) {
                this->data[j + 4 * idxRecv.idx[i]] = communicator.dataRecv[j + 4 * i];
            }
        }
    }
#endif
}

template class GaugeField<float, 1, KUQCD_GAUGE_RECON_NO>;
template class GaugeField<double, 1, KUQCD_GAUGE_RECON_NO>;
template class GaugeField<float, 16, KUQCD_GAUGE_RECON_NO>;
template class GaugeField<double, 8, KUQCD_GAUGE_RECON_NO>;

template class GaugeField<float, 1, KUQCD_GAUGE_RECON_12>;
template class GaugeField<double, 1, KUQCD_GAUGE_RECON_12>;
template class GaugeField<float, 16, KUQCD_GAUGE_RECON_12>;
template class GaugeField<double, 8, KUQCD_GAUGE_RECON_12>;

template class GaugeField<float, 1, KUQCD_GAUGE_RECON_U14>;
template class GaugeField<double, 1, KUQCD_GAUGE_RECON_U14>;
template class GaugeField<float, 16, KUQCD_GAUGE_RECON_U14>;
template class GaugeField<double, 8, KUQCD_GAUGE_RECON_U14>;

} // namespace kuqcd
