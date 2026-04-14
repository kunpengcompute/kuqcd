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
#include "inverter_cgm.h"
#include "dirac.h"

namespace kuqcd {

template<typename T, typename VecT>
int CGM<T, VecT>::Invert(int noShifts, T shifts[], const VecT & src, VecT *sol[], const double tol, int maxIters)
{
    int maxTerm = noShifts;
    T a[noShifts];
    T B[noShifts];
    T Z[noShifts];
    T Zm1[noShifts];

    std::vector<VecT> pi(noShifts, src);
    VecT s(src.latBase);
    VecT r(src);
    VecT pi0(src.latBase);
    for (int i = 0; i < noShifts; ++i) {
        a[i] = 0.0;
        B[i] = 1.0;
        Z[i] = 1.0;
        Zm1[i] = 1.0;
        sol[i] -> Clear();
    }

    T rNorm2 = Norm2(r);
    T pAp, lambda, lambda2, rr1, Bm1 = 1.0;

    int iter = 0;
    do {
        iter++;
        pi0 = pi[0];
        this->dirac->MdaggM(pi0, s, true);
        Axmy(s, pi0, shifts[0]);
        pAp = RealDot(pi0, s);
        B[0] = -rNorm2 / pAp;
        
        Axpy(r, s, B[0]);
        for (int j = 1; j < maxTerm; ++j) {
            rr1 = Bm1 * Zm1[j] / ( B[0] * a[0] * (Zm1[j] - Z[j]) 
                    + Zm1[j] * Bm1 * (1.0 - shifts[j] * B[0]) );
            Zm1[j] = Z[j];
            Z[j]   = Z[j] * rr1;
            B[j]   = B[0] * rr1;
        }
        Bm1 = B[0];
        lambda2 = Norm2(r);
        a[0] = lambda2 / rNorm2;
        rNorm2 = lambda2;
        for (int i = 0; i < maxTerm; ++i) {
            Ymax(*sol[i], pi[i], B[i]);
        }
        for (int j = 1; j < maxTerm; ++j) {
            a[j] = a[0] * Z[j] * B[j] / (Zm1[j] * B[0]);
        }
        for (int i = 0; i < maxTerm; ++i) {
            Axpby(pi[i], r, Z[i], a[i]);
        }

        do {
            lambda = Z[maxTerm - 1] * Z[maxTerm - 1] * lambda2;
            if (lambda < tol) {
                maxTerm--;
            }
        } while ( (maxTerm > 0) && (lambda < tol) );

    } while ( (maxTerm > 0) && (iter < maxIters) );

    if (iter >= maxIters - 1) {
        Rank0Printf("CGM: Warning max iteration reached %d\n", iter);
    } else {
        Rank0Printf("CGM: # iterations %d\n", iter);
    }

    return iter;
}


template class CGM<float, SpinorField<float, 1, 1>>;

template class CGM<float, SpinorField<float, 16, 1>>;
template class CGM<double, SpinorField<double, 1, 1>>;
template class CGM<double, SpinorField<double, 8, 1>>;
}