/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* This material model implementation is based on the following paper: 
Peirlinck, M., Hurtado, J.A., Rausch, M.K. et al. A universal material model subroutine 
for soft matter systems. Engineering with Computers 41, 905–927 (2025). 
https://doi.org/10.1007/s00366-024-02031-w */

#ifndef UANISOHYPER_INV_H
#define UANISOHYPER_INV_H

#include "mat_fun.h"
#include "utils.h"
#include "Parameters.h"
#include <vector>

namespace UAnisoHyper_inv {
    void uCANN_h0(const double x, const int kf, double &f, double &df, double &ddf);
    void uCANN_h1(const double x, const int kf, const double W, double &f, double &df, double &ddf);
    void uCANN_h2(const double x, const int kf, const double W, double &f, double &df, double &ddf);
    void uCANN(const double xInv, const int kf0, const int kf1, const int kf2, const double W0, const double W1, const double W2, double &psi, double (&dpsi)[9], double (&ddpsi)[9]);
    void uanisohyper_inv(const double aInv[9],const constArtificialNeuralNetworkModel paramTable, double &psi, double (&dpsi)[9], double (&ddpsi)[9]);
}

#endif // UANISOHYPER_INV_H