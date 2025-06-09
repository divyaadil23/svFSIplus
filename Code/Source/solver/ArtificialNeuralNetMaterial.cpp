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

// The functions required for CANN material model implementations are defined here

#include "mat_fun.h"
#include "utils.h"
#include "Parameters.h"
#include "ComMod.h"

/// @brief 0th layer output of CANN for activation func kf, input x
void ArtificialNeuralNetMaterial::uCANN_h0(const double x, const int kf, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = x;
        df = 1;
        ddf = 0;
    } else if (kf == 2) {
        f = (std::abs(x) + x) / 2;
        df = 0.5 * (std::abs(x) / x + 1);
        ddf = 0;
    } else if (kf == 3) {
        f = std::abs(x);
        df = std::abs(x) / x;
        ddf = 0;
    }
}

/// @brief 1st layer output of CANN for activation func kf, input x, weight W
void ArtificialNeuralNetMaterial::uCANN_h1(const double x, const int kf, const double W, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = W * x;
        df = W;
        ddf = 0;
    } else if (kf == 2) {
        f = W * W * x * x;
        df = 2 * W * W * x;
        ddf = 2 * W * W;
    }
}

/// @brief 2nd layer output of CANN for activation func kf, input x, weight W
void ArtificialNeuralNetMaterial::uCANN_h2(const double x, const int kf, const double W, double &f, double &df, double &ddf) const {
    if (kf == 1) {
        f = W * x;
        df = W;
        ddf = 0;
    } else if (kf == 2) {
        f = std::exp(W * x) - 1;
        df = W * std::exp(W * x);
        ddf = W * W * std::exp(W * x);
    } else if (kf == 3) {
        f = -std::log(1 - W * x);
        df = W / (1 - W * x);
        ddf = -W * W / ((1 - W * x) * (1 - W * x));
    }
}

/// @brief Updates psi and its derivatives
void ArtificialNeuralNetMaterial::uCANN(
    const double xInv, const int kInv,
    const int kf0, const int kf1, const int kf2,
    const double W0, const double W1, const double W2,
    double &psi, double (&dpsi)[9], double (&ddpsi)[9]
) const {
    double f0, df0, ddf0;
    uCANN_h0(xInv, kf0, f0, df0, ddf0);
    double f1, df1, ddf1;
    uCANN_h1(f0, kf1, W0, f1, df1, ddf1);
    double f2, df2, ddf2;
    uCANN_h2(f1, kf2, W1, f2, df2, ddf2);

    psi += W2 * f2;
    dpsi[kInv - 1] += W2 * df2 * df1 * df0;
    ddpsi[kInv - 1] += W2 * ((ddf2 * df1 * df1 + df2 * ddf1) * df0 * df0 + df2 * df1 * ddf0);
}

/// @brief function to build psi and dpsidI1 to 9
void ArtificialNeuralNetMaterial::evaluate(const double aInv[9], double &psi, double (&dpsi)[9], double (&ddpsi)[9]) const {
    // Initializing
    psi = 0;
    for (int i = 0; i < 9; ++i) {
        dpsi[i] = 0;
        ddpsi[i] = 0;
    }

    double ref[9] = {3, 3, 1, 1, 1, 0, 0, 1, 1};

    for (int i = 0; i < num_rows; ++i) {
        int kInv = this->invariant_indices(i);
        int kf0 = this->activation_functions(i, 0);
        int kf1 = this->activation_functions(i, 1);
        int kf2 = this->activation_functions(i, 2);
        double W0 = this->weights(i, 0);
        double W1 = this->weights(i, 1);
        double W2 = this->weights(i, 2);

        double xInv = aInv[kInv - 1] - ref[kInv - 1];
        uCANN(xInv, kInv, kf0, kf1, kf2, W0, W1, W2, psi, dpsi, ddpsi);
    }
}