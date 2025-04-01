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
    void uanisohyper_inv(const double aInv[9], const std::vector<CANNRow> CANNTable, double &psi, double (&dpsi)[9], double (&ddpsi)[9]);
}

#endif // UANISOHYPER_INV_H