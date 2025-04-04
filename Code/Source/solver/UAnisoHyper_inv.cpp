#include "mat_fun.h"
#include "utils.h"
#include "Parameters.h"

namespace UAnisoHyper_inv {

//void ten_init(const int nd);
/// @brief 0th layer output of CANN for activation func kf, input x
void uCANN_h0(const double x, const int kf, double &f, double &df, double &ddf){
    
    if (kf==1)
    {
        f = x;
        df = 1;
        ddf = 0;
    }
    else if (kf==2){
        f = (abs(x)+x)/2;
        df = 1/2*(abs(x)/x + 1);
        ddf = 0;
    }
    else if(kf==3){
        f = abs(x);
        df = abs(x)/x;
        ddf = 0;
    }
}

/// @brief 1st layer output of CANN for activation func kf, input x, weight W0
void uCANN_h1(const double x, const int kf, const double W, double &f, double &df, double &ddf){
    
    if (kf==1){
        f = W*x;
        df = W*1;
        ddf = 0;
    }
    else if (kf==2){
        f = W*W*x*x;
        df = 2*W*W*x;
        ddf = 2*W*W;
    }
}

/// @brief 2nd layer output of CANN for activation func kf, input x, weight W1
void uCANN_h2(const double x, const int kf, const double W, double &f, double &df, double &ddf){
    
    if (kf==1){
        f = W*x;
        df = W*1;
        ddf = 0;
    }
    else if(kf==2){
        f = exp(W*x) - 1;
        df = W*exp(W*x);
        ddf = W*W*exp(W*x);
    }
    else if (kf==3){
        f = -log(1-W*x);
        df = W/(1-W*x);
        ddf = -W*W/((1-W*x)*(1-W*x));
    }
}

/// @brief Updates psi and its derivatives
void uCANN(const double xInv,const int kInv,const int kf0, const int kf1, const int kf2, const double W0, const double W1,const double W2, double &psi, double (&dpsi)[9], double (&ddpsi)[9]){
    
    double f0,df0,ddf0;
    uCANN_h0(xInv, kf0, f0,df0,ddf0);
    double f1,df1,ddf1;
    uCANN_h1(f0, kf1, W0, f1, df1, ddf1);
    double f2,df2,ddf2;
    uCANN_h2(f1, kf2, W1,f2,df2,ddf2);

    //updating
    psi = psi + W2*f2;
    dpsi[kInv-1] = dpsi[kInv-1] + W2*df2*df1*df0;
    ddpsi[kInv-1] = ddpsi[kInv-1] + W2*((ddf2*df1*df1+df2*ddf1)*df0*df0+df2*df1*ddf0);
}

/// @brief function to build psi and dpsidI1 to 5
void uanisohyper_inv(const double aInv[9],const std::vector<CANNRow> CANNTable, double &psi, double (&dpsi)[9], double (&ddpsi)[9]){
    //initialising
    for (int i = 0; i < 9; i++)
    {
        dpsi[i] = 0;
        ddpsi[i] = 0;
    }
    int kInv=0,kf0=0,kf1=0,kf2=0;//activation function for layer
    double W0=0,W1=0,W2=0;//weight for layers

    //reference config
    double ref[9] = {3, 3, 1, 1, 1, 0, 0, 1, 1};
    int nRows = CANNTable.size();

    for (int i = 0; i < nRows; i++) //each row of param table
    {
        //extract invariant, activation function and weight
        kInv = CANNTable[i].invariant_index.value_;
        kf0 = CANNTable[i].activation_functions.value_[0];
        kf1 = CANNTable[i].activation_functions.value_[1];
        kf2 = CANNTable[i].activation_functions.value_[2];
        W0 = CANNTable[i].weights.value_[0];
        W1 = CANNTable[i].weights.value_[1];
        W2 = CANNTable[i].weights.value_[2];

        //invariants in reference configuration
        double xInv = aInv[kInv-1] - ref[kInv-1];

        //psi and 1st and 2nd derivative
        uCANN(xInv,kInv,kf0,kf1,kf2,W0,W1,W2,psi,dpsi,ddpsi);
    }

}
}