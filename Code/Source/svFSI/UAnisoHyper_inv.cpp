/*#include "UAnisoHyper_inv.h"
namespace UAnisoHyper_inv {
Array<int> t_ind;

//----------
// ten_init
//----------
// Initialize tensor index pointer
//
void ten_init(const int nd)
{
  if (t_ind.size() != 0) {
    return;
  }

  int nn = pow(nd, 4);
  t_ind.resize(4, nn);

  int ii = 0;
  for (int l = 0; l < nd; l++) {
    for (int k = 0; k < nd; k++) {
      for (int j = 0; j < nd; j++) {
        for (int i = 0; i < nd; i++) {
          t_ind(0,ii) = i;
          t_ind(1,ii) = j;
          t_ind(2,ii) = k;
          t_ind(3,ii) = l;
          ii = ii + 1;
        }
      }
    }
  }
}
};*/

#include "mat_fun.h"
#include "utils.h"

namespace UAnisoHyper_inv {

// extern Array<int> t_ind;

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
void uCANN(const double xInv,const int kf0, const int kf1, const int kf2, const double W0, const double W1,const double W2, double &psi, double dpsi[9], double ddpsi[9]){
    double f0,df0,ddf0;
    uCANN_h0(xInv, kf0, f0,df0,ddf0);
    double f1,df1,ddf1;
    uCANN_h1(f0, kf1, W0, f1, df1, ddf1);
    double f2,df2,ddf2;
    uCANN_h2(f1, kf2, W1,f2,df2,ddf2);
    //updating
    psi = psi + W2*f2;
    for (int i = 0; i < 9; i++)
    {
        dpsi[i] = dpsi[i] + W2*df2*df1*df0;
        ddpsi[i] = ddpsi[i] + W2*((ddf2*df1*df1+df2*ddf1)*df0*df0+df2*df1*ddf0);
    }
    
}

/// @brief function to build psi and dpsidI1 to 5
void uanisohyper_inv(const double aInv[9],const std::vector<std::vector<double>> w, double &psi, double dpsi[9], double ddpsi[9]){
    //initialising
    for (int i = 0; i < 9; i++)
    {
        dpsi[i] = 0;
        ddpsi[i] = 0;
    }
    int kInv=0,kf0=0,kf1=0,kf2=0;//activation function for layer
    double W0=0,W1=0,W2=0;//weight for layers
    //reference config
    double ref[9] = {3, 3, 1, 1, 1, 1, 1, 1, 1};
    for (int i = 0; i < 4; i++) //each row of param table
    {
        //extract invariant, activation function and weight
        kInv = w[i][0];
        kf0 = w[i][1];
        kf1 = w[i][2];
        kf2 = w[i][3];
        W0 = w[i][4];
        W1 = w[i][5];
        W2 = w[i][6];

        //invariants in reference configuration
        double xInv = aInv[kInv-1] - ref[kInv-1];

        //psi and 1st and 2nd derivatives
        uCANN(xInv,kf0,kf1,kf2,W0,W1,W2,psi,dpsi,ddpsi);
    }

}
}