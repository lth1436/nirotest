
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1111255803602454333) {
   out_1111255803602454333[0] = delta_x[0] + nom_x[0];
   out_1111255803602454333[1] = delta_x[1] + nom_x[1];
   out_1111255803602454333[2] = delta_x[2] + nom_x[2];
   out_1111255803602454333[3] = delta_x[3] + nom_x[3];
   out_1111255803602454333[4] = delta_x[4] + nom_x[4];
   out_1111255803602454333[5] = delta_x[5] + nom_x[5];
   out_1111255803602454333[6] = delta_x[6] + nom_x[6];
   out_1111255803602454333[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_9093210972024974323) {
   out_9093210972024974323[0] = -nom_x[0] + true_x[0];
   out_9093210972024974323[1] = -nom_x[1] + true_x[1];
   out_9093210972024974323[2] = -nom_x[2] + true_x[2];
   out_9093210972024974323[3] = -nom_x[3] + true_x[3];
   out_9093210972024974323[4] = -nom_x[4] + true_x[4];
   out_9093210972024974323[5] = -nom_x[5] + true_x[5];
   out_9093210972024974323[6] = -nom_x[6] + true_x[6];
   out_9093210972024974323[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_7102727360988229312) {
   out_7102727360988229312[0] = 1.0;
   out_7102727360988229312[1] = 0.0;
   out_7102727360988229312[2] = 0.0;
   out_7102727360988229312[3] = 0.0;
   out_7102727360988229312[4] = 0.0;
   out_7102727360988229312[5] = 0.0;
   out_7102727360988229312[6] = 0.0;
   out_7102727360988229312[7] = 0.0;
   out_7102727360988229312[8] = 0.0;
   out_7102727360988229312[9] = 1.0;
   out_7102727360988229312[10] = 0.0;
   out_7102727360988229312[11] = 0.0;
   out_7102727360988229312[12] = 0.0;
   out_7102727360988229312[13] = 0.0;
   out_7102727360988229312[14] = 0.0;
   out_7102727360988229312[15] = 0.0;
   out_7102727360988229312[16] = 0.0;
   out_7102727360988229312[17] = 0.0;
   out_7102727360988229312[18] = 1.0;
   out_7102727360988229312[19] = 0.0;
   out_7102727360988229312[20] = 0.0;
   out_7102727360988229312[21] = 0.0;
   out_7102727360988229312[22] = 0.0;
   out_7102727360988229312[23] = 0.0;
   out_7102727360988229312[24] = 0.0;
   out_7102727360988229312[25] = 0.0;
   out_7102727360988229312[26] = 0.0;
   out_7102727360988229312[27] = 1.0;
   out_7102727360988229312[28] = 0.0;
   out_7102727360988229312[29] = 0.0;
   out_7102727360988229312[30] = 0.0;
   out_7102727360988229312[31] = 0.0;
   out_7102727360988229312[32] = 0.0;
   out_7102727360988229312[33] = 0.0;
   out_7102727360988229312[34] = 0.0;
   out_7102727360988229312[35] = 0.0;
   out_7102727360988229312[36] = 1.0;
   out_7102727360988229312[37] = 0.0;
   out_7102727360988229312[38] = 0.0;
   out_7102727360988229312[39] = 0.0;
   out_7102727360988229312[40] = 0.0;
   out_7102727360988229312[41] = 0.0;
   out_7102727360988229312[42] = 0.0;
   out_7102727360988229312[43] = 0.0;
   out_7102727360988229312[44] = 0.0;
   out_7102727360988229312[45] = 1.0;
   out_7102727360988229312[46] = 0.0;
   out_7102727360988229312[47] = 0.0;
   out_7102727360988229312[48] = 0.0;
   out_7102727360988229312[49] = 0.0;
   out_7102727360988229312[50] = 0.0;
   out_7102727360988229312[51] = 0.0;
   out_7102727360988229312[52] = 0.0;
   out_7102727360988229312[53] = 0.0;
   out_7102727360988229312[54] = 1.0;
   out_7102727360988229312[55] = 0.0;
   out_7102727360988229312[56] = 0.0;
   out_7102727360988229312[57] = 0.0;
   out_7102727360988229312[58] = 0.0;
   out_7102727360988229312[59] = 0.0;
   out_7102727360988229312[60] = 0.0;
   out_7102727360988229312[61] = 0.0;
   out_7102727360988229312[62] = 0.0;
   out_7102727360988229312[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_26888974827175776) {
   out_26888974827175776[0] = state[0];
   out_26888974827175776[1] = state[1];
   out_26888974827175776[2] = state[2];
   out_26888974827175776[3] = state[3];
   out_26888974827175776[4] = state[4];
   out_26888974827175776[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_26888974827175776[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_26888974827175776[7] = state[7];
}
void F_fun(double *state, double dt, double *out_981366099399393920) {
   out_981366099399393920[0] = 1;
   out_981366099399393920[1] = 0;
   out_981366099399393920[2] = 0;
   out_981366099399393920[3] = 0;
   out_981366099399393920[4] = 0;
   out_981366099399393920[5] = 0;
   out_981366099399393920[6] = 0;
   out_981366099399393920[7] = 0;
   out_981366099399393920[8] = 0;
   out_981366099399393920[9] = 1;
   out_981366099399393920[10] = 0;
   out_981366099399393920[11] = 0;
   out_981366099399393920[12] = 0;
   out_981366099399393920[13] = 0;
   out_981366099399393920[14] = 0;
   out_981366099399393920[15] = 0;
   out_981366099399393920[16] = 0;
   out_981366099399393920[17] = 0;
   out_981366099399393920[18] = 1;
   out_981366099399393920[19] = 0;
   out_981366099399393920[20] = 0;
   out_981366099399393920[21] = 0;
   out_981366099399393920[22] = 0;
   out_981366099399393920[23] = 0;
   out_981366099399393920[24] = 0;
   out_981366099399393920[25] = 0;
   out_981366099399393920[26] = 0;
   out_981366099399393920[27] = 1;
   out_981366099399393920[28] = 0;
   out_981366099399393920[29] = 0;
   out_981366099399393920[30] = 0;
   out_981366099399393920[31] = 0;
   out_981366099399393920[32] = 0;
   out_981366099399393920[33] = 0;
   out_981366099399393920[34] = 0;
   out_981366099399393920[35] = 0;
   out_981366099399393920[36] = 1;
   out_981366099399393920[37] = 0;
   out_981366099399393920[38] = 0;
   out_981366099399393920[39] = 0;
   out_981366099399393920[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_981366099399393920[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_981366099399393920[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_981366099399393920[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_981366099399393920[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_981366099399393920[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_981366099399393920[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_981366099399393920[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_981366099399393920[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_981366099399393920[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_981366099399393920[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_981366099399393920[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_981366099399393920[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_981366099399393920[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_981366099399393920[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_981366099399393920[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_981366099399393920[56] = 0;
   out_981366099399393920[57] = 0;
   out_981366099399393920[58] = 0;
   out_981366099399393920[59] = 0;
   out_981366099399393920[60] = 0;
   out_981366099399393920[61] = 0;
   out_981366099399393920[62] = 0;
   out_981366099399393920[63] = 1;
}
void h_25(double *state, double *unused, double *out_7106248415861724463) {
   out_7106248415861724463[0] = state[6];
}
void H_25(double *state, double *unused, double *out_4538587127572630088) {
   out_4538587127572630088[0] = 0;
   out_4538587127572630088[1] = 0;
   out_4538587127572630088[2] = 0;
   out_4538587127572630088[3] = 0;
   out_4538587127572630088[4] = 0;
   out_4538587127572630088[5] = 0;
   out_4538587127572630088[6] = 1;
   out_4538587127572630088[7] = 0;
}
void h_24(double *state, double *unused, double *out_2025308003415535098) {
   out_2025308003415535098[0] = state[4];
   out_2025308003415535098[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8538186852106219476) {
   out_8538186852106219476[0] = 0;
   out_8538186852106219476[1] = 0;
   out_8538186852106219476[2] = 0;
   out_8538186852106219476[3] = 0;
   out_8538186852106219476[4] = 1;
   out_8538186852106219476[5] = 0;
   out_8538186852106219476[6] = 0;
   out_8538186852106219476[7] = 0;
   out_8538186852106219476[8] = 0;
   out_8538186852106219476[9] = 0;
   out_8538186852106219476[10] = 0;
   out_8538186852106219476[11] = 0;
   out_8538186852106219476[12] = 0;
   out_8538186852106219476[13] = 1;
   out_8538186852106219476[14] = 0;
   out_8538186852106219476[15] = 0;
}
void h_30(double *state, double *unused, double *out_1778133451039258899) {
   out_1778133451039258899[0] = state[4];
}
void H_30(double *state, double *unused, double *out_8546541607025698754) {
   out_8546541607025698754[0] = 0;
   out_8546541607025698754[1] = 0;
   out_8546541607025698754[2] = 0;
   out_8546541607025698754[3] = 0;
   out_8546541607025698754[4] = 1;
   out_8546541607025698754[5] = 0;
   out_8546541607025698754[6] = 0;
   out_8546541607025698754[7] = 0;
}
void h_26(double *state, double *unused, double *out_4976961111696060968) {
   out_4976961111696060968[0] = state[7];
}
void H_26(double *state, double *unused, double *out_295421804121230272) {
   out_295421804121230272[0] = 0;
   out_295421804121230272[1] = 0;
   out_295421804121230272[2] = 0;
   out_295421804121230272[3] = 0;
   out_295421804121230272[4] = 0;
   out_295421804121230272[5] = 0;
   out_295421804121230272[6] = 0;
   out_295421804121230272[7] = 1;
}
void h_27(double *state, double *unused, double *out_5521590686661771044) {
   out_5521590686661771044[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1231388240910778596) {
   out_1231388240910778596[0] = 0;
   out_1231388240910778596[1] = 0;
   out_1231388240910778596[2] = 0;
   out_1231388240910778596[3] = 1;
   out_1231388240910778596[4] = 0;
   out_1231388240910778596[5] = 0;
   out_1231388240910778596[6] = 0;
   out_1231388240910778596[7] = 0;
}
void h_29(double *state, double *unused, double *out_446118961994981662) {
   out_446118961994981662[0] = state[1];
}
void H_29(double *state, double *unused, double *out_4470828859922145384) {
   out_4470828859922145384[0] = 0;
   out_4470828859922145384[1] = 1;
   out_4470828859922145384[2] = 0;
   out_4470828859922145384[3] = 0;
   out_4470828859922145384[4] = 0;
   out_4470828859922145384[5] = 0;
   out_4470828859922145384[6] = 0;
   out_4470828859922145384[7] = 0;
}
void h_28(double *state, double *unused, double *out_7615952005031418958) {
   out_7615952005031418958[0] = state[5];
   out_7615952005031418958[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8422112739034565298) {
   out_8422112739034565298[0] = 0;
   out_8422112739034565298[1] = 0;
   out_8422112739034565298[2] = 0;
   out_8422112739034565298[3] = 0;
   out_8422112739034565298[4] = 0;
   out_8422112739034565298[5] = 1;
   out_8422112739034565298[6] = 0;
   out_8422112739034565298[7] = 0;
   out_8422112739034565298[8] = 0;
   out_8422112739034565298[9] = 0;
   out_8422112739034565298[10] = 0;
   out_8422112739034565298[11] = 0;
   out_8422112739034565298[12] = 0;
   out_8422112739034565298[13] = 0;
   out_8422112739034565298[14] = 1;
   out_8422112739034565298[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
