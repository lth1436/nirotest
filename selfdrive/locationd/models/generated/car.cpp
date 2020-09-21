
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
void err_fun(double *nom_x, double *delta_x, double *out_4490379982002751748) {
   out_4490379982002751748[0] = delta_x[0] + nom_x[0];
   out_4490379982002751748[1] = delta_x[1] + nom_x[1];
   out_4490379982002751748[2] = delta_x[2] + nom_x[2];
   out_4490379982002751748[3] = delta_x[3] + nom_x[3];
   out_4490379982002751748[4] = delta_x[4] + nom_x[4];
   out_4490379982002751748[5] = delta_x[5] + nom_x[5];
   out_4490379982002751748[6] = delta_x[6] + nom_x[6];
   out_4490379982002751748[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_454337454911244336) {
   out_454337454911244336[0] = -nom_x[0] + true_x[0];
   out_454337454911244336[1] = -nom_x[1] + true_x[1];
   out_454337454911244336[2] = -nom_x[2] + true_x[2];
   out_454337454911244336[3] = -nom_x[3] + true_x[3];
   out_454337454911244336[4] = -nom_x[4] + true_x[4];
   out_454337454911244336[5] = -nom_x[5] + true_x[5];
   out_454337454911244336[6] = -nom_x[6] + true_x[6];
   out_454337454911244336[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5974883223479611163) {
   out_5974883223479611163[0] = 1.0;
   out_5974883223479611163[1] = 0.0;
   out_5974883223479611163[2] = 0.0;
   out_5974883223479611163[3] = 0.0;
   out_5974883223479611163[4] = 0.0;
   out_5974883223479611163[5] = 0.0;
   out_5974883223479611163[6] = 0.0;
   out_5974883223479611163[7] = 0.0;
   out_5974883223479611163[8] = 0.0;
   out_5974883223479611163[9] = 1.0;
   out_5974883223479611163[10] = 0.0;
   out_5974883223479611163[11] = 0.0;
   out_5974883223479611163[12] = 0.0;
   out_5974883223479611163[13] = 0.0;
   out_5974883223479611163[14] = 0.0;
   out_5974883223479611163[15] = 0.0;
   out_5974883223479611163[16] = 0.0;
   out_5974883223479611163[17] = 0.0;
   out_5974883223479611163[18] = 1.0;
   out_5974883223479611163[19] = 0.0;
   out_5974883223479611163[20] = 0.0;
   out_5974883223479611163[21] = 0.0;
   out_5974883223479611163[22] = 0.0;
   out_5974883223479611163[23] = 0.0;
   out_5974883223479611163[24] = 0.0;
   out_5974883223479611163[25] = 0.0;
   out_5974883223479611163[26] = 0.0;
   out_5974883223479611163[27] = 1.0;
   out_5974883223479611163[28] = 0.0;
   out_5974883223479611163[29] = 0.0;
   out_5974883223479611163[30] = 0.0;
   out_5974883223479611163[31] = 0.0;
   out_5974883223479611163[32] = 0.0;
   out_5974883223479611163[33] = 0.0;
   out_5974883223479611163[34] = 0.0;
   out_5974883223479611163[35] = 0.0;
   out_5974883223479611163[36] = 1.0;
   out_5974883223479611163[37] = 0.0;
   out_5974883223479611163[38] = 0.0;
   out_5974883223479611163[39] = 0.0;
   out_5974883223479611163[40] = 0.0;
   out_5974883223479611163[41] = 0.0;
   out_5974883223479611163[42] = 0.0;
   out_5974883223479611163[43] = 0.0;
   out_5974883223479611163[44] = 0.0;
   out_5974883223479611163[45] = 1.0;
   out_5974883223479611163[46] = 0.0;
   out_5974883223479611163[47] = 0.0;
   out_5974883223479611163[48] = 0.0;
   out_5974883223479611163[49] = 0.0;
   out_5974883223479611163[50] = 0.0;
   out_5974883223479611163[51] = 0.0;
   out_5974883223479611163[52] = 0.0;
   out_5974883223479611163[53] = 0.0;
   out_5974883223479611163[54] = 1.0;
   out_5974883223479611163[55] = 0.0;
   out_5974883223479611163[56] = 0.0;
   out_5974883223479611163[57] = 0.0;
   out_5974883223479611163[58] = 0.0;
   out_5974883223479611163[59] = 0.0;
   out_5974883223479611163[60] = 0.0;
   out_5974883223479611163[61] = 0.0;
   out_5974883223479611163[62] = 0.0;
   out_5974883223479611163[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_4049466981791189978) {
   out_4049466981791189978[0] = state[0];
   out_4049466981791189978[1] = state[1];
   out_4049466981791189978[2] = state[2];
   out_4049466981791189978[3] = state[3];
   out_4049466981791189978[4] = state[4];
   out_4049466981791189978[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_4049466981791189978[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_4049466981791189978[7] = state[7];
}
void F_fun(double *state, double dt, double *out_476212456923785358) {
   out_476212456923785358[0] = 1;
   out_476212456923785358[1] = 0;
   out_476212456923785358[2] = 0;
   out_476212456923785358[3] = 0;
   out_476212456923785358[4] = 0;
   out_476212456923785358[5] = 0;
   out_476212456923785358[6] = 0;
   out_476212456923785358[7] = 0;
   out_476212456923785358[8] = 0;
   out_476212456923785358[9] = 1;
   out_476212456923785358[10] = 0;
   out_476212456923785358[11] = 0;
   out_476212456923785358[12] = 0;
   out_476212456923785358[13] = 0;
   out_476212456923785358[14] = 0;
   out_476212456923785358[15] = 0;
   out_476212456923785358[16] = 0;
   out_476212456923785358[17] = 0;
   out_476212456923785358[18] = 1;
   out_476212456923785358[19] = 0;
   out_476212456923785358[20] = 0;
   out_476212456923785358[21] = 0;
   out_476212456923785358[22] = 0;
   out_476212456923785358[23] = 0;
   out_476212456923785358[24] = 0;
   out_476212456923785358[25] = 0;
   out_476212456923785358[26] = 0;
   out_476212456923785358[27] = 1;
   out_476212456923785358[28] = 0;
   out_476212456923785358[29] = 0;
   out_476212456923785358[30] = 0;
   out_476212456923785358[31] = 0;
   out_476212456923785358[32] = 0;
   out_476212456923785358[33] = 0;
   out_476212456923785358[34] = 0;
   out_476212456923785358[35] = 0;
   out_476212456923785358[36] = 1;
   out_476212456923785358[37] = 0;
   out_476212456923785358[38] = 0;
   out_476212456923785358[39] = 0;
   out_476212456923785358[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_476212456923785358[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_476212456923785358[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_476212456923785358[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_476212456923785358[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_476212456923785358[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_476212456923785358[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_476212456923785358[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_476212456923785358[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_476212456923785358[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_476212456923785358[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_476212456923785358[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_476212456923785358[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_476212456923785358[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_476212456923785358[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_476212456923785358[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_476212456923785358[56] = 0;
   out_476212456923785358[57] = 0;
   out_476212456923785358[58] = 0;
   out_476212456923785358[59] = 0;
   out_476212456923785358[60] = 0;
   out_476212456923785358[61] = 0;
   out_476212456923785358[62] = 0;
   out_476212456923785358[63] = 1;
}
void h_25(double *state, double *unused, double *out_7028759246035361135) {
   out_7028759246035361135[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7445302580211337171) {
   out_7445302580211337171[0] = 0;
   out_7445302580211337171[1] = 0;
   out_7445302580211337171[2] = 0;
   out_7445302580211337171[3] = 0;
   out_7445302580211337171[4] = 0;
   out_7445302580211337171[5] = 0;
   out_7445302580211337171[6] = 1;
   out_7445302580211337171[7] = 0;
}
void h_24(double *state, double *unused, double *out_5664921550296873259) {
   out_5664921550296873259[0] = state[4];
   out_5664921550296873259[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4141898990419740337) {
   out_4141898990419740337[0] = 0;
   out_4141898990419740337[1] = 0;
   out_4141898990419740337[2] = 0;
   out_4141898990419740337[3] = 0;
   out_4141898990419740337[4] = 1;
   out_4141898990419740337[5] = 0;
   out_4141898990419740337[6] = 0;
   out_4141898990419740337[7] = 0;
   out_4141898990419740337[8] = 0;
   out_4141898990419740337[9] = 0;
   out_4141898990419740337[10] = 0;
   out_4141898990419740337[11] = 0;
   out_4141898990419740337[12] = 0;
   out_4141898990419740337[13] = 1;
   out_4141898990419740337[14] = 0;
   out_4141898990419740337[15] = 0;
}
void h_30(double *state, double *unused, double *out_6295385278944071789) {
   out_6295385278944071789[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5729499810405912089) {
   out_5729499810405912089[0] = 0;
   out_5729499810405912089[1] = 0;
   out_5729499810405912089[2] = 0;
   out_5729499810405912089[3] = 0;
   out_5729499810405912089[4] = 1;
   out_5729499810405912089[5] = 0;
   out_5729499810405912089[6] = 0;
   out_5729499810405912089[7] = 0;
}
void h_26(double *state, double *unused, double *out_4167870889059920426) {
   out_4167870889059920426[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8098662222629267045) {
   out_8098662222629267045[0] = 0;
   out_8098662222629267045[1] = 0;
   out_8098662222629267045[2] = 0;
   out_8098662222629267045[3] = 0;
   out_8098662222629267045[4] = 0;
   out_8098662222629267045[5] = 0;
   out_8098662222629267045[6] = 0;
   out_8098662222629267045[7] = 1;
}
void h_27(double *state, double *unused, double *out_2506394768562793394) {
   out_2506394768562793394[0] = state[3];
}
void H_27(double *state, double *unused, double *out_5833910161381330167) {
   out_5833910161381330167[0] = 0;
   out_5833910161381330167[1] = 0;
   out_5833910161381330167[2] = 0;
   out_5833910161381330167[3] = 1;
   out_5833910161381330167[4] = 0;
   out_5833910161381330167[5] = 0;
   out_5833910161381330167[6] = 0;
   out_5833910161381330167[7] = 0;
}
void h_29(double *state, double *unused, double *out_2400015747567763120) {
   out_2400015747567763120[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5684338417918773251) {
   out_5684338417918773251[0] = 0;
   out_5684338417918773251[1] = 1;
   out_5684338417918773251[2] = 0;
   out_5684338417918773251[3] = 0;
   out_5684338417918773251[4] = 0;
   out_5684338417918773251[5] = 0;
   out_5684338417918773251[6] = 0;
   out_5684338417918773251[7] = 0;
}
void h_28(double *state, double *unused, double *out_9119707630279785157) {
   out_9119707630279785157[0] = state[5];
   out_9119707630279785157[1] = state[6];
}
void H_28(double *state, double *unused, double *out_826261990081590891) {
   out_826261990081590891[0] = 0;
   out_826261990081590891[1] = 0;
   out_826261990081590891[2] = 0;
   out_826261990081590891[3] = 0;
   out_826261990081590891[4] = 0;
   out_826261990081590891[5] = 1;
   out_826261990081590891[6] = 0;
   out_826261990081590891[7] = 0;
   out_826261990081590891[8] = 0;
   out_826261990081590891[9] = 0;
   out_826261990081590891[10] = 0;
   out_826261990081590891[11] = 0;
   out_826261990081590891[12] = 0;
   out_826261990081590891[13] = 0;
   out_826261990081590891[14] = 1;
   out_826261990081590891[15] = 0;
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
