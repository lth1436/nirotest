
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
void err_fun(double *nom_x, double *delta_x, double *out_1201342623059153752) {
   out_1201342623059153752[0] = delta_x[0] + nom_x[0];
   out_1201342623059153752[1] = delta_x[1] + nom_x[1];
   out_1201342623059153752[2] = delta_x[2] + nom_x[2];
   out_1201342623059153752[3] = delta_x[3] + nom_x[3];
   out_1201342623059153752[4] = delta_x[4] + nom_x[4];
   out_1201342623059153752[5] = delta_x[5] + nom_x[5];
   out_1201342623059153752[6] = delta_x[6] + nom_x[6];
   out_1201342623059153752[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2764648982842188044) {
   out_2764648982842188044[0] = -nom_x[0] + true_x[0];
   out_2764648982842188044[1] = -nom_x[1] + true_x[1];
   out_2764648982842188044[2] = -nom_x[2] + true_x[2];
   out_2764648982842188044[3] = -nom_x[3] + true_x[3];
   out_2764648982842188044[4] = -nom_x[4] + true_x[4];
   out_2764648982842188044[5] = -nom_x[5] + true_x[5];
   out_2764648982842188044[6] = -nom_x[6] + true_x[6];
   out_2764648982842188044[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_8098358221935223241) {
   out_8098358221935223241[0] = 1.0;
   out_8098358221935223241[1] = 0.0;
   out_8098358221935223241[2] = 0.0;
   out_8098358221935223241[3] = 0.0;
   out_8098358221935223241[4] = 0.0;
   out_8098358221935223241[5] = 0.0;
   out_8098358221935223241[6] = 0.0;
   out_8098358221935223241[7] = 0.0;
   out_8098358221935223241[8] = 0.0;
   out_8098358221935223241[9] = 1.0;
   out_8098358221935223241[10] = 0.0;
   out_8098358221935223241[11] = 0.0;
   out_8098358221935223241[12] = 0.0;
   out_8098358221935223241[13] = 0.0;
   out_8098358221935223241[14] = 0.0;
   out_8098358221935223241[15] = 0.0;
   out_8098358221935223241[16] = 0.0;
   out_8098358221935223241[17] = 0.0;
   out_8098358221935223241[18] = 1.0;
   out_8098358221935223241[19] = 0.0;
   out_8098358221935223241[20] = 0.0;
   out_8098358221935223241[21] = 0.0;
   out_8098358221935223241[22] = 0.0;
   out_8098358221935223241[23] = 0.0;
   out_8098358221935223241[24] = 0.0;
   out_8098358221935223241[25] = 0.0;
   out_8098358221935223241[26] = 0.0;
   out_8098358221935223241[27] = 1.0;
   out_8098358221935223241[28] = 0.0;
   out_8098358221935223241[29] = 0.0;
   out_8098358221935223241[30] = 0.0;
   out_8098358221935223241[31] = 0.0;
   out_8098358221935223241[32] = 0.0;
   out_8098358221935223241[33] = 0.0;
   out_8098358221935223241[34] = 0.0;
   out_8098358221935223241[35] = 0.0;
   out_8098358221935223241[36] = 1.0;
   out_8098358221935223241[37] = 0.0;
   out_8098358221935223241[38] = 0.0;
   out_8098358221935223241[39] = 0.0;
   out_8098358221935223241[40] = 0.0;
   out_8098358221935223241[41] = 0.0;
   out_8098358221935223241[42] = 0.0;
   out_8098358221935223241[43] = 0.0;
   out_8098358221935223241[44] = 0.0;
   out_8098358221935223241[45] = 1.0;
   out_8098358221935223241[46] = 0.0;
   out_8098358221935223241[47] = 0.0;
   out_8098358221935223241[48] = 0.0;
   out_8098358221935223241[49] = 0.0;
   out_8098358221935223241[50] = 0.0;
   out_8098358221935223241[51] = 0.0;
   out_8098358221935223241[52] = 0.0;
   out_8098358221935223241[53] = 0.0;
   out_8098358221935223241[54] = 1.0;
   out_8098358221935223241[55] = 0.0;
   out_8098358221935223241[56] = 0.0;
   out_8098358221935223241[57] = 0.0;
   out_8098358221935223241[58] = 0.0;
   out_8098358221935223241[59] = 0.0;
   out_8098358221935223241[60] = 0.0;
   out_8098358221935223241[61] = 0.0;
   out_8098358221935223241[62] = 0.0;
   out_8098358221935223241[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_3111364141878358374) {
   out_3111364141878358374[0] = state[0];
   out_3111364141878358374[1] = state[1];
   out_3111364141878358374[2] = state[2];
   out_3111364141878358374[3] = state[3];
   out_3111364141878358374[4] = state[4];
   out_3111364141878358374[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3111364141878358374[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3111364141878358374[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2286017523026274910) {
   out_2286017523026274910[0] = 1;
   out_2286017523026274910[1] = 0;
   out_2286017523026274910[2] = 0;
   out_2286017523026274910[3] = 0;
   out_2286017523026274910[4] = 0;
   out_2286017523026274910[5] = 0;
   out_2286017523026274910[6] = 0;
   out_2286017523026274910[7] = 0;
   out_2286017523026274910[8] = 0;
   out_2286017523026274910[9] = 1;
   out_2286017523026274910[10] = 0;
   out_2286017523026274910[11] = 0;
   out_2286017523026274910[12] = 0;
   out_2286017523026274910[13] = 0;
   out_2286017523026274910[14] = 0;
   out_2286017523026274910[15] = 0;
   out_2286017523026274910[16] = 0;
   out_2286017523026274910[17] = 0;
   out_2286017523026274910[18] = 1;
   out_2286017523026274910[19] = 0;
   out_2286017523026274910[20] = 0;
   out_2286017523026274910[21] = 0;
   out_2286017523026274910[22] = 0;
   out_2286017523026274910[23] = 0;
   out_2286017523026274910[24] = 0;
   out_2286017523026274910[25] = 0;
   out_2286017523026274910[26] = 0;
   out_2286017523026274910[27] = 1;
   out_2286017523026274910[28] = 0;
   out_2286017523026274910[29] = 0;
   out_2286017523026274910[30] = 0;
   out_2286017523026274910[31] = 0;
   out_2286017523026274910[32] = 0;
   out_2286017523026274910[33] = 0;
   out_2286017523026274910[34] = 0;
   out_2286017523026274910[35] = 0;
   out_2286017523026274910[36] = 1;
   out_2286017523026274910[37] = 0;
   out_2286017523026274910[38] = 0;
   out_2286017523026274910[39] = 0;
   out_2286017523026274910[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2286017523026274910[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2286017523026274910[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2286017523026274910[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2286017523026274910[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2286017523026274910[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2286017523026274910[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2286017523026274910[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2286017523026274910[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2286017523026274910[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2286017523026274910[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2286017523026274910[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2286017523026274910[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2286017523026274910[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2286017523026274910[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2286017523026274910[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2286017523026274910[56] = 0;
   out_2286017523026274910[57] = 0;
   out_2286017523026274910[58] = 0;
   out_2286017523026274910[59] = 0;
   out_2286017523026274910[60] = 0;
   out_2286017523026274910[61] = 0;
   out_2286017523026274910[62] = 0;
   out_2286017523026274910[63] = 1;
}
void h_25(double *state, double *unused, double *out_8729249890204695209) {
   out_8729249890204695209[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8737919755792925159) {
   out_8737919755792925159[0] = 0;
   out_8737919755792925159[1] = 0;
   out_8737919755792925159[2] = 0;
   out_8737919755792925159[3] = 0;
   out_8737919755792925159[4] = 0;
   out_8737919755792925159[5] = 0;
   out_8737919755792925159[6] = 1;
   out_8737919755792925159[7] = 0;
}
void h_24(double *state, double *unused, double *out_7185823008761736623) {
   out_7185823008761736623[0] = state[4];
   out_7185823008761736623[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2327549382950519651) {
   out_2327549382950519651[0] = 0;
   out_2327549382950519651[1] = 0;
   out_2327549382950519651[2] = 0;
   out_2327549382950519651[3] = 0;
   out_2327549382950519651[4] = 1;
   out_2327549382950519651[5] = 0;
   out_2327549382950519651[6] = 0;
   out_2327549382950519651[7] = 0;
   out_2327549382950519651[8] = 0;
   out_2327549382950519651[9] = 0;
   out_2327549382950519651[10] = 0;
   out_2327549382950519651[11] = 0;
   out_2327549382950519651[12] = 0;
   out_2327549382950519651[13] = 1;
   out_2327549382950519651[14] = 0;
   out_2327549382950519651[15] = 0;
}
void h_30(double *state, double *unused, double *out_6588648484514385109) {
   out_6588648484514385109[0] = state[4];
}
void H_30(double *state, double *unused, double *out_5939412702814310723) {
   out_5939412702814310723[0] = 0;
   out_5939412702814310723[1] = 0;
   out_5939412702814310723[2] = 0;
   out_5939412702814310723[3] = 0;
   out_5939412702814310723[4] = 1;
   out_5939412702814310723[5] = 0;
   out_5939412702814310723[6] = 0;
   out_5939412702814310723[7] = 0;
}
void h_26(double *state, double *unused, double *out_772904261541737074) {
   out_772904261541737074[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1670109373613969857) {
   out_1670109373613969857[0] = 0;
   out_1670109373613969857[1] = 0;
   out_1670109373613969857[2] = 0;
   out_1670109373613969857[3] = 0;
   out_1670109373613969857[4] = 0;
   out_1670109373613969857[5] = 0;
   out_1670109373613969857[6] = 0;
   out_1670109373613969857[7] = 1;
}
void h_27(double *state, double *unused, double *out_4630137711968706042) {
   out_4630137711968706042[0] = state[3];
}
void H_27(double *state, double *unused, double *out_526244294419090717) {
   out_526244294419090717[0] = 0;
   out_526244294419090717[1] = 0;
   out_526244294419090717[2] = 0;
   out_526244294419090717[3] = 1;
   out_526244294419090717[4] = 0;
   out_526244294419090717[5] = 0;
   out_526244294419090717[6] = 0;
   out_526244294419090717[7] = 0;
}
void h_29(double *state, double *unused, double *out_388780479785955400) {
   out_388780479785955400[0] = state[1];
}
void H_29(double *state, double *unused, double *out_7492351629742088727) {
   out_7492351629742088727[0] = 0;
   out_7492351629742088727[1] = 1;
   out_7492351629742088727[2] = 0;
   out_7492351629742088727[3] = 0;
   out_7492351629742088727[4] = 0;
   out_7492351629742088727[5] = 0;
   out_7492351629742088727[6] = 0;
   out_7492351629742088727[7] = 0;
}
void h_28(double *state, double *unused, double *out_6748720444976415439) {
   out_6748720444976415439[0] = state[5];
   out_6748720444976415439[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3227202666917231703) {
   out_3227202666917231703[0] = 0;
   out_3227202666917231703[1] = 0;
   out_3227202666917231703[2] = 0;
   out_3227202666917231703[3] = 0;
   out_3227202666917231703[4] = 0;
   out_3227202666917231703[5] = 1;
   out_3227202666917231703[6] = 0;
   out_3227202666917231703[7] = 0;
   out_3227202666917231703[8] = 0;
   out_3227202666917231703[9] = 0;
   out_3227202666917231703[10] = 0;
   out_3227202666917231703[11] = 0;
   out_3227202666917231703[12] = 0;
   out_3227202666917231703[13] = 0;
   out_3227202666917231703[14] = 1;
   out_3227202666917231703[15] = 0;
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
