/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1201342623059153752);
void inv_err_fun(double *nom_x, double *true_x, double *out_2764648982842188044);
void H_mod_fun(double *state, double *out_8098358221935223241);
void f_fun(double *state, double dt, double *out_3111364141878358374);
void F_fun(double *state, double dt, double *out_2286017523026274910);
void h_25(double *state, double *unused, double *out_8729249890204695209);
void H_25(double *state, double *unused, double *out_8737919755792925159);
void h_24(double *state, double *unused, double *out_7185823008761736623);
void H_24(double *state, double *unused, double *out_2327549382950519651);
void h_30(double *state, double *unused, double *out_6588648484514385109);
void H_30(double *state, double *unused, double *out_5939412702814310723);
void h_26(double *state, double *unused, double *out_772904261541737074);
void H_26(double *state, double *unused, double *out_1670109373613969857);
void h_27(double *state, double *unused, double *out_4630137711968706042);
void H_27(double *state, double *unused, double *out_526244294419090717);
void h_29(double *state, double *unused, double *out_388780479785955400);
void H_29(double *state, double *unused, double *out_7492351629742088727);
void h_28(double *state, double *unused, double *out_6748720444976415439);
void H_28(double *state, double *unused, double *out_3227202666917231703);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
