/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1111255803602454333);
void inv_err_fun(double *nom_x, double *true_x, double *out_9093210972024974323);
void H_mod_fun(double *state, double *out_7102727360988229312);
void f_fun(double *state, double dt, double *out_26888974827175776);
void F_fun(double *state, double dt, double *out_981366099399393920);
void h_25(double *state, double *unused, double *out_7106248415861724463);
void H_25(double *state, double *unused, double *out_4538587127572630088);
void h_24(double *state, double *unused, double *out_2025308003415535098);
void H_24(double *state, double *unused, double *out_8538186852106219476);
void h_30(double *state, double *unused, double *out_1778133451039258899);
void H_30(double *state, double *unused, double *out_8546541607025698754);
void h_26(double *state, double *unused, double *out_4976961111696060968);
void H_26(double *state, double *unused, double *out_295421804121230272);
void h_27(double *state, double *unused, double *out_5521590686661771044);
void H_27(double *state, double *unused, double *out_1231388240910778596);
void h_29(double *state, double *unused, double *out_446118961994981662);
void H_29(double *state, double *unused, double *out_4470828859922145384);
void h_28(double *state, double *unused, double *out_7615952005031418958);
void H_28(double *state, double *unused, double *out_8422112739034565298);
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
