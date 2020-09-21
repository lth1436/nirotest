/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4490379982002751748);
void inv_err_fun(double *nom_x, double *true_x, double *out_454337454911244336);
void H_mod_fun(double *state, double *out_5974883223479611163);
void f_fun(double *state, double dt, double *out_4049466981791189978);
void F_fun(double *state, double dt, double *out_476212456923785358);
void h_25(double *state, double *unused, double *out_7028759246035361135);
void H_25(double *state, double *unused, double *out_7445302580211337171);
void h_24(double *state, double *unused, double *out_5664921550296873259);
void H_24(double *state, double *unused, double *out_4141898990419740337);
void h_30(double *state, double *unused, double *out_6295385278944071789);
void H_30(double *state, double *unused, double *out_5729499810405912089);
void h_26(double *state, double *unused, double *out_4167870889059920426);
void H_26(double *state, double *unused, double *out_8098662222629267045);
void h_27(double *state, double *unused, double *out_2506394768562793394);
void H_27(double *state, double *unused, double *out_5833910161381330167);
void h_29(double *state, double *unused, double *out_2400015747567763120);
void H_29(double *state, double *unused, double *out_5684338417918773251);
void h_28(double *state, double *unused, double *out_9119707630279785157);
void H_28(double *state, double *unused, double *out_826261990081590891);
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
