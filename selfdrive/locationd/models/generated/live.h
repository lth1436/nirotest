/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8433551513585401999);
void inv_err_fun(double *nom_x, double *true_x, double *out_5615633631998976527);
void H_mod_fun(double *state, double *out_6020067654100738459);
void f_fun(double *state, double dt, double *out_7133067072523873278);
void F_fun(double *state, double dt, double *out_8909620253679392972);
void h_3(double *state, double *unused, double *out_1373248960769480060);
void H_3(double *state, double *unused, double *out_5682869864149212887);
void h_4(double *state, double *unused, double *out_4653414708174628449);
void H_4(double *state, double *unused, double *out_8837565147588353398);
void h_9(double *state, double *unused, double *out_6108531424031299673);
void H_9(double *state, double *unused, double *out_1579408731260042745);
void h_10(double *state, double *unused, double *out_5117097304373056523);
void H_10(double *state, double *unused, double *out_8697786474636878091);
void h_12(double *state, double *unused, double *out_4616037805987861237);
void H_12(double *state, double *unused, double *out_8020348470470144233);
void h_31(double *state, double *unused, double *out_8201654042927603878);
void H_31(double *state, double *unused, double *out_8595444509552093617);
void h_32(double *state, double *unused, double *out_6769430688286271057);
void H_32(double *state, double *unused, double *out_6621198096819605329);
void h_13(double *state, double *unused, double *out_3080285956011812149);
void H_13(double *state, double *unused, double *out_1923905294336314226);
void h_14(double *state, double *unused, double *out_6108531424031299673);
void H_14(double *state, double *unused, double *out_1579408731260042745);
void h_19(double *state, double *unused, double *out_2670845100652305735);
void H_19(double *state, double *unused, double *out_8896645812202096099);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);