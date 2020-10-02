/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_1870969132376137805);
void inv_err_fun(double *nom_x, double *true_x, double *out_3802553487131977071);
void H_mod_fun(double *state, double *out_3762111176001179456);
void f_fun(double *state, double dt, double *out_2480536562585758481);
void F_fun(double *state, double dt, double *out_1425722008230317098);
void h_3(double *state, double *unused, double *out_580480245261766277);
void H_3(double *state, double *unused, double *out_8039239072597041883);
void h_4(double *state, double *unused, double *out_4245608367930427789);
void H_4(double *state, double *unused, double *out_7691622073353942749);
void h_9(double *state, double *unused, double *out_6350337263133992527);
void H_9(double *state, double *unused, double *out_2203357074039413582);
void h_10(double *state, double *unused, double *out_5883173185422846693);
void H_10(double *state, double *unused, double *out_6275558625866197697);
void h_12(double *state, double *unused, double *out_1029741073182827027);
void H_12(double *state, double *unused, double *out_3005095031555557250);
void h_31(double *state, double *unused, double *out_6461370520649494816);
void H_31(double *state, double *unused, double *out_5474853234853266326);
void h_32(double *state, double *unused, double *out_2101340024716927954);
void H_32(double *state, double *unused, double *out_4491552036411005924);
void h_13(double *state, double *unused, double *out_2424738802151988478);
void H_13(double *state, double *unused, double *out_8550083838302739916);
void h_14(double *state, double *unused, double *out_6350337263133992527);
void H_14(double *state, double *unused, double *out_2203357074039413582);
void h_19(double *state, double *unused, double *out_5874493945137936551);
void H_19(double *state, double *unused, double *out_6160069650627068086);
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