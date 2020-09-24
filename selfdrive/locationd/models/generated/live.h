/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_2486298187852117759);
void inv_err_fun(double *nom_x, double *true_x, double *out_1907148703811358813);
void H_mod_fun(double *state, double *out_8500660281032028337);
void f_fun(double *state, double dt, double *out_6299326005496379330);
void F_fun(double *state, double dt, double *out_8138688164867987612);
void h_3(double *state, double *unused, double *out_790521826163922144);
void H_3(double *state, double *unused, double *out_7694198342491843571);
void h_4(double *state, double *unused, double *out_8607356167012547265);
void H_4(double *state, double *unused, double *out_3391498359091349094);
void h_9(double *state, double *unused, double *out_6950665766866551305);
void H_9(double *state, double *unused, double *out_8002594424114930475);
void h_10(double *state, double *unused, double *out_6265902820605023831);
void H_10(double *state, double *unused, double *out_2639958239788819617);
void h_12(double *state, double *unused, double *out_1785371751328674611);
void H_12(double *state, double *unused, double *out_563268921265657029);
void h_31(double *state, double *unused, double *out_5205794172231818250);
void H_31(double *state, double *unused, double *out_1770534533509061491);
void h_32(double *state, double *unused, double *out_1728208116770433709);
void H_32(double *state, double *unused, double *out_2466189320526617685);
void h_13(double *state, double *unused, double *out_6182317564822547637);
void H_13(double *state, double *unused, double *out_7323804728896736870);
void h_14(double *state, double *unused, double *out_6950665766866551305);
void H_14(double *state, double *unused, double *out_8002594424114930475);
void h_19(double *state, double *unused, double *out_1912131698598172753);
void H_19(double *state, double *unused, double *out_3155539806108080961);
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