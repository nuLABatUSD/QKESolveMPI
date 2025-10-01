#ifndef _RUN_PARAMS_HH_
#define _RUN_PARAMS_HH_

#define PARAM_DELTA_M_SQUARED 1.e-18
#define PARAM_SIN_2THETA 0.8

#define EPS_MAX_LINSPACE 20.
#define EPS_LINSPACE_POINTS 201

#define IC_TCM 32.
#define IC_TEMP 32.
#define IC_NU_E 0.9
#define IC_NU_MU (2*0.9)
#define IC_NUBAR_E 0.9
#define IC_NUBAR_MU (2*0.9)
#define IC_MAX_DISTFUN 0.9

#define PARAM_N_STEPS 5
#define PARAM_DN 1

#define PARAM_DT_INIT 1.e11
#define PARAM_T_FINAL 5.e20


#define _COMPUTE_R_ERROR_ -999.
#define TOLERANCE_MIN_RATE 100.

#define ODE_SOLVER_TOLERANCE 1.e-8
#define ODE_SOLVER_TINY (ODE_SOLVER_TOLERANCE * 1.e-16)
#define ODE_SOLVER_SAFETY 0.9

#endif
