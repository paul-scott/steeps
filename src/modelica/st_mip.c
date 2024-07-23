#include "linprog/linprog.h"

#include "gurobi_c.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/*#define ST_LINPROG_DEBUG*/

#ifdef ST_LINPROG_DEBUG
# define MSG(FMT,...) fprintf(stdout,"%s:%d: " FMT "\n",__FILE__,__LINE__,##__VA_ARGS__)
# define MSG1(FMT,...) fprintf(stdout,"%s:%d: " FMT,__FILE__,__LINE__,##__VA_ARGS__)
# define MSG2(FMT,...) fprintf(stdout,FMT,##__VA_ARGS__)
# define MSGL fprintf(stdout,"\n")
#else
# define MSG(...) ((void)0)
# define MSG1(...) ((void)0)
# define MSG2(...) ((void)0)
# define MSGL ((void)0)
#endif

#define ERR(FMT,...) fprintf(stderr,"%s:%d: " FMT "\n",__FILE__,__LINE__,##__VA_ARGS__)

#define ST_ERRVAL (-999999.)

double gurobi_fun(
         MotabData *pvd
        ,MotabData *wnd
        ,double P_elec_max_pv
        ,double P_elec_max_wind
        ,double P_elec_pv_ref_size
        ,double P_elec_wind_ref_size
        ,int horizon // Time horizon for the optimization
        ,double dt // Time step duration
        ,double t0 // Initial time
        ,double eff_heater // Efficiency of the heater
        ,double eff_process // Efficiency of the process
        ,double DEmax // Maximum dispatched energy
        ,double SLmax // Maximum stored energy (storage limit)
        ,double SLinit // Initial stored energy
        ,double SLmin // Minimum stored energy
        ,double ramp_up_fraction // Fraction of DEmax allowed for ramping up
        ,double ramp_dw_fraction // Fraction of DEmax allowed for ramping up
        ,double P_elec_max // Maximum electrical power
        ,double upper_threshold // Minimum energy delivery to the process
        ,double pre_dispatched_heat // Initially pre-dispatched heat
        ,double * optimalDispatch
        ,double * runtime
){
    MSG("\n\nt = %f",t0);

    double pvdstep, wndstep;
    assert(0 == motab_check_timestep(pvd,&pvdstep));
    assert(0 == motab_check_timestep(wnd,&wndstep));
    
    static MotabData *pvdcache, *wndcache;
    if(pvdcache != pvd){
        pvdcache = pvd;
        if(pvdstep != dt)ERR("Warning: pv file timestep is %f s, different"
            " from forecasting timestep %f s (message is only shown once)",pvdstep, dt);
    }
    if(wndcache != wnd){
        wndcache = wnd;
        if(wndstep != dt)ERR("Warning: wind file timestep is %fs, different"
            " from forecasting timestep %fs (message is only shown once)",wndstep, dt);
    }

    // check that pvd and wnd can cover the required time range

    double LCOH = 3.5/3.6*100.0; // Levelized cost of heat
    double MaxRUP = ramp_up_fraction * DEmax; // Maximum ramp-up rate
    double MaxRDW = ramp_dw_fraction * DEmax; // Maximum ramp-down rate
    double MinRUP = 0.0; // Minimum ramp-up rate
    double MinRDW = 0.0; // Minimum ramp-down rate
    double UPT = upper_threshold * DEmax;
    double LPT = 0.0;
    double M = 10*DEmax; // Big M value

    if(NULL==pvd)return ST_ERRVAL;
    if(NULL==wnd)return ST_ERRVAL;

    MSG("nrows pvd = %d",pvd->nrows);
    
    //assert(wnd->nrows == 8760);
    //assert(pvd->nrows == 8760);
    
    int wind_col = motab_find_col_by_label(wnd,"power");
    assert(wind_col != -1);

    int pv_col = motab_find_col_by_label(pvd,"power");
    assert(pv_col != -1);

#define WND(I) (\
    /*MSG("WND(%d) at t=%f",I,t0+((I)-1)*dt),*/\
    1e-6*P_elec_max_wind / P_elec_wind_ref_size*motab_get_value_wraparound(wnd,    t0+(I)*dt,wind_col)\
    )
/** as noted above, PV for the ith period (counting from 1) is at t0+i*dt */
#define PV(I) (\
    /*MSG("PV(%d) at t=%f",(I),t0+(I)*dt),*/\
    1e-6*P_elec_max_pv / P_elec_pv_ref_size*motab_get_value_wraparound(pvd,    t0+(I)*dt,pv_col)\
    )

    GRBenv   *env   = NULL;
    GRBmodel *P = NULL;
    int       error = 0;

    /* Create environment */
    error = GRBemptyenv(&env);
    if (error) goto ERROR;

    error = GRBsetstrparam(env, "LogFile", "mip1.log");
    if (error) goto ERROR;

    error = GRBsetdblparam(env, "MIPGap", 0.001);
    if (error) goto ERROR;

    error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
    if (error) goto ERROR;
#ifdef ST_LINPROG_DEBUG
    error = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 1);
    if (error) goto ERROR;
#endif

    error = GRBstartenv(env);
    if (error) goto ERROR;

    // Define the number of steps in the horizon
    #define N (horizon)
    // Define macros for variable indexing
    #define SL(I) (I) // Stored energy
    #define DE(I) (SL(N) + I) // Dispatched energy
    #define SE(I) (DE(N) + I) // State of energy
    #define XE(I) (SE(N) + I) // Auxiliary variable for energy
    #define YON(I) (XE(N) + I) // Disjunct mode for ON state
    #define YOFF(I) (YON(N) + I)
    #define YPAR(I) (YOFF(N) + I)
    #define ZONOFF(I) (YPAR(N) + I)
    #define DEPAR(I) (ZONOFF(N) + I)

    /* Create an empty model */
    error = GRBnewmodel(env, &P, "mip1", 0, NULL, NULL, NULL, NULL, NULL);
    if (error) goto ERROR;

    // Set column names for readability
    #define NAMEMAX 255
    char vname[NAMEMAX];

    /* VARIABLE BOUNDS AND OBJECTIVES*/
    for(int i = 1; i <= N; ++i) {
        snprintf(vname, NAMEMAX, "SL%02d", i);
        error = GRBaddvar(P, 0, NULL, NULL, eff_process * LCOH * 0.01, SLmin, SLmax, GRB_CONTINUOUS, vname);
        if (error) goto ERROR;
    }

    for(int i = 1; i <= N; ++i) {
        snprintf(vname, NAMEMAX, "DE%02d", i);
        error = GRBaddvar(P, 0, NULL, NULL, eff_process * LCOH, 0.0, DEmax, GRB_CONTINUOUS, vname);
        if (error) goto ERROR;
    }

    for(int i = 1; i <= N; ++i) {
        double pvd_i = PV(i);
        double wnd_i = WND(i);
        double p_heat_i = (pvd_i + wnd_i) * eff_heater;
        p_heat_i = fmin(p_heat_i, P_elec_max); 
        if (pvd_i <= 0 && wnd_i <= 0) {
            snprintf(vname, NAMEMAX, "SE%02d", i);
            error = GRBaddvar(P, 0, NULL, NULL, 0.0, 0.0, 0.0, GRB_CONTINUOUS, vname);
            if (error) goto ERROR;
        } else {
            snprintf(vname, NAMEMAX, "SE%02d", i);
            error = GRBaddvar(P, 0, NULL, NULL, 0.0, 0.0, p_heat_i, GRB_CONTINUOUS, vname);
            if (error) goto ERROR;
        }
    }

    for(int i = 1; i <= N; ++i) {
        double pvd_i = PV(i);
        double wnd_i = WND(i);
        double p_heat_i = (pvd_i + wnd_i) * eff_heater;
        p_heat_i = fmin(p_heat_i, P_elec_max); 
        if (pvd_i <= 0 && wnd_i <= 0) {
            snprintf(vname, NAMEMAX, "XE%02d", i);
            error = GRBaddvar(P, 0, NULL, NULL, 0.0, 0.0, 0.0, GRB_CONTINUOUS, vname);
            if (error) goto ERROR;
        } else {
            snprintf(vname, NAMEMAX, "XE%02d", i);
            error = GRBaddvar(P, 0, NULL, NULL, 0.0, 0.0, p_heat_i, GRB_CONTINUOUS, vname);
            if (error) goto ERROR;
        }
    }

    for(int i = 1; i <= N; ++i) {
        snprintf(vname, NAMEMAX, "YON%02d", i);
        error = GRBaddvar(P, 0, NULL, NULL, 0.0, 0.0, 1.0, GRB_BINARY, vname);
        if (error) goto ERROR;
    }

    for(int i = 1; i <= N; ++i) {
        snprintf(vname, NAMEMAX, "YOFF%02d", i);
        error = GRBaddvar(P, 0, NULL, NULL, 0.0, 0.0, 1.0, GRB_BINARY, vname);
        if (error) goto ERROR;
    }

    for(int i = 1; i <= N; ++i) {
        snprintf(vname, NAMEMAX, "YPAR%02d", i);
        error = GRBaddvar(P, 0, NULL, NULL, 0.0, 0.0, 1.0, GRB_BINARY, vname);
        if (error) goto ERROR;
    }

    for(int i = 1; i <= N; ++i) {
        snprintf(vname, NAMEMAX, "ZONOFF%02d", i);
        error = GRBaddvar(P, 0, NULL, NULL, -4000.0, 0.0, 1.0, GRB_BINARY, vname);
        if (error) goto ERROR;
    }

    for(int i = 1; i <= N; ++i) {
        snprintf(vname, NAMEMAX, "DEPAR%02d", i);
        error = GRBaddvar(P, 0, NULL, NULL, -1.0*eff_process * LCOH, 0.0, DEmax, GRB_CONTINUOUS, vname);
        if (error) goto ERROR;
    }

    /* Change objective sense to maximisation */
    error = GRBsetintattr(P, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
    if (error) goto ERROR;

    /* CONSTRAINTS */
    for (unsigned i = 1; i <= N; ++i) {
        double pvd_i = PV(i);
        double wnd_i = WND(i);
        double p_heat_i = (pvd_i + wnd_i) * eff_heater;
        p_heat_i = fmin(p_heat_i, P_elec_max);  

        /* 1. Dispatched energy balance (N equations) SLi - SLi-1 = SEi - DEi
        Change in storage level equals net (stored minus dispatched) energy. */
        if (i == 1) {
            error = GRBaddconstr(P, 3,
                (int[]){SL(i)-1, SE(i)-1, DE(i)-1},
                (double[]){1.0, -1.0, 1.0},
            GRB_EQUAL, SLinit, NULL);
            if (error) goto ERROR;
        } else {
            error = GRBaddconstr(P, 4,
                (int[]){SL(i)-1, SL(i - 1)-1, SE(i)-1, DE(i)-1},
                (double[]){1.0, -1.0, -1.0, 1.0},
            GRB_EQUAL, 0.0, NULL);
            if (error) goto ERROR;
        }

        /* 2. Stored energy balance (N equations) SEi + XEi = ηH·(PVi + WNDi)·Δt
        Collected energy is either stored (SE) or dumped (XE). */
        error = GRBaddconstr(P, 2,
            (int[]){SE(i)-1, XE(i)-1},
            (double[]){1.0, 1.0},
        GRB_EQUAL, p_heat_i, NULL);
        if (error) goto ERROR;
    }

    /* 3. Long-term energy balance (1 equation)  ∑(DEi - SEi) = 0       
        No net change in storage level over the forecast interval.
        -- or --
        SL(N) = SLinit. */
/*    error = GRBaddconstr(P, 1,*/
/*        (int[]){SL(N)-1},*/
/*        (double[]){1.0},*/
/*    GRB_EQUAL, SLinit, NULL);*/
/*    if (error) goto ERROR;*/

    // MODE DETECTION (YON,YOFF,YPAR)
    for (unsigned i = 1; i <= N; ++i) {
        // DE(i) >= UPT - M*[1-YON(i)]
        error = GRBaddconstr(P, 2,
            (int[]){DE(i)-1, YON(i)-1},
            (double[]){1.0, -M},
        GRB_GREATER_EQUAL, UPT-M, NULL);
        if (error) goto ERROR;
    }

    for (unsigned i = 1; i <= N; ++i) {
        // DE(i) <= UPT + M*[1-YPAR(i)]
        error = GRBaddconstr(P, 2,
            (int[]){DE(i)-1, YPAR(i)-1},
            (double[]){1.0, +M},
        GRB_LESS_EQUAL, UPT+M, NULL);
        if (error) goto ERROR;
    }

    for (unsigned i = 1; i <= N; ++i) {
        // DE(i) >= LPT - M*[1-YPAR(i)]
        error = GRBaddconstr(P, 2,
            (int[]){DE(i)-1, YPAR(i)-1},
            (double[]){1.0, -M},
        GRB_GREATER_EQUAL, LPT-M, NULL);
        if (error) goto ERROR;
    }

    for (unsigned i = 1; i <= N; ++i) {
        // DE(i) <= LPT + M*[1-YOFF(i)]
        error = GRBaddconstr(P, 2,
            (int[]){DE(i)-1, YOFF(i)-1},
            (double[]){1.0, +M},
        GRB_LESS_EQUAL, LPT+M, NULL);
        if (error) goto ERROR;
    }

    for (unsigned i = 1; i <= N; ++i) {
        error = GRBaddconstr(P, 3,
            (int[]){YON(i)-1, YOFF(i)-1, YPAR(i)-1},
            (double[]){1.0, 1.0, 1.0},
        GRB_EQUAL, 1.0, NULL);
        if (error) goto ERROR;
    }

    // ************************************
    // RESTRICTION OF RAMPING RATE
    // ************************************
    // DE(1) - pre_dispatched_heat <= MaxRUP
    error = GRBaddconstr(P, 1,
        (int[]){DE(1)-1},
        (double[]){+1.0},
    GRB_LESS_EQUAL, MaxRUP+pre_dispatched_heat, NULL);

    error = GRBaddconstr(P, 2,
        (int[]){DE(1)-1, YPAR(1)-1},
        (double[]){+1.0, -M},
    GRB_GREATER_EQUAL, -M+pre_dispatched_heat, NULL);

    for(int i = 2; i <= N; ++i) {
        // DE(i) - DE(i-1) <= MaxRUP
        error = GRBaddconstr(P, 2,
            (int[]){DE(i)-1, DE(i-1)-1},
            (double[]){+1.0, -1.0},
        GRB_LESS_EQUAL, MaxRUP, NULL);

        // M*[1-YPAR(i)] + DE(i) - DE(i-1) >= 0
        error = GRBaddconstr(P, 2,
            (int[]){DE(i)-1, DE(i-1)-1, YPAR(i)-1},
            (double[]){+1.0, -1.0, -M},
        GRB_GREATER_EQUAL, -M, NULL);
    }

    // ************************************
    // SWITCH DETECTION: ON-->OFF
    // ************************************
    //YON(0) >= ZONOFF(1)
    double YON0 = (pre_dispatched_heat >= UPT) ? 1 : 0;
    error = GRBaddconstr(P, 1,
        (int[]){ZONOFF(1)-1},
        (double[]){+1.0},
    GRB_LESS_EQUAL, YON0, NULL);
    //YON(i-1) >= ZONOFF(i)
    for(int i = 2; i <= N; ++i) {
        error = GRBaddconstr(P, 2,
            (int[]){YON(i-1)-1, ZONOFF(i)-1},
            (double[]){+1.0, -1.0},
        GRB_GREATER_EQUAL, 0.0, NULL);
    }
    //YOFF(i) >= ZONOFF(i)
    for(int i = 1; i <= N; ++i) {
        error = GRBaddconstr(P, 2,
            (int[]){YOFF(i)-1, ZONOFF(i)-1},
            (double[]){+1.0, -1.0},
        GRB_GREATER_EQUAL, 0.0, NULL);
    }
    //YON(0)+YOFF(1)-1 <= ZONOFF(1)
    error = GRBaddconstr(P, 2,
        (int[]){YOFF(1)-1, ZONOFF(1)-1},
        (double[]){+1.0, -1.0},
    GRB_LESS_EQUAL, 1.0-YON0, NULL);
    //YON(i-1)+YOFF(i)-1 <= ZONOFF(i)
    for(int i = 2; i <= N; ++i) {
        error = GRBaddconstr(P, 3,
            (int[]){YON(i-1)-1, YOFF(i)-1, ZONOFF(i)-1},
            (double[]){+1.0, +1.0, -1.0},
        GRB_LESS_EQUAL, 1.0, NULL);
    }

    // ************************************
    // FORBIDEN TRANSITION: ON-->PAR
    // ************************************
    // YON0 + YPAR(1) <= 1
    error = GRBaddconstr(P, 1,
        (int[]){YPAR(1)-1},
        (double[]){+1.0},
    GRB_LESS_EQUAL, 1.0-YON0, NULL);
    // YON(i-1) + YPAR(i) <= 1
    for(int i = 2; i <= N; ++i) {
        error = GRBaddconstr(P, 2,
            (int[]){YON(i-1)-1, YPAR(i)-1},
            (double[]){+1.0, +1.0},
        GRB_LESS_EQUAL, 1.0, NULL);
    }

    // ************************************
    // FORBIDEN TRANSITION: PAR-->OFF
    // ************************************
    // YPAR0 + YOFF(1) <= 1
    double YPAR0 = (pre_dispatched_heat < UPT && pre_dispatched_heat > LPT) ? 1 : 0;
    error = GRBaddconstr(P, 1,
        (int[]){YOFF(1)-1},
        (double[]){+1.0},
    GRB_LESS_EQUAL, 1.0-YPAR0, NULL);
    // YPAR(i-1) + YOFF(i) <= 1
    for(int i = 2; i <= N; ++i) {
        error = GRBaddconstr(P, 2,
            (int[]){YPAR(i-1)-1, YOFF(i)-1},
            (double[]){+1.0, +1.0},
        GRB_LESS_EQUAL, 1.0, NULL);
    }

    // ************************************
    // PRODUCT OF YPAR(i) AND DE(i)
    // ************************************
    // DEPAR(i) ≤ M⋅YPAR(i)
    for(int i = 1; i <= N; ++i) {
        error = GRBaddconstr(P, 2,
            (int[]){DEPAR(i)-1, YPAR(i)-1},
            (double[]){+1.0, -M},
        GRB_LESS_EQUAL, 0.0, NULL);
    }
    // DEPAR(i) ≥ DE(i) - M⋅[1-YPAR(i)]
    for(int i = 1; i <= N; ++i) {
        error = GRBaddconstr(P, 3,
            (int[]){DEPAR(i)-1, DE(i)-1, YPAR(i)-1},
            (double[]){+1.0, -1.0, -M},
        GRB_GREATER_EQUAL, -M, NULL);
    }
    // DEPAR(i) ≤ DE(i)
    for(int i = 1; i <= N; ++i) {
        error = GRBaddconstr(P, 2,
            (int[]){DEPAR(i)-1, DE(i)-1},
            (double[]){+1.0, -1.0},
        GRB_LESS_EQUAL, 0.0, NULL);
    }
    // DEPAR(i) ≥ 0.0 (Already covered)

    /* Optimize model */
    error = GRBoptimize(P);
    if (error) goto ERROR;

    /* Write model to 'mip1.lp' */
    error = GRBwrite(P, "mip1.lp");
    if (error) goto ERROR;

    // Get the optimization result
    int optimstatus;
    double objval;
    error = GRBgetintattr(P, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto ERROR;

    double pre_runtime;
    error = GRBgetdblattr(P, GRB_DBL_ATTR_RUNTIME, &pre_runtime);
    if (error) goto ERROR;
#ifdef ST_LINPROG_DEBUG
    printf("Runtime: %f seconds\n", pre_runtime);
#endif
    runtime[0] = pre_runtime;

    if (optimstatus == GRB_OPTIMAL) {
        error = GRBgetdblattr(P, GRB_DBL_ATTR_OBJVAL, &objval);
        if (error) goto ERROR;
        MSG("Optimal objective: %.4e\n", objval);

        error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, DE(1)-1, 1, optimalDispatch);
        if (error) goto ERROR;

#ifdef ST_LINPROG_DEBUG
        MSG("      \t SL \t SE \t XE \t DE \t YON \t YOFF \t YPAR \t ONOFF \t DEPAR");
        for(int i = 1; i <= N; ++i){
            double SL_val, SE_val, XE_val, DE_val, YON_val, YOFF_val, YPAR_val, ZONOFF_val, DEPAR_val;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, SL(i)-1, 1, &SL_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, SE(i)-1, 1, &SE_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, XE(i)-1, 1, &XE_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, DE(i)-1, 1, &DE_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, YON(i)-1, 1, &YON_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, YOFF(i)-1, 1, &YOFF_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, YPAR(i)-1, 1, &YPAR_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, ZONOFF(i)-1, 1, &ZONOFF_val);
            if (error) goto ERROR;
            error = GRBgetdblattrarray(P, GRB_DBL_ATTR_X, DEPAR(i)-1, 1, &DEPAR_val);
            if (error) goto ERROR;
            MSG("%3d: \t %.1f \t %.1f \t %.1f \t %.1f \t %.1f \t %.1f \t %.1f \t %.1f \t %.1f", i,
                SL_val, SE_val, XE_val, DE_val, YON_val, YOFF_val, YPAR_val, ZONOFF_val, DEPAR_val);
        }
#endif
    } else {
        MSG("Optimization was stopped with status %d\n", optimstatus);
    }

    MSG("OPTIMAL DISPATCH FOR THE NEXT HOUR: %f", optimalDispatch[0]);

ERROR:
    if (error) {
        ERR("Gurobi error: %s\n", GRBgeterrormsg(env));
        return ST_ERRVAL;
    }

    GRBfreemodel(P);
    GRBfreeenv(env);

    /* End of the code */
}

double st_mip(
         MotabData *pvd
        ,MotabData *wnd
        ,double P_elec_max_pv
        ,double P_elec_max_wind
        ,double P_elec_pv_ref_size
        ,double P_elec_wind_ref_size
        ,int horizon // Time horizon for the optimization
        ,double dt // Time step duration
        ,double t0 // Initial time
        ,double eff_heater // Efficiency of the heater
        ,double eff_process // Efficiency of the process
        ,double DEmax // Maximum dispatched energy
        ,double SLmax // Maximum stored energy (storage limit)
        ,double SLinit // Initial stored energy
        ,double SLmin // Minimum stored energy
        ,double ramp_up_fraction // Fraction of DEmax allowed for ramping up
        ,double ramp_dw_fraction // Fraction of DEmax allowed for ramping up
        ,double P_elec_max // Maximum electrical power
        ,double upper_threshold // Minimum energy delivery to the process
        ,double pre_dispatched_heat // Initially pre-dispatched heat
        ,double pre_Q_flow_dis
        ,double pre_blk_state
        ,double pre_startup_next
        ,double t_shutdown_min
        ,double * optimalDispatch
        ,double * t_start_up_next
        ,double * runtime
){
    double pvdstep, wndstep;
    assert(0 == motab_check_timestep(pvd,&pvdstep));
    assert(0 == motab_check_timestep(wnd,&wndstep));
    
    static MotabData *pvdcache, *wndcache;
    if(pvdcache != pvd){
        pvdcache = pvd;
        if(pvdstep != dt)ERR("Warning: pv file timestep is %f s, different"
            " from forecasting timestep %f s (message is only shown once)",pvdstep, dt);
    }
    if(wndcache != wnd){
        wndcache = wnd;
        if(wndstep != dt)ERR("Warning: wind file timestep is %fs, different"
            " from forecasting timestep %fs (message is only shown once)",wndstep, dt);
    }

    // check that pvd and wnd can cover the required time range

    double LCOH = 3.5/3.6*100.0; // Levelized cost of heat
    double MaxRUP = ramp_up_fraction * DEmax; // Maximum ramp-up rate
    double MaxRDW = ramp_dw_fraction * DEmax; // Maximum ramp-down rate
    double MinRUP = 0.0; // Minimum ramp-up rate
    double MinRDW = 0.0; // Minimum ramp-down rate
    double UPT = upper_threshold * DEmax;
    double LPT = 0.0;
    double M = 10*DEmax; // Big M value

    // Print the initial values for debugging purposes
    MSG("t                = %.2f s", t0);
    MSG("dt               = %.2f s", dt);
    MSG("P_elec_max       = %.2f MWt", P_elec_max);
    MSG("eff_heater       = %.2f MWt", eff_heater);
    MSG("eff_process      = %.2f MWt", eff_process);
    MSG("DEmax            = %.2f MWth", DEmax);
    MSG("SLmax            = %.2f MWhth", SLmax);
    MSG("SLinit           = %.2f MWhth", SLinit);
    MSG("SLmin            = %.2f MWhth", SLmin);
    MSG("ramp_up_fraction = %.2f", ramp_up_fraction);
    MSG("Max ramp-up rate = %.2f", MaxRUP);
    MSG("Max ramp-dw rate = %.2f", MaxRDW);
    MSG("DEinit           = %.2f", pre_dispatched_heat);
    MSG("pre_Q_flow_dis   = %.2f", pre_Q_flow_dis);
    MSG("pre_blk_state    = %.2f", pre_blk_state);

    if(NULL==pvd)return ST_ERRVAL;
    if(NULL==wnd)return ST_ERRVAL;

    int wind_col = motab_find_col_by_label(wnd,"power");
    assert(wind_col != -1);

    int pv_col = motab_find_col_by_label(pvd,"power");
    assert(pv_col != -1);

#define WND1(I) (\
    /*MSG("WND(%d) at t=%f",I,t0+((I)-1)*dt),*/\
    1e-6*P_elec_max_wind / P_elec_wind_ref_size*motab_get_value_wraparound(wnd,    t0+(I)*dt,wind_col)\
    )
/** as noted above, PV for the ith period (counting from 1) is at t0+i*dt */
#define PV1(I) (\
    /*MSG("PV(%d) at t=%f",(I),t0+(I)*dt),*/\
    1e-6*P_elec_max_pv / P_elec_pv_ref_size*motab_get_value_wraparound(pvd,    t0+(I)*dt,pv_col)\
    )

    // Define the number of steps in the horizon
    #define N (horizon)

    // This piece of code determines the next startup time
    if (pre_blk_state == 4){
        t_start_up_next[0] = t0 + t_shutdown_min;
    } else {
        t_start_up_next[0] = pre_startup_next;
    }
    MSG("pre_blk_state: %.f",pre_blk_state);
    // While in state 1, we check if we can start the plant
    if (pre_blk_state == 1 && t0 >= t_start_up_next[0]){
        double SL[N+1];
        SL[0] = SLinit;
        int i = 1;
        int start_ramp_up = 1;
        while (i <= N) {
            double pvd_i = PV1(i);
            double wnd_i = WND1(i);
            double p_heat_i = (pvd_i + wnd_i) * eff_heater;
            p_heat_i = fmin(p_heat_i, P_elec_max);
            SL[i] = SL[i-1] + p_heat_i - MaxRUP * i;
            if (SL[i] <= SLmin && i<16) {
                // The plant is not ready to start
                start_ramp_up = 0;
                optimalDispatch[0] = 0.0;
                break;
            }
            ++i;
        }
        MSG("start_ramp_up: %d",start_ramp_up);
        if (start_ramp_up == 1){
            optimalDispatch[0] = pre_dispatched_heat + MaxRUP;
        }
    }
    if (pre_blk_state == 2){
        optimalDispatch[0] = fmin(pre_dispatched_heat + MaxRUP, DEmax);
    }

    if ((int)t0 % 86400 == 0) {
        printf("Day %.1f\n",t0/86400);
    }
    runtime[0] = 0;

    if (pre_blk_state == 3){
    gurobi_fun(
         pvd
        ,wnd
        ,P_elec_max_pv
        ,P_elec_max_wind
        ,P_elec_pv_ref_size
        ,P_elec_wind_ref_size
        ,horizon // Time horizon for the optimization
        ,dt // Time step duration
        ,t0 // Initial time
        ,eff_heater // Efficiency of the heater
        ,eff_process // Efficiency of the process
        ,DEmax // Maximum dispatched energy
        ,SLmax // Maximum stored energy (storage limit)
        ,SLinit // Initial stored energy
        ,SLmin // Minimum stored energy
        ,ramp_up_fraction // Fraction of DEmax allowed for ramping up
        ,ramp_dw_fraction // Fraction of DEmax allowed for ramping up
        ,P_elec_max // Maximum electrical power
        ,upper_threshold // Minimum energy delivery to the process
        ,pre_dispatched_heat // Initially pre-dispatched heat
        ,optimalDispatch
        ,runtime
    );
    }

}

// vim: ts=4:sw=4:noet:tw=80
