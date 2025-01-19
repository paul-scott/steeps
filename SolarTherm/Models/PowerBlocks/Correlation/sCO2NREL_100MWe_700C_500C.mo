within SolarTherm.Models.PowerBlocks.Correlation;
model sCO2NREL_100MWe_700C_500C "SCO2 cycle based on regression of NREL PB, Sodium HTF"
	import Modelica.SIunits.Conversions.*;
	input SI.Temperature T_amb;
	extends Cycle(T_in_ref = 700.0 + 273.15);
	parameter SI.Temperature T_amb_ref = 35.0 + 273.15 "Design ambient temperature";
    Real t_in "Non-dimensional inlet temperature, called theta_in";
    Real t_amb "Non-dimensional ambient temperature, called theta_amb";
    parameter Real[14] a_array = {a_0,a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8,a_9,a_10,a_11,a_12,a_13} "Array of a";
    parameter Real[10] b_array = {b_0,b_1,b_2,b_3,b_4,b_5,b_6,b_7,b_8,b_9} "Array of b";
    Real[14] w_array = {1.0,t_in,load,load^3,t_in*t_amb,(t_in*t_amb)^2,(t_in*t_amb)^3,t_amb*load,(t_amb*load)^2,(t_amb*load)^3,(t_in*load)^3,t_in/t_amb,load/t_in,load/t_amb} "Array of w variables";
    Real[10] x_array = {1.0,t_amb,t_amb^2,t_amb^3,load^3,(t_amb*load)^2,(t_in*load)^2,(t_in*load)^3,t_amb/t_in,load/t_amb} "Array of x variables";
    
    
    //Constants of a refer to regression constants of k_q
    constant Real a_0 = 25.3559;
    constant Real a_1 = -20.8812;
    constant Real a_2 = 30.5759;
    constant Real a_3 = 1.49268;
    constant Real a_4 = -69.4431;
    constant Real a_5 = 91.1203;
    constant Real a_6 = -34.6784;
    constant Real a_7 = -10.7473;
    constant Real a_8 = 1.34033;
    constant Real a_9 = -0.862459;
    constant Real a_10 = -1.43521;
    constant Real a_11 = 8.62727;
    constant Real a_12 = -10.5658;
    constant Real a_13 = -8.93598;
    //Obtained by Zeb using StatGraphics Software
    //Adjusted R2 = 99.0044%
    //Standard Error = 0.0216871
    //Mean Absolute Error = 0.0163272
    
    //Constants of b refer to regression constants of k_w
    constant Real b_0 = 8.73027;
    constant Real b_1 = -29.8164;
    constant Real b_2 = 36.5692;
    constant Real b_3 = -14.236;
    constant Real b_4 = -0.156959;
    constant Real b_5 = 0.281185;
    constant Real b_6 = -0.977885;
    constant Real b_7 = 0.38809;
    constant Real b_8 = -1.01536;
    constant Real b_9 = 0.706115;
    //Obtained by Zeb using StatGraphics Software
    //Adjusted R2 = 99.0054%
    //Standard Error = 0.00274429
    //Mean Absolute Error = 0.00206392
    
    parameter Real eta_q_des = sum(a_array) "Design HX efficiency";
    parameter Real eta_w_des = sum(b_array) "Design Cycle efficiency";
    parameter Real eta_gross_des = eta_q_des*eta_w_des "Design gross efficiency, product of eta_q and eta_w";
    /*
    //These constants are for eta_Q or k_q
    constant Real a_0 = -7.252764;
    constant Real a_1 = 3.652907;
    constant Real a_2 = 7.853152;
    constant Real a_3 = -4.364979;
    constant Real a_4 = 1.433866;
    constant Real a_5 = -0.374598;
    //These constants are for eta_gross or k_w
    constant Real b_0 = 8.102387;
    constant Real b_1 = 5.380067;
    constant Real b_2 = -2.580075;
    constant Real b_3 = -35.960902;
    constant Real b_4 = 40.725961;
    constant Real b_5 = -15.354071;
    constant Real b_6 = 0.307790;
    constant Real b_7 = -0.145251;
    //Regression design values (when all variables are 1)
    parameter Real eta_q_des = 0.947584 "Design HX efficiency";
    parameter Real eta_w_des = 0.475905 "Design Cycle efficiency";
    parameter Real eta_gross_des = 0.450960 "Design gross efficiency, product of eta_q and eta_w";
    */
equation
    //Note that k_q and k_w get passed as eta_Q and eta_gross
    t_in = T_in/T_in_ref;
    t_amb = T_amb/T_amb_ref;
    k_q = sum(a_array.*w_array);
    k_w = sum(b_array.*x_array);
    /*
    k_q = a_0 + a_1*(T_ND) + a_2*(T_amb_ND) + a_3*(T_amb_ND^2) + a_4*(load) + a_5*(load^2);
    k_w = b_0 + b_1*(T_ND) + b_2*(T_ND^2) + b_3*(T_amb_ND) + b_4*(T_amb_ND^2) + b_5*(T_amb_ND^3) + b_6*(load) + b_7*(load^2);
    */
end sCO2NREL_100MWe_700C_500C;