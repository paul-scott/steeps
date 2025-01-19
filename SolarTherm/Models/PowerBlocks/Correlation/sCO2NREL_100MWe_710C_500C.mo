within SolarTherm.Models.PowerBlocks.Correlation;
model sCO2NREL_100MWe_710C_500C "SCO2 cycle based on regression of NREL PB, Sodium HTF"
	extends Cycle(T_in_ref = 710.0 + 273.15);
	import Modelica.SIunits.Conversions.*;
	input SI.Temperature T_amb;
	parameter SI.Temperature T_amb_ref = 35.0 + 273.15 "Design ambient temperature";
    Real t_in "Non-dimensional inlet temperature, called theta_in";
    Real t_amb "Non-dimensional ambient temperature, called theta_amb";
    parameter Real[14] a_array = {a_0,a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8,a_9,a_10,a_11,a_12,a_13} "Array of a";
    parameter Real[10] b_array = {b_0,b_1,b_2,b_3,b_4,b_5,b_6,b_7,b_8,b_9} "Array of b";
    Real[14] w_array = {1.0,t_in,load,load^3,t_in*t_amb,(t_in*t_amb)^2,(t_in*t_amb)^3,t_amb*load,(t_amb*load)^2,(t_amb*load)^3,(t_in*load)^3,t_in/t_amb,load/t_in,load/t_amb} "Array of w variables";
    Real[10] x_array = {1.0,t_amb,t_amb^2,t_amb^3,load^3,(t_amb*load)^2,(t_in*load)^2,(t_in*load)^3,t_amb/t_in,load/t_amb} "Array of x variables";
    
    
    //Constants of a refer to regression constants of k_q
    constant Real a_0 = 25.58108;
    constant Real a_1 = -20.89183;
    constant Real a_2 = 30.69272;
    constant Real a_3 = 1.45557;
    constant Real a_4 = -70.22424;
    constant Real a_5 = 91.88505;
    constant Real a_6 = -34.90353;
    constant Real a_7 = -11.03358;
    constant Real a_8 = 1.38411;
    constant Real a_9 = -0.88088;
    constant Real a_10 = -1.39702;
    constant Real a_11 = 8.66269;
    constant Real a_12 = -10.22232;
    constant Real a_13 = -9.14393;
    //Obtained by Zeb using StatGraphics Software
    //Adjusted R2 = 99.0553%
    //Standard Error = 0.0210561
    //Mean Absolute Error = 0.0158249
    
    //Constants of b refer to regression constants of k_w
    constant Real b_0 = 8.6325;
    constant Real b_1 = -29.5362;
    constant Real b_2 = 36.2477;
    constant Real b_3 = -14.1233;
    constant Real b_4 = -0.1658;
    constant Real b_5 = 0.2867;
    constant Real b_6 = -0.9962;
    constant Real b_7 = 0.4034;
    constant Real b_8 = -0.9911;
    constant Real b_9 = 0.7146;
    //Obtained by Zeb using StatGraphics Software
    //Adjusted R2 = 98.9686%
    //Standard Error = 0.00275459
    //Mean Absolute Error = 0.0020361
    parameter Real eta_q_des = sum(a_array) "Design HX efficiency";
    parameter Real eta_w_des = sum(b_array) "Design Cycle efficiency";
    parameter Real eta_gross_des = eta_q_des*eta_w_des "Design gross efficiency, product of eta_q and eta_w";
    /*
    //These constants are for eta_Q or k_q
    constant Real a_0 = -7.063737;
    constant Real a_1 = 3.363395;
    constant Real a_2 = 8.051343;
    constant Real a_3 = 1.454204;
    constant Real a_4 = -4.473867;
    constant Real a_5 = -0.382406;
    //These constants are for eta_gross or k_w
    constant Real b_0 = 8.237821;
    constant Real b_1 = 4.357904;
    constant Real b_2 = -2.074309;
    constant Real b_3 = -34.690258;
    constant Real b_4 = 39.327013;
    constant Real b_5 = -14.842294;
    constant Real b_6 = 0.301619;
    constant Real b_7 = -0.141392;
    //Regression design values (when all variables are 1)
    parameter Real eta_q_des = 0.948933 "Design HX efficiency";
    parameter Real eta_w_des = 0.476106 "Design Cycle efficiency";
    parameter Real eta_gross_des = 0.451793 "Design gross efficiency, product of eta_q and eta_w";
    */
//algorithm
    //k_q := sum(a_array.*w_array);
    //k_w := sum(b_array.*x_array);
equation
    //Note that k_q and k_w get passed as eta_Q and eta_gross
    t_in = T_in/T_in_ref;
    t_amb = T_amb/T_amb_ref;
    //k_q = 0.98;
    //k_w = 0.45;
    k_q = sum(a_array.*w_array);
    k_w = sum(b_array.*x_array);
    /*
    k_q = a_0 + a_1*(T_ND) + a_2*(T_amb_ND) + a_3*(T_amb_ND^2) + a_4*(load) + a_5*(load^2);
    k_w = b_0 + b_1*(T_ND) + b_2*(T_ND^2) + b_3*(T_amb_ND) + b_4*(T_amb_ND^2) + b_5*(T_amb_ND^3) + b_6*(load) + b_7*(load^2);
    */
end sCO2NREL_100MWe_710C_500C;