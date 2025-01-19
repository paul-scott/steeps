within SolarTherm.Models.Chemistry.H2DRI;

model Reaction_Extent_Regression  
  //Note that this model is really bad, pwease do not use.
  //Inputs
  SI.Temperature T_rxn "Reaction temperature (K)"; //horizontal table axis u2
  Real t = T_rxn/1000.0 "Reaction temperature divided by 1000";
  Real r_input "Ratio of H2 to Fe2O3 (-)"; //vertical table axis u1
  
  //Raw Outputs
  Real xi_1_raw "Extent of reaction one";
  Real xi_2_raw "Extent of reaction two";
  Real xi_3_raw "Extent of reaction three";
  
  //Real Outputs
  Real xi_1 "Extent of reaction one after checking bounds";
  Real xi_2 "Extent of reaction two after checking bounds";
  Real xi_3 "Extent of reaction three after checking bounds";
  
  Real x1 = t;
  Real x2 = t^2.0;
  Real x3 = t^3.0;
  Real x4 = 1.0/t;
  Real x5 = 1.0/(t^2.0);
  Real x6 = r_input^2.0;
  Real x7 = r_input*t;
  Real x8 = (r_input^2.0)*t;
  
protected
  constant Real a0 = -2419.14785473075;
  constant Real a1 = 2337.47162756575;
  constant Real a2 = -1119.80638229692;
  constant Real a3 = 212.9770219969;
  constant Real a4 = 1241.8324126233;
  constant Real a5 = -252.4600235929;
  constant Real a6 = -0.0019612343;
  constant Real a7 = 0.0098496379;
  constant Real a8 = 0.0018203787;
  
  constant Real b0 = 13894.9914367608;
  constant Real b1 = -13424.5805400223;
  constant Real b2 = 6434.5016436538;
  constant Real b3 = -1224.6745302312;
  constant Real b4 = -7128.3585451157;
  constant Real b5 = 1448.8197053121;
  constant Real b6 = 0.0085022139;
  constant Real b7 = 0.0448779496;
  constant Real b8 = -0.0095800269;
  
  constant Real c0 = -5369.3437707716;
  constant Real c1 = 5188.373480934;
  constant Real c2 = -2486.8468533033;
  constant Real c3 = 473.1578865384;
  constant Real c4 = 2754.31318468;
  constant Real c5 = -560.0172188843;
  constant Real c6 = -0.0091673729;
  constant Real c7 = 0.1136891087;
  constant Real c8 = 0.0079547496;
  
algorithm    
    xi_1 := min(max(0.0,xi_1_raw),1.0);// 0.0 <= xi_1 <= 1.0
    xi_2 := min(max(0.0,xi_2_raw),xi_1);// 0.0 <= xi_2 <= xi_1
    xi_3 := min(max(0.0,xi_3_raw),xi_2);// 0.0 <= xi_3 <= xi_2
equation
    xi_1_raw = a0 + a1*x1 + a2*x2 + a3*x3 + a4*x4 + a5*x5 + a6*x6 + a7*x7 + a8*x8;
    xi_2_raw = b0 + b1*x1 + b2*x2 + b3*x3 + b4*x4 + b5*x5 + b6*x6 + b7*x7 + b8*x8;
    xi_3_raw = c0 + c1*x1 + c2*x2 + c3*x3 + c4*x4 + c5*x5 + c6*x6 + c7*x7 + c8*x8;
    
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Reaction_Extent_Regression;