within SolarTherm.Calculators;

model Planar_Internal_Firebrick_Insulation "Sweep U with time, solves for t1 for each U"
  import CN = Modelica.Constants;
  import MA = Modelica.Blocks.Math;
  import SI = Modelica.SIunits;
  //Polynomial coefficients of firebrick insulation
  // k(W/mK) = A + BT + CT^2 where T is in Kelvins
  parameter Real A = 0.072685;
  parameter Real B = 0.0001;
  parameter Real C = 0.0;
  
  parameter Real c1 = 102.375 "USD/m3 Firebrick";
  parameter SI.CoefficientOfHeatTransfer hc = 10.0 "W/m2K ambient convection";
  
  parameter SI.Temperature T_amb = 25.0 + 273.15 "(K) ambient";
  parameter SI.Temperature T0_fixed = 1200.0 + 273.15 "(K)";
  
  SI.Temperature T0(start = T0_fixed) "Storage temperature (K)";
  
  //Sweep this
  Real y(start=-3);
  SI.CoefficientOfHeatTransfer U= 10.0^y;
  Real R = 1.0/U;
  SI.Temperature T1;
  //Vary this manually
  //parameter Real U = 3.0 "W/m2K Target U value";

  //Calculation
  Real t1;

  Real CpA "Cost per Area";
  Real q "Heat loss rate per cross-section area";

equation
  y = -2.0+time*0.1;

  der(T0) = 0.0;
  q = (A*(T0-T1)/(t1))+(B*(T0^2.0-T1^2.0)/(2.0*t1))+(C*(T0^3.0-T1^3.0)/(3.0*t1));

  q = hc*(T1-T_amb);
  
  CpA = max(0.0,c1*t1);

  q = U*(T0-T_amb);
annotation(experiment(StopTime = 20, StartTime = 0, Tolerance = 1e-6, Interval = 1));
end Planar_Internal_Firebrick_Insulation;