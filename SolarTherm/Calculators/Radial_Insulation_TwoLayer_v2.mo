within SolarTherm.Calculators;

model Radial_Insulation_TwoLayer_v2
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import MA = Modelica.Blocks.Math;
  
  //Design parameters
  parameter SI.Temperature T1 = 1100.0 + 273.15 "Inner Pipe temperature (K)";
  parameter SI.MassFlowRate m_flow_gas = 26.6 "Maximum gas (kg/s)";
  
  parameter SI.Temperature T2_Max = 700.0 + 273.15 "Maximum tolerable insulation temperature (K)";
  parameter SI.Temperature T3_Max = 50.0+273.15;
  
  //Cubic coefficients of Refractory (High-Temp Material)
  parameter Real A = 4.6752E-2;
  parameter Real B = -2.4434E-5;
  parameter Real C = 1.8063E-7;
  parameter Real D = -6.3000E-11;
  //Cubic coefficients of Insulation (Lower-Temp Material)
  parameter Real E = -5.4630E-2;
  parameter Real F = 2.000E-4;
  parameter Real G = 0.0;
  parameter Real H = 0.0;
  parameter Real c1 = 6457.0 "Cost per unit volume of  Refractory (USD)";
  parameter Real c2 = 108.0 "Cost per unit volume of Insulation (USD)";
  parameter SI.Temperature T4 = 25.0 + 273.15 "Ambient temperature (K)";
  parameter SI.CoefficientOfHeatTransfer h_amb = 10.0 "Convective heat transfer coefficient of the outer surface (K)";
  
  
  //Receiver Thermal Output is 4.0485 MWth
  //Downcomer length is  86.7315 m, Outer diameter is 0.864m
  //Inputs
  //parameter SI.MassFlowRate m_flow_air =  "Mass flow rate of air (kg/s)";
  parameter SI.Velocity u_air = 61.0 "Maximum speed of air in chamber (m/s)";
  //Turton pipe heuristics
  parameter SI.Density rho_air = SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.rho_T(T1) "Density of air (kg/m3)";
  //parameter SI.Area
  //System is ~18.5 MW tentative
  //0.1% per meter?
  parameter SI.Length d1 = (4.0 * m_flow_gas / (CN.pi * rho_air * u_air)) ^ 0.5 "Minimum Pipe diameter (m)";
  SI.Length d3 = 2.0 * r3 "Total pipe diameter, including the insulating layers (m)";
  SI.Length t_sum "Total thickness of refractory and insulation (m)";
  //parameter Real QpL = 0.0001*18.5e6/1.0 "Acceptable Heat Rate per metre of piping (W/m)";
  //parameter SI.Temperature T3 = 100.0 + 273.15;
  //Calculation
  SI.Temperature T2;
  SI.Temperature T3(start = 200.0 + 273.15);
  
  parameter Real r1 = d1 / 2.0;
  SI.Length r2(min = 0.0,start = d1 / 2.0 + 2e-1);
  SI.Length r3(min = 0.0);
  SI.Length t1(min = 0.0) "Thickness of Refractory";
  SI.Length t2(min = 0.0, start = 1e-1) "Thickness of Insulation";
  Real CpL "Cost per length of piping (USD/m)";
  Real Penalty "Very large cost penalty for invalid design (USD/m)";
  //parameter Real QpL = CN.pi * d1 * h_amb * (373.15 - T4) "Heat loss rate per unit length (W/m)";
  Real QpL;
algorithm
  if T2 < T4 or T3 < T4 or r3 < r2 or t2 < 1e-6 then
    terminate("Invalid");
  end if;
equation
  if T3 > T3_Max then
    Penalty = 100000.0;
  else
    Penalty = 0.0;
  end if;
//Sweep t1
  //t1 = 1e-3 + time * 1e-3;
  t2 = 1e-3 + time*1e-3;
//der(t_sum) = 0.01;
  r2 = r1 + t1;
  r3 = r2 + t2;
  t1 + t2 = t_sum;
  T2 = T2_Max;
  QpL=(2.0*CN.pi/log(r1/r2))*(((A/1)*(T2^1-T1^1))+((B/2)*(T2^2-T1^2))+((C/3)*(T2^3-T1^3))+((D/4)*(T2^4-T1^4)));
  QpL=(2.0*CN.pi/log(r2/r3))*(((E/1)*(T3^1-T2^1))+((F/2)*(T3^2-T2^2))+((G/3)*(T3^3-T2^3))+((H/4)*(T3^4-T2^4)));
  //QpL = 2.0 * CN.pi / log(r1 / r2) * (A / 1 * (T2 ^ 1 - T1 ^ 1) + B / 2 * (T2 ^ 2 - T1 ^ 2) + C / 3 * (T2 ^ 3 - T1 ^ 3));
  //QpL = 2.0 * CN.pi / log(r2 / r3) * (E / 1 * (T3 ^ 1 - T2 ^ 1) + F / 2 * (T3 ^ 2 - T2 ^ 2) + G / 3 * (T3 ^ 3 - T2 ^ 3));
  QpL = CN.pi * d3 * h_amb * (T3 - T4);
  CpL = CN.pi * (r2 ^ 2 - r1 ^ 2) * c1 + CN.pi * (r3 ^ 2 - r2 ^ 2) * c2 + Penalty;
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    experiment(StopTime = 1000, StartTime = 0, Tolerance = 1.0e-5, Interval = 1, maxStepSize = 1, initialStepSize = 1));
end Radial_Insulation_TwoLayer_v2;