within SolarTherm.Calculators;

model Piping_Cost_TwoLayer_CS_TCF
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import MA = Modelica.Blocks.Math;
  
  //Input Parameters
  parameter SI.Temperature T1 = 1100.0 + 273.15 "Inner Pipe temperature (K)";
  parameter SI.MassFlowRate m_flow_gas = 54.8*2/2 "Maximum mass flow rate of air (kg/s)";  
  
  parameter SI.Temperature T2_Max = 1260.0 + 273.15 "Maximum tolerable insulation temperature (K)";
  
  parameter Real F_m_duct = 1.0 "Duct material factor (galvanised steel = 1.0, stainless steel = 4.3)";
  
  
  //Properties and Costs
  //Cubic coefficients of Refractory (High-Temp Material)
  parameter Real A = 4.675E-2;
  parameter Real B = -2.441E-5;
  parameter Real C = 1.806E-7;
  parameter Real D = -6.303E-11;
  
  //Cubic coefficients of Insulation (Lower-Temp Material)
  //Thermal Ceramic Blanket
  parameter Real E = -5.463E-2;
  parameter Real F = 2.0E-4;
  parameter Real G = 0.0;
  parameter Real H = 0.0;

/*
  //Cubic coefficient of Insulation (RockWool)
  parameter Real E = -5.463E-2;
  parameter Real F = 2.0E-4;
  parameter Real G = 0.0;
  parameter Real H = 0.0;
*/
  
  parameter Real c1 = 7433.57 "Cost per unit volume of  Refractory (USD_2022)";
  parameter Real c2 = 124.33 "Cost per unit volume of Thermal ceramic fibre Insulation (USD_2024)";
  //parameter Real c2 = 66.28 "Cost per unit volume of RockWool Insulation (USD_2022)";
  
  parameter SI.Length d1 = (4.0 * m_flow_gas / (CN.pi * rho_air * u_air)) ^ 0.5 "Minimum Pipe diameter (m)";
  parameter SI.Temperature T4 = 25.0 + 273.15 "Ambient temperature (K)";
  parameter SI.CoefficientOfHeatTransfer h_amb = 10.0 "Convective heat transfer coefficient of the outer surface (K)";

  
  
  //Receiver Thermal Output is 4.0485 MWth
  //Downcomer length is  86.7315 m, Outer diameter is 0.864m
  //Inputs
  //parameter SI.MassFlowRate m_flow_air =  "Mass flow rate of air (kg/s)";
  parameter SI.Velocity u_air = 20.0 "Maximum speed of air in chamber (m/s)"; //Turton pipe heuristics
  parameter SI.Density rho_air = SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.rho_T(T1) "Density of air (kg/m3)";
    parameter SI.DynamicViscosity vis_air=SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.mu_T(T1) "Dynamic viscosity (Ns/m2) of air at ambient pressure as a function of temperature";
  
  //parameter SI.Area 
  //System is ~18.5 MW tentative
  //0.1% per meter?
  
  // Pressure drop alone a 50 m pipe
  parameter SI.Length L=50 "Length of the pipe";  
  parameter Real Re = rho_air*u_air*d1/vis_air "Reynolds Number";
  parameter SI.Length roughness= 0.26e-3 "Roughness of cast iron 0.26 (mm), concrete 0.3--3, since the refractory is made of calcium silica, the roughness of 0.3 is probably approperate. https://www.pipeflow.com/pipe-pressure-drop-calculations/pipe-roughness";
  parameter Real relative_roughness=roughness/d1 "Relative roughness";
  parameter Real f= if Re<2300 then 64/Re else 0.25/(log(relative_roughness/3.7+5.74/(Re^0.9)))^2 "friction factor, from Moody chart https://www.pipeflow.com/pipe-pressure-drop-calculations/pipe-friction-factors";
  parameter SI.Pressure dP = f*(L/d1)*(rho_air*u_air^2/2) "Pressure drop alone the pipe (Pa)";
  
  SI.Length d3(min = 0.0) "Total pipe diameter, including the insulating layers (m)";
  SI.Length d2(min = 0.0) "Total pipe diameter, including the 1st insulating layer (m)";
  SI.Length t_sum(min = 0.0) "Total thickness of refractory and insulation (m)";
  //parameter Real QpL = 0.0001*18.5e6/1.0 "Acceptable Heat Rate per metre of piping (W/m)";
  parameter SI.Temperature T3 = 50.0 + 273.15;

  //Calculation
  SI.Temperature T2(min=273.15, start = 1080.0+273.5);
  //SI.Temperature T3 (start = 200.0 + 273.15);
  parameter Real r1 = d1/2.0;
  
  SI.Length r2(min = 1e-6,start = r1 + 1e-2);
  SI.Length r3(min = 1e-6,start = r1 + 1e-2+0.8);
  SI.Length t1(min = 0.0,start = 1e-2) "Thickness of Refractory";
  SI.Length t2(min = 1e-6, start = 0.8) "Thickness of Insulation";

  
  Real CpL_1 "Cost per length of piping due to inner layer (USD/m)";
  Real CpL_2 "Cost per length of piping due to outer layer (USD/m)";
  Real CpL_duct "FCI Cost per length of piping due to outer duct layer (USD/m)";
  Real Penalty "Very large cost penalty for invalid design (USD/m)";
  
  Real CpL_total "Cost per length of total piping (USD/m)";
  //parameter Real QpL = CN.pi*d1*h_amb*(373.15-T4) "Heat loss rate per unit length (W/m)"; 
  Real QpL (min = 1e-12);
  
algorithm
  if (T2 < T4) or (T3 < T4) or (r3 < r2) or (t2 < 1e-6) then
    terminate("Invalid");
  end if;
  
equation
  if T2 > T2_Max then
    Penalty = 100000.0;
  else
    Penalty = 0.0;
  end if;
  //Sweep t1
  t1 = 1e-2 + time*1e-3;
  //der(t_sum) = 0.01;
  r2 = r1 + t1;
  r3 = r2 + t2;
  t1 + t2 = t_sum;
  
  d3 = 2.0*r3;
  d2 = 2.0*r2;
  
  
  //QpL=(2.0*CN.pi/log(r1/r2))*(((A/1)*(T2^1-T1^1))+((B/2)*(T2^2-T1^2))+((C/3)*(T2^3-T1^3))+((D/4)*(T2^4-T1^4)));
  QpL=(2.0*CN.pi/log(r2/r3))*(((E/1)*(T3^1-T2^1))+((F/2)*(T3^2-T2^2))+((G/3)*(T3^3-T2^3))+((H/4)*(T3^4-T2^4)));
  (2.0*CN.pi/log(r1/r2))*(((A/1)*(T2^1-T1^1))+((B/2)*(T2^2-T1^2))+((C/3)*(T2^3-T1^3))+((D/4)*(T2^4-T1^4))) = (2.0*CN.pi/log(r2/r3))*(((E/1)*(T3^1-T2^1))+((F/2)*(T3^2-T2^2))+((G/3)*(T3^3-T2^3))+((H/4)*(T3^4-T2^4)));
  QpL=(CN.pi*d3*h_amb)*(T3-T4);


  CpL_1 = CN.pi*(r2^2-r1^2)*c1;
  CpL_2 = CN.pi*(r3^2-r2^2)*c2;
  CpL_duct = 1.18*(816.0/1000.0)*F_m_duct*120.1539*(d3^2.14);
  CpL_total = CpL_1 + CpL_2 + CpL_duct + Penalty;

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),experiment(StopTime = 1000, StartTime = 0, Tolerance = 1.0e-5, Interval = 1, maxStepSize = 1, initialStepSize = 1));
end Piping_Cost_TwoLayer_CS_TCF;