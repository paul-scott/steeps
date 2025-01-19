within SolarTherm.Validation.Gan_HBS;

model HBS_Validation_Charging_Heatloss
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  
  package Medium = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_107kPa;
  package Fluid_Package = SolarTherm.Validation.Gan_HBS.Materials.Flue_Gas_107kPa;
  
  //Validation Parameters
  parameter SI.Length H_tank = 37.34 "Tank height (m)";
  parameter SI.Length D_tank = 7.748 "Tank diameter (m)";
  parameter Real ar = H_tank/D_tank;
  
  parameter SI.Length d_p = 19.6e-3 "Filler pore diameter (m)";
  parameter SI.Length s_p = 31.7e-3 "Filler pore separation (m)";
  
  parameter SI.Temperature T_max = 1326.85 + 273.15 "Maximum design system temperature (K)";
  parameter SI.Temperature T_min = 224.85 + 273.15 "Minimum design system temperature (K)";
  
  parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.576 "Heat loss coeff of surfaces (W/m2K)";
  
  parameter SI.Time t_inlet_data[15] = {0.00, 159.71, 407.05, 655.74, 892.23, 1203.97, 1788.11, 2392.57, 2978.07, 3581.18, 4184.97, 4779.95, 5383.73, 5978.71, 6162.00} "Inlet temperature signal, time axis (s)";
  parameter SI.Temperature T_inlet_data[15] = {1437.59, 1468.91, 1517.41, 1546.29, 1565.01, 1577.83, 1590.09, 1594.78, 1596.91, 1595.58, 1594.08, 1595.58, 1595.77, 1596.82, 1598.21} "Inlet temperature signal, temperature axis (K)";
  
  parameter SI.SpecificEnthalpy Delta_h1 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_1.h_Tf(T_max,0.0) - SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_1.h_Tf(T_min,0.0);
  parameter SI.SpecificEnthalpy Delta_h2 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_2.h_Tf(T_max,0.0) - SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_2.h_Tf(T_min,0.0);
  parameter SI.SpecificEnthalpy Delta_h3 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_3.h_Tf(T_max,0.0) - SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_3.h_Tf(T_min,0.0);
  parameter SI.SpecificEnthalpy Delta_h4 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_4.h_Tf(T_max,0.0) - SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_4.h_Tf(T_min,0.0);
  parameter SI.SpecificEnthalpy Delta_h5 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_5.h_Tf(T_max,0.0) - SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_5.h_Tf(T_min,0.0);
  parameter SI.SpecificEnthalpy Delta_h6 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_6.h_Tf(T_max,0.0) - SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_6.h_Tf(T_min,0.0);
  
  parameter SI.Mass m_p1[5] = HBS.Tank_A.m_p[1:5];
  parameter SI.Mass m_p2[37] = HBS.Tank_A.m_p[6:42];
  parameter SI.Mass m_p3[16] = HBS.Tank_A.m_p[43:58];
  parameter SI.Mass m_p4[8] = HBS.Tank_A.m_p[59:66];
  parameter SI.Mass m_p5[6] = HBS.Tank_A.m_p[67:72];
  parameter SI.Mass m_p6[28] = HBS.Tank_A.m_p[73:100];
  
  parameter SI.Energy E_max = sum(m_p1)*Delta_h1 + sum(m_p2)*Delta_h2 + sum(m_p3)*Delta_h3 + sum(m_p4)*Delta_h4 + sum(m_p5)*Delta_h5 + sum(m_p6)*Delta_h6;
  
  parameter SI.MassFlowRate m_flow_charging_des = 65.0 "Fixed inlet mass flow rate (kg/s)";
  
  //Initial temperature profile
  parameter SI.Length z_f_start_data[27] = {1.38, 3.20, 4.55, 5.81, 7.17, 8.56, 9.90, 11.03, 12.48, 13.94, 15.59, 17.36, 19.07, 20.16, 20.96, 22.42, 24.13, 25.72, 27.24, 28.77, 30.36, 31.45, 32.13, 33.53, 35.06, 36.28, 37.11};
  parameter SI.Length z_p_start_data[29] = {1.41, 3.20, 4.55, 5.89, 7.23, 8.57, 9.91, 11.25, 12.60, 13.94, 15.28, 16.62, 17.97, 19.31, 20.59, 21.87, 23.21, 24.50, 25.72, 26.94, 28.16, 29.38, 30.60, 31.82, 33.04, 34.26, 35.42, 36.53, 37.06};
  parameter SI.Temperature T_f_start_data[27] = {509.05, 526.04, 543.44, 567.53, 595.70, 626.36, 661.10, 688.91, 724.56, 762.99, 808.56, 855.62, 901.21, 934.41, 959.62, 1005.08, 1058.40, 1114.38, 1167.36, 1220.43, 1273.46, 1307.91, 1333.20, 1381.43, 1434.84, 1486.60, 1536.83};
  parameter SI.Temperature T_p_start_data[29] = {505.53, 516.34, 533.76, 556.09, 582.43, 612.82, 644.83, 677.92, 710.18, 744.53, 778.81, 814.02, 851.06, 887.40, 923.84, 963.37, 1002.57, 1043.46, 1084.85, 1125.67, 1168.40, 1209.98, 1251.57, 1293.07, 1332.45, 1373.27, 1412.84, 1454.80, 1481.26};
  
  
  parameter SI.Length z_start[100] = HBS.Tank_A.z_f;
  parameter SI.SpecificEnthalpy h_p_start[100] = Initialise_h_p_start(z_p_start_data,T_p_start_data,z_start);
  parameter SI.SpecificEnthalpy h_f_start[100] = Initialise_h_f_start(z_f_start_data,T_f_start_data,z_start);
  parameter Integer Correlation = 1;
  
  SI.MassFlowRate m_flow_charging_signal;
  SI.Temperature T_charging_signal;
  
  
  Modelica.Fluid.Sources.Boundary_pT Charging_Inlet(redeclare package Medium = Medium, nPorts = 1, p = 107209, use_T_in = true) annotation(
    Placement(visible = true, transformation(origin = {-76, 56}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Charging_Pump(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-36, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow_charging(y = m_flow_charging_signal) annotation(
    Placement(visible = true, transformation(origin = {-75, 90}, extent = {{-23, -10}, {23, 10}}, rotation = 0)));
  SolarTherm.Validation.Gan_HBS.HBS_6Layer_Tank HBS(redeclare package Medium = Medium, redeclare replaceable package Fluid_Package = Fluid_Package, N_f = 100, T_max = T_max, T_min = T_min, Correlation = Correlation, ar = ar, d_p= d_p, s_p = s_p,Tank_A.H_tank=H_tank,Tank_A.D_tank=D_tank,Tank_A.h_p_start=h_p_start,Tank_A.h_f_start=h_f_start,U_loss_tank=U_loss_tank) annotation(
    Placement(visible = true, transformation(origin = {0, -6}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink2 Fluid_Sink(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = { -59, -71}, extent = {{21, -21}, {-21, 21}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression T_charging(y = T_charging_signal) annotation(
    Placement(visible = true, transformation(origin = {-133, 56}, extent = {{-23, -10}, {23, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Tamb(y = 298.15) annotation(
    Placement(visible = true, transformation(origin = {-62, 4}, extent = {{-12, -18}, {12, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression p_amb(y = 107209) annotation(
    Placement(visible = true, transformation(origin = {-61, -24}, extent = {{-13, -16}, {13, 16}}, rotation = 0)));
equation
  m_flow_charging_signal = m_flow_charging_des;
  T_charging_signal = Modelica.Math.Vectors.interpolate(t_inlet_data,T_inlet_data,time);



  connect(m_flow_charging.y, Charging_Pump.m_flow) annotation(
    Line(points = {{-50, 90}, {-36, 90}, {-36, 64}}, color = {0, 0, 127}));
  connect(Charging_Inlet.ports[1], Charging_Pump.fluid_a) annotation(
    Line(points = {{-64, 56}, {-46, 56}, {-46, 56}, {-46, 56}}, color = {0, 127, 255}));
  connect(Charging_Pump.fluid_b, HBS.fluid_a) annotation(
    Line(points = {{-26, 56}, {0, 56}, {0, 30}, {0, 30}}));
  connect(HBS.fluid_b, Fluid_Sink.port_a) annotation(
    Line(points = {{0, -42}, {0, -42}, {0, -72}, {-38, -72}, {-38, -70}}, color = {0, 127, 255}));
  connect(T_charging.y, Charging_Inlet.T_in) annotation(
    Line(points = {{-108, 56}, {-100, 56}, {-100, 60}, {-90, 60}, {-90, 60}}, color = {0, 0, 127}));
  connect(Tamb.y, HBS.T_amb) annotation(
    Line(points = {{-48, 4}, {-22, 4}, {-22, 2}, {-20, 2}}, color = {0, 0, 127}));
  connect(p_amb.y, HBS.p_amb) annotation(
    Line(points = {{-46, -24}, {-38, -24}, {-38, -16}, {-20, -16}, {-20, -14}}, color = {0, 0, 127}));

annotation(experiment(StopTime = 6600, StartTime = 600, Tolerance = 1e-5, Interval = 10),
    Diagram(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1)));
end HBS_Validation_Charging_Heatloss;