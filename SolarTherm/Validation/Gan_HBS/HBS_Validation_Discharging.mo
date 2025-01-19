within SolarTherm.Validation.Gan_HBS;

model HBS_Validation_Discharging
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  
  package Medium = SolarTherm.Validation.Gan_HBS.Media.Air_533kPa;
  package Fluid_Package = SolarTherm.Validation.Gan_HBS.Materials.Air_533kPa;
  
  //Validation Parameters
  parameter SI.Length H_tank = 37.34 "Tank height (m)";
  parameter SI.Length D_tank = 7.748 "Tank diameter (m)";
  parameter Real ar = H_tank/D_tank;
  
  parameter SI.Length d_p = 16.0e-3 "Filler pore diameter (m)";
  parameter Real eta = 0.4034;
  //parameter SI.Length s_p = 31.7e-3 "Filler pore separation (m)";
  
  parameter SI.Temperature T_max = 1326.85 + 273.15 "Maximum design system temperature (K)";
  parameter SI.Temperature T_min = 224.85 + 273.15 "Minimum design system temperature (K)";
  
  parameter SI.Time t_inlet_data[2] = {0.0,13200.0} "Inlet temperature signal, time axis (s)";
  parameter SI.Temperature T_inlet_data[2] = {498.0,498.0} "Inlet temperature signal, temperature axis (K)";
  
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
  
  parameter SI.MassFlowRate m_flow_discharging_des = 55.0 "Fixed inlet mass flow rate (kg/s)";
  
  //Initial temperature profile
  parameter SI.Length z_f_start_data[28] = {1.3851956751, 2.4441220088, 3.7493709459, 5.6575568977, 6.8733081674, 8.2755396306, 9.7820648213, 10.7949415777, 12.3969157084, 13.9374085464, 15.3135653814, 16.3670003068, 17.3671087609, 19.0643203353, 20.034203208, 21.0040634445, 22.1776798841, 23.2381678886, 24.1964333816, 25.7320755418, 27.6718579691, 29.1065339853, 30.035975209, 31.7938610628, 33.2283947074, 34.4001558291, 36.1576555515, 37.1699422914};
  parameter SI.Length z_p_start_data[29] = {1.3172325925, 2.8220471581, 4.1552093909, 5.5228122862, 6.8216016343, 8.1548699414, 9.4276119759, 10.700398331, 12.0337629337, 13.3188267009, 14.5067902795, 15.7918090009, 17.1253428376, 18.580103123, 19.9136442888, 21.2471679983, 22.5201295745, 23.7930893644, 25.1266432812, 26.278426519, 27.3695802886, 28.5819526223, 29.7943090246, 31.0673043662, 32.4008010924, 33.7342225567, 35.0674539953, 36.4004977659, 37.1881033991};
  parameter SI.Temperature T_f_start_data[28] = {528.7889427532, 577.236312062, 608.2654238346, 649.0871622958, 675.430505304, 709.3969714766, 749.1924883593, 776.713553414, 820.9545401097, 865.2203938094, 905.7043761469, 936.411172441, 965.872081714, 1017.8807141091, 1049.2608820006, 1079.8526328088, 1115.1866057658, 1149.0033305593, 1177.5315952943, 1227.0175806864, 1290.358946754, 1338.7816632763, 1367.8142090274, 1423.7554785486, 1467.2193926859, 1498.8439127612, 1541.3362222881, 1558.9956049221};
  parameter SI.Temperature T_p_start_data[29] = {561.6347783633, 605.9807747773, 634.2659473316, 662.1193947003, 693.1969919267, 725.1767340624, 758.3789160979, 793.1247857627, 828.4585051164, 868.9561356056, 902.4182242472, 941.3469087291, 982.5750447054, 1027.4153546793, 1068.8987648235, 1109.7741704, 1150.6229900301, 1191.4095944027, 1233.3371221604, 1272.2494512893, 1308.7899236655, 1348.6679479429, 1387.9910793827, 1430.0159511716, 1469.9515270454, 1507.265732844, 1537.9613458072, 1562.1204775615, 1572.8419496123};
  
  
  parameter SI.Length z_start[100] = HBS.Tank_A.z_f;
  parameter SI.SpecificEnthalpy h_p_start[100] = Initialise_h_p_start(z_p_start_data,T_p_start_data,z_start);
  parameter SI.SpecificEnthalpy h_f_start[100] = Initialise_h_f_start_discharging(z_f_start_data,T_f_start_data,z_start);
  parameter Integer Correlation = 1;
  
  SI.MassFlowRate m_flow_discharging_signal;
  SI.Temperature T_discharging_signal;
  
  
  Modelica.Fluid.Sources.Boundary_pT Discharging_Inlet(redeclare package Medium = Medium, nPorts = 1, p = 107209, use_T_in = true) annotation(
    Placement(visible = true, transformation(origin = {-70, -66}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Discharging_Pump(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-24, -66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow_discharging(y = m_flow_discharging_signal) annotation(
    Placement(visible = true, transformation(origin = {-121, -42}, extent = {{-23, -10}, {23, 10}}, rotation = 0)));
  SolarTherm.Validation.Gan_HBS.HBS_6Layer_Tank HBS(redeclare package Medium = Medium, redeclare replaceable package Fluid_Package = Fluid_Package, N_f = 100, T_max = T_max, T_min = T_min, Correlation = Correlation, ar = ar, d_p= d_p, eta = eta,Tank_A.H_tank=H_tank,Tank_A.D_tank=D_tank,Tank_A.h_p_start=h_p_start,Tank_A.h_f_start=h_f_start) annotation(
    Placement(visible = true, transformation(origin = {0, -6}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink2 Fluid_Sink(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-59, 65}, extent = {{21, -21}, {-21, 21}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression T_discharging(y = T_discharging_signal) annotation(
    Placement(visible = true, transformation(origin = {-133, -62}, extent = {{-23, -10}, {23, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Tamb(y = 298.15) annotation(
    Placement(visible = true, transformation(origin = {-62, 4}, extent = {{-12, -18}, {12, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression p_amb(y = 532818) annotation(
    Placement(visible = true, transformation(origin = {-61, -24}, extent = {{-13, -16}, {13, 16}}, rotation = 0)));
equation
  m_flow_discharging_signal = m_flow_discharging_des;
  T_discharging_signal = Modelica.Math.Vectors.interpolate(t_inlet_data,T_inlet_data,time);
  connect(Discharging_Inlet.ports[1], Discharging_Pump.fluid_a) annotation(
    Line(points = {{-58, -66}, {-34, -66}}, color = {0, 127, 255}));
  connect(Tamb.y, HBS.T_amb) annotation(
    Line(points = {{-48, 4}, {-22, 4}, {-22, 2}, {-20, 2}}, color = {0, 0, 127}));
  connect(p_amb.y, HBS.p_amb) annotation(
    Line(points = {{-46, -24}, {-38, -24}, {-38, -16}, {-20, -16}, {-20, -14}}, color = {0, 0, 127}));
  connect(Discharging_Inlet.T_in, T_discharging.y) annotation(
    Line(points = {{-84, -62}, {-106, -62}, {-106, -62}, {-108, -62}}, color = {0, 0, 127}));
  connect(m_flow_discharging.y, Discharging_Pump.m_flow) annotation(
    Line(points = {{-96, -42}, {-24, -42}, {-24, -58}, {-24, -58}}, color = {0, 0, 127}));
  connect(Discharging_Pump.fluid_b, HBS.fluid_b) annotation(
    Line(points = {{-14, -66}, {0, -66}, {0, -42}, {0, -42}}, color = {0, 127, 255}));
  connect(HBS.fluid_a, Fluid_Sink.port_a) annotation(
    Line(points = {{0, 30}, {0, 30}, {0, 66}, {-38, 66}, {-38, 66}}, color = {0, 127, 255}));

annotation(experiment(StopTime = 2820, StartTime = 0, Tolerance = 1e-5, Interval = 60),
    Diagram(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1)));
end HBS_Validation_Discharging;