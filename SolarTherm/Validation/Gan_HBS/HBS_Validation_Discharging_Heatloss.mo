within SolarTherm.Validation.Gan_HBS;

model HBS_Validation_Discharging_Heatloss
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  
  package Medium = SolarTherm.Validation.Gan_HBS.Media.Air_533kPa;
  package Fluid_Package = SolarTherm.Validation.Gan_HBS.Materials.Air_533kPa;
  
  //Validation Parameters
  parameter SI.Length H_tank = 37.34 "Tank height (m)";
  parameter SI.Length D_tank = 7.748 "Tank diameter (m)";
  parameter Real ar = H_tank/D_tank;
  
  parameter SI.Length d_p = 19.6e-3 "Filler pore diameter (m)";
  parameter SI.Length s_p = 31.7e-3 "Filler pore separation (m)";
  
  parameter SI.Temperature T_max = 1600.0 + 273.15 "Maximum design system temperature (K)";
  parameter SI.Temperature T_min = 224.85 + 273.15 "Minimum design system temperature (K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.576 "Heat loss coeff of surfaces (W/m2K)";
  
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
  parameter SI.Length z_f_start_data[100] = {0.1867, 0.5601, 0.9335, 1.3069, 1.6803, 2.0537, 2.4271, 2.8005, 3.1739, 3.5473, 3.9207, 4.2941, 4.6675, 5.0409, 5.4143, 5.7877, 6.1611, 6.5345, 6.9079, 7.2813, 7.6547, 8.0281, 8.4015, 8.7749, 9.1483, 9.5217, 9.8951, 10.2685, 10.6419, 11.0153, 11.3887, 11.7621, 12.1355, 12.5089, 12.8823, 13.2557, 13.6291, 14.0025, 14.3759, 14.7493, 15.1227, 15.4961, 15.8695, 16.2429, 16.6163, 16.9897, 17.3631, 17.7365, 18.1099, 18.4833, 18.8567, 19.2301, 19.6035, 19.9769, 20.3503, 20.7237, 21.0971, 21.4705, 21.8439, 22.2173, 22.5907, 22.9641, 23.3375, 23.7109, 24.0843, 24.4577, 24.8311, 25.2045, 25.5779, 25.9513, 26.3247, 26.6981, 27.0715, 27.4449, 27.8183, 28.1917, 28.5651, 28.9385, 29.3119, 29.6853, 30.0587, 30.4321, 30.8055, 31.1789, 31.5523, 31.9257, 32.2991, 32.6725, 33.0459, 33.4193, 33.7927, 34.1661, 34.5395, 34.9129, 35.2863, 35.6597, 36.0331, 36.4065, 36.7799, 37.1533};
  parameter SI.Length z_p_start_data[100] = {0.1867, 0.5601, 0.9335, 1.3069, 1.6803, 2.0537, 2.4271, 2.8005, 3.1739, 3.5473, 3.9207, 4.2941, 4.6675, 5.0409, 5.4143, 5.7877, 6.1611, 6.5345, 6.9079, 7.2813, 7.6547, 8.0281, 8.4015, 8.7749, 9.1483, 9.5217, 9.8951, 10.2685, 10.6419, 11.0153, 11.3887, 11.7621, 12.1355, 12.5089, 12.8823, 13.2557, 13.6291, 14.0025, 14.3759, 14.7493, 15.1227, 15.4961, 15.8695, 16.2429, 16.6163, 16.9897, 17.3631, 17.7365, 18.1099, 18.4833, 18.8567, 19.2301, 19.6035, 19.9769, 20.3503, 20.7237, 21.0971, 21.4705, 21.8439, 22.2173, 22.5907, 22.9641, 23.3375, 23.7109, 24.0843, 24.4577, 24.8311, 25.2045, 25.5779, 25.9513, 26.3247, 26.6981, 27.0715, 27.4449, 27.8183, 28.1917, 28.5651, 28.9385, 29.3119, 29.6853, 30.0587, 30.4321, 30.8055, 31.1789, 31.5523, 31.9257, 32.2991, 32.6725, 33.0459, 33.4193, 33.7927, 34.1661, 34.5395, 34.9129, 35.2863, 35.6597, 36.0331, 36.4065, 36.7799, 37.1533};
  parameter SI.Temperature T_f_start_data[100] = {580.92, 587.70, 594.38, 601.10, 607.70, 614.28, 620.77, 627.52, 634.53, 641.82, 649.39, 657.23, 665.35, 673.75, 682.42, 691.36, 700.54, 709.68, 718.93, 728.26, 737.68, 747.15, 756.68, 766.27, 775.91, 785.62, 795.39, 805.11, 814.72, 824.30, 833.96, 843.76, 853.70, 863.79, 874.03, 884.43, 895.00, 905.60, 916.20, 926.87, 937.61, 948.43, 959.33, 970.53, 981.84, 993.25, 1004.69, 1016.07, 1027.47, 1038.92, 1050.42, 1061.99, 1073.62, 1085.31, 1097.08, 1108.77, 1120.42, 1132.06, 1143.69, 1155.62, 1167.55, 1179.49, 1191.42, 1203.32, 1215.08, 1226.77, 1238.52, 1251.31, 1264.22, 1277.25, 1290.40, 1303.63, 1316.83, 1329.86, 1342.94, 1356.08, 1369.29, 1382.57, 1395.93, 1409.26, 1422.54, 1435.85, 1449.24, 1462.68, 1476.10, 1489.41, 1502.49, 1515.13, 1527.24, 1538.67, 1549.29, 1558.94, 1567.49, 1574.86, 1580.99, 1585.86, 1589.55, 1592.17, 1593.92, 1595.04};
  parameter SI.Temperature T_p_start_data[100] = {568.35, 574.71, 580.83, 587.10, 593.45, 600.05, 605.77, 611.79, 618.12, 624.78, 631.75, 639.06, 646.69, 654.67, 662.97, 671.58, 680.50, 689.63, 698.91, 708.31, 717.81, 727.38, 737.01, 746.69, 756.42, 766.20, 776.03, 785.92, 795.77, 805.39, 814.96, 824.67, 834.52, 844.52, 854.66, 864.95, 875.40, 885.99, 896.65, 907.37, 918.16, 929.01, 939.58, 950.79, 962.09, 973.49, 984.99, 996.51, 1008.02, 1019.55, 1031.13, 1042.76, 1054.46, 1066.21, 1078.02, 1089.87, 1101.69, 1113.49, 1124.82, 1136.90, 1148.98, 1161.07, 1173.15, 1185.22, 1197.22, 1208.96, 1219.26, 1232.04, 1244.94, 1257.95, 1271.08, 1284.32, 1297.95, 1311.05, 1324.18, 1337.37, 1350.62, 1363.95, 1377.35, 1390.83, 1404.22, 1417.53, 1431.00, 1444.61, 1458.30, 1472.02, 1485.69, 1499.14, 1512.24, 1524.84, 1536.79, 1547.92, 1558.06, 1567.07, 1574.83, 1581.24, 1586.29, 1590.03, 1592.60, 1593.55};
  
  
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
  SolarTherm.Validation.Gan_HBS.HBS_6Layer_Tank HBS(redeclare package Medium = Medium, redeclare replaceable package Fluid_Package = Fluid_Package, N_f = 100, T_max = T_max, T_min = T_min, Correlation = Correlation, ar = ar, d_p= d_p, s_p = s_p,Tank_A.H_tank=H_tank,Tank_A.D_tank=D_tank,Tank_A.h_p_start=h_p_start,Tank_A.h_f_start=h_f_start,U_loss_tank=U_loss_tank) annotation(
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

annotation(experiment(StopTime = 7038, StartTime = 0, Tolerance = 1e-5, Interval = 60),
    Diagram(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1)));
end HBS_Validation_Discharging_Heatloss;