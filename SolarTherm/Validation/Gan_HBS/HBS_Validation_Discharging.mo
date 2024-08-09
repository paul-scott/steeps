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
  
  parameter SI.Length d_p = 19.6e-3 "Filler pore diameter (m)";
  parameter SI.Length s_p = 31.7e-3 "Filler pore separation (m)";
  
  parameter SI.Temperature T_max = 1600.0 + 273.15 "Maximum design system temperature (K)";
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
  parameter SI.Length z_f_start_data[100] = {0.1867, 0.5601, 0.9335, 1.3069, 1.6803, 2.0537, 2.4271, 2.8005, 3.1739, 3.5473, 3.9207, 4.2941, 4.6675, 5.0409, 5.4143, 5.7877, 6.1611, 6.5345, 6.9079, 7.2813, 7.6547, 8.0281, 8.4015, 8.7749, 9.1483, 9.5217, 9.8951, 10.2685, 10.6419, 11.0153, 11.3887, 11.7621, 12.1355, 12.5089, 12.8823, 13.2557, 13.6291, 14.0025, 14.3759, 14.7493, 15.1227, 15.4961, 15.8695, 16.2429, 16.6163, 16.9897, 17.3631, 17.7365, 18.1099, 18.4833, 18.8567, 19.2301, 19.6035, 19.9769, 20.3503, 20.7237, 21.0971, 21.4705, 21.8439, 22.2173, 22.5907, 22.9641, 23.3375, 23.7109, 24.0843, 24.4577, 24.8311, 25.2045, 25.5779, 25.9513, 26.3247, 26.6981, 27.0715, 27.4449, 27.8183, 28.1917, 28.5651, 28.9385, 29.3119, 29.6853, 30.0587, 30.4321, 30.8055, 31.1789, 31.5523, 31.9257, 32.2991, 32.6725, 33.0459, 33.4193, 33.7927, 34.1661, 34.5395, 34.9129, 35.2863, 35.6597, 36.0331, 36.4065, 36.7799, 37.1533};
  parameter SI.Length z_p_start_data[100] = {0.1867, 0.5601, 0.9335, 1.3069, 1.6803, 2.0537, 2.4271, 2.8005, 3.1739, 3.5473, 3.9207, 4.2941, 4.6675, 5.0409, 5.4143, 5.7877, 6.1611, 6.5345, 6.9079, 7.2813, 7.6547, 8.0281, 8.4015, 8.7749, 9.1483, 9.5217, 9.8951, 10.2685, 10.6419, 11.0153, 11.3887, 11.7621, 12.1355, 12.5089, 12.8823, 13.2557, 13.6291, 14.0025, 14.3759, 14.7493, 15.1227, 15.4961, 15.8695, 16.2429, 16.6163, 16.9897, 17.3631, 17.7365, 18.1099, 18.4833, 18.8567, 19.2301, 19.6035, 19.9769, 20.3503, 20.7237, 21.0971, 21.4705, 21.8439, 22.2173, 22.5907, 22.9641, 23.3375, 23.7109, 24.0843, 24.4577, 24.8311, 25.2045, 25.5779, 25.9513, 26.3247, 26.6981, 27.0715, 27.4449, 27.8183, 28.1917, 28.5651, 28.9385, 29.3119, 29.6853, 30.0587, 30.4321, 30.8055, 31.1789, 31.5523, 31.9257, 32.2991, 32.6725, 33.0459, 33.4193, 33.7927, 34.1661, 34.5395, 34.9129, 35.2863, 35.6597, 36.0331, 36.4065, 36.7799, 37.1533};
  parameter SI.Temperature T_f_start_data[100] = {581.41, 587.999, 594.68, 601.398, 607.996, 614.581, 621.064, 627.819, 634.843, 642.137, 649.706, 657.554, 665.685, 674.095, 682.781, 691.734, 700.911, 710.062, 719.319, 728.664, 738.084, 747.566, 757.105, 766.698, 776.35, 786.065, 795.85, 805.56, 815.172, 824.758, 834.423, 844.228, 854.179, 864.279, 874.531, 884.941, 895.515, 906.118, 916.726, 927.403, 938.154, 948.984, 959.898, 971.111, 982.426, 993.852, 1005.29, 1016.68, 1028.09, 1039.55, 1051.07, 1062.64, 1074.28, 1085.99, 1097.76, 1109.46, 1121.11, 1132.76, 1144.4, 1156.34, 1168.28, 1180.23, 1192.18, 1204.07, 1215.84, 1227.54, 1239.31, 1252.13, 1265.07, 1278.14, 1291.33, 1304.59, 1317.85, 1330.93, 1344.07, 1357.29, 1370.58, 1383.96, 1397.42, 1410.84, 1424.23, 1437.67, 1451.21, 1464.79, 1478.34, 1491.78, 1504.95, 1517.66, 1529.82, 1541.26, 1551.83, 1561.4, 1569.84, 1577.05, 1582.99, 1587.67, 1591.19, 1593.68, 1595.38, 1596.54};
  parameter SI.Temperature T_p_start_data[100] = {568.961, 574.937, 581.071, 587.342, 593.699, 600.281, 606.006, 612.034, 618.373, 625.03, 632.009, 639.318, 646.963, 654.944, 663.253, 671.879, 680.807, 689.948, 699.235, 708.644, 718.152, 727.733, 737.374, 747.062, 756.796, 766.582, 776.426, 786.314, 796.17, 805.784, 815.361, 825.076, 834.933, 844.935, 855.085, 865.386, 875.843, 886.441, 897.108, 907.838, 918.633, 929.496, 940.068, 951.29, 962.605, 974.018, 985.523, 997.052, 1008.56, 1020.1, 1031.69, 1043.34, 1055.04, 1066.8, 1078.62, 1090.48, 1102.31, 1114.12, 1125.45, 1137.54, 1149.63, 1161.73, 1173.81, 1185.89, 1197.9, 1209.64, 1219.94, 1232.75, 1245.67, 1258.71, 1271.87, 1285.14, 1298.82, 1311.96, 1325.15, 1338.4, 1351.72, 1365.12, 1378.62, 1392.19, 1405.66, 1419.11, 1432.72, 1446.47, 1460.31, 1474.18, 1487.98, 1501.55, 1514.74, 1527.4, 1539.37, 1550.47, 1560.54, 1569.42, 1576.99, 1583.19, 1588.02, 1591.54, 1593.91, 1595.36};
  
  
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
  SolarTherm.Validation.Gan_HBS.HBS_6Layer_Tank HBS(redeclare package Medium = Medium, redeclare replaceable package Fluid_Package = Fluid_Package, N_f = 100, T_max = T_max, T_min = T_min, Correlation = Correlation, ar = ar, d_p= d_p, s_p = s_p,Tank_A.H_tank=H_tank,Tank_A.D_tank=D_tank,Tank_A.h_p_start=h_p_start,Tank_A.h_f_start=h_f_start) annotation(
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

annotation(experiment(StopTime = 438, StartTime = 0, Tolerance = 1e-5, Interval = 6),
    Diagram(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1)));
end HBS_Validation_Discharging;