within SolarTherm.Systems.HBSTES_Reference_Cases;

model HBSTES_Reference_1_SystemLevel
  //Total simulation time = 8448s, Capacity Factor = 98.58%, Total annual heat delivered = 774.4 GWh
  extends Modelica.Icons.Example;
  import Modelica.SIunits.Conversions.*;
  import Modelica.Constants.*;
  parameter String PV_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_pv.motab");
  parameter String Wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_wind.motab");
  parameter String schd_input = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Schedules/schedule_Qflow_Calciner.motab");
  replaceable package Medium = SolarTherm.Media.Air.Air_amb_p_curvefit;
  replaceable package Fluid = SolarTherm.Materials.Air_amb_p_curvefit;
  replaceable package Filler = SolarTherm.Materials.Mullite_20pct_porosity_averaged;
  //Parameter Inputs
  parameter Real RM = 2.0 "Renewable Multiple (pre-transmission oversizing)";
  parameter Real HM = 2.0 "Heater Multiple";
  parameter Real PV_fraction = 0.5 "PV_fraction";
  parameter Real t_storage = 55.0 "Hours of storage (hours)";
  parameter Real util_storage_des = 0.325446;
  //Utilisation determined via component-level analysis 0.5767
  parameter Real level_storage_mid = 0.441777;
  //Midpoint of minimum and maximum storage levels determine via component-level analysis 0.4637
  parameter SI.SpecificEnthalpy h_tol = 0.05 * (TES.Tank_A.h_f_max - TES.Tank_A.h_f_min);
  //Heater Parameters
  parameter Real eff_heater = 1.00 "Electrical-to-heat conversion efficiency of the heater";
  //Renewable Parameters
  parameter SI.Power P_renewable_des = RM * P_heater_des;
  parameter SI.HeatFlowRate Q_heater_des = HM * Q_process_des;
  parameter SI.Power P_heater_des = Q_heater_des / eff_heater;
  parameter SI.Power PV_ref_size = 50.0e6;
  parameter SI.Power Wind_ref_size = 50.0e6;
  parameter Real LOF_PV = 1.300;
  parameter Real LOF_Wind = 1.131;
  parameter SI.Power P_wind_net = (1.0 - PV_fraction)*P_renewable_des;
  parameter SI.Power P_PV_net = PV_fraction*P_renewable_des;
  parameter SI.Power P_wind_gross = P_wind_net*LOF_Wind;
  parameter SI.Power P_PV_gross = P_PV_net*LOF_PV;
  //Results
  SI.Energy E_supplied(start = 0) "Energy supplied by the boiler to the industrial process (J)";
  SI.Energy E_demand(start = 0) "Energy demanded by the industrial process (J)";
  SI.Energy E_renewable(start = 0) "Electrical energy produced by the renewable source (J)";
  SI.Energy E_wind(start = 0) "Electrical energy produced by the wind renewable source (J)";
  SI.Energy E_pv(start = 0) "Electrical energy produced by the pv renewable source (J)";
  //SI.Energy E_heater_in(start = 0) "Electrical energy transmitted to the heater (J)";
  SI.Energy E_heater_raw(start = 0) "Heat energy produced by the heater (after efficiency losses)";
  SI.Energy E_heater_out(start = 0) "Heat energy transferred to the air (after curtailment)";
  Real Capacity_Factor(start = 0) "Capacity factor of the system";
  //Discretisation and geometry
  parameter Integer N_f = 50;
  //50
  parameter SI.Length d_p = 0.03 "Hole diameter in the filler (m)";
  parameter Real ar = 4.8 "Tank H/D ratio (m)";
  parameter Real eta = 0.51 "Packed-bed porosity";
  //parameter SI.Length L_pipe = 50.0;
  //parameter SI.Length D_pipe = 0.050;
  //parameter SI.Length D_solid = 0.075;
  //Misc Parameters
  parameter SI.CoefficientOfHeatTransfer U_loss_top = 0.588 "Heat loss coefficient at the top of the tank (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_bot = 2.0 "Heat loss coefficient at the bottom of the tank (W/m2K)";
  parameter Integer Correlation = 1;
  //1=Liq 2=Air
  
  //Costs
  parameter SI.Mass m_filler_pertank = TES.Tank_A.rho_p * (1.0 - TES.Tank_A.eta) * (3.142 * TES.Tank_A.D_tank * TES.Tank_A.D_tank * TES.Tank_A.H_tank/4.0);
  
  /*
    parameter SI.Temperature T_max = 1100.0 + 273.15 "Maximum temperature (K)";
    parameter SI.Temperature T_process_des = 1000.0 + 273.15 "Design process inlet temperature (K)";
    parameter SI.Temperature T_high_set = 1000.0 + 273.15 "TES hot blend temperature temperature (K)";
    parameter SI.Temperature T_process_min = 1000.0 + 273.15 "Minimum tolerated outlet temperature to process (K)";
    parameter SI.Temperature T_heater_max = 740.0 + 273.15 "Maximum tolerated outlet temperature to heater (K)";
    parameter SI.Temperature T_low_set = 740.0 + 273.15 "TES cold blend temperature (K)";
    parameter SI.Temperature T_heater_des = 640.0 + 273.15 "Design receiver inlet temperature (K)";
    
    parameter SI.Temperature T_max = 500.0 + 273.15 "Maximum system temperature (K)";
    parameter SI.Temperature T_boiler_start = 400.0 + 273.15 "Temperature above-which TES can start discharge (K)"; //50K buffer
    parameter SI.Temperature T_boiler_min = 350.0 + 273.15 "Temperature below-which TES stops discharge (K)";
    parameter SI.Temperature T_heater_max = 300.0 + 273.15 "Temperature above-which TES stops charging (K)";
    parameter SI.Temperature T_heater_start = 250.0 + 273.15 "Temperature below-which TES can start charging (K)"; //50k buffer
    parameter SI.Temperature T_min = 150.0 + 273.15 "Minimum system temperature (K)";
      //Legacy terms from the CSP field. PB (power block) refers to the boiler, Recv (receiver) refers to the electric heater.
    */
  //Temperature Controls
  parameter SI.Temperature T_max = 1300.0 + 273.15 "Maximum system temperature (K)";
  parameter SI.Temperature T_process_start = 1100.0 + 273.15 "Temperature above-which TES can start discharge (K)";
  //50K buffer
  parameter SI.Temperature T_process_min = 1000.0 + 273.15 "Temperature below-which TES stops discharge (K)";
  parameter SI.Temperature T_process_des = 1000.0 + 273.15 "Design process inlet temperature (K)";
  parameter SI.Temperature T_heater_max = 400.0 + 273.15 "Temperature above-which TES stops charging (K)";
  parameter SI.Temperature T_heater_des = 400.0 + 273.15 "Design receiver inlet temperature (K)";
  parameter SI.Temperature T_heater_start = 300.0 + 273.15 "Temperature below-which TES can start charging (K)";
  //50k buffer
  parameter SI.Temperature T_min = 200.0 + 273.15 "Minimum system temperature (K)";
  //Legacy terms from the CSP field. PB (power block) refers to the boiler, Recv (receiver) refers to the electric heater.
  //Level-Controls
  parameter SI.Time t_stor_start_dis = 1.0 * 3600.0 "Number of effective storage seconds stored before TES can start discharging (1 hour)";
  //Calculated Parameters
  parameter Modelica.SIunits.Energy E_max = t_storage * 3600.0 * Q_process_des "Maximum tank stored energy (J)";
  parameter Modelica.SIunits.HeatFlowRate Q_process_des = 15.03e6 "Heat-rate to process at design (W)";
  parameter Modelica.SIunits.MassFlowRate m_process_des = Q_process_des / (h_air_process_des - h_air_min_des) "Design process input mass flow rate (kg/s)";
  
parameter Medium.ThermodynamicState state_air_min_des = Medium.setState_pTX(Medium.p_default, T_min) "Thermodynamic state of air at the minimum system temperature T_min.";
  parameter Medium.ThermodynamicState state_air_max_des = Medium.setState_pTX(Medium.p_default, T_max) "Thermodynamic state of air at the maximum system temperature T_min.";
  parameter Medium.ThermodynamicState state_air_process_des = Medium.setState_pTX(Medium.p_default, T_process_des) "Thermodynamic state of air at the design process inlet T_process_des.";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_min_des = Medium.specificEnthalpy(state_air_min_des) "Specific enthalpy of air at minimum system temperature T_min (J/kg)";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_max_des = Medium.specificEnthalpy(state_air_max_des) "Specific enthalpy of air at maximum system temperature T_max (J/kg)";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_process_des = Medium.specificEnthalpy(state_air_process_des) "Specific enthalpy of air at design process inlet temperature T_process_des (J/kg)";
  SolarTherm.Models.Storage.Thermocline.Thermocline_HBS_LC_SingleTank_Final TES(redeclare package Medium = Medium, redeclare package Fluid_Package = Fluid, redeclare package Filler_Package = Filler, N_f = N_f, T_max = T_max, T_min = T_min, Correlation = Correlation, E_max = E_max, ar = ar, d_p = d_p, eta = eta, U_loss_top = U_loss_top, U_loss_bot = U_loss_bot) annotation(
    Placement(visible = true, transformation(origin = {32, 0}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpCold(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-18, -78}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction_LoopBreaker Splitter_Top(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {32, 67.6249}, extent = {{-18, -13.9366}, {18, 13.9366}}, rotation = 0)));
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction Splitter_Bot(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {31, -43.8998}, extent = {{17, 0}, {-17, -22.039}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpHot(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {91, 79}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Tamb(y = 298.15) annotation(
    Placement(visible = true, transformation(origin = {-5, 0}, extent = {{-9, -12}, {9, 12}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression p_amb(y = 101325) annotation(
    Placement(visible = true, transformation(origin = {107, 2.22045e-16}, extent = {{11, -12}, {-11, 12}}, rotation = 0)));
  SolarTherm.Models.Control.WindPV_Thermocline_Control Control(redeclare package HTF = Medium, E_max = E_max, Q_boiler_des = Q_process_des, T_boiler_min = T_process_min, T_boiler_start = T_process_start, T_heater_max = T_heater_max, T_heater_start = T_heater_start, T_target = T_max, util_storage_des = util_storage_des, h_target = h_air_max_des, level_mid = level_storage_mid, m_0 = 1e-8, m_boiler_des = m_process_des*(h_air_process_des-h_air_min_des)/(h_air_max_des-h_air_min_des), m_min = 1e-8, m_tol = 0.01 * m_process_des, t_stor_start_dis = t_stor_start_dis, t_wait = 1.0 * 3600.0) annotation(
    Placement(visible = true, transformation(origin = {114, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0))); //m_boiler_des needs to be scaled down to the mass flow rate at the overdesigned temperature T_max.
  SolarTherm.Models.Fluid.HeatExchangers.Boiler_Basic Process(redeclare package Medium = Medium, T_cold_set = T_min, T_hot_set = T_process_des) annotation(
    Placement(visible = true, transformation(origin = {158, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.CSP.CRS.Receivers.Basic_Heater basic_Heater(redeclare package Medium = Medium, P_heater_des = P_heater_des, Q_heater_des = Q_heater_des, eff_heater = eff_heater, T_cold_set = T_min, T_hot_set = T_max) annotation(
    Placement(visible = true, transformation(origin = {-46, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable PV_input(fileName = PV_file, tableName = "Power", tableOnFile = true, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(visible = true, transformation(origin = {-124, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add Grid_Sum(k1 = P_PV_gross / PV_ref_size, k2 = P_wind_gross / Wind_ref_size) annotation(
    Placement(visible = true, transformation(origin = {-84, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Wind_input(fileName = Wind_file, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "Power", tableOnFile = true) annotation(
    Placement(visible = true, transformation(origin = {-124, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression ConstantDemand(y = Q_process_des) annotation(
    Placement(visible = true, transformation(origin = {182, 42}, extent = {{16, -10}, {-16, 10}}, rotation = 0)));
equation
  der(E_pv) = Grid_Sum.k1*Grid_Sum.u1;
  der(E_wind) = Grid_Sum.k2*Grid_Sum.u2;
  der(E_renewable) = Grid_Sum.y;
  //der(E_heater_in) = basic_Heater.P_heater_out;
  der(E_heater_raw) = basic_Heater.Q_heater_raw;
  der(E_heater_out) = basic_Heater.Q_out;
  der(E_supplied) = Process.Q_flow;
  
  der(E_demand) = Control.Q_demand;
  if time > 86400.0 then
    Capacity_Factor = E_supplied / E_demand;
  else
    Capacity_Factor = 0.0;
  end if;
//grid_input.Q_defocus_y = min(gridInput.grid_input.y[1], scheduler.y[1] * (h_salt_hot_set - h_salt_cold_set));
  connect(TES.fluid_b, Splitter_Bot.fluid_c) annotation(
    Line(points = {{32, -30}, {32, -56.5}, {31, -56.5}, {31, -60}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Tamb.y, TES.T_amb) annotation(
    Line(points = {{5, 0}, {15, 0}}, color = {0, 0, 127}));
  connect(p_amb.y, TES.p_amb) annotation(
    Line(points = {{95, 0}, {49, 0}}, color = {0, 0, 127}));
  connect(Splitter_Bot.fluid_b, pumpCold.fluid_a) annotation(
    Line(points = {{17, -78}, {-8, -78}}, color = {0, 127, 255}, thickness = 0.5));
  connect(TES.T_top_measured, Control.T_top_tank) annotation(
    Line(points = {{49, 21}, {78, 21}, {78, 24}, {103, 24}}, color = {0, 0, 127}));
  connect(TES.T_bot_measured, Control.T_bot_tank) annotation(
    Line(points = {{49, -21}, {82, -21}, {82, 20}, {103, 20}}, color = {0, 0, 127}));
  connect(TES.Level, Control.Level) annotation(
    Line(points = {{49, 8}, {90, 8}, {90, 28}, {103, 28}, {103, 29}}, color = {0, 0, 127}));
  connect(Splitter_Top.fluid_c, TES.fluid_a) annotation(
    Line(points = {{32, 67}, {32, 30}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Splitter_Top.fluid_b, pumpHot.fluid_a) annotation(
    Line(points = {{47, 79}, {82, 79}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Process.fluid_b, Splitter_Bot.fluid_a) annotation(
    Line(points = {{158, -10}, {158, -78}, {45, -78}}, color = {0, 127, 255}, thickness = 0.5));
  connect(pumpHot.fluid_b, Process.fluid_a) annotation(
    Line(points = {{100, 79}, {158, 79}, {158, 10}}, color = {0, 127, 255}, thickness = 0.5));
  connect(pumpCold.fluid_b, basic_Heater.fluid_a) annotation(
    Line(points = {{-28, -78}, {-46, -78}, {-46, 1}}, color = {0, 127, 255}, thickness = 0.5));
  connect(basic_Heater.fluid_b, Splitter_Top.fluid_a) annotation(
    Line(points = {{-46, 19}, {-46, 79}, {18, 79}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Control.curtail, basic_Heater.curtail) annotation(
    Line(points = {{125, 22}, {128, 22}, {128, -36}, {-68, -36}, {-68, 1}, {-57, 1}}, color = {255, 0, 255}));
  connect(Control.Q_curtail, basic_Heater.Q_curtail) annotation(
    Line(points = {{103, 17}, {-10, 17}, {-10, -12}, {-72, -12}, {-72, 2}, {-57, 2}, {-57, 5}}, color = {0, 0, 127}));
  connect(TES.h_top_outlet, Control.h_tank_top) annotation(
    Line(points = {{22, 25}, {22, 46}, {105, 46}, {105, 37}}, color = {0, 0, 127}));
  connect(PV_input.y[1], Grid_Sum.u1) annotation(
    Line(points = {{-113, 34}, {-104, 34}, {-104, 24}, {-96, 24}}, color = {0, 0, 127}));
  connect(Wind_input.y[1], Grid_Sum.u2) annotation(
    Line(points = {{-113, 4}, {-104, 4}, {-104, 12}, {-96, 12}}, color = {0, 0, 127}));
  connect(Grid_Sum.y, basic_Heater.P_supply) annotation(
    Line(points = {{-73, 18}, {-66, 18}, {-66, 10}, {-57, 10}}, color = {0, 0, 127}));
  connect(Process.h_out_signal, Control.h_boiler_outlet) annotation(
    Line(points = {{148, -8}, {136, -8}, {136, 52}, {116, 52}, {116, 38}, {116, 38}}, color = {0, 0, 127}));
  connect(TES.h_bot_outlet, Control.h_tank_bot) annotation(
    Line(points = {{22, -25}, {22, -30}, {132, -30}, {132, 46}, {110, 46}, {110, 38}}, color = {0, 0, 127}));
  connect(basic_Heater.Q_heater_raw, Control.Q_heater_raw) annotation(
    Line(points = {{-34, 10}, {-20, 10}, {-20, 52}, {86, 52}, {86, 32}, {104, 32}, {104, 32}}, color = {0, 0, 127}));
  connect(Control.m_heater_signal, pumpCold.m_flow) annotation(
    Line(points = {{126, 32}, {140, 32}, {140, -58}, {-18, -58}, {-18, -70}, {-18, -70}}, color = {0, 0, 127}));
  connect(Control.m_boiler_signal, pumpHot.m_flow) annotation(
    Line(points = {{126, 26}, {130, 26}, {130, 94}, {90, 94}, {90, 86}, {92, 86}}, color = {0, 0, 127}));
  connect(ConstantDemand.y, Control.Q_demand) annotation(
    Line(points = {{164, 42}, {120, 42}, {120, 38}, {122, 38}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}, initialScale = 0.1), graphics = {Text(origin = {85, 68}, extent = {{-11, 4}, {23, -10}}, textString = "Hot Pump"), Text(origin = {-21, -90}, extent = {{-11, 4}, {23, -10}}, textString = "Cold Pump"), Text(origin = {-51, -4}, extent = {{-11, 4}, {13, -6}}, textString = "Heater")}),
    Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}}, preserveAspectRatio = false)),
    experiment(StopTime = 3.1536e+07, StartTime = 0, Tolerance = 1.0e-5, Interval = 300, maxStepSize = 60, initialStepSize = 60));
end HBSTES_Reference_1_SystemLevel;