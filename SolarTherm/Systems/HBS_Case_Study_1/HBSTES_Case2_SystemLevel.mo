within SolarTherm.Systems.HBS_Case_Study_1;

model HBSTES_Case2_SystemLevel
  extends Modelica.Icons.Example;
  import Modelica.SIunits.Conversions.*;
  import Modelica.Constants.*;
  parameter String PV_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_pv.motab");
  parameter String Wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_wind.motab");
  parameter String schd_input = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Schedules/schedule_Qflow_Calciner.motab");
  replaceable package Medium = SolarTherm.Media.Air.Air_amb_p_curvefit;
  replaceable package Fluid = SolarTherm.Materials.Air_amb_p_curvefit;
  replaceable package Filler = SolarTherm.Materials.Mullite_20pct_porosity;
  //Parameter Inputs
  parameter Real RM = 2.0 "Renewable Multiple (pre-transmission oversizing)";
  parameter Real HM = 2.0 "Heater Multiple";
  parameter Real PV_fraction = 0.5 "PV_fraction";
  
  parameter Real t_storage = 20.0 "Hours of storage (hours)";
  parameter Real util_storage_des = 0.2779;
  //UtilisaSystemstion determined via component-level analysis 0.5767
  parameter Real level_storage_mid = 0.6105;
  //Midpoint of minimum and maximum storage levels determine via component-level analysis 0.4637
  //parameter SI.SpecificEnthalpy h_tol = 0.05 * (TES.Tank_A.h_f_max - TES.Tank_A.h_f_min);
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
  parameter Integer N_f = 100;
  //50
  parameter SI.Length d_p = 0.03 "Hole diameter in the filler (m)";
  parameter Real ar = 4.8/(sqrt(0.5)) "Tank H/D ratio (m)";
  parameter Real eta = 0.53 "Packed-bed porosity";
  //parameter SI.Length L_pipe = 50.0;
  //parameter SI.Length D_pipe = 0.050;
  //parameter SI.Length D_solid = 0.075;
  //Misc Parameters
  parameter SI.CoefficientOfHeatTransfer U_loss_top = 10.0*(323.15-298.15)/(T_max-298.15) "Heat loss coefficient at the top of the tank (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_bot = 10.0*(323.15-298.15)/(T_max-298.15) "Heat loss coefficient at the bottom of the tank (W/m2K)";
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
  parameter SI.Temperature T_max = 1100.0 + 273.15 "Maximum system temperature (K)";
  parameter SI.Temperature T_process_start = 1000.0 + 273.15 "Temperature above-which TES can start discharge (K)";
  //50K buffer
  parameter SI.Temperature T_process_min = 1000.0 + 273.15 "Temperature below-which TES stops discharge (K)";
  parameter SI.Temperature T_process_des = 1000.0 + 273.15 "Design process inlet temperature (K)";
  parameter SI.Temperature T_heater_max = 450.0 + 273.15 "Temperature above-which TES stops charging (K)";
  parameter SI.Temperature T_heater_des = 450.0 + 273.15 "Design receiver inlet temperature (K)";
  parameter SI.Temperature T_heater_start = 317.3 + 273.15 "Temperature below-which TES can start charging (K)";
  //50k buffer
  parameter SI.Temperature T_min = 317.3 + 273.15 "Minimum system temperature (K)";
  //Legacy terms from the CSP field. PB (power block) refers to the boiler, Recv (receiver) refers to the electric heater.
  //Level-Controls
  parameter SI.Time t_stor_start_dis = 1.0 * 3600.0 "Number of effective storage seconds stored before TES can start discharging (1 hour)";
  //Calculated Parameters
  parameter Modelica.SIunits.Energy E_max = t_storage * 3600.0 * Q_process_des "Maximum tank stored energy (J)";
  parameter Modelica.SIunits.HeatFlowRate Q_process_des = 20.398e6 "Heat-rate to process at design (W)";
  parameter Modelica.SIunits.MassFlowRate m_process_des = Q_process_des / (h_air_process_des - h_air_min_des) "Design process input mass flow rate (kg/s)";
  
parameter Medium.ThermodynamicState state_air_min_des = Medium.setState_pTX(Medium.p_default, T_min) "Thermodynamic state of air at the minimum system temperature T_min.";
  parameter Medium.ThermodynamicState state_air_max_des = Medium.setState_pTX(Medium.p_default, T_max) "Thermodynamic state of air at the maximum system temperature T_min.";
  parameter Medium.ThermodynamicState state_air_process_des = Medium.setState_pTX(Medium.p_default, T_process_des) "Thermodynamic state of air at the design process inlet T_process_des.";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_min_des = Medium.specificEnthalpy(state_air_min_des) "Specific enthalpy of air at minimum system temperature T_min (J/kg)";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_max_des = Medium.specificEnthalpy(state_air_max_des) "Specific enthalpy of air at maximum system temperature T_max (J/kg)";
  parameter Modelica.SIunits.SpecificEnthalpy h_air_process_des = Medium.specificEnthalpy(state_air_process_des) "Specific enthalpy of air at design process inlet temperature T_process_des (J/kg)";
  SolarTherm.Models.Storage.Thermocline.Parallel.Thermocline_HBS_LC_2P_MixedOutlet TES(redeclare package Medium = Medium, redeclare package Fluid_Package = Fluid, redeclare package Filler_Package_A = Filler, redeclare package Filler_Package_B = Filler, N_f_A = N_f, T_max = T_max, T_min = T_min, Correlation = Correlation, E_max = E_max, ar_A = ar, d_p_A = d_p, eta_A = eta, U_loss_top_A = U_loss_top, U_loss_bot_A = U_loss_bot, T_recv_set = T_heater_max, T_PB_set = T_process_des) annotation(
    Placement(visible = true, transformation(origin = {34, -8}, extent = {{-38, -38}, {38, 38}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpCold(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-20, -76}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction_LoopBreaker Splitter_Top(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {38, 62.6249}, extent = {{-22, -20.9366}, {22, 20.9366}}, rotation = 0)));
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction Splitter_Bot(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {32, -43.0917}, extent = {{20, 0}, {-20, -21.8345}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpHot(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {89, 79}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Tamb(y = 298.15) annotation(
    Placement(visible = true, transformation(origin = {-11, 2}, extent = {{-9, -12}, {9, 12}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression p_amb(y = 101325) annotation(
    Placement(visible = true, transformation(origin = {-11, -14}, extent = {{-9, -12}, {9, 12}}, rotation = 0)));
  SolarTherm.Models.Control.WindPV_Thermocline_Control Control(redeclare package HTF = Medium, E_max = E_max, Q_boiler_des = Q_process_des, T_boiler_min = T_process_min, T_boiler_start = T_process_start, T_heater_max = T_heater_max, T_heater_start = T_heater_start, T_target = T_max, util_storage_des = util_storage_des, h_target = h_air_max_des, level_mid = level_storage_mid, m_0 = 1e-8, m_boiler_des = m_process_des*(h_air_process_des-h_air_min_des)/(h_air_max_des-h_air_min_des), m_min = 1e-8, m_tol = 0.01 * m_process_des, t_stor_start_dis = t_stor_start_dis, t_wait = 1.0 * 3600.0) annotation(
    Placement(visible = true, transformation(origin = {114, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0))); //m_boiler_des needs to be scaled down to the mass flow rate at the overdesigned temperature T_max.
  SolarTherm.Models.Fluid.HeatExchangers.Boiler_Basic Process(redeclare package Medium = Medium, T_cold_set = T_min, T_hot_set = T_process_des) annotation(
    Placement(visible = true, transformation(origin = {158, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.CSP.CRS.Receivers.Basic_Heater basic_Heater(redeclare package Medium = Medium, P_heater_des = P_heater_des, Q_heater_des = Q_heater_des, eff_heater = eff_heater, T_cold_set = T_min, T_hot_set = T_max) annotation(
    Placement(visible = true, transformation(origin = {-41, 11}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable PV_input(fileName = PV_file, tableName = "Power", tableOnFile = true, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
    Placement(visible = true, transformation(origin = {-124, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add Grid_Sum(k1 = P_PV_gross / PV_ref_size, k2 = P_wind_gross / Wind_ref_size) annotation(
    Placement(visible = true, transformation(origin = {-84, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Wind_input(fileName = Wind_file, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "Power", tableOnFile = true) annotation(
    Placement(visible = true, transformation(origin = {-124, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression ConstantDemand(y = Q_process_des) annotation(
    Placement(visible = true, transformation(origin = {182, 42}, extent = {{16, -10}, {-16, 10}}, rotation = 0)));
    
  //Cost-base parameters
  parameter SI.MassFlowRate m_flow1_max = 26.6 "Maximum mass flow rate of blower 1 (kg/s) at 317.3C";
  parameter SI.MassFlowRate m_flow2_max = HM*m_flow1_max*(h_air_max_des-h_air_min_des)/(h_air_process_des-Fluid.h_Tf(T_heater_des,0.0)) "Maximum mass flow rate of blower 2 (kg/s) at 450C";
  parameter SI.VolumeFlowRate V_flow1_max = m_flow1_max/Fluid.rho_Tf(25.0+273.15,0.0);
  parameter SI.VolumeFlowRate V_flow2_max = m_flow2_max/Fluid.rho_Tf(T_heater_des,0.0);
  parameter SI.Pressure p_high = 122000.0 "Blower outlet pressure assuming pressure drop of ~3psig (Pa)";
  parameter SI.Pressure p_amb_des = 100000.0 "Design ambient pressure (Pa)";
  parameter SI.Power P_C1 = (0.9855*(1.4/0.4)*(V_flow1_max*p_high/0.75)*(((p_high/p_amb_des)^(0.4/1.4))-1) )/0.90 "Sizing power of blower 1 (W)";
  parameter SI.Power P_C2 = (0.9855*(1.4/0.4)*(V_flow2_max*p_high/0.75)*(((p_high/p_amb_des)^(0.4/1.4))-1) )/0.90 "Sizing power of blower 2 (W)";
  
  parameter SI.Length L_piping1 = max(1.5*TES.Tank_A.D_tank, 50) "Length of piping 1 (m)";
  parameter SI.Length L_piping2 = max(1.5*TES.Tank_A.D_tank + 0.75*TES.Tank_A.H_tank, 50) "Length of piping 2 (m)";
  parameter SI.Length L_piping3 = max(1.5*TES.Tank_A.D_tank + 0.75*TES.Tank_A.H_tank, 50) "Length of piping 3 (m)";
  parameter SI.Length L_piping4 = max(1.5*TES.Tank_A.D_tank, 50) "Length of piping 4 (m)";
  parameter SI.Length L_piping5 = max(1.5*TES.Tank_A.H_tank, 50) "Length of piping 5 (m)";
  parameter Real FCIpL_piping1 = (I_year/816.0)*479.77 "FCI cost per metre of piping 1 (USD/m)";
  parameter Real FCIpL_piping2 = (I_year/816.0)*1299.45 "FCI cost per metre of piping 2 (USD/m)";
  parameter Real FCIpL_piping3 = (I_year/816.0)*2060.85 "FCI cost per metre of piping 3 (USD/m)";
  parameter Real FCIpL_piping4 = (I_year/816.0)*1291.48"FCI cost per metre of piping 4 (USD/m)";
  parameter Real FCIpL_piping5 = (I_year/816.0)*60.21 "FCI cost per metre of piping 5 (USD/m)";
  //Costs
  //TES
  parameter Real I_year = 816.0 "CEPCI for 2022";
  
  parameter Real FOB_filler = (I_year/816.0)*TES.C_filler "FOB cost of checkerbrick (USD)";
  parameter Real FOB_insulation = (I_year/816.0)*TES.C_insulation "FOB cost of TES insulation (USD)";
  parameter Real FOB_tank = (I_year/816.0)*TES.C_tank "FOB cost of TES tank shell (USD)";
  parameter Real FOB_heater = P_heater_des*0.206/(4.0*1.05) "FOB cost of the air heaters (USD)";
  parameter Real FOB_HX = (I_year/500.0)*(div(730.406,185.8)*6200.0*((10.764*185.8)^0.42) +  6200.0*((10.764*rem(730.406,185.8))^0.42)) "FOB cost of gas-gas HXs at the return air stream";
  parameter Real FOB_blower1 = (I_year/500.0)*1.0*(div(P_C1,745700.0)*exp(6.8929+0.79*log(745700.0/745.7)) + exp(6.8929+0.79*log(rem(P_C1,745700.0)/745.7)))"FOB cost of blower 1, cast iron (USD)";
  parameter Real FOB_blower2 = (I_year/500.0)*1.0*(div(P_C2,745700.0)*exp(6.8929+0.79*log(745700.0/745.7)) + exp(6.8929+0.79*log(rem(P_C2,745700.0)/745.7)))"FOB cost of blower 2, cast iron (USD)";
  
  parameter Real FCI_piping1 = L_piping1*FCIpL_piping1;
  parameter Real FCI_piping2 = L_piping2*FCIpL_piping2;
  parameter Real FCI_piping3 = L_piping3*FCIpL_piping3;
  parameter Real FCI_piping4 = L_piping4*FCIpL_piping4;
  parameter Real FCI_piping5 = L_piping5*FCIpL_piping5;
  
  parameter Real FCI_filler = FOB_filler*1.05*4.0; // 1.05 is for delivery cost, 4 is a hand factor (from Sieder)
  parameter Real FCI_insulation = FOB_insulation*1.05*4.0;
  parameter Real FCI_tank = FOB_tank*1.05*4.0;
  parameter Real FCI_heater = FOB_heater*1.05*4.0;
  parameter Real FCI_HX = FOB_HX*1.05*3.5;
  parameter Real FCI_blower1 = FOB_blower1*1.05*4.0; // fresh air blower
  parameter Real FCI_blower2 = FOB_blower2*1.05*4.0; // charging loop blower (before heater)

  
  parameter Real FCI_PV = 1.075*P_PV_gross;
  parameter Real FCI_wind = 1.4622*P_wind_gross;
  
  parameter Real FCI_total = FCI_filler + FCI_insulation + FCI_tank + FCI_heater + FCI_HX + FCI_blower1 +  FCI_blower2 + FCI_piping1 + FCI_piping2 + FCI_piping3 + FCI_piping4 + FCI_piping5 + FCI_PV + FCI_wind;

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
    Line(points = {{34, -38}, {34, -47.5}, {32, -47.5}, {32, -59}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Splitter_Bot.fluid_b, pumpCold.fluid_a) annotation(
    Line(points = {{16, -76}, {-10, -76}}, color = {0, 127, 255}, thickness = 0.5));
  connect(TES.T_top_measured, Control.T_top_tank) annotation(
    Line(points = {{51, 13}, {78, 13}, {78, 28}, {103, 28}}, color = {0, 0, 127}));
  connect(TES.T_bot_measured, Control.T_bot_tank) annotation(
    Line(points = {{51, -29}, {82, -29}, {82, 24}, {103, 24}}, color = {0, 0, 127}));
  connect(TES.Level, Control.Level) annotation(
    Line(points = {{51, 0}, {90, 0}, {90, 28}, {103, 28}, {103, 33}}, color = {0, 0, 127}));
  connect(Splitter_Top.fluid_c, TES.fluid_a) annotation(
    Line(points = {{38, 62}, {38, 42}, {34, 42}, {34, 22}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Splitter_Top.fluid_b, pumpHot.fluid_a) annotation(
    Line(points = {{56, 79}, {80, 79}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Process.fluid_b, Splitter_Bot.fluid_a) annotation(
    Line(points = {{158, -10}, {158, -76}, {48, -76}}, color = {0, 127, 255}, thickness = 0.5));
  connect(pumpHot.fluid_b, Process.fluid_a) annotation(
    Line(points = {{98, 79}, {158, 79}, {158, 10}}, color = {0, 127, 255}, thickness = 0.5));
  connect(Control.curtail, basic_Heater.curtail) annotation(
    Line(points = {{125, 26}, {128, 26}, {128, -52}, {-68, -52}, {-68, 20}, {-51, 20}}, color = {255, 0, 255}));
  connect(TES.h_top_outlet, Control.h_tank_top) annotation(
    Line(points = {{24, 17}, {24, 46}, {105, 46}, {105, 41}}, color = {0, 0, 127}));
  connect(PV_input.y[1], Grid_Sum.u1) annotation(
    Line(points = {{-113, 34}, {-104, 34}, {-104, 24}, {-96, 24}}, color = {0, 0, 127}));
  connect(Wind_input.y[1], Grid_Sum.u2) annotation(
    Line(points = {{-113, 4}, {-104, 4}, {-104, 12}, {-96, 12}}, color = {0, 0, 127}));
  connect(Grid_Sum.y, basic_Heater.P_supply) annotation(
    Line(points = {{-73, 18}, {-66, 18}, {-66, 17}, {-51, 17}}, color = {0, 0, 127}));
  connect(Process.h_out_signal, Control.h_boiler_outlet) annotation(
    Line(points = {{148, -8}, {136, -8}, {136, 52}, {115, 52}, {115, 41}}, color = {0, 0, 127}));
  connect(TES.h_bot_outlet, Control.h_tank_bot) annotation(
    Line(points = {{24, -33}, {132, -33}, {132, 46}, {110, 46}, {110, 41}}, color = {0, 0, 127}));
  connect(basic_Heater.Q_heater_raw, Control.Q_heater_raw) annotation(
    Line(points = {{-31, 18}, {-20, 18}, {-20, 52}, {86, 52}, {86, 37}, {103, 37}}, color = {0, 0, 127}));
  connect(Control.m_heater_signal, pumpCold.m_flow) annotation(
    Line(points = {{125, 36}, {140, 36}, {140, -56}, {-20, -56}, {-20, -67}}, color = {0, 0, 127}));
  connect(Control.m_boiler_signal, pumpHot.m_flow) annotation(
    Line(points = {{125, 31}, {130, 31}, {130, 94}, {90, 94}, {90, 87}, {89, 87}}, color = {0, 0, 127}));
  connect(ConstantDemand.y, Control.Q_demand) annotation(
    Line(points = {{164, 42}, {120, 42}, {120, 41}, {121, 41}}, color = {0, 0, 127}));
  connect(Control.Q_curtail, basic_Heater.Q_curtail) annotation(
    Line(points = {{104, 22}, {-60, 22}, {-60, 14}, {-52, 14}, {-52, 14}, {-52, 14}}, color = {0, 0, 127}));
  connect(basic_Heater.fluid_b, Splitter_Top.fluid_a) annotation(
    Line(points = {{-32, 12}, {-26, 12}, {-26, 79}, {20, 79}}, color = {0, 127, 255}, thickness = 0.5));
  connect(pumpCold.fluid_b, basic_Heater.fluid_a) annotation(
    Line(points = {{-30, -76}, {-58, -76}, {-58, 12}, {-50, 12}, {-50, 12}}, color = {11, 133, 255}, thickness = 0.55));
  connect(Tamb.y, TES.T_amb) annotation(
    Line(points = {{-2, 2}, {9, 2}, {9, -8}, {17, -8}}, color = {0, 0, 127}));
  connect(p_amb.y, TES.p_amb) annotation(
    Line(points = {{-2, -14}, {9, -14}, {9, -8}, {51, -8}}, color = {0, 0, 127}));
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-200, -100}, {200, 100}}, initialScale = 0.1)),
    Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}}, preserveAspectRatio = false)),
    experiment(StopTime = 3.1536e+07, StartTime = 0, Tolerance = 1.0e-5, Interval = 300, maxStepSize = 60, initialStepSize = 60));
end HBSTES_Case2_SystemLevel;