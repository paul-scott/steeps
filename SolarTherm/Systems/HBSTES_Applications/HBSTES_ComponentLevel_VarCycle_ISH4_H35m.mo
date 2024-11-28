within SolarTherm.Systems.HBSTES_Applications;

model HBSTES_ComponentLevel_VarCycle_ISH4_H35m "This component analysis allows a variable cycle length depending on t_storage and t_standby. Tank heights are fixed to 35m for all cases."
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  package Medium = SolarTherm.Media.Air.Air_CoolProp_1bar;
  package Fluid_Package = SolarTherm.Materials.Air_CoolProp_Table_1bar;
  package Filler_Package = SolarTherm.Materials.Mullite_20pct_porosity;
  //These parameters are varied
  parameter Real HM = 2.0 "Heater Multiple";
  parameter SI.Time t_discharge = 20.0 * 3600.0 "Rated discharging period (s)";
  parameter SI.Time t_standby = 2.0*3600.0 "Additional standby time after both charging and discharging cycle (s)";
  //Numerical Discretisation Settings
  parameter Integer N_f = 100;
  parameter Integer N_tanks = 1 "Number of tanks in the storage array, needed to adjust aspect ratio such that all tanks are 35m high";
  parameter Real f_TES_oversize = 1.0;
  parameter Real t_storage_hours_ideal = t_discharge*f_TES_oversize/3600.0;
  //parameter Integer N_p = 5; //Not used
  //TES Geometric Parameters
  parameter SI.Length d_p = 0.03 "Hole diameter in the filler (m)";
  parameter Real ar =  11.336*(N_tanks^0.5)*(2.93731e11/E_max)^0.5  "Tank H/D ratio (m)";
  parameter Real eta = 0.53 "Packed-bed porosity";
  //parameter SI.Length s_p = 0.04 "Separation of holes in the filler (m)";
  parameter SI.Temperature T_ext_max = 50.0 + 273.15 "Maximum external temperature of the HBS walls under worst-case conditions (K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_top = 10.0*(T_ext_max-298.15)/(T_max-298.15) "Heat loss coefficient at the top of the tank (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_bot = 10.0*(T_ext_max-298.15)/(T_max-298.15) "Heat loss coefficient at the bottom of the tank (W/m2K)";
  //Temperature Controls
  parameter SI.Temperature T_max = 1200.0 + 273.15 "Maximum temperature (K)";
  parameter SI.Temperature T_process_des = 1000.0 + 273.15 "Design process inlet temperature (K)";
  parameter SI.Temperature T_high_set = 1000.0 + 273.15 "TES hot blend temperature temperature (K)";
  parameter SI.Temperature T_process_min = 1000.0 + 273.15 "Minimum tolerated outlet temperature to process (K)";
  parameter SI.Temperature T_heater_max = 450.0 + 273.15 "Maximum tolerated outlet temperature to heater (K)";
  parameter SI.Temperature T_low_set = 450.0 + 273.15 "TES cold blend temperature (K)";
  parameter SI.Temperature T_heater_des = 317.3 + 273.15 "Design receiver inlet temperature (K)";
  parameter SI.Temperature T_min = 317.3 + 273.15 "Minimum temperature (K)";
  parameter Integer Correlation = 1; //Gas Internal flow
  //parameter SI.SpecificEnthalpy h_tol = 0.05*(TES.Tank_A.h_f_max - TES.Tank_A.h_f_min);
  //Tank Geometry
  parameter SI.Power Q_process_des = 26.6*(SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.h_T(T_process_des)-SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.h_T(T_min)) ;//E_max/t_discharge "Design required process heat-rate (W_th)";//E_max/t_discharge
  parameter SI.Power Q_heater_des = HM * Q_process_des "Design heater output heat-rate (W_th)";
  parameter SI.Energy E_max = Q_process_des*t_discharge*f_TES_oversize "Ideal storage capacity (J_thermal)"; //Note 3 tanks
  parameter SI.Time t_charge = t_discharge / (HM - 1.0) "Charging period (s)";
  parameter SI.MassFlowRate m_charge_des = (Q_heater_des - Q_process_des) / (h_f_max - h_f_min) "Design charging mass flow rate assuming design temperature outlet (kg/s)";
  parameter SI.MassFlowRate m_discharge_des = Q_process_des / (h_f_max - h_f_min)  "Design discharging mass flow rate assuming design temperature outlet (kg/s)";
  //Derived Parameters
  parameter SI.Time t_cycle = t_charge + t_discharge + 2.0*t_standby "Total duration of one cycle (s)";
  //this is 24 hours
  parameter SI.SpecificEnthalpy h_f_min = Fluid_Package.h_Tf(T_min, 0.0);
  parameter SI.SpecificEnthalpy h_f_max = Fluid_Package.h_Tf(T_max, 1.0);
  //Output signals
  Modelica.Fluid.Sources.Boundary_pT Heater_outlet(redeclare package Medium = Medium, T = T_max, nPorts = 1, p = 100000) annotation(
    Placement(visible = true, transformation(origin = {-112, 48}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT Process_outlet(redeclare package Medium = Medium, T = T_min, nPorts = 1, p = 100000) annotation(
    Placement(visible = true, transformation(origin = {107, -59}, extent = {{15, -15}, {-15, 15}}, rotation = 0)));
  //Storage Model
  //Componenets and Connectors
  SolarTherm.Models.Fluid.Sources.FluidSink2 Recv_Sink(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-125, -59}, extent = {{21, -21}, {-21, 21}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Tamb(y = 298.15) annotation(
    Placement(visible = true, transformation(origin = {-40, -4}, extent = {{-12, -18}, {12, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow_Heater(y = m_Heater_signal) annotation(
    Placement(visible = true, transformation(origin = {-119, 3}, extent = {{-19, -17}, {19, 17}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure heater_sink_pump(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = { -57, -59}, extent = {{15, -15}, {-15, 15}}, rotation = 0)));
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction thermocline_Splitter1(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-2, 67.5547}, extent = {{-16, 0}, {16, 22.4453}}, rotation = 0)));
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction thermocline_Splitter2(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-2, -36.3493}, extent = {{-14, 0}, {14, 21.6507}}, rotation = 180)));
  Modelica.Blocks.Sources.RealExpression m_flow_Process(y = m_Process_signal) annotation(
    Placement(visible = true, transformation(origin = {120, 3}, extent = {{20, -19}, {-20, 19}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple process_source_pump(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {45, -59}, extent = {{13, -13}, {-13, 13}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink2 Process_Sink(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {112, 44}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression p_amb(y = 100000) annotation(
    Placement(visible = true, transformation(origin = {43, -4}, extent = {{13, -16}, {-13, 16}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple heater_source_pump(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-52, 48}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure process_sink_pump(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {46, 44}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  SolarTherm.Models.Fluid.HeatExchangers.mass_loop_breaker mass_loop_breaker annotation(
    Placement(visible = true, transformation(origin = {-2, 50}, extent = {{-24, -24}, {24, 24}}, rotation = -90)));
  SolarTherm.Models.Storage.Thermocline.Thermocline_HBS_LC_SingleTank_Final TES(redeclare package Medium = Medium, redeclare package Fluid_Package = Fluid_Package, redeclare package Filler_Package = Filler_Package, N_f = N_f, T_max = T_max, T_min = T_min, Correlation = Correlation, E_max = E_max, ar = ar, d_p = d_p, eta = eta, U_loss_top = U_loss_top, U_loss_bot = U_loss_bot) annotation(
    Placement(visible = true, transformation(origin = {-2, -4}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  //Mass flow Signals starts in charging state //,h_tol=h_tol
  SI.MassFlowRate m_Heater_signal(start = m_charge_des);
  SI.MassFlowRate m_Process_signal(start = 0.0);
  //Measured Temperatures
  Real T_top_degC;
  Real T_bot_degC;
  Real T_outlet_degC;
  Real T_inlet_degC;
  //Energies and Exergies for final cycle calculation
  SI.Energy E_chg(start = 0) "Charged energy into storage";
  SI.Energy E_dis(start = 0) "Discharged energy into storage";
  //SI.Energy E_lost(start = 0) "Externally lost energy from storage";
  //SI.Energy E_pump(start = 0) "Pumping energy consumed";
  parameter SI.Energy X_max = m_discharge_des * t_discharge * f_TES_oversize * (h_f_max * (1.0 - 298.15 / T_max) - h_f_min * (1.0 - 298.15 / T_min)) "Theoretical Exergetic storage capacity (J)";
  //SI.Energy X_chg(start = 0.0) "Total exergy charged i3n the last cycle";
  //SI.Energy X_dis(start = 0.0) "Discharged exergy";
  //Utilisation and efficiencies
  Real util_energy(start = 0.0) "Energetic utilisation (1st Law)";
  Real eff_energy(start = 0.0) "Energetic efficiency (1st Law)";
  //Real util_exergy(start = 0.0) "Exergetic utilisation (2nd law)";
  //Real eff_exergy(start = 0.0) "Exergetic efficiency (2nd law)";
  //Individual Tank Utilization Calculations
  parameter SI.Energy E_max_A = TES.Tank_A.E_max "Maximum energy capacity of Tank A";
  //parameter SI.Energy E_max_B = thermocline_Tank.Tank_B.E_max "Maximum energy capacity of Tank B";
  //parameter SI.Energy E_max_C = thermocline_Tank.Tank_C.E_max "Maximum energy capacity of Tank C";
  //SI.Energy E_dis_A (start = 0.0) "Discharged energy from Tank A";
  //SI.Energy E_dis_B (start = 0.0) "Discharged energy from Tank B";
  //SI.Energy E_dis_C (start = 0.0) "Discharged energy from Tank C";
  //SI.Energy E_pump_A (start = 0.0) "Pumping energy of Tank A";
  //SI.Energy E_pump_B (start = 0.0) "Pumping energy of Tank B";
  //SI.Energy E_pump_C (start = 0.0) "Pumping energy of Tank C";
  //Real util_energy_A(start = 0.0) "Storage utilisation of Tank A";
  //Real util_energy_B(start = 0.0) "Storage utilisation of Tank B";
  //Real util_energy_C(start = 0.0) "Storage utilisation of Tank C";
  //Stats to determine if stable cycle has been reached
  //SI.Energy E_dis_total(start = 0.0) "Total discharged energy for entire simulation";
  //SI.Energy E_pump_total(start = 0.0) "Total pump energy losses for entire simulation";
  //SI.Energy E_dis_prev(start = 0.001) "Total discharged energy for previous cycle";
  //SI.Energy E_pump_prev(start = 0.001) "Total pump energy losses for previous cycle";
  //Real util_prev(start = 0.001) "Previous day's level, starts at 0.1% to prevent divbyzero";
  //Real time_stop(start = 3600.0 * 24.0 * 10.0) "The maximum time this simulation can last";
  Boolean Chg(start = false) "Can storage be charged?";
  Boolean Dis(start = false) "Can storage be discharged?";
  //Integer Control_State(start = 1) "Control state 1= Charge 2=Standby 3=Discharge";
  Real Level_high(start = 1.0);
  Real Level_low(start = 0.0);
  Real Level_mid(start = 0.5);
  
  SI.Energy E_stored(start = 0.0);
  
  //Exergy Analysis
  SI.Energy X_chg(start=0.0) "Total exergy charged in the last cycle";
  SI.Energy X_dis(start=0.0) "Discharged exergy";
  Real util_exergy(start=0.0) "Exergetic utilisation (2nd law)";
  Real eff_exergy(start=0.0) "Exergetic efficiency (2nd law)";
algorithm
//Mass flow controls
  when rem(time, t_cycle) > 1e-6 then //Start of charging cycle
    Dis := false;
    if TES.T_bot_measured < T_heater_max then
      Chg := true;
    end if;
  end when;
  
  when rem(time, t_cycle) > 1e-6 + t_charge then //End of charging cycle, begin standby
    Chg := false;
  end when;
  
  when rem(time, t_cycle) > 1e-6 + t_charge + t_standby then //begin discharging cycle
    if TES.T_top_measured > T_process_min then
      Dis := true;
    end if;
  end when;
  
   when rem(time, t_cycle) > 1e-6 + t_charge + t_standby + t_discharge then //end discharging cycle
    Dis := false;
   end when;
  
   
  
  when TES.T_bot_measured > T_heater_max then
    Chg := false;
    Level_high := TES.Level;
  end when;
  when TES.T_top_measured < T_process_min then
    Dis := false;
    Level_low := TES.Level;
    Level_mid := 0.5 * (Level_high + Level_low);
  end when;
//Force stop discharging, unacceptable top outlet Temperature
//Utilisation convergence
  when time > 10.0 * t_cycle then
    terminate("Ten cycles have been reached");
  end when;
//Replace time + t_cycle - 1.0e-6 with 3600.0*24.0*10.0 to bypass the 1% check and simulate the full 10 days
//if abs((((E_dis_total - E_dis_prev) - (E_pump_total - E_pump_prev))/E_max) - util_prev)/(util_prev) < 0.01 and div(time,t_cycle) > 0.0 then //stable cycle
//stable cycle
//time + t_cycle - 1.0e-6;
//calculate utilisation
//util_prev := ((E_dis_total - E_dis_prev) - (E_pump_total - E_pump_prev))/E_max;
//E_pump_prev := E_pump_total;
equation
  der(E_stored) = TES.fluid_a.m_flow * (TES.fluid_a.h_outflow - TES.fluid_b.h_outflow);
  if Chg then
    m_Heater_signal = m_charge_des * (h_f_max - h_f_min) / (h_f_max - TES.h_bot_outlet);
  else
    m_Heater_signal = 1.0e-8;
  end if;
  if Dis then
    m_Process_signal = m_discharge_des * (h_f_max - h_f_min) / (TES.h_top_outlet - h_f_min);
  else
    m_Process_signal = 1.0e-8;
  end if;
//if Control_State == 1 then
//Charging
//m_Recv_signal = m_charge_des * (h_f_max - h_f_min) / (h_f_max - thermocline_Tank.h_bot_outlet);
//m_PB_signal = 1.0e-8;
//elseif Control_State == 2 then
//Standby
//m_Recv_signal = 1.0e-8;
//m_PB_signal = 1.0e-8;
//else
//m_Recv_signal = 1.0e-8;
//m_PB_signal = m_discharge_des * (h_f_max - h_f_min) / (thermocline_Tank.h_top_outlet - h_f_min);
//end if;
//Measured temperature sensors
  T_top_degC = TES.T_top_measured - 273.15;
  T_bot_degC = TES.T_bot_measured - 273.15;
//Combined storage outlet temperatures (25degC if no flow)
  if TES.fluid_a.m_flow > 1.0e-3 then
//there is mass flow
    T_outlet_degC = TES.fluid_bot.T - 273.15;
    T_inlet_degC = TES.fluid_top.T - 273.15;
  elseif TES.fluid_a.m_flow < (-1.0e-3) then
    T_outlet_degC = TES.fluid_top.T - 273.15;
    T_inlet_degC = TES.fluid_bot.T - 273.15;
  else
    T_outlet_degC = 25.0;
    T_inlet_degC = 25.0;
  end if;
  if time > 9.0 * t_cycle and time < 9.0 * t_cycle + t_charge then
//last charging phase
    der(E_chg) = TES.fluid_a.m_flow * (inStream(TES.fluid_a.h_outflow) - TES.fluid_b.h_outflow);
    der(E_dis) = 0.0;
    der(X_chg) = TES.fluid_a.m_flow * (inStream(TES.fluid_a.h_outflow)*(1.0 - (298.15/TES.fluid_top.T)) - (TES.fluid_b.h_outflow)*(1.0 - (298.15/TES.fluid_bot.T)));
    der(X_dis) = 0.0;
  elseif time >= 9.0 * t_cycle + t_charge + t_standby and time < 9.0 * t_cycle + t_charge + t_standby + t_discharge then
//last discharging phase
    der(E_chg) = 0.0;
    der(E_dis) = TES.fluid_b.m_flow * (TES.fluid_a.h_outflow - inStream(TES.fluid_b.h_outflow));
    der(X_chg) = 0.0;
    der(X_dis) = TES.fluid_b.m_flow * ((TES.fluid_a.h_outflow)*(1.0-(298.15/TES.fluid_top.T)) - (inStream(TES.fluid_b.h_outflow)) * (1.0-(298.15/TES.fluid_bot.T)));
  else
    der(E_chg) = 0.0;
    der(E_dis) = 0.0;
    der(X_chg) = 0.0;
    der(X_dis) = 0.0;  
  end if;
  if time > 9.0 * t_cycle + t_charge + t_standby then
//we can calculate utilisation and efficiency because denominator is not zero;
    util_energy = E_dis / E_max;
    eff_energy = E_dis / E_chg;
    util_exergy = (X_dis) / X_max;
    eff_exergy = (X_dis) / X_chg;
  else
    util_energy = 0.0;
    eff_energy = 0.0;
    util_exergy = 0.0;
    eff_exergy = 0.0;
  end if;
//Connectors
  connect(thermocline_Splitter2.fluid_b, heater_sink_pump.fluid_a) annotation(
    Line(points = {{-13, -59}, {-42, -59}}, color = {0, 127, 255}));
  connect(m_flow_Heater.y, heater_sink_pump.m_flow) annotation(
    Line(points = {{-98, 3}, {-57, 3}, {-57, -46}}, color = {0, 0, 127}));
  connect(process_source_pump.fluid_b, thermocline_Splitter2.fluid_a) annotation(
    Line(points = {{32, -59}, {9, -59}}, color = {0, 127, 255}));
  connect(heater_source_pump.fluid_b, thermocline_Splitter1.fluid_a) annotation(
    Line(points = {{-40, 48}, {-30, 48}, {-30, 92}, {-15, 92}}, color = {0, 127, 255}));
  connect(thermocline_Splitter1.fluid_b, process_sink_pump.fluid_a) annotation(
    Line(points = {{11, 92}, {22, 92}, {22, 44}, {34, 44}}, color = {0, 127, 255}));
  connect(Heater_outlet.ports[1], heater_source_pump.fluid_a) annotation(
    Line(points = {{-96, 48}, {-64, 48}}, color = {0, 127, 255}));
  connect(Process_outlet.ports[1], process_source_pump.fluid_a) annotation(
    Line(points = {{92, -59}, {58, -59}}, color = {0, 127, 255}));
  connect(m_flow_Heater.y, heater_source_pump.m_flow) annotation(
    Line(points = {{-98, 3}, {-74, 3}, {-74, 72}, {-52, 72}, {-52, 58}}, color = {0, 0, 127}));
  connect(m_flow_Process.y, process_sink_pump.m_flow) annotation(
    Line(points = {{98, 3}, {70, 3}, {70, 70}, {46, 70}, {46, 54}}, color = {0, 0, 127}));
  connect(m_flow_Process.y, process_source_pump.m_flow) annotation(
    Line(points = {{98, 3}, {70, 3}, {70, -30}, {45, -30}, {45, -48}}, color = {0, 0, 127}));
  connect(Recv_Sink.port_a, heater_sink_pump.fluid_b) annotation(
    Line(points = {{-104, -59}, {-72, -59}}, color = {0, 127, 255}));
  connect(Process_Sink.port_a, process_sink_pump.fluid_b) annotation(
    Line(points = {{94, 44}, {58, 44}}, color = {0, 127, 255}));
  connect(thermocline_Splitter1.fluid_c, mass_loop_breaker.port_a) annotation(
    Line(points = {{-2, 74}, {-2, 64}}, color = {0, 127, 255}));
  connect(mass_loop_breaker.port_b, TES.fluid_a) annotation(
    Line(points = {{-2, 36}, {-2, 20}}, color = {0, 127, 255}));
  connect(TES.fluid_b, thermocline_Splitter2.fluid_c) annotation(
    Line(points = {{-2, -28}, {-2, -42}}, color = {0, 127, 255}));
  connect(Tamb.y, TES.T_amb) annotation(
    Line(points = {{-27, -4}, {-16, -4}}, color = {0, 0, 127}));
  connect(p_amb.y, TES.p_amb) annotation(
    Line(points = {{29, -4}, {12, -4}}, color = {0, 0, 127}));
  annotation(
    experiment(StopTime = 2304000, StartTime = 0, Tolerance = 1e-4, Interval = 60),
    Diagram(coordinateSystem(extent = {{-150, -100}, {150, 100}}, preserveAspectRatio = false)),
    Icon(coordinateSystem(extent = {{-150, -100}, {150, 100}}, preserveAspectRatio = false)));
end HBSTES_ComponentLevel_VarCycle_ISH4_H35m;