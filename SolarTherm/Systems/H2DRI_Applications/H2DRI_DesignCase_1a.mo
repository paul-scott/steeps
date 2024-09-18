within SolarTherm.Systems.H2DRI_Applications;

model H2DRI_DesignCase_1a
  //0D Model
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  
  //Medium Packages
  replaceable package Medium_H2 = Modelica.Media.IdealGases.SingleGases.H2;
  replaceable package Medium_H2O = Modelica.Media.Water.WaterIF97_ph;
  replaceable package Medium_Fe2O3 = SolarTherm.Media.SolidParticles.Fe2O3_ph;
  
  //Chemical Constants
  parameter SI.MolarMass M_Fe2O3 = SolarTherm.Models.Chemistry.ChemTable.Fe2O3.M;
  parameter SI.MolarMass M_H2O = SolarTherm.Models.Chemistry.ChemTable.H2O.M;
  parameter SI.MolarMass M_H2 = SolarTherm.Models.Chemistry.ChemTable.H2.M;

  //Plant Constants e.g. reference/feedstock conditions
  parameter SI.Pressure p_des = 100000.0 "Design pressure of the plant (Pa)";
  parameter SI.Temperature T_H2_feedstock_des = 25.0+273.15 "Design feedstock H2 temperature (K)";
  parameter SI.Temperature T_Fe2O3_feedstock_des = 25.0+273.15 "Design feedstock Fe2O3 temperature (K)";
  parameter SI.SpecificEnthalpy h_H2_feedstock_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_H2_feedstock_des) "Specific Enthalpy of feedstock H2 (J/kg)";
  
  //Design Inputs
  parameter SI.MassFlowRate m_flow_Fe_des = 47.6 "Design maximum mass flow rate output of Iron (kg/s)";
  //parameter SI.Temperature T_reactants_des = 757.0 + 273.15 "Design reactant temperature (K)";
  parameter SI.Temperature T_products_des = 600.0 + 273.15 "Design product temperature (K)";
  parameter Real eff_GGHX_des = 0.80 "Effectiveness of the GGHX (-)";
  
  //Design Reactor Inlet Mixed Enthalpies
  parameter SI.SpecificEnthalpy h_H2_reactant_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_reactants_des);
  parameter SI.SpecificEnthalpy h_Fe2O3_reactant_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_reactants_des);
  
  //Calculated Parameters
    //Rectant and product temperature relationship is estimated using 2' polynomial
  parameter SI.Temperature T_reactants_des = -304.2 + 1.851*T_products_des - ((4.127E-4)*(T_products_des^2.0)) "Design estimated reactant temperature (K)";
  //parameter SI.Temperature T_products_des = 263.0 + 0.3380*T_reactants_des - ((2.744E-4)*(T_reactants_des^2.0)) "Design estimated product temperature (K)";
  parameter SI.Temperature T_condenser_out_des = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p_H2O_offgas_des) "Design saturation temperature of the condenser (K)";
  parameter SI.SpecificEnthalpy h_H2_condenser_out_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_condenser_out_des) "Design specific enthalpy of H2 leaving the condenser (J/kg)";
  
  parameter SI.SpecificEnthalpy h_Fe2O3_hot_des = ( m_flow_H2_des*(h_H2_reactant_des-h_H2_hot_des) + m_flow_Fe2O3_des*h_Fe2O3_reactant_des )/m_flow_Fe2O3_des;

  parameter SI.SpecificEnthalpy h_H2_mix_des = (m_flow_H2_stoi*h_H2_feedstock_des + m_flow_H2_excess*h_H2_condenser_out_des)/(m_flow_H2_stoi+m_flow_H2_excess) "Design specific enthalpy of H2 entering the GGXH cold stream after being mixed with topup stream (J/kg)";
  parameter SI.Temperature T_H2_mix_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_mix_des) "Design temperature of H2 entering the GGXH cold stream after being mixed with topup stream (J/kg)";
  
  
  parameter SI.Temperature T_Fe2O3_hot_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_hot_des) "Design temperature of Fe2O3 in the hot TES (K)";

  //GGHX Outlet design temperature
  parameter SI.SpecificEnthalpy h_H2_GGHX_max=Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_products_des) "Design specific enthalpy of H2 leaving the GGHX cold stream assuming it had been heated to the maximum possible temperature i.e. reactor product temperature (J/kg)";
  parameter SI.SpecificEnthalpy h_H2_hot_des=h_H2_mix_des + eff_GGHX_des*(h_H2_GGHX_max-h_H2_mix_des) "Design of H2 leaving the GGHX cold stream after considering HX effectiveness and assuming C_C is C_min (J/kg)";



  //Storage Temperatures
  parameter SI.Temperature T_H2_hot_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_hot_des) "Design temperature of the hot H2 reservoir";
  
  //Flow Rates
  parameter SI.MassFlowRate m_flow_Fe2O3_des = 1.4297448*m_flow_Fe_des "Design required Fe2O3 flow rate (kg/s)";
  
  parameter SI.MolarFlowRate n_flow_Fe2O3_des = m_flow_Fe2O3_des/M_Fe2O3 "Design required molar Fe2O3 flow rate (mol/s)";
  parameter SI.MolarFlowRate n_flow_H2_stoi = 3.0*n_flow_Fe2O3_des "Design stoichiometric molar flow rate of H2 consumed assuming 100% yield of iron, also required top-up moles (mol/s)";
  parameter SI.MolarFlowRate n_flow_H2O_des = n_flow_H2_stoi "Design molar flow rate of H2O produced in the reactor assuming 100% yield of iron (mol/s)";
  parameter SI.MassFlowRate m_flow_H2O_des = n_flow_H2O_des*M_H2O "Design mass flow rate of H2O produced in the reactor assuming 100% yield of iron (kg/s)";
  parameter SI.MassFlowRate m_flow_H2_stoi = n_flow_H2_stoi*M_H2 "Design stoichiometric mass flow rate of H2 consumed assuming 100% yield of iron, also required top-up mass flow (kg/s)";
  
  parameter SI.MolarFlowRate n_flow_H2_des = n_flow_Fe2O3_des*min(1.05*(389.6 - (0.9702 * T_reactants_des) + (8.343e-4 * (T_reactants_des ^ 2)) - (2.413e-7 * (T_reactants_des ^ 3))), 20.0) "Design molar flow rate of H2 required to ensure 100% yield of iron (mol/s)";
  parameter SI.MassFlowRate m_flow_H2_des = n_flow_H2_des*M_H2 "Design mass flow rate of H2 required to ensure 100% yield of iron (kg/s)";
  parameter SI.MassFlowRate m_flow_H2_excess = m_flow_H2_des - m_flow_H2_stoi "Excess mass flow rate of H2 exiting the reactor (kg/s)";
  
  //Reactor Outlet Pressures
  parameter SI.Pressure p_H2_offgas_des = p_des*(n_flow_H2_des - n_flow_H2_stoi)/(n_flow_H2_des - n_flow_H2_stoi + n_flow_H2O_des) "Design partial pressure of H2 in the offgas (Pa)";
  parameter SI.Pressure p_H2O_offgas_des = p_des*(n_flow_H2O_des)/(n_flow_H2_des - n_flow_H2_stoi + n_flow_H2O_des) "Design partial pressure of H2O in the offgas (Pa)";
  
  //Reference enthalpies
  parameter SI.SpecificEnthalpy h_H2_ref = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(1.0e5,298.15);
  parameter SI.SpecificEnthalpy h_H2O_ref = Modelica.Media.Water.IF97_Utilities.h_pT(1.0e5,298.15);
  //Analytics
  Real yield = Reactor.m_flow_Fe_out/m_flow_Fe_des;
  
  Real mflowcp_H2;
  Real mflowcp_Fe2O3;

  SI.SpecificHeatCapacity cp_H2;
  SI.SpecificHeatCapacity cp_Fe2O3;
  
  

  //1kg/s Fe maximum. This is scaled till 10x
  //parameter SI.Temperature T_inlet_base = 975.0+273.15;
  //parameter SI.Temperature T_Fe2O3_hot_des = 1150.0 + 273.15 "Design reactor inlet temperature of Fe2O3 (K)";
  //parameter SI.Temperature T_H2_hot_des = 700.0 + 273.15 "Design reactor inlet temperature of H2 (K)";
  
  SI.Temperature T_Fe2O3_in(start=T_Fe2O3_hot_des);
  //SI.Temperature T_H2_in(start = 700.0 + 273.15);
  SI.HeatFlowRate Q_flow_cooling_signal;
  SI.SpecificEnthalpy h_sat_l;
  SI.Temperature T_sat;
  SI.HeatFlowRate Q_flow_H2_req "Required heating rate of preheated H2 gas stream";
  SI.HeatFlowRate Q_flow_Fe2O3_req "Required heating rate of preheated Fe2O3";
  Real r_min "Required molar ratio of H2 to Fe2O3 to ensure 100% yield";
  //parameter Real frac_1 = 0.7138 "Fraction of off-gas that is sent to the gas-gas HX";
  Modelica.Fluid.Sources.Boundary_pT Fe2O3_source(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph, nPorts = 1, p = 100000, use_T_in = true) annotation(
    Placement(visible = true, transformation(origin = {-132, -58}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  //Storage Model
  //Componenets and Connectors
  Modelica.Blocks.Sources.RealExpression m_flow_Fe2O3(y = m_Fe2O3_signal) annotation(
    Placement(visible = true, transformation(origin = {-175, -9}, extent = {{-19, -17}, {19, 17}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow_H2(y = m_H2_signal) annotation(
    Placement(visible = true, transformation(origin = {-178, 83}, extent = {{-20, -19}, {20, 19}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression p_amb(y = 100000) annotation(
    Placement(visible = true, transformation(origin = {31, 38}, extent = {{19, -14}, {-19, 14}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Q_flow_cooling(y = Q_flow_cooling_signal) annotation(
    Placement(visible = true, transformation(origin = {195, 72}, extent = {{13, -16}, {-13, 16}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Fe2O3_pump(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph) annotation(
    Placement(visible = true, transformation(origin = {-47, -55}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple H2_pump(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(origin = {-56, 38}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  //Mass flow Signals starts in charging state
  SI.MassFlowRate m_Fe2O3_signal(start = m_flow_Fe2O3_des);
  SI.MassFlowRate m_H2_signal(start = m_flow_H2_des);
  SolarTherm.Models.Reactors.Reactor_H2DRI_Overall Reactor(T_rxn(start = T_reactants_des))  annotation(
    Placement(visible = true, transformation(origin = {7, -19}, extent = {{-37, -37}, {37, 37}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink H2O_Sink(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) annotation(
    Placement(visible = true, transformation(origin = {206, 28}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression T_Fe2O3_signal(y = T_Fe2O3_in) annotation(
    Placement(visible = true, transformation(origin = {-175, -53}, extent = {{-19, -17}, {19, 17}}, rotation = 0)));
  SolarTherm.Models.Fluid.Condenser.Condenser_H2_H2O_Isobaric Condenser(p_total = p_des)  annotation(
    Placement(visible = true, transformation(origin = {161, 37}, extent = {{-31, -31}, {31, 31}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink H2O_Liq_Sink(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) annotation(
    Placement(visible = true, transformation(origin = {206, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Storage.Tank.Silo TES(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph, redeclare package Filler_Package = SolarTherm.Materials.Fe2O3, D_tank = 0.1, H_tank = 0.5, L_start = 0.50, P_max = 1e6, T_max = T_Fe2O3_hot_des, T_min = T_Fe2O3_hot_des - 100.0, T_set = T_Fe2O3_hot_des - 50.0, T_start = T_Fe2O3_hot_des, U_loss_tank = 0.0, eff_heater = 0.99, epsilon = 0.26) annotation(
    Placement(visible = true, transformation(origin = {-74, -45}, extent = {{-16, -15}, {16, 15}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Fe2O3_source_pump(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph) annotation(
    Placement(visible = true, transformation(origin = {-107, -57}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT H2_Topup_Source(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2, T = 25 + 273.15, nPorts = 1, p = 100000, use_T_in = false) annotation(
    Placement(visible = true, transformation(origin = {196, 138}, extent = {{12, -12}, {-12, 12}}, rotation = 0)));
  SolarTherm.Models.Fluid.HeatExchangers.H2DRI.HX_H2H2O_H2_Simple GG_HX(T_C_in_des = T_H2_mix_des, T_H_in_des = T_products_des, eff_des = eff_GGHX_des, m_flow_C_des = m_flow_H2_des, m_flow_HA_des = m_flow_H2_excess, m_flow_HB_des = m_flow_H2O_des, p_C_des = p_des, p_HA_des = p_H2_offgas_des, p_HB_des = p_H2O_offgas_des) annotation(
    Placement(visible = true, transformation(origin = {87, 71}, extent = {{-23, 23}, {23, -23}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure H2_Topup_Pump(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(origin = {162, 120}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow_H2_Topup(y = m_H2_signal - Reactor.m_flow_H2_out) annotation(
    Placement(visible = true, transformation(origin = {75, 162}, extent = {{-53, -16}, {53, 16}}, rotation = 0)));
  SolarTherm.Models.Fluid.Valves.Mixer H2_Topup_Mixer(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2, p = p_des) annotation(
    Placement(visible = true, transformation(origin = {131, 101}, extent = {{15, 15}, {-15, -15}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow_Fe2O3_cold(y = m_Fe2O3_signal) annotation(
    Placement(visible = true, transformation(origin = {191, -27}, extent = {{19, -17}, {-19, 17}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT Fe2O3_source_cold(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph, T = 25.0 + 273.15, nPorts = 1, p = 100000, use_T_in = false) annotation(
    Placement(visible = true, transformation(origin = {192, -64}, extent = {{12, -12}, {-12, 12}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure Fe2O3_pump_cold(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph) annotation(
    Placement(visible = true, transformation(origin = {145, -81}, extent = {{7, -7}, {-7, 7}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink Fe2O3_Preheated_Sink(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph) annotation(
    Placement(visible = true, transformation(origin = {-6, -90}, extent = {{22, -22}, {-22, 22}}, rotation = 0)));
  SolarTherm.Models.Storage.Tank.Reservoir H2_Reservoir(T_start = T_H2_hot_des, m_start = 10.0) annotation(
    Placement(visible = true, transformation(origin = {-104, 50}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
equation
  mflowcp_H2 = H2_Reservoir.fluid_a.m_flow*cp_H2;
  mflowcp_Fe2O3 = TES.fluid_a.m_flow*cp_Fe2O3;
  
  cp_Fe2O3 = (SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(TES.medium.T) -SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(Reactor.T_reactants))/(TES.medium.T - Reactor.T_reactants);
  cp_H2 = (Modelica.Media.IdealGases.Common.Functions.h_T(Modelica.Media.IdealGases.Common.SingleGasesData.H2,Reactor.T_reactants) - Modelica.Media.IdealGases.Common.Functions.h_T(Modelica.Media.IdealGases.Common.SingleGasesData.H2,H2_Reservoir.T))/(Reactor.T_reactants - H2_Reservoir.T);
//Vary molar ratios
//m_H2_signal = m_H2_base * 12.0 / 3.0;
//m_H2_signal = m_H2_base * (time + 3.0) / 3.0;
  r_min = min(1.05*(389.6 - 0.9702 * Reactor.T_reactants + 8.343e-4 * Reactor.T_reactants ^ 2 - 2.413e-7 * Reactor.T_reactants ^ 3), 20.0);
  m_H2_signal = r_min * m_Fe2O3_signal * SolarTherm.Models.Chemistry.ChemTable.H2.M / SolarTherm.Models.Chemistry.ChemTable.Fe2O3.M;
  m_Fe2O3_signal = m_flow_Fe2O3_des;
  T_Fe2O3_in = T_Fe2O3_hot_des;
  //T_H2_in = 700.0 + 273.15;
//Dewaterer control
  h_sat_l = Modelica.Media.Water.IF97_Utilities.hl_p(Reactor.fluid_H2O_out.p);
  T_sat = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(Reactor.fluid_H2O_out.p);
  Q_flow_cooling_signal = (Condenser.fluid_H2O_in.m_flow * (Condenser.fluid_H2O_in.h_outflow - h_sat_l) + Condenser.fluid_H2_in.m_flow * (Condenser.fluid_H2_in.h_outflow - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(Reactor.fluid_H2_out.p, T_sat))) * 1.00;
//Vary temperature
//T_Fe2O3_in = T_inlet_base + (time*100);
//T_H2_in = T_inlet_base + (time*100);
//m_Fe2O3_signal = m_Fe2O3_base;
//m_H2_signal = m_H2_base*(7.5)/3.0;
//Missing Connectors
  TES.T_amb = 25.0 + 273.15;
  TES.p = 100000;
  Q_flow_H2_req = GG_HX.m_flow_C * (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(1.0e5, T_H2_hot_des) - GG_HX.fluid_C_out.h_outflow);
//H2_Source.medium.h - inStream(H2_Topup_Sink.port_a.h_outflow));
  Q_flow_Fe2O3_req = Fe2O3_Preheated_Sink.port_a.m_flow * (Fe2O3_source.medium.h - inStream(Fe2O3_Preheated_Sink.port_a.h_outflow));
//End Missing connectors
//Connectors
  connect(T_Fe2O3_signal.y, Fe2O3_source.T_in) annotation(
    Line(points = {{-154, -53}, {-146, -53}}, color = {0, 0, 127}));
  connect(Q_flow_cooling.y, Condenser.Q_flow_removed) annotation(
    Line(points = {{181, 72}, {161, 72}, {161, 54}}, color = {0, 0, 127}));
  connect(m_flow_Fe2O3.y, Fe2O3_pump.m_flow) annotation(
    Line(points = {{-154, -9}, {-47, -9}, {-47, -49}}, color = {0, 0, 127}));
  connect(m_flow_H2.y, H2_pump.m_flow) annotation(
    Line(points = {{-156, 83}, {-56, 83}, {-56, 45}}, color = {0, 0, 127}));
  connect(Fe2O3_source.ports[1], Fe2O3_source_pump.fluid_a) annotation(
    Line(points = {{-120, -58}, {-120, -57.5}, {-114, -57.5}, {-114, -57}}, color = {0, 127, 255}));
  connect(Fe2O3_source_pump.fluid_b, TES.fluid_a) annotation(
    Line(points = {{-100, -57}, {-94, -57}, {-94, -37.5}, {-90, -37.5}}, color = {0, 127, 255}));
  connect(TES.fluid_b, Fe2O3_pump.fluid_a) annotation(
    Line(points = {{-58, -55.5}, {-54, -55.5}, {-54, -55}}, color = {0, 127, 255}));
  connect(Fe2O3_pump.fluid_b, Reactor.fluid_Fe2O3_in) annotation(
    Line(points = {{-40, -55}, {-36, -55}, {-36, -24}, {-18, -24}}, color = {0, 127, 255}));
  connect(m_flow_Fe2O3.y, Fe2O3_source_pump.m_flow) annotation(
    Line(points = {{-154, -9}, {-107, -9}, {-107, -51}}, color = {0, 0, 127}));
  connect(H2_pump.fluid_b, Reactor.fluid_H2_in) annotation(
    Line(points = {{-48, 38}, {-36, 38}, {-36, -15}, {-18, -15}}, color = {0, 170, 0}));
  connect(H2_Topup_Source.ports[1], H2_Topup_Pump.fluid_a) annotation(
    Line(points = {{184, 138}, {176, 138}, {176, 120}, {170, 120}}, color = {0, 170, 0}));
  connect(m_flow_H2_Topup.y, H2_Topup_Pump.m_flow) annotation(
    Line(points = {{133, 162}, {162, 162}, {162, 127}}, color = {0, 0, 127}));
  connect(Condenser.fluid_H2_out, H2_Topup_Mixer.fluid_A_in) annotation(
    Line(points = {{190, 50}, {214, 50}, {214, 98}, {145, 98}}, color = {0, 170, 0}));
  connect(H2_Topup_Pump.fluid_b, H2_Topup_Mixer.fluid_B_in) annotation(
    Line(points = {{154, 120}, {148, 120}, {148, 103}, {145, 103}}, color = {0, 170, 0}));
  connect(H2_Topup_Mixer.fluid_out, GG_HX.fluid_C_in) annotation(
    Line(points = {{117, 101}, {108, 101}, {108, 76}}, color = {0, 170, 0}));
  connect(p_amb.y, Reactor.p_reactor) annotation(
    Line(points = {{10, 38}, {10, 25.5}, {7, 25.5}, {7, 1}}, color = {0, 0, 127}));
  connect(Fe2O3_source_cold.ports[1], Fe2O3_pump_cold.fluid_a) annotation(
    Line(points = {{180, -64}, {167, -64}, {167, -81}, {152, -81}}, color = {0, 127, 255}));
  connect(m_flow_Fe2O3_cold.y, Fe2O3_pump_cold.m_flow) annotation(
    Line(points = {{170, -27}, {146, -27}, {146, -75}, {145, -75}}, color = {0, 0, 127}));
  connect(Condenser.fluid_H2O_gas_out, H2O_Sink.port_a) annotation(
    Line(points = {{190, 42}, {190, 34}, {194, 34}, {194, 28}}, color = {0, 127, 255}));
  connect(Condenser.fluid_H2O_liq_out, H2O_Liq_Sink.port_a) annotation(
    Line(points = {{190, 18}, {190, -4}, {196, -4}}, color = {0, 127, 255}));
  connect(GG_HX.fluid_C_out, H2_Reservoir.fluid_a) annotation(
    Line(points = {{66, 76}, {52, 76}, {52, 104}, {-136, 104}, {-136, 58}, {-120, 58}}, color = {0, 170, 0}));
  connect(H2_Reservoir.fluid_b, H2_pump.fluid_a) annotation(
    Line(points = {{-88, 38}, {-64, 38}, {-64, 38}, {-64, 38}}, color = {0, 127, 255}));
  connect(Fe2O3_pump_cold.fluid_b, Fe2O3_Preheated_Sink.port_a) annotation(
    Line(points = {{138, -80}, {70, -80}, {70, -90}, {16, -90}}, color = {0, 127, 255}));
  connect(GG_HX.fluid_HB_out, Condenser.fluid_H2O_in) annotation(
    Line(points = {{108, 64}, {128, 64}, {128, 42}, {132, 42}}, color = {0, 127, 255}));
  connect(GG_HX.fluid_HA_out, Condenser.fluid_H2_in) annotation(
    Line(points = {{108, 68}, {132, 68}, {132, 50}}, color = {0, 170, 0}));
  connect(GG_HX.fluid_HA_in, Reactor.fluid_H2_out) annotation(
    Line(points = {{66, 68}, {54, 68}, {54, -15}, {32, -15}}, color = {0, 170, 0}));
  connect(GG_HX.fluid_HB_in, Reactor.fluid_H2O_out) annotation(
    Line(points = {{66, 64}, {58, 64}, {58, -24}, {32, -24}}, color = {0, 127, 255}));
  annotation(
    experiment(StopTime = 1000.0, StartTime = 0.0, Tolerance = 1e-3, Interval = 1.0),
    Diagram(coordinateSystem(extent = {{-200, -100}, {220, 200}}, preserveAspectRatio = false, initialScale = 0.1), graphics = {Text(origin = {-132, 173}, extent = {{-50, 17}, {144, 5}}, textString = "Case 1: Hot Fe2O3 TES sustains the process", fontSize = 12), Text(origin = {-120, 163}, extent = {{-36, 19}, {154, 1}}, textString = "H2 is only heated using gas-gas HX.", fontSize = 12), Line(origin = {-102, -92}, points = {{72, 0}, {-72, 0}}, pattern = LinePattern.Dash), Line(origin = {-174, -78}, points = {{0, -14}, {0, 14}}, pattern = LinePattern.Dash, arrow = {Arrow.None, Arrow.Filled}), Text(origin = {-102, -95}, extent = {{-56, 3}, {42, -5}}, textString = "Heated to target temperature"), Text(origin = {-62, -21}, extent = {{-38, 7}, {12, -5}}, textString = "T_Fe2O3_hot_des"), Text(origin = {-96, 73}, extent = {{-28, 9}, {12, -5}}, textString = "T_H2_hot_des"), Text(origin = {56, -33}, extent = {{-20, 7}, {42, -19}}, textString = "T_products_des", textStyle = {TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold, TextStyle.Bold}), Text(origin = {84, 45}, extent = {{-20, 7}, {30, -13}}, textString = "eff_GGHX_des", textStyle = {TextStyle.Bold, TextStyle.Bold})}),
    Icon(coordinateSystem(extent = {{-200, -100}, {220, 200}}, preserveAspectRatio = false)));
end H2DRI_DesignCase_1a;