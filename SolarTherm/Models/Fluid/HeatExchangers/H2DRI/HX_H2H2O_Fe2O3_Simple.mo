within SolarTherm.Models.Fluid.HeatExchangers.H2DRI;

model HX_H2H2O_Fe2O3_Simple
  //In this model, HA refers to hot H2 inlet, HB is hot H2O inlet, C is cold H2 inlet.
  //Cost model is based on an approximation by Hewitt & Pugh (2007) for a standard double-pipe HX. Effectiveness-NTU relationship is for a counterflow arrangement.
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  //Ports
  Modelica.Fluid.Interfaces.FluidPort_a fluid_HA_in annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-100, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_HB_in annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-100, 14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_HA_out annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {100, 34}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_HB_out annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {100, 14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_C_in annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {100, -22}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_C_out annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {-100, -24}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  //Design Parameters
  parameter SI.Pressure p_HA_des = 0.75 * 100000.0 "Partial pressure of hot H2 gas (Pa)";
  parameter SI.Pressure p_HB_des = 0.25 * 100000.0 "Partial pressure of hot H2O gas (Pa)";
  parameter SI.Pressure p_C_des = 1.00 * 100000.0 "Partial pressure of cold Fe2O3 (Pa)";
  parameter SI.Temperature T_H_in_des = 650.0 + 273.15 "HX Inlet Hot Temperature (K)";
  parameter SI.Temperature T_C_in_des = 100.0 + 273.15 "HX Inlet Cold Temperature (K)";
  parameter SI.MassFlowRate m_flow_HA_des = 0.5*0.162 "Mass flow rate of hot H2 gas (kg/s)";
  parameter SI.MassFlowRate m_flow_HB_des = 0.5*0.484 "Mass flow rate of hot H2O gas (kg/s)";
  parameter SI.MassFlowRate m_flow_C_des = 1.43 "Mass flow rate of cold Fe2O3 (kg/s)";
  parameter Real eff_des = 0.8 "Design Effectiveness of HX";
  parameter SI.SpecificHeatCapacity cp_HA_des = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_HA_des, T_H_in_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_HA_des, T_C_in_des)) / (T_H_in_des - T_C_in_des) "Approximate cp of hot hydrogen stream (J/kgK)";
  //parameter SI.SpecificHeatCapacity cp_HB_des = (Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_HB_des, T_H_in_des) - Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_HB_des, T_C_in_des)) / (T_H_in_des - T_C_in_des) "Approximate cp of hot water vapour stream (J/kgK)";
  parameter SI.SpecificHeatCapacity cp_HB_des = 0.5*(Modelica.Media.Water.IF97_Utilities.cp_pT(p_HB_des, T_H_in_des) + Modelica.Media.Water.IF97_Utilities.cp_pT(p_HB_des, T_C_in_des)) "Approximate cp of hot water vapour stream (J/kgK)";
  parameter SI.SpecificHeatCapacity cp_C_des = (SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_H_in_des) - SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_C_in_des)) / (T_H_in_des - T_C_in_des) "Approximate cp of cold Fe2O3 stream (J/kgK)";
  parameter SI.ThermalConductance C_H_des = m_flow_HA_des * cp_HA_des + m_flow_HB_des * cp_HB_des "Heat Capacity Rate of the hot stream";
  parameter SI.ThermalConductance C_C_des = m_flow_C_des * cp_C_des "Heat Capacity Rate of the cold stream (W/K)";
  parameter SI.ThermalConductance C_min_des = min(C_H_des, C_C_des);
  parameter SI.ThermalConductance C_max_des = max(C_H_des, C_C_des);
  parameter Real C_r_des = C_min_des / C_max_des;
  //parameter SI.ThermalConductance UA = C_min / C_r * log((eff - 1.0) / (eff * C_r - 1.0)) "Assuming counterflow HX (J/K)";
  parameter Real NTU = 1.0 / (C_r_des - 1.0) * log((eff_des - 1.0) / (eff_des * C_r_des - 1.0)) "Assuming counterflow HX";
  parameter SI.ThermalConductance UA_des = C_min_des * NTU;
  //HX Variables
  SI.ThermalConductance C_H "Instantaneous Heat Capacity Rate of the hot stream (W/K)";
  SI.ThermalConductance C_C "Instantaneous Heat Capacity Rate of the hot stream (W/K)";
  SI.ThermalConductance C_min "Instantaneous Cmin (W/K)";
  SI.ThermalConductance C_max "Instantaneous Cmax (W/K)";
  Real C_r "Instantaneous Cr";
  Real eff "Instantaneous effectiveness";
  //Single Variables
  SI.MassFlowRate m_flow_HA "Mass flow rate of hot hydrogen (kg/s)";
  SI.MassFlowRate m_flow_HB "Mass flow rate of hot water vapour (kg/s)";
  SI.MassFlowRate m_flow_C "Mass flow rate of cold hydrogen (kg/s)";
  SI.Pressure p_HA "Pressure of hot hydrogen (Pa)";
  SI.Pressure p_HB "Pressure of hot water vapour (Pa)";
  SI.Pressure p_C "Pressure of cold hydrogen (Pa)";
  SI.HeatFlowRate Q_flow "Total heat flow rate from hot stream to cold stream (J/s)";
  //Temperatures and enthalpies
  SI.SpecificEnthalpy h_HA_in "Specific Enthalpy of hot inlet hydrogen stream (J/kg)";
  SI.SpecificEnthalpy h_HB_in "Specific Enthalpy of hot inlet water stream (J/kg)";
  SI.SpecificEnthalpy h_C_in "Specific Enthalpy of cold inlet hydrogen stream (J/kg)";
  SI.SpecificEnthalpy h_HA_out "Specific Enthalpy of hot outlet hydrogen stream (J/kg)";
  SI.SpecificEnthalpy h_HB_out "Specific Enthalpy of hot outlet water stream (J/kg)";
  SI.SpecificEnthalpy h_C_out "Specific Enthalpy of cold outlet hydrogen stream (J/kg)";
  SI.Temperature T_H_in "Temperature of hot inlet gas mixture (K)";
  SI.Temperature T_H_out "Temperature of hot outlet gas mixture (K)";
  SI.Temperature T_C_in "Temperature of cold inlet H2 (K)";
  SI.Temperature T_C_out "Temperature of cold outlet H2 (K)";
  //Component cost
  //Based on Seider double-pipe, SS/SS, USD_2022
  parameter SI.CoefficientOfHeatTransfer U_approx = 10.0 "Approximate overall (W.m2K)";
  parameter SI.Area A_approx = UA_des / U_approx "Approximate HX area (m2)";
  parameter Real C_component = div(A_approx, 18.58) * 14503.48 + 4.896 * exp(7.146 + 0.16 * log(10.7639 * rem(A_approx, 18.58))) "Total component cost USD_2022";
  //(if A_approx < 1115.0 then 4.896*exp(7.146+0.16*log(10.7639*A_approx)) else
algorithm
  if m_flow_HA + m_flow_HB > 1.0e-9 and m_flow_C > 1.0e-9 then
    C_H := m_flow_HA * cp_HA_des + m_flow_HB * cp_HB_des;
    C_C := m_flow_C * cp_C_des;
    C_min := min(C_H, C_C);
    C_max := max(C_H, C_C);
    C_r := C_min / C_max;
    eff := (1.0 - exp(-1.0 * NTU * (1.0 - C_r))) / (1.0 - C_r * exp(-1.0 * NTU * (1.0 - C_r)));
  else
    C_H := 0.0;
    C_C := 0.0;
    C_min := 0.0;
    C_max := 0.0;
    C_r := 0.0;
    eff := 0.0;
  end if;
//flowing
  Q_flow := eff * C_min * (T_H_in - T_C_in);
//for i in 1:N loop
//T_H_avg[i] := Modelica.Media.IdealGases.SingleGases.H2.T_h(0.5 * (h_HA[i] + h_HA[i + 1]));
//T_C_avg[i] := Modelica.Media.IdealGases.SingleGases.H2.T_h(0.5 * (h_C[i] + h_C[i + 1]));
//Q_flow[i] := UA / N * (T_H_avg[i] - T_C_avg[i]);
//end for;
/*
  Q_flow[1] := (UA/N)*(T_H[1]-T_C[2]);
  for i in 2:N-1 loop
    Q_flow[i] := (UA/N)*(T_H[i-1]-T_C[i+1]);
  end for;
  Q_flow[N] := (UA/N)*(T_H[N]-T_C[N+1]);
  */

equation
//Ports
  m_flow_HA = fluid_HA_in.m_flow;
  m_flow_HB = fluid_HB_in.m_flow;
  m_flow_C = fluid_C_in.m_flow;
  fluid_HA_out.m_flow = -1.0 * m_flow_HA;
  fluid_HB_out.m_flow = -1.0 * m_flow_HB;
  fluid_C_out.m_flow = -1.0 * m_flow_C;
  p_HA = fluid_HA_in.p;
  p_HB = fluid_HB_in.p;
  p_C = fluid_C_in.p;
  fluid_HA_out.p = p_HA;
  fluid_HB_out.p = p_HB;
  fluid_C_out.p = p_C;
//Dont forget the outlet enthalpies
  fluid_HA_in.h_outflow = inStream(fluid_HA_in.h_outflow);
  fluid_HB_in.h_outflow = inStream(fluid_HB_in.h_outflow);
  fluid_C_in.h_outflow = inStream(fluid_C_in.h_outflow);
  fluid_HA_out.h_outflow = h_HA_out;
  fluid_HB_out.h_outflow = h_HB_out;
  fluid_C_out.h_outflow = h_C_out;
  T_H_in = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_HA, h_HA_in);
  T_C_in = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_C_in); 
  Q_flow = m_flow_HA * (h_HA_in - h_HA_out) + m_flow_HB * (h_HB_in - h_HB_out);
  Q_flow = m_flow_C * (h_C_out - h_C_in);
  h_HA_in = inStream(fluid_HA_in.h_outflow);
  h_HB_in = inStream(fluid_HB_in.h_outflow);
  h_C_in = inStream(fluid_C_in.h_outflow);
  h_HA_out = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_HA, T_H_out);
  h_HB_out = Modelica.Media.Water.IF97_Utilities.h_pT(p_HB, T_H_out);
  T_C_out = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_C_out);
//T_H_out = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_HA,h_HA_out);
//Calculate instantaneous variables
/*
  h_HA[1] = inStream(fluid_HA_in.h_outflow);
  h_HB[1] = inStream(fluid_HB_in.h_outflow);
  
  
//Left-most node (i=1)
//--Properties
//----Port Flows
  h_HA[1] = inStream(fluid_HA_in.h_outflow);
  h_HB[1] = inStream(fluid_HB_in.h_outflow);
  fluid_C_out.h_outflow = h_C[1];
//----Properties
  h_HA[1] = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_HA, T_H[1]);
  h_HB[1] = Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p_HB, T_H[1]);
//h_C[1] = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_C,T_C[1]);
//Mid nodes (i=2 to i=N)
//--Properties
  for i in 2:N loop
    h_HA[i] = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_HA, T_H[i]);
    h_HB[i] = Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p_HB, T_H[i]);
//h_C[i] = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_C,T_C[i]);
  end for;
//Right-most node (i=N+1)
//--Properties
  h_HA[N + 1] = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_HA, T_H[N + 1]);
  h_HB[N + 1] = Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p_HB, T_H[N + 1]);
//h_C[N+1] = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_C,T_C[N+1]);
//----Port Flows
  fluid_HA_out.h_outflow = h_HA[N + 1];
  fluid_HB_out.h_outflow = h_HB[N + 1];
  h_C[N + 1] = inStream(fluid_C_in.h_outflow);
//For each element
  for i in 1:N loop
//--Average temperature across each HX section
//T_H_plus[i] = 0.5*(T_H[i] + T_H[i+1]);
//T_C_plus[i] = 0.5*(T_C[i] + T_C[i+1]);
//--Energy Flow Balance across each HX section
    Q_flow[i] = m_flow_HA * (h_HA[i] - h_HA[i + 1]) + m_flow_HB * (h_HB[i] - h_HB[i + 1]);
    Q_flow[i] = m_flow_C * (h_C[i] - h_C[i + 1]);
//Q_flow[i] = (UA/N)*(T_H_avg[i]-T_C_avg[i]);
  end for;
  */
  annotation(
    Diagram(coordinateSystem(extent = {{-110, -100}, {110, 100}})),
    Icon(graphics = {Text(origin = {22, 204}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {103, -138}}, textString = "%name"), Rectangle(origin = {0, 23}, fillColor = {170, 0, 0}, fillPattern = FillPattern.Sphere, extent = {{-68, 23}, {70, -23}}), Rectangle(origin = {0, -23}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Sphere, extent = {{-68, 23}, {70, -23}}), Line(origin = {-83.539, 34.461}, points = {{-16, 0}, {16, 0}}), Line(origin = {-84, 13.539}, points = {{-16, 0.460992}, {16, 0.460992}}), Line(origin = {-84, -24}, points = {{-16, 0}, {16, 0}}), Line(origin = {85, 34}, points = {{-15, 0}, {15, 0}}), Line(origin = {85, 14}, points = {{-15, 0}, {15, 0}}), Line(origin = {85, -22}, points = {{-15, 0}, {15, 0}}), Line(origin = {3, 24}, points = {{-37, -2}, {33, -2}}, color = {255, 255, 255}, thickness = 1, arrow = {Arrow.None, Arrow.Open}, arrowSize = 10), Line(origin = {72.15, -21.18}, points = {{-37, -2}, {-107, -2}}, color = {255, 255, 255}, thickness = 1, arrow = {Arrow.None, Arrow.Open}, arrowSize = 10), Text(origin = {-79, 24}, extent = {{-17, 28}, {9, 14}}, textString = "H2"), Text(origin = {89, 24}, extent = {{-17, 28}, {9, 14}}, textString = "H2"), Text(origin = {-79, -16}, extent = {{-17, 28}, {9, 14}}, textString = "H2O"), Text(origin = {89, -16}, extent = {{-17, 28}, {9, 14}}, textString = "H2O"), Text(origin = {-87, -56}, extent = {{-17, 28}, {17, 4}}, textString = "Fe2O3"), Text(origin = {89, -54}, extent = {{-17, 28}, {17, 4}}, textString = "Fe2O3")}, coordinateSystem(extent = {{-110, -100}, {110, 100}}, initialScale = 0.1)));
end HX_H2H2O_Fe2O3_Simple;