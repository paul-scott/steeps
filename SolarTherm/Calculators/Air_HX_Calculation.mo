within SolarTherm.Calculators;

model Air_HX_Calculation
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import MA = Modelica.Blocks.Math;

  parameter SI.MassFlowRate m_flow_C = 26.6 "Mass flow rate of cold stream (fresh air) (kg/s)";
  parameter SI.MassFlowRate m_flow_H = 26.6 "Mass flow rate of hot stream (dusty return air) (kg/s)";
  parameter SI.Temperature T_C_in = 25.0 + 273.15 "Inlet temperature of cold stream (fresh air) (K)";
  parameter SI.Temperature T_H_in = 388.0 + 273.15 "Inlet temperature of hot stream (dusty return air) (K)";
  Real eff_HX(start=0.80) "Heat exchanger effectiveness (-)";
  //SI.Temperature T_C_out "Outlet temperature of heated fresh air stream (K)";
  //SI.Temperature T_H_out "Outlet temperature of cooled dusty return air stream (K)";
  
  parameter SI.Pressure p_des = 100000.0 "Design pressure of both air streams (Pa)";
  
  parameter SI.SpecificHeatCapacity cp_H = (Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(p_des,T_H_in) - Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(p_des,T_C_in))/(T_H_in-T_C_in);
  parameter SI.SpecificHeatCapacity cp_C = (Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(p_des,T_H_in) - Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(p_des,T_C_in))/(T_H_in-T_C_in);
  //parameter SI.SpecificHeatCapacity cp_C 
  
  
  parameter SI.ThermalConductance C_H = m_flow_H*cp_H "Heat capacity rate of dusty return air (W/K)";
  parameter SI.ThermalConductance C_C = m_flow_C*cp_C "Heat capacity rate of fresh air (W/K)";
  
  parameter SI.ThermalConductance C_min = min(C_H,C_C) "Minimum heat capacity rate of the HX streams (W/K)";
  
  parameter SI.CoefficientOfHeatTransfer U_HX = 30.0 "Overall heat exchanger heat transfer coefficient (W/m2K)";
  SI.Area A_HX "Heat transfer area required by HX (m2)";
  parameter SI.HeatFlowRate Q_flow_max = C_min*(T_H_in-T_C_in) "Maximum possible heat transfer rate between the two fluid streams (W)";
  SI.HeatFlowRate Q_flow_real = eff_HX*Q_flow_max "Actual heat transfer rate between the two fluid streams (W)";
  
  parameter SI.SpecificEnthalpy h_H_in = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(p_des,T_H_in) "Specific Enthalpy of hot inlet dusty return air (J/kg)";
  parameter SI.SpecificEnthalpy h_C_in = Modelica.Media.Air.ReferenceAir.Air_Utilities.h_pT(p_des,T_C_in) "Specific Enthalpy of cold inlet fresh air (J/kg)";
  SI.SpecificEnthalpy h_H_out "Specific Enthalpy of cooled dusty return air (J/kg)";
  SI.SpecificEnthalpy h_C_out "Specific Enthalpy of heated fresh air (J/kg)";
  
  SI.Temperature T_H_out = Modelica.Media.Air.ReferenceAir.Air_Utilities.T_ph(p_des,h_H_out);
  SI.Temperature T_C_out = Modelica.Media.Air.ReferenceAir.Air_Utilities.T_ph(p_des,h_C_out);
equation
  Q_flow_real = m_flow_H*(h_H_in-h_H_out);
  Q_flow_real = m_flow_C*(h_C_out-h_C_in);
  A_HX = Q_flow_real/(U_HX*(T_H_in-T_C_in));
  der(eff_HX) = 0.0;//No change in HX eff for the moment


annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Air_HX_Calculation;