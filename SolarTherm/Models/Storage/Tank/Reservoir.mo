within SolarTherm.Models.Storage.Tank;

model Reservoir
  extends Interfaces.Models.StorageFluid;
  redeclare replaceable package Medium = Modelica.Media.IdealGases.SingleGases.H2;
  parameter SI.Temperature T_start=from_degC(700) "Starting temperature of fluid (kg)";
  parameter SI.Mass m_start = 10.0 "Starting mass of fluid (kg)";
  parameter SI.Pressure p_start = 100000.0 "Pressure of fluid (Pa)";
  SI.Mass m(start = m_start) "Mass of fluid (kg)";
  Medium.BaseProperties medium;
  SI.SpecificHeatCapacity h(start=Medium.specificEnthalpy_pT(p_start,T_start)) "Specific Enthalpy of fluid (J/kg)";
  SI.Temperature T "Temperature of fluid (K)";
  SI.HeatFlowRate der_m_h "Intermediate calculation variable (J/s)";
  SI.SpecificEnthalpy h_in;
algorithm
  der_m_h := fluid_a.m_flow*h_in + fluid_b.m_flow*h;
equation
//Pressure
  fluid_a.p = fluid_b.p;
  
  //Energy and Mass Balance
  der(m) = fluid_a.m_flow + fluid_b.m_flow;
  der(h) = (der_m_h - h*der(m))/m;
  
  //fluid_a.m_flow + fluid_b.m_flow = 0.0;
  //State
  medium.h = h;
  medium.p = p_start;
  T = medium.T;
  
  h_in = inStream(fluid_a.h_outflow);
  fluid_a.h_outflow = inStream(fluid_a.h_outflow);
  fluid_b.h_outflow = h;
  //der(fluid_state.h) = (fluid_a.m_flow*inStream(fluid_a.h_outflow)+fluid_b.m_flow*fluid_state.h)/m;
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Text(origin = {12, 0}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {125, -140}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
end Reservoir;