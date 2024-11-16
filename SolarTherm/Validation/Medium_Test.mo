within SolarTherm.Validation;

model Medium_Test
  replaceable package Medium = SolarTherm.Media.SolidParticles.Fe2O3H2O_ph;
  replaceable package Utilities = SolarTherm.Media.SolidParticles.Fe2O3H2O_utilities;
  Modelica.SIunits.Temperature T_signal(start=298.15);
  Modelica.SIunits.MassFlowRate m_flow_signal = 1.0;
  
  Modelica.SIunits.SpecificEnthalpy h_fluid = Utilities.h_T(T_signal);
  Modelica.SIunits.SpecificHeatCapacityAtConstantPressure cp_fluid = Utilities.cp_T(T_signal);
  Modelica.SIunits.Density rho_fluid = Utilities.rho_T(T_signal);
  Modelica.SIunits.ThermalConductivity k_fluid = Utilities.k_T(T_signal);
  Modelica.SIunits.SpecificEntropy s_fluid = Utilities.s_T(T_signal);

  Medium.BaseProperties State;
  
  Modelica.SIunits.SpecificEnthalpy h_state = State.h;
  Modelica.SIunits.Density rho_state = State.d;

  
  
  Modelica.Fluid.Sources.Boundary_pT Medium_source(redeclare package Medium = Medium, nPorts = 1, p = 101325, use_T_in = true, use_p_in = false) annotation(
    Placement(visible = true, transformation(origin = {-19, 7}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink Sink(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {82, 2}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression T_input(y = T_signal) annotation(
    Placement(visible = true, transformation(origin = {-81, 11}, extent = {{-19, -17}, {19, 17}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Pump(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {28, 4}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow(y = m_flow_signal) annotation(
    Placement(visible = true, transformation(origin = {-5, 55}, extent = {{-19, -17}, {19, 17}}, rotation = 0)));
equation
  State.p = 101325.0;
  State.T = T_signal;
  der(T_signal) = 1.0;
  connect(T_input.y, Medium_source.T_in) annotation(
    Line(points = {{-60, 11}, {-18, 11}, {-18, 14}, {-39, 14}}, color = {0, 0, 127}));
  connect(m_flow.y, Pump.m_flow) annotation(
    Line(points = {{16, 56}, {28, 56}, {28, 11}}, color = {0, 0, 127}));
  connect(Medium_source.ports[1], Pump.fluid_a) annotation(
    Line(points = {{-2, 8}, {8, 8}, {8, 4}, {20, 4}, {20, 4}}, color = {0, 127, 255}));
  connect(Pump.fluid_b, Sink.port_a) annotation(
    Line(points = {{36, 4}, {44, 4}, {44, 0}, {56, 0}, {56, 2}}, color = {0, 127, 255}));

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)), experiment(StopTime = 1200, StartTime = 0, Tolerance = 1.0e-5, Interval = 10, maxStepSize = 10, initialStepSize = 10));
end Medium_Test;