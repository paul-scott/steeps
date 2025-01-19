within SolarTherm.Models.Fluid.HeatExchangers;

model Boiler_Basic//Simplify Boiler
import SI = Modelica.SIunits;
import Modelica.SIunits.Conversions.*;
  replaceable package Medium = SolarTherm.Media.Air.Air_amb_p;
  //parameter Modelica.SIunits.HeatFlowRate Q_flow_ref=50e6 "Design thermal power";
  parameter Modelica.SIunits.Temperature T_cold_set = from_degC(290);
  parameter Modelica.SIunits.Temperature T_hot_set = from_degC(565);
  final parameter Medium.ThermodynamicState state_cold_set = Medium.setState_pTX(Medium.p_default, T_cold_set) "Cold fluid thermodynamic state at design";
  final parameter Medium.ThermodynamicState state_hot_set = Medium.setState_pTX(Medium.p_default, T_hot_set) "Cold fluid thermodynamic state at design";
  final parameter Modelica.SIunits.SpecificEnthalpy h_cold_set = Medium.specificEnthalpy(state_cold_set) "Cold fluid specific enthalpy at design";  
  final parameter Modelica.SIunits.SpecificEnthalpy h_hot_set = Medium.specificEnthalpy(state_hot_set) "Cold fluid specific enthalpy at design";
  //final parameter SI.MassFlowRate m_flow_ref = Q_flow_ref / (h_hot_set - h_cold_set) "Mass flow rate at design";
  //final parameter Real nu_eps=1e-3 "Minimum load";
  //parameter Real nu_min=1e-3 "Minimum turbine operation";

  Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = Medium)
    "Fluid connector a (positive design flow direction is from port_a to port_b)"
    annotation (Placement(visible = true,transformation(extent={{-110,-12},{-90,8}},  rotation=
           0), iconTransformation(extent = {{-10, 90}, {10, 110}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium =
        Medium)
    "Fluid connector b (positive design flow direction is from port_a to port_b)"
    annotation (Placement(visible = true,transformation(extent={{110,-10},{90,10}}, rotation=0), iconTransformation(extent = {{10, -112}, {-10, -92}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.RealOutput h_out_signal "outlet specific enthalpy" annotation (Placement(visible = true, transformation(origin = {-80, -72}, extent = {{20, -20}, {-20, 20}}, rotation = 0), iconTransformation(origin = {-92, -72}, extent = {{8, -8}, {-8, 8}}, rotation = 90))) ;

  //Real load;
  //Boolean logic;
  Medium.BaseProperties state_in;
  SI.HeatFlowRate Q_flow;
  SI.Energy E_thermal(start=0.0);
  SI.SpecificEnthalpy h_in;
  SI.SpecificEnthalpy h_out;

equation
  state_in.h = h_in;

  state_in.p = 1e5;
  //load=max(nu_eps,fluid_a.m_flow/m_flow_ref);
  //logic=load>nu_min;
  h_in=inStream(fluid_a.h_outflow);
  h_out=h_cold_set;
  h_out_signal=h_cold_set;
  
  if fluid_a.m_flow > 1.0e-6 then
    Q_flow = (h_in-h_out)*fluid_a.m_flow;
  else
    Q_flow = 0.0;
  end if;
  
  fluid_b.h_outflow = h_cold_set;
  fluid_a.h_outflow = h_in;
  
  fluid_a.m_flow + fluid_b.m_flow = 0.0;
  fluid_a.p=fluid_b.p;
  
  der(E_thermal)=Q_flow;
  annotation(
    Icon(graphics = {Text(origin = {0, -6}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {151, -154}}, textString = "%name"), Rectangle(fillColor = {135, 135, 135}, fillPattern = FillPattern.Backward, extent = {{-100, 64}, {100, -64}}), Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, -44}, {100, 44}}), Rectangle(fillColor = {0, 127, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, -44}, {100, 44}}), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{18, 0}, {48, -29}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-1, 29}, {29, 0}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{43, 31}, {73, 2}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-31, 1}, {-1, -28}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{50, 15}, {80, -14}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{-72, 25}, {-42, -4}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{71, -11}, {101, -40}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{72, 28}, {102, -1}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{71, 40}, {101, 11}}, endAngle = 360), Line(points = {{0, -64}, {0, -100}}, color = {0, 127, 255}), Line(origin = {-40, 0}, points = {{40, 100}, {40, 64}}, color = {0, 127, 255}), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{58, -11}, {88, -40}}, endAngle = 360), Ellipse(lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, pattern = LinePattern.None, fillPattern = FillPattern.Solid, extent = {{71, 1}, {101, -28}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)));
end Boiler_Basic;