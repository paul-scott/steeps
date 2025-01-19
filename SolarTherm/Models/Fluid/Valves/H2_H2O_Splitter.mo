within SolarTherm.Models.Fluid.Valves;

model H2_H2O_Splitter
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  
  Modelica.Fluid.Interfaces.FluidPort_a fluid_A_in(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-94, 18}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_B_in(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-94, -16}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_A1_out(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {4, 94}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_A2_out(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {4, -94}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_B1_out(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {38, 94}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_B2_out(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {36, -94}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  parameter SI.Pressure p = 100000.0;
  Modelica.Blocks.Interfaces.RealInput frac_1 annotation(
		Placement(visible = true, 
		transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), 
		iconTransformation(origin = {49, 5}, extent = {{-7, -7}, {7, 7}}, rotation = 180)));
		
  SI.MolarFlowRate n_flow_H2_in "Inlet molar flow rate of Hydrogen gas (mol/s)";
  SI.MolarFlowRate n_flow_H2O_in "Inlet molar flow rate of Water vapour (mol/s)";
  SI.Pressure p_H2 "Partial pressure of H2 (Pa)";
  SI.Pressure p_H2O "Partial pressure of H2O vapour (Pa)";
equation
  n_flow_H2_in = fluid_A_in.m_flow/2.01588e-3;
  n_flow_H2O_in = fluid_B_in.m_flow/18.01528e-3;
  
  p_H2 = (p*n_flow_H2_in)/(n_flow_H2O_in+n_flow_H2_in);
  p_H2O = (p*n_flow_H2O_in)/(n_flow_H2O_in+n_flow_H2_in);
  
  
  fluid_A1_out.m_flow = -1.0*frac_1*fluid_A_in.m_flow;
  fluid_A2_out.m_flow = -1.0*(1.0-frac_1)*fluid_A_in.m_flow;
  fluid_B1_out.m_flow = -1.0*frac_1*fluid_B_in.m_flow;
  fluid_B2_out.m_flow = -1.0*(1.0-frac_1)*fluid_B_in.m_flow;
  
  fluid_A_in.h_outflow = inStream(fluid_A_in.h_outflow);
  fluid_B_in.h_outflow = inStream(fluid_B_in.h_outflow);
  fluid_A1_out.h_outflow = inStream(fluid_A_in.h_outflow);
  fluid_A2_out.h_outflow = inStream(fluid_A_in.h_outflow);
  fluid_B1_out.h_outflow = inStream(fluid_B_in.h_outflow);
  fluid_B2_out.h_outflow = inStream(fluid_B_in.h_outflow);
  
  fluid_A1_out.p = p_H2;
  fluid_A2_out.p = p_H2;
  fluid_B1_out.p = p_H2O;
  fluid_B2_out.p = p_H2O;
  
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Rectangle(origin = {-2, 1}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Sphere, extent = {{-68, 21}, {42, -19}}), Line(origin = {-85, 18}, points = {{15, 0}, {-15, 0}}), Line(origin = {-85, -15.7969}, points = {{15, 0}, {-15, 0}}), Text(origin = {-82, 32}, extent = {{-10, 10}, {6, -8}}, textString = "H2"), Text(origin = {-84, -32}, extent = {{-10, 10}, {12, -12}}, textString = "H2O"), Text(origin = {60, 32}, extent = {{-10, 10}, {30, -28}}, textString = "frac_1"), Line(origin = {4, -80}, points = {{0, 10}, {0, -8}}), Line(origin = {36.2694, -80.461}, points = {{0, 10}, {0, -8}}), Line(origin = {4, 77.5956}, points = {{0, 10}, {0, -8}}), Line(origin = {37.6524, 78.2126}, points = {{0, 10}, {0, -8}}), Text(origin = {-18, 92}, extent = {{-10, 10}, {6, -8}}, textString = "H2"), Text(origin = {-16, -94}, extent = {{-10, 10}, {6, -8}}, textString = "H2"), Text(origin = {60, -92}, extent = {{-10, 10}, {12, -12}}, textString = "H2O"), Text(origin = {58, 94}, extent = {{-10, 10}, {14, -12}}, textString = "H2O"), Text(origin = {88, 106}, extent = {{-10, 10}, {42, -36}}, textString = "flow 1"), Text(origin = {88, -84}, extent = {{-10, 10}, {46, -26}}, textString = "flow 2"), Rectangle(origin = {66, 13}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Sphere, extent = {{-66, 57}, {-26, -83}}), Line(origin = {-29.05, 1.92}, points = {{-23, 0}, {13, 0}}, color = {255, 255, 255}, thickness = 1, arrow = {Arrow.None, Arrow.Open}, arrowSize = 5), Line(origin = {104.91, 34.1629}, points = {{-23, 0}, {-23, 34}}, color = {255, 255, 255}, thickness = 1, arrow = {Arrow.None, Arrow.Open}, arrowSize = 5), Line(origin = {42.2146, 26.787}, points = {{-23, 0}, {-23, 34}}, color = {255, 255, 255}, thickness = 1, arrow = {Arrow.None, Arrow.Open}, arrowSize = 5), Line(origin = {42.3476, -26.5021}, points = {{-23, 0}, {-23, -32}}, color = {255, 255, 255}, thickness = 1, arrow = {Arrow.None, Arrow.Open}, arrowSize = 5)}, coordinateSystem(initialScale = 0.1)));
end H2_H2O_Splitter;