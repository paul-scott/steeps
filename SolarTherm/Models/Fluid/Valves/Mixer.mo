within SolarTherm.Models.Fluid.Valves;

model Mixer
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  replaceable package Medium = Modelica.Media.IdealGases.SingleGases.H2
    constrainedby Modelica.Media.Interfaces.PartialMedium
    "Medium in the component"
      annotation (choicesAllMatching = true,Dialog(group="Working fluid"));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_A_in(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-94, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_B_in(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-94, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_out(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {94, 2}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  parameter SI.Pressure p = 100000.0;
equation
  fluid_out.p = p;
  fluid_A_in.m_flow + fluid_B_in.m_flow + fluid_out.m_flow = 0.0;
  fluid_A_in.h_outflow = inStream(fluid_A_in.h_outflow);
  fluid_B_in.h_outflow = inStream(fluid_B_in.h_outflow);
  if (fluid_A_in.m_flow+fluid_B_in.m_flow) > 1.0e-8 then
    fluid_out.h_outflow = (fluid_A_in.m_flow*inStream(fluid_A_in.h_outflow) + fluid_B_in.m_flow*inStream(fluid_B_in.h_outflow))/(fluid_A_in.m_flow+fluid_B_in.m_flow);
  else
    fluid_out.h_outflow = 0.5*(inStream(fluid_A_in.h_outflow)+inStream(fluid_B_in.h_outflow));
  end if;
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Rectangle(origin = {-2, 3}, fillColor = {0, 0, 127}, fillPattern = FillPattern.Sphere, extent = {{-68, 31}, {70, -31}}), Text(origin = {-4, 220}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {151, -154}}, textString = "%name"), Line(origin = {-85, 20}, points = {{15, 0}, {-15, 0}}), Line(origin = {-0.67, 4.02}, points = {{-45, 0}, {45, 0}}, color = {255, 255, 255}, thickness = 2, arrow = {Arrow.None, Arrow.Open}, arrowSize = 10), Line(origin = {-85, -13.7969}, points = {{15, 0}, {-15, 0}}), Text(origin = {-78, 28}, extent = {{-10, 10}, {6, -8}}, textString = "A"), Text(origin = {-78, -26}, extent = {{-10, 10}, {6, -8}}, textString = "B"), Line(origin = {82, 2}, points = {{-14, 0}, {14, 0}})}, coordinateSystem(initialScale = 0.1)));
end Mixer;