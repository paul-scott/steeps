within SolarTherm.Models.Fluid.HeatExchangers;

model mass_loop_breaker
  replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
  Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium=Medium) annotation (
  Placement(visible = true,transformation(extent={{-110, -10}, {-90, 10}},  rotation=0), iconTransformation(origin = {-60, 0}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium=Medium) annotation (
  Placement(visible = true,transformation(extent={{88, -10}, {110, 10}},  rotation=0), iconTransformation(origin = {60, -2.22045e-16}, extent = {{4, -4}, {-4, 4}}, rotation = 0)));

equation
  port_a.h_outflow = inStream(port_b.h_outflow); // shouldn't flow backwards anyway
  port_b.p = port_a.p;
  port_b.h_outflow= inStream(port_a.h_outflow);
  port_b.Xi_outflow = inStream(port_a.Xi_outflow);
  port_a.Xi_outflow = inStream(port_b.Xi_outflow);
  port_b.C_outflow = inStream(port_a.C_outflow);
  port_a.C_outflow = inStream(port_b.C_outflow);

annotation(
    Diagram(graphics = {Rectangle(extent = {{-80, 40}, {80, -40}}), Line(origin = {-89, 0}, points = {{-9, 0}, {9, 0}, {9, 0}}), Line(origin = {90, 0}, points = {{-10, 0}, {10, 0}, {10, 0}}), Line(origin = {-10, 0}, points = {{10, 40}, {-10, -40}, {-10, -40}}), Line(origin = {10, 0}, points = {{10, 40}, {-10, -40}, {-10, -40}}), Text(origin = {-30, 0}, extent = {{-32, 16}, {32, -16}}, textString = "M"), Text(origin = {30, 0}, extent = {{-32, 16}, {32, -16}}, textString = "M")}, coordinateSystem(initialScale = 0.1)),
    Icon(graphics = {Line(origin = {-10, 0}, points = {{10, 20}, {-10, -20}, {-10, -20}}), Line(origin = {10, 0}, points = {{10, 20}, {-10, -20}, {-10, -20}}), Text(origin = {-19, 11}, extent = {{-11, 15}, {11, -15}}, textString = "m"), Text(origin = {21, -9}, extent = {{-11, 15}, {11, -15}}, textString = "m"), Rectangle(extent = {{-30, 20}, {30, -20}}), Line(origin = {-45, -1}, points = {{-15, 1}, {15, 1}, {15, -1}}), Line(origin = {45, 0}, points = {{-15, 0}, {15, 0}, {15, 0}})}, coordinateSystem(initialScale = 0.1)));
end mass_loop_breaker;