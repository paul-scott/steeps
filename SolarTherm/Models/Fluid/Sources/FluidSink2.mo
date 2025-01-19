within SolarTherm.Models.Fluid.Sources;

model FluidSink2 "Fluid sink but does not fix any pressure or mass-fraction value."
 replaceable package Medium = Modelica.Media.Interfaces.PartialMedium;
 Modelica.Fluid.Interfaces.FluidPort_a port_a(
  redeclare package Medium=Medium,
  m_flow(min=0))     annotation (Placement(visible = true,transformation(extent={{-110,-12},{-90,8}},  rotation=
           0), iconTransformation(origin = {-100, 8.88178e-16}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
           
  //Medium.BaseProperties Fluid "Fluid state calculator";
equation
  //Fluid.h = inStream(port_a.h_outflow);
  //Fluid.p = port_a.p;
 port_a.h_outflow = inStream(port_a.h_outflow);//0; // shouldn't flow backwards anyway
annotation(
    Icon(graphics = {Rectangle(origin = {0, -1}, extent = {{-60, 41}, {60, -39}}), Text(origin = {3, 0}, extent = {{-23, 16}, {23, -16}}, textString = "Sink"), Line(origin = {-80, 0}, points = {{20, 0}, {-20, 0}, {-20, 0}}), Text(origin = {0, 72}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {151, -154}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
end FluidSink2;