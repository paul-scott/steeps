within SolarTherm.Models.Fluid.Sources;
model FluidSink "Infinite fluid sink"
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
    Icon(graphics = {Rectangle(origin = {0, -1}, extent = {{-68, 27}, {82, -25}}), Text(origin = {-11, 14}, extent = {{-23, 10}, {57, -38}}, textString = "Sink"), Line(origin = {-80, 0}, points = {{12, 0}, {-20, 0}, {-20, 0}}), Text(origin = {-2, 184}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {151, -154}}, textString = "%name")}, coordinateSystem(initialScale = 0.1)));
end FluidSink;