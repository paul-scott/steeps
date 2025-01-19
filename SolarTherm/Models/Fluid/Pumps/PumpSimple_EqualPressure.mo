within SolarTherm.Models.Fluid.Pumps;

model PumpSimple_EqualPressure
  extends SolarTherm.Interfaces.Models.Pump;
  Modelica.Blocks.Interfaces.RealInput m_flow annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=-90,
        origin={-10,78}), iconTransformation(
        extent={{-14,-14},{14,14}},
        rotation=-90,
        origin={0,86})));
  parameter Real k_loss(unit="J/kg")=0.55e3;
  SI.Power W_loss;
equation
  fluid_b.m_flow=-m_flow;
  fluid_a.m_flow=m_flow;
  fluid_a.h_outflow=fluid_b.h_outflow;
  fluid_b.h_outflow=inStream(fluid_a.h_outflow);
  fluid_a.Xi_outflow=fluid_b.Xi_outflow;
  fluid_b.Xi_outflow=inStream(fluid_a.Xi_outflow);
  fluid_a.p=fluid_b.p;
  W_loss=k_loss*m_flow;
  annotation (Documentation(revisions = "<html>
<ul>
<li>Zebedee Kee:<br>Modified to allow equal pressures on both sides. </li>
</ul>
</html>"),
    Icon(graphics = {Text(origin = {-10, 14}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {151, -154}}, textString = "%name")}));
end PumpSimple_EqualPressure;