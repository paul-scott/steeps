within SolarTherm.Models.Fluid.Valves;

model PBS_TeeJunction_LoopBreaker
  replaceable package Medium = SolarTherm.Media.Sodium.Sodium_pT
    constrainedby Modelica.Media.Interfaces.PartialMedium
    "Medium in the component"
      annotation (choicesAllMatching = true);
  Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = Medium) annotation (Placement(visible = true,transformation(extent = {{-10, 0}, {10, 20}}, rotation = 0), iconTransformation(origin = {-80, 80}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = Medium) annotation (Placement(visible = true,transformation(extent = {{-6, 0}, {14, 20}}, rotation = 0), iconTransformation(origin = {82, 80}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_c(redeclare package Medium = Medium) annotation (Placement(visible = true,transformation(extent = {{-68, 60}, {-48, 80}}, rotation = 0),
        iconTransformation(origin = {0, -2}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));

equation
  //This component is a loop breaker, thus no pressure equality or mass conservation is needed.




  if fluid_a.m_flow + fluid_b.m_flow < 0.0 then //fluid_a and fluid_c merge into fluid_b
    fluid_a.h_outflow = inStream(fluid_a.h_outflow);
    fluid_b.h_outflow = ((fluid_a.m_flow*(inStream(fluid_c.h_outflow)-inStream(fluid_a.h_outflow)))+(fluid_b.m_flow*inStream(fluid_c.h_outflow)))/(fluid_b.m_flow);
    fluid_c.h_outflow = inStream(fluid_c.h_outflow);
  else
    fluid_a.h_outflow = inStream(fluid_a.h_outflow);
    fluid_b.h_outflow = inStream(fluid_a.h_outflow);
    fluid_c.h_outflow = inStream(fluid_a.h_outflow);
  end if;
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,0},
            {100,100}})),              Icon(coordinateSystem( preserveAspectRatio = false, initialScale = 0.1), graphics = {Line(origin = {-40.0259, 80.013}, points = {{-40, 0}, {40, 0}, {120, 0}}, thickness = 1), Line(origin = {0, 40}, points = {{0, 40}, {0, -20}, {0, -40}}, thickness = 1), Polygon(origin = {40, 80}, fillColor = {135, 135, 135}, fillPattern = FillPattern.Solid, points = {{-10, 10}, {-10, -10}, {10, 0}, {-10, 10}, {-10, 10}}), Polygon(origin = {-40, 80}, fillColor = {135, 135, 135}, fillPattern = FillPattern.Solid, points = {{-10, 10}, {-10, -10}, {10, 0}, {-10, 10}, {-10, 10}}), Polygon(origin = {0, 43}, fillColor = {124, 124, 124}, fillPattern = FillPattern.Solid, points = {{-8, 7}, {8, 7}, {0, -7}, {-8, 7}, {-8, 7}}), Polygon(origin = {0, 59}, rotation = 180, fillColor = {124, 124, 124}, fillPattern = FillPattern.Solid, points = {{-8, 7}, {8, 7}, {0, -7}, {-8, 7}, {-8, 7}}), Line(origin = {0, 59.2254}, points = {{-12, -34}, {12, -34}, {12, 8}, {60, 8}, {60, 34}, {-60, 34}, {-60, 8}, {-12, 8}, {-12, -34}, {-12, -34}}), Text(origin = {1, 102}, extent = {{-17, 8}, {17, -8}}, textString = "LB")}),
    Documentation(revisions ="<html>
		<p>By Zebedee Kee on 03/12/2020</p>
		</html>",info="<html>
		<p>This component is a pipe teejunction where fluid_a is always an inlet and fluid_b is always an outlet. Fluid_c however can flow two-ways therefore depending on the relative flow rates of fluid_a and fluid_b, outlet enthalpies of b and c are calculated. This component is a loop breaker, it assumes that the pressure equality and mass conservation equations are already calculated elsewhere in the loop.</p>
		</html>"));
end PBS_TeeJunction_LoopBreaker;