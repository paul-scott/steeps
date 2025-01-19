within SolarTherm.Interfaces.Models;

partial model Electric_Heater
  //extends SolarTherm.Interfaces.Models.Receiver;
  Modelica.Blocks.Interfaces.RealInput P_supply "Supplied electrical power" annotation (Placement(
        visible = true,transformation(
        
        origin={44,84},extent={{-12,-12},{12,12}},
        rotation=-90), iconTransformation(
        
        origin={-112, 68},extent={{12, -12}, {-12, 12}},
        rotation= 180)));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = Medium) annotation (Placement(
        visible = true,transformation(extent = {{10, -100}, {30, -80}}, rotation = 0), iconTransformation(extent = {{-110, -8}, {-90, 12}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = Medium) annotation (Placement(
        visible = true,transformation(extent = {{24, 38}, {44, 58}}, rotation = 0), iconTransformation(extent = {{90, -10}, {110, 10}}, rotation = 0)));
  replaceable package Medium = SolarTherm.Media.MoltenSalt.MoltenSalt_ph
    constrainedby Modelica.Media.Interfaces.PartialMedium
    "Medium in the component"
      annotation (choicesAllMatching = true,Dialog(group="Working fluid"));
  annotation(
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Rectangle(lineColor = {200, 200, 200}, fillColor = {248, 248, 248}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, -100}, {100, 100}}, radius = 25), Ellipse(origin = {1, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}}, endAngle = 360), Polygon(origin = {1, -3}, fillColor = {252, 233, 79}, fillPattern = FillPattern.Solid, points = {{1, 49}, {-31, -5}, {5, -1}, {-1, -43}, {27, 17}, {-7, 11}, {1, 49}}), Line(origin = {-75, 2}, points = {{25, 0}, {-25, 0}}), Line(origin = {75.0846, 0.57827}, points = {{25, 0}, {-25, 0}})}));

end Electric_Heater;