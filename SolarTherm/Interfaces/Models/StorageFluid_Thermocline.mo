within SolarTherm.Interfaces.Models;
partial model StorageFluid_Thermocline
  extends SolarTherm.Icons.ThermoclineModel;
  Modelica.Fluid.Interfaces.FluidPort_a fluid_a(redeclare package Medium = Medium) annotation (Placement(
        visible = true,transformation(extent = {{-10, 70}, {10, 90}}, rotation = 0), iconTransformation(origin = {6.66134e-16, 80}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_b(redeclare package Medium = Medium) annotation (Placement(
        visible = true,transformation(extent = {{-10, -90}, {10, -70}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  replaceable package Medium = SolarTherm.Media.MoltenSalt.MoltenSalt_ph
    constrainedby Modelica.Media.Interfaces.PartialMedium
    "Medium in the component"
      annotation (choicesAllMatching = true,Dialog(group="Working fluid"));
  annotation(
    Icon(graphics = {Line(origin = {0, 70}, points = {{0, -10}, {0, 10}, {0, 10}}), Line(origin = {0, -70}, points = {{0, 10}, {0, -10}, {0, -10}})}));
end StorageFluid_Thermocline;