within SolarTherm.Materials;
partial package PartialMaterial
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  package Tables = Modelica.Blocks.Tables;
  import SolarTherm.Utilities.Interpolation.Interpolate1D;
  
  constant SI.MolarMass MM = 0.0 "Molar mass (kg/mol)";
  constant SI.Temperature T_melt = 0.0 "Melting point (K)";
  constant Real cost = 0.0 "USD/kg";
  
  replaceable partial function h_Tf
    "Return specific enthalpy as a function of temperature and liquid mass fraction f"
    extends Modelica.Icons.Function;
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  end h_Tf;
  
  replaceable partial function rho_Tf
    "Return specific enthalpy as a function of temperature and liquid mass fraction f"
    extends Modelica.Icons.Function;
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.Density rho "Specific Enthalpy (J/kg)";
  end rho_Tf;
  
  replaceable partial model State "A model which calculates state and properties"
    	SI.SpecificEnthalpy h "Specific Enthalpy wrt 298.15K (J/kg)";
	SI.Temperature T "Temperature (K)";
	Real f "Liquid Mass Fraction";
	SI.Density rho "Density (kg/m3)";
	SI.ThermalConductivity k "Thermal conductivity (W/mK)";
  end State;
  
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Rectangle(lineColor = {200, 200, 200}, fillColor = {248, 248, 248}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, -100}, {100, 100}}, radius = 25), Rectangle(origin = {24, -19}, lineThickness = 1, extent = {{-44, 43}, {44, -43}}), Line(origin = {-47, -43}, points = {{-27, 19}, {27, -19}}, thickness = 1), Rectangle(origin = {-30, 19}, lineThickness = 1, extent = {{-44, 43}, {44, -43}}), Line(origin = {40.9688, 42.9689}, points = {{-27, 19}, {27, -19}}, thickness = 1), Ellipse(origin = {-73, -23}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360), Ellipse(origin = {67, 23}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360), Line(origin = {40.9375, -42.6719}, points = {{-27, 19}, {27, -19}}, thickness = 1), Ellipse(origin = {-21, -59}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360), Line(origin = {-47, 42.9688}, points = {{-27, 19}, {27, -19}}, thickness = 1), Ellipse(origin = {13, 61}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360), Ellipse(origin = {-19, 23}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360), Ellipse(origin = {13, -23}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360), Ellipse(origin = {69, -61}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360), Ellipse(origin = {-73, 61}, fillColor = {200, 200, 200}, fillPattern = FillPattern.Sphere, lineThickness = 1, extent = {{-13, 13}, {13, -13}}, endAngle = 360)}));
end PartialMaterial;
