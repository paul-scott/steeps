within SolarTherm.Models.Chemistry;

partial record Chemical
  parameter Modelica.SIunits.MolarMass M "Molar mass (kg/mol)";
  parameter Modelica.SIunits.MolarEnthalpy Hf0 "Standard enthalpy of formation (J/mol)";
  parameter Modelica.SIunits.MolarEntropy S0 "Absolute entropy at standard conditions (J/molK)";
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Chemical;