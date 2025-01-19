within SolarTherm.Models.Chemistry;

package ChemTable_NIST
  //Format M (kg/mol), Hf0 (J/mol), S0 (J/molK)
  //       Molar mass, Standard enthalpy of formation, Absolute entropy at standard conditions
  final Chemical H2(M=2.01588e-3, Hf0=0.0, S0=130.68);
  final Chemical H2O(M=18.0153e-3,Hf0=-285.83e3, S0=69.95);
  final Chemical Fe(M=55.845e-3,Hf0=0.0, S0=27.32);
  final Chemical FeO(M=71.844e-3,Hf0=-272.04e3, S0=60.75);
  final Chemical Fe3O4(M=231.533e-3,Hf0=-1120.89e3, S0=145.245);
  final Chemical Fe2O3(M=159.688e-3,Hf0=-825.5e3, S0=87.33);

  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Rectangle(lineColor = {200, 200, 200}, fillColor = {248, 248, 248}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, -100}, {100, 100}}, radius = 25), Rectangle(origin = {0, -11}, lineThickness = 1, extent = {{-72, 69}, {72, -53}}), Line(origin = {-1, 34}, points = {{-71, 0}, {73, 0}}, thickness = 1), Line(origin = {-40, -3}, points = {{0, 61}, {0, -61}}, thickness = 1), Line(origin = {-2, -3}, points = {{0, 61}, {0, -61}}, thickness = 1), Line(origin = {36, -4}, points = {{0, 62}, {0, -60}}, thickness = 1)}, coordinateSystem(initialScale = 0.1)),
  Documentation(info = "<html><head></head><body>Molar mass (kg/mol), Standard enthalpy of formation (J/mol) and Standard absolute entropy (J/mol K) for chemical substances.<div><br></div><div>Values were obtained from NIST-JANAF website.</div></body></html>"));
end ChemTable_NIST;