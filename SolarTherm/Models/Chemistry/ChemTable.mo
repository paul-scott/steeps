within SolarTherm.Models.Chemistry;

package ChemTable
  //Data obtained from FactSage
  //Format M (kg/mol), Hf0 (J/mol), S0 (J/molK)
  //       Molar mass, Standard enthalpy of formation, Absolute entropy at standard conditions
  final Chemical H2(M=2.01588e-3, Hf0=0.0, S0=130.680);
  final Chemical H2O(M=18.01528e-3, Hf0=-285.8300e3, S0=69.950);
  final Chemical Fe(M=55.845e-3, Hf0=0.0, S0=27.280);
  final Chemical FeO(M=71.8444e-3, Hf0=-265.8322e3, S0=59.496);//HfO changed from -272.04e3 to -265.832e3
  final Chemical Fe3O4(M=231.5326e-3, Hf0=-1117.4025e3, S0=146.114);
  final Chemical Fe2O3(M=159.6882e-3, Hf0=-825.7869e3, S0=87.729);
  
  final Chemical Fe2O3H2O(M=177.70348e-3, Hf0=-1117.998e3, S0=118.826);
  final Chemical Al2O3(M=101.9612772e-3, Hf0=-1675.7e3, S0=50.820);
  final Chemical SiO2(M=60.0843e-3, Hf0=-910.6999e3, S0=41.460);
  
  final Chemical MgO(M=40.3044e-3, Hf0=-601.500e3, S0=26.951);
  final Chemical Al6Si2O13(M=426.0524316e-3, Hf0=-6819.210e3, S0=274.9);
  final Chemical CaO(M=56.0774e-3, Hf0=-635.0900e3, S0=37.750);
  final Chemical CaCO3(M=100.0869e-3, Hf0=-1206.6656e3, S0=87.990);
  
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Rectangle(lineColor = {200, 200, 200}, fillColor = {248, 248, 248}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-100, -100}, {100, 100}}, radius = 25), Rectangle(origin = {0, -11}, lineThickness = 1, extent = {{-72, 69}, {72, -53}}), Line(origin = {-1, 34}, points = {{-71, 0}, {73, 0}}, thickness = 1), Line(origin = {-40, -3}, points = {{0, 61}, {0, -61}}, thickness = 1), Line(origin = {-2, -3}, points = {{0, 61}, {0, -61}}, thickness = 1), Line(origin = {36, -4}, points = {{0, 62}, {0, -60}}, thickness = 1)}, coordinateSystem(initialScale = 0.1)),
  Documentation(info = "<html><head></head><body>Molar mass (kg/mol), Standard enthalpy of formation (J/mol) and Standard absolute entropy (J/mol K) for chemical substances.<div><br></div><div>Values were obtained from FactSage.</div><div><br></div><div><font face=\"Arial, Helvetica, sans-serif\">C. W. Bale, E. 
Bélisle, P. Chartrand, S. A. Decterov, G. Eriksson, A.E. Gheribi, K. 
Hack, I. H. Jung, Y. B. Kang, J. Melançon, A. D. Pelton, S. Petersen, C.
 Robelin. J. Sangster, P. Spencer and M-A. Van Ende, <em>FactSage Thermochemical Software and Databases</em> - 2010 - 2016, Calphad, vol. 54, pp 35-53, 2016 &lt;www.factsage.com&gt;</font></div></body></html>"));
end ChemTable;