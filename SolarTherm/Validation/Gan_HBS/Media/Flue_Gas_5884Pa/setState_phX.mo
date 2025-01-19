within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_5884Pa;

function setState_phX "Return thermodynamic state from p, h, and X or Xi"
  extends Modelica.Icons.Function;
  input AbsolutePressure p "Pressure";
  input SpecificEnthalpy h "Specific enthalpy";
  input MassFraction X[:] = reference_X "Mass fractions";
  output ThermodynamicState state "Thermodynamic state record";
algorithm
  state := ThermodynamicState(p = p, h = h);
end setState_phX;