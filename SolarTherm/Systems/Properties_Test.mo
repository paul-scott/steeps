within SolarTherm.Systems;

model Properties_Test
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  
  parameter SI.Pressure p = 0.20*1e5 "Partial pressure of gas (Pa)";
  SI.Temperature T(start = 298.15) "Temperature sweep (K)";
  
  SI.Temperature T_sat "Saturation temperature at pressure p (K)";
  
  SI.SpecificHeatCapacity cp "cp of substance (J/kgK)";
  
  SI.Density rho_Fe2O3 "Density of pure Fe2O3 (kg/m3)";
  
equation
  T = 765.15 + time;
  cp = Modelica.Media.Water.IF97_Utilities.cp_pT(p,T);
  T_sat = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
  rho_Fe2O3 = SolarTherm.Media.SolidParticles.Fe2O3_utilities.rho_T(T);
  
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Properties_Test;