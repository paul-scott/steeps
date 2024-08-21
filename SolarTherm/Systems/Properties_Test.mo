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
  
  
  
equation
  T = 298.15 + time;
  cp = Modelica.Media.Water.IF97_Utilities.cp_pT(p,T);
  T_sat = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p);
  
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Properties_Test;