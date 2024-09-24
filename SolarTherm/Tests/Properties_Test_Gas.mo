within SolarTherm.Tests;

model Properties_Test_Gas
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  
  //package Utility_Package = SolarTherm.Media.SolidParticles.Fe2O3_utilities;
  
  SI.Temperature T;
  Real T_C;
  SI.Pressure p;
  parameter SI.Temperature T_start = 25.0+273.15;
  parameter SI.Temperature T_end = 1500.0;
  parameter SI.Pressure p_fixed = 100000.0; //1bar
  
  SI.SpecificEnthalpy h_H2;
  SI.SpecificEnthalpy h_H2O;
  
  SI.SpecificHeatCapacity cp_H2;
  SI.SpecificHeatCapacity cp_H2O;
  //SI.SpecificEntropy s;
  SI.Density rho_H2;
  //SI.ThermalConductivity k;
equation
  p = p_fixed;
  T = 54.59 + 273.15;//T_start + (time/100)*(T_end-T_start);
  T_C = T - 273.15;
  //h = Utility_Package.h_T(T);
  //s = Utility_Package.s_T(T);
  rho_H2 = Modelica.Media.IdealGases.SingleGases.H2.density(Modelica.Media.IdealGases.SingleGases.H2.setState_pT(p,T));
  //k = Utility_Package.k_T(T);
  h_H2 = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_fixed,T);
  h_H2O = Modelica.Media.Water.IF97_Utilities.h_pT(p_fixed,T);
  
  cp_H2 = Modelica.Media.IdealGases.Common.Functions.cp_T(Modelica.Media.IdealGases.Common.SingleGasesData.H2,T);
  cp_H2O = Modelica.Media.Water.IF97_Utilities.cp_pT(p_fixed,T);
  
annotation(experiment(StopTime = 100, StartTime = 0, Tolerance = 1e-3, Interval = 0.1),
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Properties_Test_Gas;