within SolarTherm.Validation.Gan_HBS.Materials;

package Flue_Gas_5884Pa
//Composition {CO2 = 30.7 mol%, H2O = 1.8 mol%, N2 = 67.2 mol%, O2 = 0.3 mol%} Pressure = 5884 Pa 
  extends SolarTherm.Materials.PartialMaterial(MM = 56.077e-3, T_melt = 3200.0, cost = 0.2);
  redeclare model State "A model which calculates state and properties of Flue Gas"
    parameter SI.SpecificEnthalpy h_start;
    SI.Density rho;
    SI.SpecificEnthalpy h (start=h_start);
    SI.Temperature T;
    SI.DynamicViscosity mu;
    SI.SpecificHeatCapacity cp;
    SI.ThermalConductivity k;
  equation
    T = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.T_h(h);
    rho = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.rho_T(T);
    cp = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.cp_T(T);
    mu = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.mu_T(T);
    k = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.k_T(T);
  end State;
  
  redeclare function h_Tf "Specific enthalpy of air vs Temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  algorithm
    h := SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.h_T(T);
  end h_Tf;
  
  redeclare function rho_Tf "Density of air vs Temperature"
    input SI.Temperature T;
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.Density rho;
  algorithm
    rho := SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.rho_T(T);
  end rho_Tf;


  function T_h "Temperature interpolated from enthalpy"
    input SI.SpecificEnthalpy h;
    output SI.Temperature T;
  algorithm
    T := SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities.T_h(h);
  end T_h;
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Flue_Gas_5884Pa;