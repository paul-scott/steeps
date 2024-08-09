within SolarTherm.Validation.Gan_HBS.Materials;

package Air_533kPa
  //Composition {N2 = 78 mol%, O2 = 21 mol%, Ar = 1 mol%} Pressure = 532818 Pa
  extends SolarTherm.Materials.PartialMaterial(MM = 56.077e-3, T_melt = 3200.0, cost = 0.2);

  model State "A model which calculates state and properties"
    parameter SI.SpecificEnthalpy h_start;
    SI.Density rho;
    SI.SpecificEnthalpy h(start = h_start);
    SI.Temperature T;
    SI.DynamicViscosity mu;
    SI.SpecificHeatCapacity cp;
    SI.ThermalConductivity k;
  equation
    T = SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.T_h(h);
    rho = SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.rho_T(T);
    cp = SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.cp_T(T);
    mu = SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.mu_T(T);
    k = SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.k_T(T);
  end State;

  function h_Tf "Specific enthalpy of air vs Temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  algorithm
    h := SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.h_T(T);
  end h_Tf;

  function rho_Tf "Density of air vs Temperature"
    input SI.Temperature T;
    input Real f = 0 "Liquid mass melt fraction (No effect on result)";
    output SI.Density rho;
  algorithm
    rho := SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.rho_T(T);
  end rho_Tf;

  function T_h "Temperature interpolated from enthalpy"
    input SI.SpecificEnthalpy h;
    output SI.Temperature T;
  algorithm
    T := SolarTherm.Validation.Gan_HBS.Media.Air_533kPa_Utilities.T_h(h);
  end T_h;
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Air_533kPa;