within SolarTherm.Validation.Gan_HBS.Materials;
//Experimental Data Checkerbrick_5
package Checkerbrick_5
  extends SolarTherm.Materials.PartialMaterial(MM = 426.0524e-3, T_melt = 2000.0+273.15, cost = 0.65);
  import SolarTherm.Utilities.Interpolation.Interpolate1D;
  
  redeclare model State "A model which calculates state and properties"
    SI.SpecificEnthalpy h "Specific Enthalpy wrt 298.15K (J/kg)";
    SI.Temperature T "Temperature (K)";
    Real f "Liquid Mass Fraction";
    SI.Density rho "Density (kg/m3)";
    SI.ThermalConductivity k "Thermal conductivity (W/mK)";
  equation
    f = 0.0;
    h = h_Tf(T, 0);
    rho = rho_Tf(T, 0);
    k = k_Tf(T, 0);
  end State;

  redeclare function h_Tf "find specific enthalpy from Temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  protected
    Real T_data[8] = {298.15, 400.00, 600.00, 800.00, 1000.00, 1200.00, 1400.00, 1600.00};
    Real h_data[8] = {0.0, 88606.3, 288106.3, 516506.3, 754506.3, 993406.3, 1238606.3, 1489606.3};
  algorithm
    h := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,h_data,T);
  end h_Tf;

  redeclare function rho_Tf "find density from temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.Density rho "Density (kg/m3)";
  algorithm
    rho := 2500.0;
  end rho_Tf;

  function k_Tf "find thermal conductivity from temperature"
    input SI.Temperature T;
    input Real f;
    output SI.ThermalConductivity k;
  protected
    Real T_data[3] = {473.0, 973.0, 1573.0};
    Real k_data[3] = {1.65, 1.65, 1.95};
  algorithm
    k := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
  end k_Tf;

  function Tf_h "Find temperature and liquid fraction from temperature"
    input SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
    output SI.Temperature T "Absoulte temperature (K)";
    output Real f "mass liquid fraction";
protected
    Real T_data[8] = {298.15, 400.00, 600.00, 800.00, 1000.00, 1200.00, 1400.00, 1600.00};
    Real h_data[8] = {0.0, 88606.3, 288106.3, 516506.3, 754506.3, 993406.3, 1238606.3, 1489606.3};
  algorithm
    T := SolarTherm.Utilities.Interpolation.Interpolate1D(h_data,T_data,h);
    f := 0.0;
  end Tf_h;
end Checkerbrick_5;