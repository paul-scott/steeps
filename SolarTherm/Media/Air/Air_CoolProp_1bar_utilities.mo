within SolarTherm.Media.Air;
package Air_CoolProp_1bar_utilities
 extends Modelica.Icons.UtilitiesPackage;

  function cp_T "Specific heat capacity (J/kgK) of air at ambient pressure as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real cp_data[16] = {1005.93, 1005.99, 1013.76, 1029.49, 1050.81, 1074.56, 1098.26, 1120.46, 1140.53, 1158.33, 1173.98, 1187.71, 1199.79, 1210.47, 1219.97, 1228.48};
  algorithm
    cp := Modelica.Math.Vectors.interpolate(T_data, cp_data, T);
  end cp_T;

  function mu_T "Dynamic viscosity (Ns/m2) of air at ambient pressure as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.DynamicViscosity mu "Dynamic Viscosity (Ns/m2)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real mu_data[16] = {0.00001839, 0.00001848, 0.00002298, 0.00002701, 0.00003068, 0.00003407, 0.00003726, 0.00004028, 0.00004316, 0.00004592, 0.00004859, 0.00005119, 0.00005371, 0.00005618, 0.00005860, 0.00006097};
  algorithm
    mu := Modelica.Math.Vectors.interpolate(T_data, mu_data, T);
  end mu_T;

  function h_T "Specific enthalpy (J/kg) of air at ambient pressure as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real h_data[16] = {298266.93, 300127.95, 401045.95, 503148.93, 607130.53, 713390.74, 822039.63, 932991.58, 1046059.73, 1161021.03, 1277653.20, 1395752.52, 1515140.51, 1635664.50, 1757195.57, 1879625.36};
  algorithm
    h := Modelica.Math.Vectors.interpolate(T_data, h_data, T);
  end h_T;

  function k_T "Thermal conductivity (W/mK) of air at ambient pressure as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.ThermalConductivity k "Thermal Conductivity (W/mK)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real k_data[16] = {0.02586, 0.02600, 0.03296, 0.03934, 0.04529, 0.05091, 0.05628, 0.06144, 0.06644, 0.07131, 0.07607, 0.08073, 0.08532, 0.08984, 0.09431, 0.09873};
  algorithm
    k := Modelica.Math.Vectors.interpolate(T_data, k_data, T);
  end k_T;

  function rho_T "Density (kg/m3) of air at ambient pressure as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.Density rho "Density (kg/m3)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real rho_data[16] = {1.1690, 1.1617, 0.8709, 0.6966, 0.5805, 0.4976, 0.4354, 0.3870, 0.3483, 0.3167, 0.2903, 0.2679, 0.2488, 0.2322, 0.2177, 0.2049};
  algorithm
    rho := Modelica.Math.Vectors.interpolate(T_data, rho_data, T);
  end rho_T;

  function T_h "Temperature (K) of air at ambient pressure as a function of specific enthalpy"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
    output Modelica.SIunits.Temperature T "Temperature (K)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real h_data[16] = {298266.93, 300127.95, 401045.95, 503148.93, 607130.53, 713390.74, 822039.63, 932991.58, 1046059.73, 1161021.03, 1277653.20, 1395752.52, 1515140.51, 1635664.50, 1757195.57, 1879625.36};
  algorithm
    T := Modelica.Math.Vectors.interpolate(h_data, T_data, h);
  end T_h;
  annotation();
end Air_CoolProp_1bar_utilities;