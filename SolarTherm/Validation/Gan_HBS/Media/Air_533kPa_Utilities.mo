within SolarTherm.Validation.Gan_HBS.Media;

package Air_533kPa_Utilities
  function cp_T "Specific heat capacity (J/kgK) of air at 532.818kPa at ambient pressure as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real cp_data[16] = {1012.96, 1012.92, 1017.26, 1031.56, 1052.17, 1075.51, 1098.95, 1120.98, 1140.94, 1158.65, 1174.24, 1187.93, 1199.97, 1210.62, 1220.10, 1228.58};
  algorithm
    cp := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data, cp_data, T);
  end cp_T;

  function mu_T "Dynamic viscosity (Ns/m2) of air at 532.818kPa as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.DynamicViscosity mu "Dynamic Viscosity (Ns/m2)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real mu_data[16] = {0.00001846, 0.00001855, 0.00002303, 0.00002704, 0.00003071, 0.00003410, 0.00003728, 0.00004030, 0.00004317, 0.00004594, 0.00004861, 0.00005120, 0.00005372, 0.00005619, 0.00005861, 0.00006098};
  algorithm
    mu := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data, mu_data, T);
  end mu_T;

  function h_T "Specific enthalpy (J/kg) of air at 532.818kPa as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real h_data[16] = {297277.26, 299151.20, 400561.88, 502934.41, 607084.09, 713457.97, 822188.13, 933200.60, 1046315.22, 1161313.08, 1277974.57, 1396097.79, 1515505.52, 1636045.99, 1757590.95, 1880032.55};
  algorithm
    h := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data, h_data, T);
  end h_T;

  function k_T "Thermal conductivity (W/mK) of air at 532.818kPa as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.ThermalConductivity k "Thermal Conductivity (W/mK)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real k_data[16] = {0.02602, 0.02615, 0.03308, 0.03943, 0.04537, 0.05098, 0.05634, 0.06149, 0.06649, 0.07135, 0.07610, 0.08076, 0.08535, 0.08987, 0.09434, 0.09876};
  algorithm
    k := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data, k_data, T);
  end k_T;

  function rho_T "Density (kg/m3) of air at 532.818kPa as a function of temperature"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Temperature T "Temperature (K)";
    output Modelica.SIunits.Density rho "Density (kg/m3)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real rho_data[16] = {6.236, 6.197, 4.636, 3.706, 3.088, 2.647, 2.316, 2.059, 1.853, 1.685, 1.545, 1.426, 1.324, 1.236, 1.159, 1.091};
  algorithm
    rho := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data, rho_data, T);
  end rho_T;

  function T_h "Temperature (K) of air at 532.818kPa as a function of specific enthalpy"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
    output Modelica.SIunits.Temperature T "Temperature (K)";
  protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real h_data[16] = {297277.26, 299151.20, 400561.88, 502934.41, 607084.09, 713457.97, 822188.13, 933200.60, 1046315.22, 1161313.08, 1277974.57, 1396097.79, 1515505.52, 1636045.99, 1757590.95, 1880032.55};
  algorithm
    T := SolarTherm.Utilities.Interpolation.Interpolate1D(h_data, T_data, h);
  end T_h;
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Air_533kPa_Utilities;