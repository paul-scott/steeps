within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities;
function T_h "Temperature (K) of air at ambient pressure as a function of specific enthalpy"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
	output Modelica.SIunits.Temperature T "Temperature (K)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real h_data[16] = {170.2, 2010.3, 102844.1, 206363.8, 312569.6, 421461.3, 533039.0, 647302.8, 764252.5, 883888.3, 1006210.0, 1131217.7, 1258911.5, 1389291.2, 1522357.0, 1658108.7};
algorithm
    T := SolarTherm.Utilities.Interpolation.Interpolate1D(h_data,T_data,h);
end T_h;
