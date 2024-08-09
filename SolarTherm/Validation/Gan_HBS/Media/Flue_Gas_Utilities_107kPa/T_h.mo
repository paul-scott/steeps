within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa;
function T_h "Temperature (K) of air at ambient pressure as a function of specific enthalpy"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
	output Modelica.SIunits.Temperature T "Temperature (K)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real h_data[16] = {412434.46, 414229.86, 513319.41, 616262.56, 722915.80, 833144.11, 946724.94, 1063356.42, 1182701.78, 1304427.18, 1428223.77, 1553817.20, 1680969.34, 1809476.26, 1939164.38, 2069886.40};
algorithm
    T := SolarTherm.Utilities.Interpolation.Interpolate1D(h_data,T_data,h);
end T_h;