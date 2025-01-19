within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa_Utilities;
function T_h "Temperature (K) of air at ambient pressure as a function of specific enthalpy"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
	output Modelica.SIunits.Temperature T "Temperature (K)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real h_data[16] = {-599.8, 1244.3, 101943.1, 204641.9, 309340.7, 416039.5, 524738.3, 635437.1, 748135.9, 862834.7, 979533.5, 1098232.3, 1218931.1, 1341629.9, 1466328.7, 1593027.5};
algorithm
    T := SolarTherm.Utilities.Interpolation.Interpolate1D(h_data,T_data,h);
end T_h;