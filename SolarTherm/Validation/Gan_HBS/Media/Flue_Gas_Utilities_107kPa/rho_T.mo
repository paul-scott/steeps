within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa;
function rho_T "Density (kg/m3) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.Density rho "Density (kg/m3)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real rho_data[16] = {1.419, 1.410, 1.056, 0.845, 0.704, 0.603, 0.528, 0.469, 0.422, 0.384, 0.352, 0.325, 0.302, 0.282, 0.264, 0.248};
algorithm
    rho := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,rho_data,T);
end rho_T;