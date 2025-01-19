within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function rho_T "Density (kg/m3) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.Density rho "Density (kg/m3)";
algorithm
    rho := 351.9/T; //R^2 = 0.9999 from 100 to 1300 K (force intercept of 1/T vs rho to 0)
end rho_T;