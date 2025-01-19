within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function alpha_T "Thermal diffusivity (m2/s) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalDiffusivity alpha "Thermal Diffusivity (m2/s)";
algorithm
    alpha := k_T(T)/(rho_T(T)*cp_T(T));
end alpha_T;