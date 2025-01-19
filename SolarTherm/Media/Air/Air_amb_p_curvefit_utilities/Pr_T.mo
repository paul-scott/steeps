within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function Pr_T "Prandtl number of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalDiffusivity Pr "Prandtl Number";
algorithm
    Pr := cp_T(T)*mu_T(T)/k_T(T);
end Pr_T;