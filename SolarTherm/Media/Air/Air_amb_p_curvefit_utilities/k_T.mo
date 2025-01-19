within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function k_T "Thermal conductivity (W/mK) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity k "Thermal Conductivity (W/mK)";
algorithm
    k := (2.138E-11)*T*T*T - (6.383E-8)*T*T + (1.109E-4)*T - 1.538E-3; //R^2=0.9999 between 100 and 1300K 
end k_T;