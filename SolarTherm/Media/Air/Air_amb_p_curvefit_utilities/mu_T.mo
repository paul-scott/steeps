within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function mu_T "Dynamic viscosity (Ns/m2) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.DynamicViscosity mu "Dynamic Viscosity (Ns/m2)";
algorithm
    mu := (1.112E-14)*T*T*T - (3.854E-11)*T*T + (6.892E-8)*T + 8.271E-7; //R^2=0.9999 between 100 and 1300K 
end mu_T;