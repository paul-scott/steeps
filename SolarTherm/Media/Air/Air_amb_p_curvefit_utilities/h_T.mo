within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function h_T "Specific enthalpy (J/kg) of air at ambient pressure as a function of temperature (reference 0 J/kg at 100K)"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
algorithm
    h := (9.392E-2)*T*T + (9.489E2)*T - 9.277E4; //R^2=0.9999 between 100 and 1300K 
end h_T;