within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function T_h "Temperature (K) of air at ambient pressure as a function of specific enthalpy"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
	output Modelica.SIunits.Temperature T "Temperature (K)";
algorithm
    T := (-7.486E-11)*h*h + (1.021E-3)*h + 97.99; //R^2 = 0.9999 from 100K to 1300K
end T_h;