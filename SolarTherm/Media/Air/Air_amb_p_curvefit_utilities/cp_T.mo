within SolarTherm.Media.Air.Air_amb_p_curvefit_utilities;
function cp_T "Specific heat capacity (J/kgK) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
algorithm
    cp := 0.18784*T+948.9; //Derivative of h-T
end cp_T;