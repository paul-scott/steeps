within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function h_T "Specific enthalpy (J/kg) of DownThermA saturated liquid as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real h_data[9] = {0.0, 65180.0, 135480.0, 229680.0, 330855.0, 438955.0, 554055.0, 676555.0, 807855.0};
algorithm
    h := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,h_data,T);
end h_T;