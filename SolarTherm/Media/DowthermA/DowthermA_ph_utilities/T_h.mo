within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function T_h "Temperature (K) of DowthermA saturated temperature as a function of specific enthalpy"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.SpecificEnthalpy h1 "Specific Enthalpy (J/kg)";
	output Modelica.SIunits.Temperature T "Temperature (K)";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real h_data[9] = {0.0, 65180.0, 135480.0, 229680.0, 330855.0, 438955.0, 554055.0, 676555.0, 807855.0};
algorithm
    T := SolarTherm.Utilities.Interpolation.Interpolate1D(h_data,T_data,h1);
end T_h;