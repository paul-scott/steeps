within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function alpha_T "Thermal diffusivity (m2/s) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalDiffusivity alpha "Thermal Diffusivity (m2/s)";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real alpha_data[9] = {82.6741, 75.5177, 69.6110, 63.2289, 57.7575, 53.0632, 48.9617, 46.4950, 42.0722}; //times 10^-9
algorithm
    alpha := (1.0e-9)*SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,alpha_data,T);
end alpha_T;