within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function rho_T "Density (kg/m3) of DowthermA saturated liquid as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.Density rho "Density (kg/m3)";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real rho_data[9] = {1063.5, 1023.7, 990.7, 947.8, 902.5, 854.0, 801.3, 724.3, 672.5};
algorithm
    rho := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,rho_data,T);
end rho_T;