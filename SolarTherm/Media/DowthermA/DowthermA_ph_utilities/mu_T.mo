within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function mu_T "Dynamic viscosity (Ns/m2) of DowthermA saturated liquid as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.DynamicViscosity mu "Dynamic Viscosity (Ns/m2)";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real mu_data[9] = {0.00500, 0.00158, 0.00091, 0.00056, 0.00038, 0.00027, 0.00020, 0.00016, 0.00012};
algorithm
    mu := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,mu_data,T);
end mu_T;