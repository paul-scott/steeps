within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function k_T "Thermal conductivity (W/mK) of DowthermA saturated liquid as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity k "Thermal Conductivity (W/mK)";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real k_data[9] = {0.1395, 0.1315, 0.1251, 0.1171, 0.1091, 0.1011, 0.0931, 0.0851, 0.0771};
algorithm
    k := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
end k_T;