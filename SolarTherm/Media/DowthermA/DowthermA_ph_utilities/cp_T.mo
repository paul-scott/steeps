within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function cp_T "Specific heat capacity (J/kgK) of DowthermA saturated liquid as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real cp_data[9] = {1586.6, 1701.0, 1814.0, 1954.0, 2093.0, 2231.0, 2373.0, 2527.0, 2725.0};
algorithm
    cp := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,cp_data,T);
end cp_T;