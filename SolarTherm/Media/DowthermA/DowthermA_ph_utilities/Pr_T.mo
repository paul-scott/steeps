within SolarTherm.Media.DowthermA.DowthermA_ph_utilities;
function Pr_T "Prandtl number of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalDiffusivity Pr "Prandtl Number";
protected
    Real T_data[9] = {298.15, 338.15, 378.15, 428.15, 478.15, 528.15, 578.15, 628.15, 678.15};
    Real Pr_data[9] = {56.8674, 20.4379, 13.1954, 9.3445, 7.2900, 5.9582, 5.0977, 4.7511, 4.2412};
algorithm
    Pr := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,Pr_data,T);
end Pr_T;