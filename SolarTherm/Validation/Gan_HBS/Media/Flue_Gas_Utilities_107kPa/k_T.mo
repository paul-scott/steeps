within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa;
function k_T "Thermal conductivity (W/mK) of flue gas at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity k "Thermal Conductivity (W/mK)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real k_data[16] = {0.02288, 0.02302, 0.03022, 0.03710, 0.04369, 0.05000, 0.05607, 0.06193, 0.06759, 0.07310, 0.07845, 0.08367, 0.08878, 0.09378, 0.09870, 0.10353};
algorithm
    k := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
end k_T;