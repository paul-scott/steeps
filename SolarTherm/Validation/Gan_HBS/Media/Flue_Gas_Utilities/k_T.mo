within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities;
function k_T "Thermal conductivity (W/mK) of flue gas at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity k "Thermal Conductivity (W/mK)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real k_data[16] = {0.02228, 0.02242, 0.02964, 0.03656, 0.04317, 0.04945, 0.05540, 0.06101, 0.06627, 0.07119, 0.07573, 0.07991, 0.08371, 0.08712, 0.09014, 0.09276};
algorithm
    k := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
end k_T;
