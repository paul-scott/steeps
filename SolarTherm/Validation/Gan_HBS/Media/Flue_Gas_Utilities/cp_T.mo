within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities;
function cp_T "Specific heat capacity (J/kgK) of flue gas at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real cp_data[16] = {995.92, 996.43, 1025.79, 1057.65, 1090.28, 1122.36, 1152.92, 1181.31, 1207.12, 1230.16, 1250.39, 1267.92, 1282.95, 1295.75, 1306.63, 1315.91};
algorithm
    cp := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,cp_data,T);
end cp_T;
