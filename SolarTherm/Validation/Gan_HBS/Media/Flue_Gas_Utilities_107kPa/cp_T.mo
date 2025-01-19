within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa;
function cp_T "Specific heat capacity (J/kgK) of flue gas at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real cp_data[16] = {970.11, 970.87, 1010.49, 1048.16, 1084.70, 1119.50, 1151.61, 1180.45, 1205.89, 1228.10, 1247.38, 1264.10, 1278.61, 1291.24, 1302.28, 1311.96};
algorithm
    cp := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,cp_data,T);
end cp_T;