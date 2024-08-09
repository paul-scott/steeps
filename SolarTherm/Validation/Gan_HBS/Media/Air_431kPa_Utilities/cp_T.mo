within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa_Utilities;
function cp_T "Specific heat capacity (J/kgK) of air at 431.493kPa at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real cp_data[16] = {1010.98, 1011.13, 1022.76, 1040.02, 1060.66, 1082.95, 1105.54, 1127.47, 1148.06, 1166.91, 1183.80, 1198.68, 1211.64, 1222.83, 1232.47, 1240.83};
algorithm
    cp := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,cp_data,T);
end cp_T;