within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa_Utilities;
function rho_T "Density (kg/m3) of air at 431.493kPa as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.Density rho "Density (kg/m3)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real rho_data[16] = {5.0343712, 5.0033252, 3.7524750, 3.0019710, 2.5016374, 2.1442576, 1.8762234, 1.6677527, 1.5009764, 1.3645233, 1.2508124, 1.1545956, 1.0721242, 1.0006489, 0.9381081, 0.8829251
};
algorithm
    rho := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,rho_data,T);
end rho_T;