within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa;
function mu_T "Dynamic viscosity (Ns/m2) of flue gas at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.DynamicViscosity mu "Dynamic Viscosity (Ns/m2)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real mu_data[16] = {0.00001668, 0.00001677, 0.00002120, 0.00002521, 0.00002889, 0.00003230, 0.00003549, 0.00003851, 0.00004138, 0.00004413, 0.00004678, 0.00004935, 0.00005185, 0.00005428, 0.00005666, 0.00005899};
algorithm
    mu := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,mu_data,T);
end mu_T;