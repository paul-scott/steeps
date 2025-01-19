within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa_Utilities;
function mu_T "Dynamic viscosity (Ns/m2) of air at 431.493kPa as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.DynamicViscosity mu "Dynamic Viscosity (Ns/m2)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real mu_data[16] = {0.00001763, 0.00001772, 0.00002220, 0.00002644, 0.00003043, 0.00003417, 0.00003767, 0.00004092, 0.00004393, 0.00004669, 0.00004920, 0.00005147, 0.00005349, 0.00005526, 0.00005679, 0.00005807};
algorithm
    mu := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,mu_data,T);
end mu_T;