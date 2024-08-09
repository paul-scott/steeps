within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities;
function mu_T "Dynamic viscosity (Ns/m2) of flue gas at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.DynamicViscosity mu "Dynamic Viscosity (Ns/m2)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real mu_data[16] = {0.00001612, 0.00001621, 0.00002061, 0.00002475, 0.00002862, 0.00003222, 0.00003555, 0.00003861, 0.00004140, 0.00004392, 0.00004617, 0.00004815, 0.00004986, 0.00005130, 0.00005247, 0.00005337};
algorithm
    mu := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,mu_data,T);
end mu_T;
