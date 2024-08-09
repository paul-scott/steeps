within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities;
function rho_T "Density (kg/m3) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.Density rho "Density (kg/m3)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real rho_data[16] = {0.0743300, 0.0738716, 0.0554034, 0.0443225, 0.0369353, 0.0316588, 0.0277014, 0.0246235, 0.0221611, 0.0201464, 0.0184676, 0.0170470, 0.0158293, 0.0147740, 0.0138506, 0.0130359};
algorithm
    rho := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,rho_data,T);
end rho_T;
