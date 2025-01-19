within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa_Utilities;
function k_T "Thermal conductivity (W/mK) of air at 431.493kPa as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity k "Thermal Conductivity (W/mK)";
protected
    Real T_data[16] = {298.15, 300.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1100.00, 1200.00, 1300.00, 1400.00, 1500.00, 1600.00, 1700.00};
    Real k_data[16] = {0.02562, 0.02575, 0.03243, 0.03872, 0.04460, 0.05008, 0.05517, 0.05985, 0.06413, 0.06802, 0.07150, 0.07458, 0.07727, 0.07955, 0.08143, 0.08292};
algorithm
    k := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
end k_T;