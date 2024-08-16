within SolarTherm.Media.SolidParticles.Fe_utilities;
function lamda_T "Thermal conductivity of Fe as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity lamda "Thermal conductivity (W/mK)";
protected
    constant Modelica.SIunits.Temperature T_data[20] = {100.00, 150.00, 200.00, 250.00, 273.20, 300.00, 350.00, 400.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00, 1059.00, 1183.00, 1200.00, 1300.00, 1400.00, 1500.00};
    constant Modelica.SIunits.ThermalConductivity k_data[20] = {132.00, 104.00, 94.00, 86.50, 83.50, 80.30, 74.40, 69.40, 61.30, 54.70, 48.70, 43.30, 38.00, 32.60, 29.60, 29.90, 28.20, 29.90, 30.90, 31.80};
algorithm
	lamda := Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
end lamda_T;