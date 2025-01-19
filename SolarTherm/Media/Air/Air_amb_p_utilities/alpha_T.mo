within SolarTherm.Media.Air.Air_amb_p_utilities;
function alpha_T "Thermal diffusivity (m2/s) of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalDiffusivity alpha "Thermal Diffusivity (m2/s)";
protected
    Real T_data[35] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000};
    Real alpha_data[35] = {0.00000254,0.00000584,0.0000103,0.0000159,0.0000225,0.0000299,0.0000383,0.0000472,0.0000567,0.0000667,0.0000769,0.0000873,0.000098,0.000109,0.00012,0.000131,0.000143,0.000155,0.000168,0.000195,0.000224,0.000257,0.000303,0.00035,0.00039,0.000435,0.000482,0.000534,0.000589,0.000646,0.000714,0.000783,0.000869,0.00096,0.00157};
algorithm
    alpha := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,alpha_data,T);
end alpha_T;