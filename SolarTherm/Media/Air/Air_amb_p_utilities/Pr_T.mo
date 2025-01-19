within SolarTherm.Media.Air.Air_amb_p_utilities;
function Pr_T "Prandtl number of air at ambient pressure as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalDiffusivity Pr "Prandtl Number";
protected
    Real T_data[35] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000};
    Real Pr_data[35] = {0.786,0.758,0.737,0.72,0.707,0.7,0.69,0.686,0.684,0.683,0.685,0.69,0.695,0.702,0.709,0.716,0.72,0.723,0.726,0.728,0.728,0.719,0.703,0.685,0.688,0.685,0.683,0.677,0.672,0.667,0.655,0.647,0.63,0.613,0.536};
algorithm
    Pr := SolarTherm.Utilities.Interpolation.Interpolate1D(T_data,Pr_data,T);
end Pr_T;