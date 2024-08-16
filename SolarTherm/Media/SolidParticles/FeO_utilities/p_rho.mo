within SolarTherm.Media.SolidParticles.FeO_utilities;
function p_rho "Pressure of FeO as a function of density. Just set this to 1atm"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Density rho "Density (kg/m3)";
	output Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure (Pa)";

algorithm
	p := 101325.0; // Constant pressure.
end p_rho;