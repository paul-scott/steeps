within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa;

function setState_pTX
		"Return thermodynamic state from p, T, and X or Xi"
		extends Modelica.Icons.Function;
		input AbsolutePressure p "Pressure";
		input Temperature T "Temperature";
		input MassFraction X[:]=reference_X "Mass fractions";
		output ThermodynamicState state "Thermodynamic state record";
	algorithm
		state := ThermodynamicState(p=p, h=h_T(T));
	end setState_pTX;