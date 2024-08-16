function setState_dTX
		"Return thermodynamic state from d, T, and X or Xi"
		extends Modelica.Icons.Function;
		input Density d "Pressure";
		input Temperature T "Specific entropy";
		input MassFraction X[:]=reference_X "Mass fractions";
		output ThermodynamicState state "Thermodynamic state record";
	algorithm
		state := ThermodynamicState(p=p_rho(d), h=h_T(T));
	end setState_dTX;