within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa;

function extends specificEnthalpy "Return specific enthalpy"
	algorithm
		h := state.h;
		annotation (Inline=true);
		end specificEnthalpy;