within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_107kPa;

function extends specificEnthalpy "Return specific enthalpy"
	algorithm
		h := state.h;
		annotation (Inline=true);
		end specificEnthalpy;