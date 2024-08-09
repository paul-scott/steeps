within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa;

function extends temperature "Return temperature"
	algorithm
		T := T_h(state.h);
		annotation (Inline=true);
	end temperature;