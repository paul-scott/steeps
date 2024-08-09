within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_107kPa;

function extends temperature "Return temperature"
	algorithm
		T := T_h(state.h);
		annotation (Inline=true);
	end temperature;