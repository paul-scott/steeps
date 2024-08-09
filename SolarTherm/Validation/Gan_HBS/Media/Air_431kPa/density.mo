within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa;

function extends density "Return density"
	algorithm
		d := rho_T(T_h(state.h));
		annotation (Inline=true);
	end density;