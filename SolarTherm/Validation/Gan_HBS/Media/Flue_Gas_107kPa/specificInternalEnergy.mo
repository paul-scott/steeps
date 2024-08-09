within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_107kPa;

function extends specificInternalEnergy "Return specific internal energy"
	algorithm
		u := state.h - state.p / rho_T(T_h(state.h));
		annotation (Inline=true);
	end specificInternalEnergy;