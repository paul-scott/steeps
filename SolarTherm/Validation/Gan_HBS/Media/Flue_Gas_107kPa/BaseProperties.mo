within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_107kPa;

model extends BaseProperties(final standardOrderComponents=true)
		"Base properties of medium"

	equation
        h = state.h;
		d = rho_T(T);
		//h = h_Tf(T,0);
		u = h - p/d;
		MM = 0.03276;
		R = 8.3144/MM;
		state.p = p;
		T = T_h(h);
	end BaseProperties;