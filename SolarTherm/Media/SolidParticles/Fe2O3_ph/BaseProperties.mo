model extends BaseProperties(final standardOrderComponents=true)
		"Base properties of medium"

	equation
	T = T_h(h);
	p = state.p;
	d = rho_T(T);
	h = state.h;
	u = h - p / d;
	MM = 0.10313381;
	R = 8.3144 / MM;

	end BaseProperties;