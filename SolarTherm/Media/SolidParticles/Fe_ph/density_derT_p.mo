function extends density_derT_p
		"Return density derivative w.r.t. temperature at constant pressure"
	algorithm
		ddTp := drho_dT_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end density_derT_p;