within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa;

function extends dynamicViscosity "Return dynamic viscosity"
	algorithm
		eta := mu_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end dynamicViscosity;