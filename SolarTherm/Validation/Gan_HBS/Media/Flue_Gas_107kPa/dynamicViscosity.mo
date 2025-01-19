within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_107kPa;

function extends dynamicViscosity "Return dynamic viscosity"
	algorithm
		eta := mu_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end dynamicViscosity;