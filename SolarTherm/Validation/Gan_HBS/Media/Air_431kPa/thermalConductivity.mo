within SolarTherm.Validation.Gan_HBS.Media.Air_431kPa;

function extends thermalConductivity "Return thermal conductivity"
	algorithm
		lambda := k_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end thermalConductivity;