function extends thermalConductivity "Return thermal conductivity"
	algorithm
		lambda := lamda_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end thermalConductivity;