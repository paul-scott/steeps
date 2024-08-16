function extends specificHeatCapacityCv
		"Return specific heat capacity at constant volume"
	algorithm
		cv := cp_T(T_h(state.h)); // for a solid substance: cp is almost equal to cv
		annotation (Documentation(info="<html>

			</html>"));
	end specificHeatCapacityCv;