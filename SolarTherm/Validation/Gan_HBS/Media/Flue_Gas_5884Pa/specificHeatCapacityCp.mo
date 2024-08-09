within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_5884Pa;

function extends specificHeatCapacityCp "Return specific heat capacity at constant pressure"
  algorithm
    cp := cp_T(T_h(state.h));
  annotation(
    Documentation(info = "<html>

			</html>"));
end specificHeatCapacityCp;