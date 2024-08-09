within SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_5884Pa;

function extends pressure "Return pressure"
  algorithm
    p := 101325.0;
  annotation(
    Inline = true);
end pressure;