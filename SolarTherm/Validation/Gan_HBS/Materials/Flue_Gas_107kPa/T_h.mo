within SolarTherm.Validation.Gan_HBS.Materials.Flue_Gas_107kPa;

function T_h "Temperature interpolated from enthalpy"
  input SI.SpecificEnthalpy h;
  output SI.Temperature T;
algorithm
  T := SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.T_h(h);
end T_h;
