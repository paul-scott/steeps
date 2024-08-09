within SolarTherm.Validation.Gan_HBS.Materials.Flue_Gas_107kPa;

function h_Tf "Specific enthalpy of air vs Temperature"
  input SI.Temperature T "Absolute temperature (K)";
  input Real f = 0 "Liquid mass melt fraction (No effect on result)";
  output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
algorithm
  h := SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.h_T(T);
end h_Tf;
