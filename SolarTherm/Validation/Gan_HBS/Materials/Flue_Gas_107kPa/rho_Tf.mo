within SolarTherm.Validation.Gan_HBS.Materials.Flue_Gas_107kPa;

function rho_Tf "Density of air vs Temperature"
  input SI.Temperature T;
  input Real f = 0 "Liquid mass melt fraction (No effect on result)";
  output SI.Density rho;
algorithm
  rho := SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.rho_T(T);
end rho_Tf;
