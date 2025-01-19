within SolarTherm.Validation.Gan_HBS.Materials.Flue_Gas_107kPa;

model State "A model which calculates state and properties"
    parameter SI.SpecificEnthalpy h_start;
    SI.Density rho;
    SI.SpecificEnthalpy h (start=h_start);
    SI.Temperature T;
    SI.DynamicViscosity mu;
    SI.SpecificHeatCapacity cp;
    SI.ThermalConductivity k;
  equation
    T = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.T_h(h);
    rho = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.rho_T(T);
    cp = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.cp_T(T);
    mu = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.mu_T(T);
    k = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_Utilities_107kPa.k_T(T);
end State;
