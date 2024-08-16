within SolarTherm.Materials;
package Fe2O3
  extends SolarTherm.Materials.PartialMaterial(MM = 159.7e-3, T_melt = 1539.0+273.15, cost = 0.2);
  import SolarTherm.Utilities.Interpolation.Interpolate1D;
  
  redeclare model State "A model which calculates state and properties"
	
	SI.SpecificEnthalpy h "Specific Enthalpy wrt 298.15K (J/kg)";
	SI.Temperature T "Temperature (K)";
	Real f "Liquid Mass Fraction";
	SI.Density rho "Density (kg/m3)";
	SI.ThermalConductivity k "Thermal conductivity (W/mK)";

	
  equation
	T = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h);
    rho = SolarTherm.Media.SolidParticles.Fe2O3_utilities.rho_T(T);
    k = SolarTherm.Media.SolidParticles.Fe2O3_utilities.k_T(T); //some genius spelt it like that, not Zeb
    f = 0.0; //Never melts
  end State;
  
  function Tf_h "Find temperature and liquid fraction from temperature"
    input SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
    output SI.Temperature T "Absoulte temperature (K)";
    output Real f "mass liquid fraction";
  algorithm
    T := SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h);
    f := 0.0;
    //T := Interpolate1D(h_data,T_data,h);
    //f := 0.0;
  end Tf_h;

  redeclare function h_Tf "find specific enthalpy from Temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.SpecificEnthalpy h "Specific Enthalpy (J/kg)";
  algorithm
    h := SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T);
    //h := Interpolate1D(T_data,h_data,T);
  end h_Tf;
  
  function k_Tf "find thermal conductivity from temperature"
    input SI.Temperature T;
    input Real f;
    output SI.ThermalConductivity k;
  algorithm
    k := SolarTherm.Media.SolidParticles.Fe2O3_utilities.k_T(T); //because some genius decided to spell it that way
    //k := Interpolate1D(T_data,k_data,T);
  end k_Tf;
    
  redeclare function rho_Tf "find density from temperature"
    input SI.Temperature T "Absolute temperature (K)";
    input Real f "Liquid mass fraction";
    output SI.Density rho "Density (kg/m3)";
  algorithm
    rho := SolarTherm.Media.SolidParticles.Fe2O3_utilities.rho_T(T);
    //rho := Interpolate1D(T_data,rho_data,T);
  end rho_Tf;
end Fe2O3;