within SolarTherm.Media.SolidParticles;

package Fe_utilities
  extends Modelica.Icons.UtilitiesPackage;
  import SolarTherm.Models.Chemistry.Property_Tables.Fe.*;
  
  function cp_T "Specific heat capacity of as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificHeatCapacity cp "Specific heat capacity (J/kgK)";
  protected
	constant Modelica.SIunits.Temperature T_data[122] = T_table;
	constant Modelica.SIunits.SpecificHeatCapacity cp_data[122] = cp_table;

  algorithm
	cp := Modelica.Math.Vectors.interpolate(T_data,cp_data,T);
  end cp_T;
  
  function h_T "Specific enthalpyas a function of temperature with respect to 298.15K"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy (J/kgK)";
  protected
	constant Modelica.SIunits.Temperature T_data[122] = T_table;
	constant Modelica.SIunits.SpecificEnthalpy h_data[122] = h_table;
  algorithm
	h := Modelica.Math.Vectors.interpolate(T_data,h_data,T);
  end h_T;

  function s_T "Specific absolute entropy obtained from NIST JANAF tables."
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.SpecificEntropy s "Specific absolute entropy (J/kgK)";
  protected
	constant Modelica.SIunits.Temperature T_data[122] = T_table;
	constant Modelica.SIunits.SpecificEntropy s_data[122] = s_table;
  algorithm
	s := Modelica.Math.Vectors.interpolate(T_data,s_data,T);
  end s_T;
  
  function T_h "Temperature as a function of specific enthalpy"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy (J/kg)";
	output Modelica.SIunits.Temperature T "Temperature (K)";
  protected
	constant Modelica.SIunits.Temperature T_data[122] = T_table;
	constant Modelica.SIunits.SpecificEnthalpy h_data[122] = h_table;
  algorithm
	T := Modelica.Math.Vectors.interpolate(h_data,T_data,h);
  end T_h;
  
  function T_s "Temperature as a function of specific absolute entropy."
	extends Modelica.Icons.Function;
	input Modelica.SIunits.SpecificEntropy s "Specific absolute entropy";
	output Modelica.SIunits.Temperature T "Temperature";
  protected
	constant Modelica.SIunits.Temperature T_data[122] = T_table;
	constant Modelica.SIunits.SpecificEntropy s_data[122] = s_table;
  algorithm
	T := Modelica.Math.Vectors.interpolate(s_data,T_data,s);
  end T_s;
  
  function rho_T "Density as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.Density rho "Density (kg/m3)";
  protected
	constant Modelica.SIunits.Temperature T_data[122] = T_table;
	constant Modelica.SIunits.Density rho_data[122] = rho_table;
  algorithm
	rho := Modelica.Math.Vectors.interpolate(T_data,rho_data,T);
  end rho_T;

  function k_T "Thermal conductivity as a function of temperature"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Temperature T "Temperature (K)";
	output Modelica.SIunits.ThermalConductivity k "Thermal conductivity (W/mK)";
  protected
    constant Modelica.SIunits.Temperature T_data[122] = T_table;
    constant Modelica.SIunits.ThermalConductivity k_data[122] = k_table;
  algorithm
	k := Utilities.Interpolation.Interpolate1D(T_data,k_data,T);
  end k_T;
  
  function p_rho "Pressure a function of density. Just set this to 1 atm"
	extends Modelica.Icons.Function;
	input Modelica.SIunits.Density rho "Density (kg/m3)";
	output Modelica.Media.Interfaces.Types.AbsolutePressure p "Pressure (Pa)";
  algorithm
	p := 101325.0;
  end p_rho;

  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Fe_utilities;