within SolarTherm.Media.SolidParticles;

package Fe3O4_ph
	extends Modelica.Media.Interfaces.PartialMedium(
		ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph,
		final mediumName="Fe3O4",
		final substanceNames={"Fe3O4"},
		final singleState=false,
		final reducedX=true,
		final fixedX=true,
		Temperature(
			min=298.15,
			max=1812.0,
			start=50.0+273.15));

	import SolarTherm.Media.SolidParticles.Fe3O4_utilities.*;

	redeclare record extends ThermodynamicState
		"A selection of variables that uniquely defines the thermodynamic state"
		AbsolutePressure p "Absolute pressure of medium";
		SpecificEnthalpy h "Specific enthalpy";
		 annotation (Documentation(info="<html>

			</html>"));
	end ThermodynamicState;

	redeclare model extends BaseProperties(final standardOrderComponents=true)
		"Base properties of medium"

	equation
	T = T_h(h);
	p = state.p;
	d = rho_T(T);
	h = state.h;
	u = h - p / d;
	MM = 231.5326e-3;
	R = 8.3144 / MM;

	end BaseProperties;

	/* Provide implementations of the following optional properties.
		If not available, delete the corresponding function.
		The record "ThermodynamicState" contains the input arguments
		of all the function and is defined together with the used
		type definitions in PartialMedium. The record most often contains two of the
	variables "p, T, d, h" (e.g., medium.T)
	*/

	redeclare function setState_pTX
		"Return thermodynamic state from p, T, and X or Xi"
		extends Modelica.Icons.Function;
		input AbsolutePressure p "Pressure";
		input Temperature T "Temperature";
		input MassFraction X[:]=reference_X "Mass fractions";
		output ThermodynamicState state "Thermodynamic state record";
	algorithm
		state := ThermodynamicState(p=p, h=h_T(T));
	end setState_pTX;

	redeclare function setState_phX
		"Return thermodynamic state from p, h, and X or Xi"
		extends Modelica.Icons.Function;
		input AbsolutePressure p "Pressure";
		input SpecificEnthalpy h "Specific enthalpy";
		input MassFraction X[:]=reference_X "Mass fractions";
		output ThermodynamicState state "Thermodynamic state record";
	algorithm
		state := ThermodynamicState(p=p, h=h);
	end setState_phX;

	redeclare function setState_psX
		"Return thermodynamic state from p, s, and X or Xi"
		extends Modelica.Icons.Function;
		input AbsolutePressure p "Pressure";
		input SpecificEntropy s "Specific entropy";
		input MassFraction X[:]=reference_X "Mass fractions";
	output ThermodynamicState state "Thermodynamic state record";
	algorithm
		state := ThermodynamicState(p=p, h=h_s(s));
	end setState_psX;

	redeclare function setState_dTX
		"Return thermodynamic state from d, T, and X or Xi"
		extends Modelica.Icons.Function;
		input Density d "Pressure";
		input Temperature T "Specific entropy";
		input MassFraction X[:]=reference_X "Mass fractions";
		output ThermodynamicState state "Thermodynamic state record";
	algorithm
		state := ThermodynamicState(p=p_rho(d), h=h_T(T));
	end setState_dTX;

	redeclare function extends pressure "Return pressure"
	algorithm
		p := state.p;
		annotation (Inline=true);
	end pressure;

	redeclare function extends temperature "Return temperature"
	algorithm
		T := T_h(state.h);
		annotation (Inline=true);
	end temperature;

	redeclare function extends specificEnthalpy "Return specific enthalpy"
	algorithm
		h := state.h;
		annotation (Inline=true);
	end specificEnthalpy;

	redeclare function extends density "Return density"
	algorithm
		d := rho_T(T_h(state.h));
		annotation (Inline=true);
	end density;

	redeclare function extends specificInternalEnergy "Return specific internal energy"
	algorithm
		u := state.h - state.p / rho_T(T_h(state.h));
		annotation (Inline=true);
	end specificInternalEnergy;

	redeclare function extends thermalConductivity "Return thermal conductivity"
	algorithm
		lambda := k_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end thermalConductivity;

	redeclare function extends specificEntropy "Return specific entropy"
	algorithm
		s := s_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end specificEntropy;

	redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
	algorithm
		g := 0; //TODO: to be implemented.
		annotation (Documentation(info="<html>

			</html>"));	
	end specificGibbsEnergy;

	redeclare function extends specificHelmholtzEnergy "Return specific Helmholtz energy"
	algorithm
		f := 0; //TODO: to be implemented.
		annotation (Documentation(info="<html>

			</html>"));	
	end specificHelmholtzEnergy;

	redeclare function extends specificHeatCapacityCp
		"Return specific heat capacity at constant pressure"
	algorithm
		cp := cp_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end specificHeatCapacityCp;

	redeclare function extends specificHeatCapacityCv
		"Return specific heat capacity at constant volume"
	algorithm
		cv := cp_T(T_h(state.h)); // for a solid substance: cp is almost equal to cv
		annotation (Documentation(info="<html>

			</html>"));
	end specificHeatCapacityCv;

	redeclare function extends isentropicEnthalpy
		"Return isentropic enthalpy"
	algorithm
		h_is := 0;//TODO: to be implemented.
		annotation (Documentation(info="<html>

			</html>"));
	end isentropicEnthalpy;

	redeclare function extends density_derT_p
		"Return density derivative w.r.t. temperature at constant pressure"
	algorithm
		ddTp := 0;//TODO: to be implemented. drho_dT_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end density_derT_p;

	annotation(
    Documentation(info = "<html><head></head><body><p><span style=\"font-family: Arial,sans-serif;\">Calculation of thermo-physical properties for pure solid Fe3O4 in the temperature region of 298.15 K to 1500 K. Thermodynamic properties are explicit in terms of enthalpy and pressure. Properties are calculated via linear interpolation of data-tables with 10 K temperature intervals.
		</span></p><b><span style=\"font-family: Arial, sans-serif;\"><div><b style=\"font-family: 'DejaVu Sans Mono';\"><span style=\"font-family: Arial, sans-serif;\">Restrictions</span></b><div style=\"font-family: 'DejaVu Sans Mono'; font-weight: normal;\"><span style=\"font-family: Arial, sans-serif;\"><br></span></div><div style=\"font-family: 'DejaVu Sans Mono'; font-weight: normal;\"><span style=\"font-family: Arial, sans-serif;\">The functions provided by this package shall be used inside of the restricted limits according to the referenced literature.</span><ul><li><b><span style=\"font-family: Arial, sans-serif;\">298.15 K ≤ T ≤ 1500 K</span></b></li><li><b><span style=\"font-family: Arial, sans-serif;\">Explicit for pressure and enthalpy.</span></b></li></ul><div><font face=\"Arial, sans-serif\">The following sources were used:</font></div><div><ul><li><font face=\"Arial, sans-serif\"><b>cp, h, s: FactSage</b></font></li><li><font face=\"Arial, sans-serif\"><b>rho, k: Beygelzimer</b></font></li></ul><div><font face=\"Arial, sans-serif\">The following phases were considered in this model:</font></div></div><div><ul><li><font face=\"Arial, sans-serif\"><b>S1: 298.15 K - 848.12 K</b></font></li><li><font face=\"Arial, sans-serif\"><b>S2: 848.12 - 1500.0 K</b></font></li></ul><div><font face=\"Arial, sans-serif\">Specific enthalpy,&nbsp;</font><i style=\"font-family: Arial, sans-serif;\">h</i><span style=\"font-family: Arial, sans-serif;\">(</span><i style=\"font-family: Arial, sans-serif;\">T</i><span style=\"font-family: Arial, sans-serif;\">)</span><span style=\"font-family: Arial, sans-serif;\">&nbsp;[J/kg] measured with respect to phase S1 at 298.15 K.</span></div></div><div><font face=\"Arial, sans-serif\">Specific absolute entropy,&nbsp;<i>s</i>(<i>T</i>) [J/kgK] provided as an absolute measurement.</font></div><p><b><span style=\"font-family: Arial, sans-serif;\">References</span></b></p><p style=\"margin-left: 30px;\">C. W. Bale, E. Bélisle, P. Chartrand, S. A. Decterov, G. Eriksson, A.E. Gheribi, K. Hack, I. H. Jung, Y. B. Kang, J. Melançon, A. D. Pelton, S. Petersen, C. Robelin. J. Sangster, P. Spencer and M-A. Van Ende, FactSage Thermochemical Software and Databases - 2010 - 2016, Calphad, vol. 54, pp 35-53, 2016 &lt;www.factsage.com&gt;</p><p style=\"margin-left: 30px;\">Emmanuil Beygelzimer and Yan Beygelzimer. Thermal conductivity of oxide scale and its com- ponents in the range from 0 C to 1300 C: Generalized estimates with account for movability of phase transitions, October 2021.</p><p style=\"margin-left: 30px;\">Emmanuil Beygelzimer and Yan Beygelzimer. Generalized estimates for the density of oxide scale in the range from 0 C to 1300 C, October 2021.</p></div></div></span></b><div>
		</div></body></html>"));
end Fe3O4_ph;
