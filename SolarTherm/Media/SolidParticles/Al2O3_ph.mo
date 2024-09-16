within SolarTherm.Media.SolidParticles;

package Al2O3_ph
	/* The statement below extends from PartialMedium and sets some
		package constants. Provide values for these constants
		that are appropriate for your medium model. Note that other
		constants (such as nX, nXi) are automatically defined by
		definitions given in the base class Interfaces.PartialMedium"
	*/

	extends Modelica.Media.Interfaces.PartialMedium(
		ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph,
		final mediumName="Al2O3",
		final substanceNames={"Al2O3"},
		final singleState=false,
		final reducedX=true,
		final fixedX=true,
		Temperature(
			min=293.15,
			max=1812.0,
			start=900+273.15));

	import SolarTherm.Media.SolidParticles.Al2O3_utilities.*;

	// Provide medium constants here
	//constant SpecificHeatCapacity cp_const=123456
		//"Constant specific heat capacity at constant pressure";

	/* The vector substanceNames is mandatory, as the number of
		substances is determined based on its size. Here we assume
		a single-component medium.
		singleState is true if u and d do not depend on pressure, but only
		on a thermal variable (temperature or enthalpy). Otherwise, set it
		to false.
		For a single-substance medium, just set reducedX and fixedX to true, and there's
		no need to bother about medium compositions at all. Otherwise, set
		final reducedX = true if the medium model has nS-1 independent mass
		fraction, or reducedX = false if the medium model has nS independent
		mass fractions (nS = number of substances).
		If a mixture has a fixed composition set fixedX=true, otherwise false.
		The modifiers for reducedX and fixedX should normally be final
		since the other equations are based on these values.

		It is also possible to redeclare the min, max, and start attributes of
		Medium types, defined in the base class Interfaces.PartialMedium
		(the example of Temperature is shown here). Min and max attributes
		should be set in accordance to the limits of validity of the medium
		model, while the start attribute should be a reasonable default value
		for the initialization of nonlinear solver iterations */

		/* Provide an implementation of model BaseProperties,
		that is defined in PartialMedium. Select two independent
		variables from p, T, d, u, h. The other independent
		variables are the mass fractions "Xi", if there is more
		than one substance. Provide 3 equations to obtain the remaining
		variables as functions of the independent variables.
		It is also necessary to provide two additional equations to set
		the gas constant R and the molar mass MM of the medium.
		Finally, the thermodynamic state vector, defined in the base class
		Interfaces.PartialMedium.BaseProperties, should be set, according to
		its definition (see ThermodynamicState below).
		The computation of vector X[nX] from Xi[nXi] is already included in
		the base class Interfaces.PartialMedium.BaseProperties, so it should not
		be repeated here.
		The code fragment below is for a single-substance medium with
		p and h as independent variables.
	*/

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
	MM = 101.961e-3;
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

	annotation (Documentation(info= "<html><head></head><body><p><span style=\"font-family: Arial,sans-serif;\">Calculation of thermo-physical properties for pure solid Al2O3 in the temperature region of 298.15 K to 1500 K at its most stable phase at 1 bar. Thermodynamic properties are explicit in terms of enthalpy and pressure. Properties are calculated via linear interpolation of data-tables with 10 K temperature intervals.
		</span></p><p><b><span style=\"font-family: Arial, sans-serif;\"></span></b></p><div><b style=\"font-family: 'DejaVu Sans Mono';\"><span style=\"font-family: Arial, sans-serif;\">Restrictions</span></b><div style=\"font-family: 'DejaVu Sans Mono'; font-weight: normal;\"><span style=\"font-family: Arial, sans-serif;\"><br></span></div><div style=\"font-family: 'DejaVu Sans Mono'; font-weight: normal;\"><span style=\"font-family: Arial, sans-serif;\">The functions provided by this package shall be used inside of the restricted limits according to the referenced literature.</span><ul><li><b><span style=\"font-family: Arial, sans-serif;\">298.15 K ≤ T ≤ 1500 K</span></b></li><li><b><span style=\"font-family: Arial, sans-serif;\">Explicit for pressure and enthalpy.</span></b></li></ul><div><font face=\"Arial, sans-serif\">The following sources were used:</font></div><div><ul><li><font face=\"Arial, sans-serif\"><b>cp, h, s: FactSage</b></font></li><li><font face=\"Arial, sans-serif\"><b>rho: Touloukian Vol 4</b></font></li><li><font face=\"Arial, sans-serif\"><b>k: CRC Materials Science and Engineering Handbook</b></font></li></ul><div><font face=\"Arial, sans-serif\">The following phases were considered in this model:</font></div></div><div><ul><li><font face=\"Arial, sans-serif\"><b>S4 (Corundum alpha): 298.15 K - 1500.0 K</b></font></li></ul><div><font face=\"Arial, sans-serif\">Specific enthalpy,&nbsp;</font><i style=\"font-family: Arial, sans-serif;\">h</i><span style=\"font-family: Arial, sans-serif;\">(</span><i style=\"font-family: Arial, sans-serif;\">T</i><span style=\"font-family: Arial, sans-serif;\">)</span><span style=\"font-family: Arial, sans-serif;\">&nbsp;[J/kg] measured with respect to phase S4 at 298.15 K.</span></div></div><div><font face=\"Arial, sans-serif\">Specific absolute entropy,&nbsp;<i>s</i>(<i>T</i>) [J/kgK] provided as an absolute measurement.</font></div><p><b><span style=\"font-family: Arial, sans-serif;\">References</span></b></p><p style=\"margin-left: 30px;\">C. W. Bale, E. Bélisle, P. Chartrand, S. A. Decterov, G. Eriksson, A.E. Gheribi, K. Hack, I. H. Jung, Y. B. Kang, J. Melançon, A. D. Pelton, S. Petersen, C. Robelin. J. Sangster, P. Spencer and M-A. Van Ende, FactSage Thermochemical Software and Databases - 2010 - 2016, Calphad, vol. 54, pp 35-53, 2016 &lt;www.factsage.com&gt;</p><p style=\"margin-left: 30px;\">Y. S. Touloukian 1967, Thermophysical Properties of High Temperature Solid Materials Volume 4: Oxides and their solutions ad mixtures.</p><p style=\"margin-left: 30px;\">J. F. Shackelford, W. Alexander 2001, Materials Science and Engineering Handbook Third Edition. CRC Press.</p></div></div><div>
		</div></body></html>"));
end Al2O3_ph;