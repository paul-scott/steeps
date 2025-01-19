within SolarTherm.Media.Air;
package Air_CoolProp_1bar "Air, explicit in h"
	/* For a new medium, make a copy of this package and remove
		the "partial" keyword from the package definition above.
		The statement below extends from PartialMedium and sets some
		package constants. Provide values for these constants
		that are appropriate for your medium model. Note that other
		constants (such as nX, nXi) are automatically defined by
		definitions given in the base class Interfaces.PartialMedium"
	*/
	extends Modelica.Media.Interfaces.PartialMedium(
		ThermoStates = Modelica.Media.Interfaces.Choices.IndependentVariables.ph,
		final mediumName="Air",
		final substanceNames={"Air"},
		final singleState=false,
		final reducedX=true,
		final fixedX=true,
		Temperature(
			min=298.15,
			max=1700.0,
			start=800.0));
	//import SolarTherm.Media.MoltenSalt.MoltenSalt_utilities.*;
	import SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.*;
	// Provide medium constants here
	//constant SpecificHeatCapacity cp_const=123456
	// "Constant specific heat capacity at constant pressure";

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
		p,T as independent variables.
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
        h = state.h;
		d = rho_T(T);
		//h = h_Tf(T,0);
		u = h - p/d;
		MM = 0.02897;
		R = 8.3144/MM;
		state.p = p;
		T = T_h(h);
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

	redeclare function extends pressure "Return pressure"
	algorithm
		p := 100000.0;
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


	redeclare function extends dynamicViscosity "Return dynamic viscosity"
	algorithm
		eta := mu_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end dynamicViscosity;

	redeclare function extends thermalConductivity "Return thermal conductivity"
	algorithm
		lambda := k_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end thermalConductivity;

	redeclare function extends specificHeatCapacityCp
		"Return specific heat capacity at constant pressure"
	algorithm
		cp := cp_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end specificHeatCapacityCp;
/*
	redeclare function extends specificHeatCapacityCv
		"Return specific heat capacity at constant volume"
	algorithm
		cv := 0;
		annotation (Documentation(info="<html>

			</html>"));
	end specificHeatCapacityCv;

	redeclare function extends isentropicExponent "Return isentropic exponent"
		extends Modelica.Icons.Function;
	algorithm
		gamma := 1;
		annotation (Documentation(info="<html>

			</html>"));
	end isentropicExponent;

	redeclare function extends isentropicEnthalpy
		"Return isentropic enthalpy"
	algorithm
		h_is := 0; // To be completed
		annotation (Documentation(info="<html>

			</html>"));

	end isentropicEnthalpy;

	redeclare function extends velocityOfSound "Return velocity of sound"
		extends Modelica.Icons.Function;
	algorithm
		a := 0;
		annotation (Documentation(info="<html>

			</html>"));
	end velocityOfSound;

	redeclare function extends isobaricExpansionCoefficient
		"Return overall the isobaric expansion coefficient beta"
	algorithm
		beta := beta_T(T_h(state.h));
			annotation (Documentation(info="<html>

				</html>"));
	end isobaricExpansionCoefficient;

	redeclare function extends isothermalCompressibility
		"Return overall the isothermal compressibility factor"
	algorithm
		kappa := kappa_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end isothermalCompressibility;


	function enthalpyOfVaporization
		"Return vaporization enthalpy of condensing fluid"
		extends Modelica.Icons.Function;
		input ThermodynamicState state "Thermodynamic state record";
		output SpecificEnthalpy r0 "Vaporization enthalpy";
	algorithm
		r0 := h_fg_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end enthalpyOfVaporization;

	redeclare function extends density_derT_p
		"Return density derivative w.r.t. temperature at constant pressure"
	algorithm
		ddTp := drho_dT_T(T_h(state.h));
		annotation (Documentation(info="<html>

			</html>"));
	end density_derT_p;
*/
	annotation (Documentation(info= "<html><head></head><body><p><span style=\"font-family: Arial,sans-serif;\">Calculation of air properties at 1bar using CoolProp. The composition of air by mole percentage is assumed to be 78% N2, 21% O2 and 1% Ar. Properties of air are tabulated at 298.15K, then 300.0K to 1700K at 100K intervals.</span></p><p><font face=\"Arial, sans-serif\">Thermodynamic properties are defined by specific enthalpy and pressure. Pressure in this case is set to 1e5 Pa.</font></p><p><b><span style=\"font-family: Arial,sans-serif; color: #008000;\">Restriction</span></b></p>
<p><span style=\"font-family: Arial,sans-serif;\">The functions provided by this package shall be used inside of the restricted limits according to the referenced literature. </span></p>
<ul>
<li><b><span style=\"font-family: Arial,sans-serif;\">298.15 Kelvin ≤ T ≤ 1700.0 Kelvin </span></b></li>
<li><b><span style=\"font-family: Arial,sans-serif;\">explicit for pressure and specific enthalpy </span></b></li>
</ul>
<p><b><span style=\"font-family: Arial,sans-serif;\">References</span></b> </p>
<p style=\"margin-left: 30px;\"><font face=\"Arial\" size=\"4\"><span role=\"presentation\" dir=\"ltr\" style=\"left: 168.298px; top: 420.075px;\">Ian H. Bell, Jorrit Wronski, Sylvain Quoilin, and Vincent Lemort.</span><span role=\"presentation\" dir=\"ltr\" style=\"left: 736.745px; top: 420.075px;\"> </span><span role=\"presentation\" dir=\"ltr\" style=\"left: 757.575px; top: 420.075px;\">Pure and pseudo-pure&nbsp;</span></font><span style=\"font-family: Arial; font-size: large;\">fluid thermophysical property evaluation and the open-source thermophysical property li</span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 194.758px; top: 472.543px;\">brary coolprop.</span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 319.927px; top: 472.543px;\"> </span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 333.483px; top: 472.543px;\">Industrial &amp; Engineering Chemistry Research</span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 693.583px; top: 472.543px;\">, 53(6):2498–2508, 2014.</span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 903.348px; top: 472.543px;\"> </span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 916.972px; top: 472.543px;\">doi:&nbsp;</span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 194.758px; top: 498.777px;\">10.1021/ie4033999. URL</span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 393.119px; top: 498.777px;\"> </span><span role=\"presentation\" dir=\"ltr\" style=\"font-family: Arial; font-size: large; left: 397.284px; top: 498.777px;\">http://pubs.acs.org/doi/abs/10.1021/ie4033999</span></p></body></html>"));
end Air_CoolProp_1bar;