within SolarTherm.Models.PowerBlocks.sCO2Cycle.OnDesign;
model Turbine
	extends SolarTherm.Media.CO2.PropCO2;
	replaceable package MedPB = SolarTherm.Media.CO2.CO2_ph;
	outer Modelica.Fluid.System system;
	Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = MedPB) annotation(
		Placement(visible = true, transformation(origin = {-32, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-38, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
	Modelica.Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = MedPB) annotation(
		Placement(visible = true, transformation(origin = {60, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {60, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
	parameter Real eta_turb = 0.9 "isentropic efficiency of the turbine";
	parameter Real PR = 2.313 "Pressure ratio";
	parameter Modelica.SIunits.ThermodynamicTemperature T_amb = 273.15 + 40 "Outlet temperature in Kelvin";
	parameter Boolean is_second_turbine = false;
	parameter Modelica.SIunits.MassFlowRate flowGuess=100;
	
	// This last parameter, when put to true, adds a mass_flow equality equation. Reason is, you need one and only one component in the cycle without mass flow or else you will have too much equation and circular equality issues. therefore, when adding turbines, this equation should be present (or when it's an open-cycle).
	MedPB.ThermodynamicState state_a "thermodynamic state at the entrance";
	MedPB.ThermodynamicState state_isen "thermodynamic state at the end of the isentropic decompression";
	MedPB.ThermodynamicState state_b "thermodynamic state at the end of the real decompresssion";
	Modelica.SIunits.Power W_turb "Outlet power";
	SolarTherm.Types.SpecificExergy ex_d "destroyed exergy";
	Modelica.SIunits.SpecificEntropy s_entrance " entropy at the entrance of the turbine";
	Modelica.SIunits.Area A_nozzle (start=10^(-3));
	parameter Real N_shaft = 75000 * 5 / 6 * 2 * 3.14159 / 60;
	Real d_outlet;
	Real C_spouting (start=500);
	Real diam_turb;
	Real tipSpeed (start=400);
	protected
	Modelica.SIunits.MassFlowRate mStart (start=flowGuess);
	equation
	state_a = MedPB.setState_phX(port_a.p, inStream(port_a.h_outflow));
	s_entrance = MedPB.specificEntropy(state_a);
	state_isen = MedPB.setState_psX(state_a.p / PR, s_entrance);
	state_b = MedPB.setState_phX(state_a.p / PR, state_a.h + (state_isen.h - state_a.h) * eta_turb);
	port_b.p = state_b.p;
	port_b.h_outflow = state_b.h;
	W_turb = port_a.m_flow * (state_b.h - state_a.h);
	port_a.m_flow + port_b.m_flow = 0;
	mStart= port_a.m_flow;
	ex_d = W_turb + port_a.m_flow * (state_a.h - T_amb * MedPB.specificEntropy(state_a)) + port_b.m_flow * (state_b.h - T_amb * MedPB.specificEntropy(state_b));
	d_outlet = MedPB.density(state_b);
	port_a.m_flow = C_spouting * A_nozzle * d_outlet;
	C_spouting ^ 2 = 2 * (state_a.h - state_isen.h);
	tipSpeed = N_shaft * diam_turb / 2;
	tipSpeed / C_spouting = 0.707;
	port_a.h_outflow = inStream(port_b.h_outflow);
	annotation(
		Diagram(graphics = {Text(origin = {-36, -28}, extent = {{18, 80}, {78, 16}}, textString = "TURBINE"), Polygon(origin = {15, 20}, points = {{-35, 44}, {-35, -52}, {35, -68}, {35, 68}, {-35, 44}, {35, 68}, {-35, 44}})}, coordinateSystem(initialScale = 0.1)),
		Icon(graphics = {Text(origin = {-10, 26}, extent = {{-10, 12}, {52, -34}}, textString = "TURBINE"), Ellipse(extent = {{56, 58}, {56, 58}}, endAngle = 360), Polygon(origin = {11, 17}, points = {{-37, 49}, {-37, -51}, {37, -71}, {37, 71}, {-37, 49}})}, coordinateSystem(initialScale = 0.1)),
		Documentation(info = "<html>
			<p>&nbsp;</p>
			<div>0D Model of a sCO2 turbine.</div>
			<div>The model is based on the thesis of J. Dyreby. The diameter of the turbine is calculated in order to maximize the efficiency at the design point.</div>
			<div></div>
			<div>The on-design calculation of the cycle gives us as results to be integrated in the off-design power blocks:</div>
			<div>
			<ul>
			<li>Diameter of the turbine</li>
			<li>Area of the nozzle</li>
			</ul>
			</div>
			<p>J. J. Dyreby, &laquo;&nbsp;Modeling the supercritical
			carbon dioxide Brayton cycle with recompression&nbsp;&raquo;, The University of
			Wisconsin-Madison, 2014. Available at https://sel.me.wisc.edu/publications-theses.shtml</p>
			<p>&nbsp;</p>
			</html>"));
	end Turbine;
