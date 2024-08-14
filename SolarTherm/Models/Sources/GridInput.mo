within SolarTherm.Models.Sources;
model GridInput
	extends SolarTherm.Icons.GridInput;
	parameter String pv_file "File with the reference PV farm output";
	parameter String wind_file "File with the reference Wind farm output";
	parameter Modelica.SIunits.Power P_elec_min = 25e6 "Minimum power input";
	parameter Modelica.SIunits.Power P_elec_max = 100e6 "Maximum power input";
	parameter Modelica.SIunits.Efficiency pv_fraction = 0.5 "PV fraction of renewable input at design";
	parameter Modelica.SIunits.Power P_elec_pv_ref_size = 50e6 "PV farm reference size";
	parameter Modelica.SIunits.Power P_elec_wind_ref_size = 50e6 "Wind farm reference size";
	parameter Real renewable_multiple = 2 "Oversizing factor with respect to the heat input";
	final parameter Modelica.SIunits.Power P_elec_max_pv = pv_fraction*P_elec_max "Maximum PV capacity";
	final parameter Modelica.SIunits.Power P_elec_max_wind = (1-pv_fraction)*P_elec_max "Maximum Wind farm capacity";

	Modelica.Blocks.Sources.CombiTimeTable P_elec_ref_pv(
		fileName = pv_file, 
		tableName = "Power", 
		tableOnFile = true,
		smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
	Modelica.Blocks.Sources.CombiTimeTable P_elec_rec_wind(
		fileName = wind_file, 
		tableName = "Power", 
		tableOnFile = true,
		smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative);

	Modelica.Blocks.Sources.RealExpression P_elec_net_switch1(y = P_elec_net) annotation(
		Placement(visible = true, transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Sources.RealExpression P_elec_curtail_switch1(y = min(P_schedule,P_elec_net)) annotation(
		Placement(visible = true, transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Logical.Switch switch1 annotation(
		Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.BooleanInput curtail annotation(
		Placement(visible = true, 
		transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), 
		iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealInput P_schedule annotation(
		Placement(visible = true, 
		transformation(origin = {-100, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), 
		iconTransformation(origin = {-100, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.RealOutput electricity annotation(
		Placement(visible = true, 
		transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), 
		iconTransformation(origin = {100, 2.66454e-15}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	Modelica.Blocks.Interfaces.BooleanOutput on_renewable annotation (
		Placement(transformation(origin={0,-114}, extent={{20,-20},{-20,20}},rotation=90),
		iconTransformation(origin = {0,-114}, extent = {{10,-10},{-10,10}}, rotation = 90)));

	Modelica.SIunits.Power P_elec_net;
	Modelica.SIunits.Power P_elec_raw;
	Modelica.SIunits.Power P_elec_raw_pv;
	Modelica.SIunits.Power P_elec_raw_wind;

initial equation
	on_renewable = P_elec_raw > P_elec_min;

equation
	P_elec_raw_pv = P_elec_max_pv / P_elec_pv_ref_size * P_elec_ref_pv.y[1];
	P_elec_raw_wind = P_elec_max_wind / P_elec_wind_ref_size * P_elec_rec_wind.y[1];
	
	P_elec_raw = min(P_elec_raw_pv + P_elec_raw_wind, P_elec_max);
	
	on_renewable = P_elec_raw > P_elec_min;
	
	P_elec_net = if on_renewable then P_elec_raw else 0.0;

	connect(P_elec_curtail_switch1.y, switch1.u1) annotation(
		Line(points = {{-39, -30}, {-25.5, -30}, {-25.5, -8}, {-12, -8}}, color = {0, 0, 127}));
	connect(curtail, switch1.u2) annotation(
	Line(points = {{-100, 0}, {-12, 0}}, color = {255, 0, 255}));
	connect(P_elec_net_switch1.y, switch1.u3) annotation(
		Line(points = {{-38, 30}, {-26, 30}, {-26, 8}, {-12, 8}}, color = {0, 0, 127}));
	connect(switch1.y, electricity) annotation(
		Line(points = {{12, 0}, {100, 0}}, color = {0, 0, 127}));

annotation(
	Documentation(info="<html>
	<p>
	<b>GridInput</b> models the input of electrical power from renewable sources, such as photovoltaic (PV) and wind farms, into an electrical grid. The model calculates the net electrical power input based on reference power data, curtailment, and user-defined parameters.
	</p>
	<p>
	The <b>GridInput</b> model has the following connectors and parameters:
	</p>
	<ul>
	<li> Parameters:
	<ul>
	<li> <b>pv_file</b>: File with the reference PV farm output.</li>
	<li> <b>wind_file</b>: File with the reference Wind farm output.</li>
	<li> <b>P_elec_min</b>: Minimum power input, in Watts. Default: 25e6 W.</li>
	<li> <b>P_elec_max</b>: Maximum power input, in Watts. Default: 100e6 W.</li>
	<li> <b>pv_fraction</b>: PV fraction of renewable input at design. Default: 0.5.</li>
	<li> <b>P_elec_pv_ref_size</b>: PV farm reference size, in Watts. Default: 50e6 W.</li>
	<li> <b>P_elec_wind_ref_size</b>: Wind farm reference size, in Watts. Default: 50e6 W.</li>
	</ul>
	</li>
	<li> Connectors:
	<ul>
	<li> <b>P_elec_net</b>: Net electrical power output to the grid.</li>
	<li> <b>curtail</b>: Input for curtailment control.</li>
	<li> <b>P_schedule</b>: Input for controlling power scheduling.</li>
	<li> <b>electricity</b>: Output for the electrical power input to the grid.</li>
	<li> <b>on_renewable</b>: Output indicating whether the power input is from renewable sources.</li>
	</ul>
	</li>
	</ul>
	</html>", revisions="<html>
	<ul>
	<li><i>September 2023</i> by <a href=\"mailto:armando.fontalvo@anu.edu.au\">Armando Fontalvo</a>:<br>
	Created documentation for GridInput.</li>
	</ul>
	</html>"));
end GridInput;
