within examples;
model WindPVTESsystem
	extends Modelica.Icons.Example;
	import Modelica.SIunits.Conversions.*;
	import Modelica.Constants.*;

	// Importing molten salt properties
	replaceable package Medium = SolarTherm.Media.MoltenSalt.MoltenSalt_ph;

	// Importing pv, wind and schedule files
	parameter String schedule_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Schedules/dummy_heat_schedule.motab");

	// Renewable energy input
	parameter String pv_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_pv.motab");
	parameter String wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_wind.motab");
	parameter Modelica.SIunits.Power pv_ref_size = 50e6 "PV farm reference size";
	parameter Modelica.SIunits.Power wind_ref_size = 50e6 "Wind farm reference size";
	parameter Modelica.SIunits.Power P_elec_min = 1e6;
	parameter Modelica.SIunits.Efficiency pv_fraction = 0.5 "Maximum hot salt mass flow rate";
	parameter Real renewable_multiple = 2 "Renewable energy to process heat demand factor";
	parameter Real heater_multiple = 2 "Heater energy to process heat demand factor";
	parameter Modelica.SIunits.Power P_grid_max = renewable_multiple * heater_multiple * Q_process_des "Maximum hot salt mass flow rate";
	parameter Modelica.SIunits.Power P_elec_max = heater_multiple * Q_process_des "Maximum hot salt mass flow rate";

	// Heater parameters
	parameter Modelica.SIunits.Efficiency heater_efficiency = 0.99 "Electric heater efficiency";

	// Thermal energy storage parameters
	parameter Real t_storage(unit = "h") = 8 "Hours of storage";
	parameter Modelica.SIunits.Energy E_max = t_storage * 3600 * Q_process_des "Maximum tank stored energy";
	parameter Modelica.SIunits.Temperature T_amb_des = from_degC(25) "Ambient temperature at design";
	parameter Modelica.SIunits.Efficiency off_design_factor = 1 / heater_efficiency "Off-design curtailment correction";

	// Thermal process parameters
	parameter Modelica.SIunits.HeatFlowRate Q_process_des = 50e6 "Process heat demand at design";

	// Control parameters
	parameter Real level_off = 5 "Hot tank empty trigger lower bound";
	parameter Real level_on = 10 "Hot tank empty trigger upper bound";
	parameter Real level_curtailment_off = 93 "Hot tank full trigger lower bound";
	parameter Real level_curtailment_on = 98 "Hot tank full trigger upper bound";
	parameter Real split_cold = 0.7 "Starting medium fraction in cold tank";
	parameter Modelica.SIunits.Time t_process_ramp_up = 3*3600 "Delay until process starts";
	parameter Modelica.SIunits.Time t_process_ramp_dw = 2*3600 "Delay until process shuts down";

	// System defaults
	inner Modelica.Fluid.System system(T_start = T_amb_des, allowFlowReversal = false, p_start = Medium.p_default) annotation(
	Placement(visible = true, transformation(origin = {-70, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	// Ambient temperature
	Modelica.Blocks.Sources.RealExpression T_ambient(y = system.T_ambient) annotation(
	Placement(visible = true, transformation(origin = {-20, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	// Renewable energy input
	SolarTherm.Models.Sources.GridInput renewable_input(
		P_elec_max = P_grid_max, 
		P_elec_min = P_elec_min, 
		P_elec_pv_ref_size = pv_ref_size, 
		P_elec_wind_ref_size = wind_ref_size, 
		pv_file = pv_file, 
		pv_fraction = pv_fraction, 
		wind_file = wind_file, 
		renewable_multiple = renewable_multiple) 
		annotation(
			Placement(visible = true, transformation(origin = {-44, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	// Curtailment input
	Modelica.Blocks.Sources.RealExpression P_elec_schedule(y = P_curtail) annotation(
		Placement(visible = true, transformation(origin = {-82, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	// Electrical Heater
	SolarTherm.Models.Fluid.HeatExchangers.SimpleResistiveHeater heater(heater_efficiency = heater_efficiency) annotation(
		Placement(visible = true, transformation(origin = {-8, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	// Modelica
	SolarTherm.Models.Storage.Tank.GenericST thermal_storage(
		E_max = E_max, 
		E_0 = (1-split_cold) * E_max, 
		Q_flow_loss_des = 0, 
		T_amb_des = T_amb_des, 
		cf = {1, 0}, 
		ca = {1}) 
		annotation(Placement(visible = true, transformation(origin = {24, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	SolarTherm.Models.Control.GenericControl controller(
		Q_schedule_des = Q_process_des, 
		level_curtailment_off = level_curtailment_off, 
		level_curtailment_on = level_curtailment_on, 
		level_off = level_off, 
		level_on = level_on, 
		t_process_ramp_dw = t_process_ramp_dw, 
		t_process_ramp_up = t_process_ramp_up) 
	annotation(Placement(visible = true, transformation(origin = {50, 30}, extent = {{-10, 10}, {10, -10}}, rotation = 0)));

	SolarTherm.Models.Processes.GenericThermalProcess thermal_process annotation(
		Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	// Process scheduler
	Modelica.Blocks.Sources.CombiTimeTable schedule(
		fileName = schedule_file, 
		tableName = "Q_flow", 
		tableOnFile = true,
		smoothness=Modelica.Blocks.Types.Smoothness.ContinuousDerivative) annotation(
		Placement(visible = true, transformation(origin = {10, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

	// Variables
	Modelica.SIunits.Energy E_thermal(start = 0, fixed = true, displayUnit = "MW.h") "Generated Energy";
	Modelica.SIunits.Energy E_schedule(start = 0, fixed = true, displayUnit = "MW.h") "Scheduled Energy";
	Modelica.SIunits.Energy E_renewable(start = 0, fixed = true, displayUnit = "MW.h") "Renewable Energy";
	Modelica.SIunits.Power P_thermal "Thermal Output power of boiler";
	Modelica.SIunits.Power P_curtail "Electrical input due to curtailment";
equation
// Renewable heating connections
  connect(renewable_input.electricity, heater.P_elec_in) annotation(
    Line(points = {{-34, 0}, {-18, 0}}, color = {0, 0, 127}));
 connect(P_elec_schedule.y, renewable_input.P_schedule) annotation(
    Line(points = {{-71, -20}, {-64.5, -20}, {-64.5, -6}, {-54, -6}}, color = {0, 0, 127}));
	connect(T_ambient.y, thermal_storage.T_amb) annotation(
		Line(points = {{-9, 24}, {20, 24}, {20, 10}}, color = {0, 0, 127}));
 connect(heater.P_thermal_out, thermal_storage.Q_in) annotation(
    Line(points = {{2, 0}, {14, 0}}, color = {0, 0, 127}));
 connect(controller.Q_discharging, thermal_process.Q_schedule) annotation(
    Line(points = {{61, 30}, {66, 30}, {66, 10}}, color = {0, 0, 127}));
	connect(thermal_storage.Q_out, thermal_process.Q_in) annotation(
		Line(points = {{34, 0}, {60, 0}}, color = {0, 0, 127}));
 P_curtail = schedule.y[1] / heater_efficiency;
	P_thermal = controller.Q_discharging;
//	P_thermal = if controller.process_state==3 then thermal_process.Q_in else 0.0;
	thermal_storage.Q_flow_set = Q_process_des;

	der(E_thermal) = P_thermal;
 der(E_renewable) = renewable_input.electricity;
	der(E_schedule) = Q_process_des;
 connect(schedule.y[1], controller.Q_schedule) annotation(
    Line(points = {{21, 36}, {39, 36}}, color = {0, 0, 127}));
 connect(thermal_storage.Level, controller.storage_level) annotation(
    Line(points = {{34, 6}, {36, 6}, {36, 24}, {39, 24}}, color = {0, 0, 127}));
 connect(controller.curtailment, renewable_input.curtail) annotation(
    Line(points = {{50, 41}, {50, 60}, {-80, 60}, {-80, 0}, {-54, 0}}, color = {255, 0, 255}, pattern = LinePattern.Dash));
 connect(renewable_input.on_renewable, heater.on_renewable) annotation(
    Line(points = {{-44, -12}, {-44, -20}, {-30, -20}, {-30, -6}, {-18, -6}}, color = {255, 0, 255}, pattern = LinePattern.Dash));

annotation(
	Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}), graphics = {Bitmap(extent = {{-18, 64}, {-18, 64}})}),
	Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}})),
	experiment(StartTime = 0, StopTime = 3.1536e+07, Tolerance = 1e-06, Interval = 300),
	__OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
	__OpenModelica_simulationFlags(lv = "LOG_STATS", noEventEmit = "()", s = "dassl"),
	Documentation(revisions = "<html>
	<ul>
	<li> A. Fontalvo Lascano (September 2023) :<br>Released first version. </li>
	</ul>
	
	</html>"));
end WindPVTESsystem;