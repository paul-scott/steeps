within SolarTherm.Systems;

model WindPVPHESsystem
  extends Modelica.Icons.Example;
  import Modelica.SIunits.Conversions.*;
  import Modelica.Constants.*;
  // Renewable energy input
  parameter String pv_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_pv.motab");
  parameter String wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_wind.motab");
  parameter Modelica.SIunits.Power pv_ref_size = 50e6 "PV farm reference size";
  parameter Modelica.SIunits.Power wind_ref_size = 50e6 "Wind farm reference size";
  parameter Modelica.SIunits.Power P_elec_min = 1e6;
  parameter Modelica.SIunits.Efficiency pv_fraction = 0.5 "Maximum hot salt mass flow rate";
  parameter Real renewable_multiple = 2 "Renewable energy to process heat demand factor";
  parameter Real heater_multiple = 2 "Heater energy to process heat demand factor";
  parameter Modelica.SIunits.Power P_elec_max = heater_multiple * Q_process_des "Maximum hot salt mass flow rate";
  // Thermal energy storage parameters
  parameter Real t_storage(unit = "h") = 8 "Hours of storage";
  parameter Modelica.SIunits.Energy E_max = t_storage * 3600 * Q_process_des "Maximum tank stored energy";
  // Thermal process parameters
  parameter Modelica.SIunits.HeatFlowRate Q_process_des = 50e6 "Process heat demand at design";
  // Control parameters
  parameter Real level_off = 5 "Hot tank empty trigger lower bound";
  parameter Real level_on = 10 "Hot tank empty trigger upper bound";
  parameter Real level_curtailment_off = 93 "Hot tank full trigger lower bound";
  parameter Real level_curtailment_on = 98 "Hot tank full trigger upper bound";
  parameter Real split_cold = 0.7 "Starting medium fraction in cold tank";
  parameter Modelica.SIunits.Time t_process_ramp_up = 3 * 3600 "Delay until process starts";
  parameter Modelica.SIunits.Time t_process_ramp_dw = 2 * 3600 "Delay until process shuts down";
  // Renewable energy input
  SolarTherm.Models.Sources.GridInputSplit renewable_input(P_elec_max = P_elec_max, P_elec_min = P_elec_min, P_elec_pv_ref_size = pv_ref_size, P_elec_wind_ref_size = wind_ref_size, pv_file = pv_file, pv_fraction = pv_fraction, wind_file = wind_file, renewable_multiple = renewable_multiple) annotation(
    Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Curtailment input
  Modelica.Blocks.Sources.RealExpression P_elec_schedule(y = P_curtail) annotation(
    Placement(visible = true, transformation(origin = {-82, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Electrical Heater
  // Modelica
  // Process scheduler
  // Variables
  Modelica.SIunits.Energy E_thermal(start = 0, fixed = true, displayUnit = "MW.h") "Generated Energy";
  Modelica.SIunits.Energy E_schedule(start = 0, fixed = true, displayUnit = "MW.h") "Scheduled Energy";
  Modelica.SIunits.Energy E_renewable(start = 0, fixed = true, displayUnit = "MW.h") "Renewable Energy";
  Modelica.SIunits.Power P_thermal "Thermal Output power of boiler";
  Modelica.SIunits.Power P_curtail "Electrical input due to curtailment";
  Modelica.Blocks.Sources.BooleanExpression curtail annotation(
    Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Storage.PumpHydro PHES annotation(
    Placement(visible = true, transformation(origin = {0, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression P_sch_PHES annotation(
    Placement(visible = true, transformation(origin = {30, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression P_sch_ren(y = Q_process_des) annotation(
    Placement(visible = true, transformation(origin = {28, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Processes.GenericThermalProcess2inputs process annotation(
    Placement(visible = true, transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
// Renewable heating connections
  P_curtail = Q_process_des;
  P_thermal = Q_process_des;
  der(E_thermal) = P_thermal;
  der(E_renewable) = renewable_input.P_elec_out1;
  der(E_schedule) = Q_process_des;
  connect(P_elec_schedule.y, renewable_input.P_schedule) annotation(
    Line(points = {{-70, -20}, {-64, -20}, {-64, -6}, {-60, -6}}, color = {0, 0, 127}));
  connect(curtail.y, renewable_input.curtail) annotation(
    Line(points = {{-72, 0}, {-60, 0}}, color = {255, 0, 255}));
  connect(renewable_input.P_elec_out2, PHES.P_elec_in) annotation(
    Line(points = {{-40, -4}, {-24, -4}, {-24, -30}, {-10, -30}}, color = {0, 0, 127}));
  connect(P_sch_PHES.y, process.Q_schedule2) annotation(
    Line(points = {{42, 44}, {54, 44}, {54, 10}}, color = {0, 0, 127}));
  connect(P_sch_ren.y, process.Q_schedule1) annotation(
    Line(points = {{40, 20}, {46, 20}, {46, 10}}, color = {0, 0, 127}));
  connect(renewable_input.P_elec_out1, process.Q_in1) annotation(
    Line(points = {{-40, 6}, {40, 6}}, color = {0, 0, 127}));
  connect(PHES.P_elec_out, process.Q_in2) annotation(
    Line(points = {{10, -30}, {22, -30}, {22, -4}, {40, -4}}, color = {0, 0, 127}));
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
end WindPVPHESsystem;
