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
  SolarTherm.Systems.GridInputSplit renewable_input(P_elec_max = P_elec_max, P_elec_min = P_elec_min, P_elec_pv_ref_size = pv_ref_size, P_elec_wind_ref_size = wind_ref_size, pv_file = pv_file, pv_fraction = pv_fraction, wind_file = wind_file, renewable_multiple = renewable_multiple) annotation(
    Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Curtailment input
  // Electrical Heater
  // Modelica
  // Process scheduler
  // Variables
  Modelica.SIunits.Energy E_thermal(start = 0, fixed = true, displayUnit = "MW.h") "Generated Energy";
  Modelica.SIunits.Energy E_schedule(start = 0, fixed = true, displayUnit = "MW.h") "Scheduled Energy";
  Modelica.SIunits.Energy E_renewable(start = 0, fixed = true, displayUnit = "MW.h") "Renewable Energy";
  Modelica.SIunits.Power P_thermal "Thermal Output power of boiler";
  Modelica.SIunits.Power P_curtail "Electrical input due to curtailment";
  SolarTherm.Systems.PumpHydroStorage PHES annotation(
    Placement(visible = true, transformation(origin = {-40, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.Beneficiation Beneficiation annotation(
    Placement(visible = true, transformation(origin = {34, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.HydrogenStorage H2_storage annotation(
    Placement(visible = true, transformation(origin = {32, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.GenericControl control annotation(
    Placement(visible = true, transformation(origin = {-72, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.Electrolyser electrolyser annotation(
    Placement(visible = true, transformation(origin = {-10, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergyMerger merger1 annotation(
    Placement(visible = true, transformation(origin = {-30, 70}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergySplitter split1 annotation(
    Placement(visible = true, transformation(origin = {-60, 0}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergySplitter split2 annotation(
    Placement(visible = true, transformation(origin = {-50, 4}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergySplitter split3 annotation(
    Placement(visible = true, transformation(origin = {-20, -40}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergyMerger merge2 annotation(
    Placement(visible = true, transformation(origin = {-2, -4}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergySplitter split4 annotation(
    Placement(visible = true, transformation(origin = {10, 70}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.FBReactor reactor annotation(
    Placement(visible = true, transformation(origin = {72, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergyMerger energyMerger annotation(
    Placement(visible = true, transformation(origin = {50, 70}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.HydrogenBurner burner annotation(
    Placement(visible = true, transformation(origin = {82, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergySplitter split5 annotation(
    Placement(visible = true, transformation(origin = {60, 70}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.EnergySplitter energySplitter annotation(
    Placement(visible = true, transformation(origin = {10, -4}, extent = {{-4, -10}, {4, 10}}, rotation = 0)));
  SolarTherm.Systems.SmelterBOF smelterBOF annotation(
    Placement(visible = true, transformation(origin = {78, -40}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression DRI_input(y=param.m_DRI) annotation(
    Placement(visible = true, transformation(origin = {4, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Systems.Parameters param annotation(
    Placement(visible = true, transformation(origin = {-70, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
// Renewable heating connections
  P_curtail = Q_process_des;
  P_thermal = Q_process_des;
  der(E_thermal) = P_thermal;
  der(E_renewable) = renewable_input.P_elec_out1;
  der(E_schedule) = Q_process_des;
  connect(renewable_input.P_elec_out, split1.u1) annotation(
    Line(points = {{-70, 0}, {-62, 0}}, color = {0, 0, 127}));
  connect(split1.y1, split2.u1) annotation(
    Line(points = {{-58, 4}, {-52, 4}}, color = {0, 0, 127}));
  connect(split1.y2, PHES.P_elec_in) annotation(
    Line(points = {{-58, -4}, {-55, -4}, {-55, -40}, {-50, -40}}, color = {0, 0, 127}));
  connect(PHES.P_elec_out, split3.u1) annotation(
    Line(points = {{-30, -40}, {-22, -40}}, color = {0, 0, 127}));
  connect(split3.y2, merge2.u2) annotation(
    Line(points = {{-18, -44}, {-7, -44}, {-7, -8}, {-4, -8}}, color = {0, 0, 127}));
  connect(split2.y1, merger1.u1) annotation(
    Line(points = {{-48, 8}, {-40, 8}, {-40, 74}, {-32, 74}}, color = {0, 0, 127}));
  connect(split3.y1, merger1.u2) annotation(
    Line(points = {{-18, -36}, {-12, -36}, {-12, 8}, {-36, 8}, {-36, 66}, {-32, 66}}, color = {0, 0, 127}));
  connect(merger1.y, electrolyser.P_elec_in) annotation(
    Line(points = {{-28.4, 70}, {-20.4, 70}}, color = {0, 0, 127}));
  connect(electrolyser.P_thermal_out, split4.u1) annotation(
    Line(points = {{0, 70}, {8, 70}}, color = {0, 0, 127}));
  connect(split4.y1, H2_storage.E_in) annotation(
    Line(points = {{11.6, 74}, {21.6, 74}}, color = {0, 0, 127}));
  connect(split2.y2, merge2.u1) annotation(
    Line(points = {{-48, 0}, {-4, 0}}, color = {0, 0, 127}));
  connect(H2_storage.E_out, energyMerger.u1) annotation(
    Line(points = {{42, 74}, {48, 74}}, color = {0, 0, 127}));
  connect(split4.y2, energyMerger.u2) annotation(
    Line(points = {{11.6, 66}, {19.6, 66}, {19.6, 58}, {43.6, 58}, {43.6, 66}, {47.6, 66}}, color = {0, 0, 127}));
  connect(energyMerger.y, split5.u1) annotation(
    Line(points = {{51.6, 70}, {57.6, 70}}, color = {0, 0, 127}));
  connect(split5.y1, burner.E_in) annotation(
    Line(points = {{61.6, 74}, {71.6, 74}}, color = {0, 0, 127}));
  connect(Beneficiation.y, reactor.u2) annotation(
    Line(points = {{44, 6}, {62, 6}}, color = {0, 0, 127}, pattern = LinePattern.DashDotDot));
  connect(energySplitter.y1, Beneficiation.E_in) annotation(
    Line(points = {{12, 0}, {24, 0}}, color = {0, 0, 127}));
  connect(merge2.y, energySplitter.u1) annotation(
    Line(points = {{0, -4}, {8, -4}}, color = {0, 0, 127}));
  connect(burner.y, reactor.u1) annotation(
    Line(points = {{92, 76}, {96, 76}, {96, 20}, {50, 20}, {50, 12}, {62, 12}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(split5.y2, reactor.u3) annotation(
    Line(points = {{62, 66}, {66, 66}, {66, 24}, {48, 24}, {48, 0}, {62, 0}}, color = {0, 0, 127}));
  connect(energySplitter.y2, smelterBOF.u2) annotation(
    Line(points = {{12, -8}, {14, -8}, {14, -36}, {59, -36}}, color = {0, 0, 127}));
  connect(reactor.y, smelterBOF.u1) annotation(
    Line(points = {{82, 6}, {96, 6}, {96, -12}, {54, -12}, {54, -24}, {59, -24}}, color = {0, 0, 127}, pattern = LinePattern.DashDotDot));
  connect(DRI_input.y, Beneficiation.E_schedule) annotation(
    Line(points = {{16, 30}, {20, 30}, {20, 6}, {24, 6}}, color = {0, 0, 127}, pattern = LinePattern.DashDotDot));
  annotation(
    Diagram(coordinateSystem(extent = {{-100, 100}, {100, -100}}), graphics = {Bitmap(extent = {{-18, 64}, {-18, 64}})}),
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