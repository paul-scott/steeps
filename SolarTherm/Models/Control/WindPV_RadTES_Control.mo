within SolarTherm.Models.Control;

model WindPV_RadTES_Control
  //Controller for a radiatively - charged TES system
  extends Icons.Control;
  replaceable package HTF = SolarTherm.Media.Air.Air_amb_p; //Fluid property model
  
  //Input Parameters
  parameter SI.Temperature T_rad_max = 800.0 + 273.15 "Maximum allowable heater element temperature at which charging is stopped (K)";
  parameter SI.Temperature T_rad_start = 750.0 + 273.15 "Temperature of the hottest heater element below which charging can start (K)";
  parameter SI.Temperature T_boiler_min = 450.0 + 273.15 "Minimum acceptable boiler HTF inlet temperature (K), also the design boiler temperature";
  parameter SI.Temperature T_boiler_start = 500.0 + 273.15 "Temperature of the top of the TES at which discharge can be started (K)";

  //Boiler Parameters
  parameter SI.MassFlowRate m_boiler_des = 100.0 "Reference boiler inlet mass flow rate (kg/s)";
  parameter SI.HeatFlowRate Q_boiler_des = 200e6 "Power block design heat rate (W)";
  parameter SI.SpecificEnthalpy h_boiler_des = HTF.specificEnthalpy(HTF.setState_pTX(101323.0, T_boiler_min)) "Boiler design HTF inlet enthalpy (J/kg)";
  
  //Mass-flow control parameters
  parameter SI.MassFlowRate m_0 = 1e-8 "Minimum mass flow rate through any pipe (kg/s)";
  parameter SI.MassFlowRate m_min = 1e-8 "minimum mass flow rate to start (kg/s)";
  
  //Connectors
  Modelica.Blocks.Interfaces.RealInput T_top_tank "Temperature of the top fluid element of the TES (K)" annotation(
    Placement(visible = true, transformation(extent = {{-124, -4}, {-84, 36}}, rotation = 0), iconTransformation(extent = {{-124, -50}, {-84, -10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput T_rad_measured "Temperature of the hottest heating element (K)" annotation(
    Placement(visible = true, transformation(extent = {{-124, -4}, {-84, 36}}, rotation = 0), iconTransformation(extent = {{-124, 44}, {-84, 84}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput m_flow_boiler(start = 1.0e-8) "Boiler loop mass flow rate signal (kg/s)" annotation(
    Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, 6}, {130, 46}}, rotation = 0)));
  Modelica.Blocks.Interfaces.BooleanOutput heater_on(start = true) "Heater on/off signal (Boolean)" annotation(
    Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, -60}, {130, -20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput Q_demand "The instantaneous boiler heat-rate demand (W)" annotation(
    Placement(visible = true, transformation(extent = {{-124, -102}, {-84, -62}}, rotation = 0), iconTransformation(origin = {70, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealInput h_boiler_out "Design HTF cold outlet from the boiler (J/kg)" annotation(
    Placement(visible = true, transformation(extent = {{-124, -54}, {-84, -14}}, rotation = 0), iconTransformation(origin = {14, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealInput h_tank_top "Enthalpy of the HTF coming out of top of tank (J/kg)" annotation(
    Placement(visible = true, transformation(extent = {{-124, 80}, {-84, 120}}, rotation = 0), iconTransformation(origin = {-90, 112}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    
  //Variables
  SI.MassFlowRate m_boiler_dem = m_boiler_des * Q_demand / Q_boiler_des "Instantaneous boiler input mass flow demand (kg/s) if it is at the design enthalpy";

  Boolean Chg(start = true) "Can the TES be charged (Boolean)?";
  Boolean Dis(start = false) "Can the TES be discharged (Boolean)?";
  
algorithm
  //The algorithm section is used for variables that are turned on/off by certain thresholds being crossed

  //Radiative charging is controlled by the temperature of the hottest heating element in the TES
  when T_rad_measured > T_rad_max then
    Chg := false;
  end when;
  
  when T_rad_measured < T_rad_start then
    Chg := true;
  end when;

  //Discharging is controlled by the temperature at the top-most TES fluid element
  when T_top_tank > T_boiler_start then
    Dis := true;
  end when;

  when T_top_tank < T_boiler_min then
    Dis := false;
  end when;
  
equation
  //The boiler mass flow rate calculation
  if Dis == true then
    m_flow_boiler = m_boiler_dem * (h_boiler_des - h_boiler_out) / (h_tank_top - h_boiler_out); //Mass flow after considering off-design
  else
    m_flow_boiler = m_min;
  end if;
  
  //The heater controller signal
  if Chg == true then
    heater_on = true;
  else
    heater_on = false;
  end if;

  annotation(
    Documentation(revisions = "<html>
		<p>By Zebedee Kee on 01/12/2023</p>
		</html>", info = "<html>
		<p>This component controls the steam-supply system involving a Wind-PV electrical input, boiler output and a TES that is charged via radiative heating.</p>
		</html>"));
end WindPV_RadTES_Control;