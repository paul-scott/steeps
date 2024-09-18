within SolarTherm.Systems.H2DRI_Applications;

model H2DRI_DesignCase_2b_Dynamic
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  
  //Material Packages
  replaceable package Material_IOE_OreD = SolarTherm.Materials.IOE_Dehydroxylated;
  replaceable package Medium_Ore_Dehydroxylated = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_ph;
  
  parameter String PV_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_pv.motab");
  parameter String Wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/dummy_wind.motab");
  
  //Parameter Inputs
  parameter Real RM = 1.0 "Renewable Multiple (pre-transmission oversizing)";
  parameter Real HM = 2.0 "Heater Multiple";
  parameter Real PV_fraction = 0.5 "PV_fraction";
  parameter SI.Time t_storage = 10.0*3600.0;
  
  
  
  parameter SI.HeatFlowRate Q_process_des = 7.20898e7;
  parameter SI.MassFlowRate m_flow_ore_des = 154.1 "Mass flow of ore out of the hot tank and into the med tank if running at design point (kg/s)";
  
  parameter Real eff_heater = 0.95 "Electrical-to-heat conversion efficiency of the heater";
  //Renewable Parameters
  parameter SI.Power P_renewable_des = RM * P_heater_des;
  parameter SI.HeatFlowRate Q_heater_des = HM * Q_process_des;
  parameter SI.Power P_heater_des = Q_heater_des / eff_heater;
  parameter SI.Power PV_ref_size = 50.0e6;
  parameter SI.Power Wind_ref_size = 50.0e6;
  parameter Real LOF_PV = 1.300;
  parameter Real LOF_Wind = 1.131;
  parameter SI.Power P_wind_net = (1.0 - PV_fraction)*P_renewable_des;
  parameter SI.Power P_PV_net = PV_fraction*P_renewable_des;
  parameter SI.Power P_wind_gross = P_wind_net*LOF_Wind;
  parameter SI.Power P_PV_gross = P_PV_net*LOF_PV;
  
  //Tank Parameters
  parameter SI.Temperature T_hot_set = 757.6 + 273.15;
  parameter SI.Temperature T_med_set = 289.0 + 273.15;
  
  parameter Real ar = 2.0;
  parameter Real eps_packing = 0.20;
  parameter Real eps_material = 0.089;
  parameter Real epsilon = eps_packing + eps_material - eps_packing*eps_material "Effective porosity";
  parameter SI.Density rho_ore_med = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_utilities.rho_T(T_med_set) "Density of pure ore (kg/m3)";
  parameter SI.Density rho_ore_hot = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_utilities.rho_T(T_hot_set) "Density of pure ore (kg/m3)";
  
  parameter SI.Mass m_ore_tank_hot = m_flow_ore_des*t_storage "Maximum mass capacity of the hot tank (kg)";
  parameter SI.Mass m_ore_tank_med = m_flow_ore_des*t_storage "Maximum mass capacity of the hot tank (kg)";
  
  parameter SI.Volume V_tank_hot = m_ore_tank_hot/(rho_ore_hot*(1.0-epsilon));
  parameter SI.Volume V_tank_med = m_ore_tank_med/(rho_ore_med*(1.0-epsilon));
  
  parameter SI.Length D_tank_hot = (4.0*V_tank_hot/(CN.pi*ar))^(1.0/3.0);
  
  parameter SI.Length D_tank_med = (4.0*V_tank_med/(CN.pi*ar))^(1.0/3.0);
  
  parameter SI.Length H_tank_hot = D_tank_hot*ar;
  parameter SI.Length H_tank_med = D_tank_med*ar;
  
  Modelica.Blocks.Sources.CombiTimeTable PV_input(fileName = PV_file, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "Power", tableOnFile = true) annotation(
    Placement(visible = true, transformation(origin = {-88, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Wind_input(fileName = Wind_file, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "Power", tableOnFile = true) annotation(
    Placement(visible = true, transformation(origin = {-88, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add Grid_Sum(k1 = P_PV_gross / PV_ref_size, k2 = P_wind_gross / Wind_ref_size) annotation(
    Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.CSP.CRS.Receivers.Basic_Heater Heater(redeclare package Medium = Medium_Ore_Dehydroxylated, P_heater_des = P_heater_des, Q_heater_des = Q_heater_des, T_cold_set = T_med_set, T_hot_set = T_hot_set, eff_heater = eff_heater) annotation(
    Placement(visible = true, transformation(origin = {28, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Storage.Tank.Silo Med_Tank(redeclare package Medium = Medium_Ore_Dehydroxylated, redeclare package Filler_Package = Material_IOE_OreD, T_min = T_med_set - 10.0, T_max = T_med_set + 10.0, U_loss_tank = 0.0, T_start = T_med_set, T_set = T_med_set - 10.0, L_start = 0.50, epsilon = epsilon, H_tank = H_tank_med, D_tank = D_tank_med) annotation(
    Placement(visible = true, transformation(origin = {-18, -8}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  SolarTherm.Models.Storage.Tank.Silo Hot_Tank(redeclare package Medium = Medium_Ore_Dehydroxylated, redeclare package Filler_Package = Material_IOE_OreD, T_min = T_hot_set - 10.0, T_max = T_hot_set + 10.0, U_loss_tank = 0.0, T_start = T_hot_set, T_set = T_hot_set - 10.0, L_start = 0.50, epsilon = epsilon, H_tank = H_tank_hot, D_tank = D_tank_hot) annotation(
    Placement(visible = true, transformation(origin = {54, -6}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Heater_Lift(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {4, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink2 Sink(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {92, -12}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT OreD_source(redeclare package Medium = Medium_Ore_Dehydroxylated, T = 289.0 + 273.15,nPorts = 1, p = 100000, use_T_in = false) annotation(
    Placement(visible = true, transformation(origin = {-74, -12}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Med_Lift(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {-42, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure Reactor_Lift(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {74, -12}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));

//Controller
  parameter SI.Time t_wait = 1.0*3600 "Waiting time between turning off process and being able to turn on";
  SI.Time t_threshold(start=0.0) "if time passes this value, Process := true";

  Boolean Chg(start = true);
  Boolean Dis(start = true);
  
  Boolean Process(start = true);
  
  Integer State(start = 1);
  
  SI.Mass m_Fe_produced(start=0);
  SI.Mass m_Fe_target(start=0);
  Real CapF(start=0);
  //Boolean Defocus(start = false);
  
  
  //if Heater.Q_heater_raw <= 0.0 then
    //if Dis == true then
      //State = 4;
    //end if;
  //end if;
  

algorithm


  //Boiler Timer Control
  when Reactor_Lift.m_flow < 0.01 * m_flow_ore_des then //take this as shutdown
    Process := false; //start the cooldown
    t_threshold := time + t_wait;
  end when;
  when time > t_threshold then
    Process := true;
  end when;


  when Hot_Tank.L > 0.98 then
    Chg := false;
  end when;
  when Hot_Tank.L < 0.95 then
    Chg := true;
  end when;
  
  when Hot_Tank.L < 0.02 then
    Dis := false;
  end when;
  when Hot_Tank.L > 0.05 then
    Dis := true;
  end when;

equation
  der(m_Fe_target) = m_flow_ore_des*0.6242;
  der(m_Fe_produced) = Sink.port_a.m_flow*0.6242;
  if time < 10.0 then
    CapF = 0.0;
  else
    CapF = m_Fe_produced/m_Fe_target;
  end if;


  if Heater.Q_heater_raw <= 1.0e-6 then
    if Dis == true then
      if Process == true then
        State = 4;
      else
        State = 6;
      end if;
    else
      State = 6;
    end if;
  elseif Heater.Q_heater_raw >= Q_process_des then
    if Process == true then
      if Chg == true then
        State = 1;
      else
        State = 3;
      end if;
    else
      if Chg == true then
        State = 5;
      else
        State = 6;
      end if;
    end if;
  else
    if Dis == true then
      if Process == true then
        State = 2;
      else
        State = 6;
      end if;
    else
      if Chg == true then
        State = 5;
      else
        State = 6;
      end if;
    end if;
  end if;
  if State == 1 then
    Heater_Lift.m_flow = (Heater.Q_heater_raw/Q_process_des)*m_flow_ore_des;
  
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des; //Not used anyway

    Med_Lift.m_flow = m_flow_ore_des;  
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 2 then
    Heater_Lift.m_flow = (Heater.Q_heater_raw/Q_process_des)*m_flow_ore_des;
    
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des; //Not used anyway
    
    Med_Lift.m_flow = m_flow_ore_des;  
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 3 then
    Heater_Lift.m_flow = m_flow_ore_des;
    
    Heater.curtail = true;
    Heater.Q_curtail = Q_process_des; //Curtailed
    
    Med_Lift.m_flow = m_flow_ore_des;  
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 4 then
    Heater_Lift.m_flow = 1.0e-9;
    
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des; //Not used anyway
    
    Med_Lift.m_flow = m_flow_ore_des;  
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 5 then
    Heater_Lift.m_flow = (Heater.Q_heater_raw/Q_process_des)*m_flow_ore_des;
    
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des; //Not used anyway
    
    Med_Lift.m_flow = 1.0e-9;  
    Reactor_Lift.m_flow = 1.0e-9;
  else
    Heater_Lift.m_flow = 1.0e-9;
    
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des; //Not used anyway
    Med_Lift.m_flow = 1.0e-9;  
    Reactor_Lift.m_flow = 1.0e-9;
  end if;

  Med_Tank.T_amb = 298.15;
  Hot_Tank.T_amb = 298.15;
  
  Med_Tank.p = 100000.0;
  Hot_Tank.p = 100000.0;
  connect(Wind_input.y[1], Grid_Sum.u2) annotation(
    Line(points = {{-76, 26}, {-70, 26}, {-70, 34}, {-62, 34}, {-62, 34}}, color = {0, 0, 127}));
  connect(PV_input.y[1], Grid_Sum.u1) annotation(
    Line(points = {{-76, 54}, {-70, 54}, {-70, 46}, {-62, 46}, {-62, 46}}, color = {0, 0, 127}));
  connect(Med_Lift.fluid_b, Med_Tank.fluid_a) annotation(
    Line(points = {{-36, -14}, {-32, -14}, {-32, -4}, {-26, -4}, {-26, -4}}, color = {0, 127, 255}));
  connect(OreD_source.ports[1], Med_Lift.fluid_a) annotation(
    Line(points = {{-66, -12}, {-58, -12}, {-58, -14}, {-48, -14}}, color = {0, 127, 255}));
  connect(Reactor_Lift.fluid_b, Sink.port_a) annotation(
    Line(points = {{80, -12}, {86, -12}}, color = {0, 127, 255}));
  connect(Hot_Tank.fluid_b, Reactor_Lift.fluid_a) annotation(
    Line(points = {{62, -12}, {68, -12}, {68, -12}, {68, -12}}, color = {0, 127, 255}));
  connect(Grid_Sum.y, Heater.P_supply) annotation(
    Line(points = {{-38, 40}, {-6, 40}, {-6, 1}, {17, 1}}, color = {0, 0, 127}));
  connect(Heater.fluid_b, Hot_Tank.fluid_a) annotation(
    Line(points = {{38, -6}, {42, -6}, {42, -2}, {46, -2}, {46, -2}}, color = {0, 127, 255}));
  connect(Med_Tank.fluid_b, Heater_Lift.fluid_a) annotation(
    Line(points = {{-10, -14}, {-2, -14}, {-2, -14}, {-2, -14}}, color = {0, 127, 255}));
  connect(Heater_Lift.fluid_b, Heater.fluid_a) annotation(
    Line(points = {{10, -14}, {14, -14}, {14, -6}, {18, -6}, {18, -6}}, color = {0, 127, 255}));

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)), experiment(StopTime = 3.1536e+07, StartTime = 0, Tolerance = 1.0e-5, Interval = 300, maxStepSize = 60, initialStepSize = 60));
end H2DRI_DesignCase_2b_Dynamic;