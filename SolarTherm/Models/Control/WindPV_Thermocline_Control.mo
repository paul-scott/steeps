within SolarTherm.Models.Control;

model WindPV_Thermocline_Control
  extends Icons.Control; //Gives it a nice icon
  replaceable package HTF = SolarTherm.Media.Air.Air_amb_p; //Define the material property model
  
  //Temperature threshold parameters
  parameter SI.Temperature T_target = 565.0 + 273.15 "Target receiver outlet temperature";
  parameter SI.Temperature T_boiler_start = 525.0 + 273.15 "PB input temperature at which it can start";
  parameter SI.Temperature T_boiler_min = 515.0 + 273.15 "Minimum tolerated PB input temperature";
  parameter SI.Temperature T_heater_max = 330.0 + 273.15 "Maximum tolerated receiver input temperature";
  parameter SI.Temperature T_heater_start = 320.0 + 273.15 "Receiver input temperature at which receiver can start";
  
  //Boiler design parameters
  parameter SI.MassFlowRate m_boiler_des = 100.0 "Reference HTF mass flow rate to the boiler at the target enthalpy";
  parameter SI.HeatFlowRate Q_boiler_des = 200e6 "Boiler design input heat rate";
  parameter SI.SpecificEnthalpy h_target = HTF.specificEnthalpy(HTF.setState_pTX(101323.0, T_target)) "Specific enthalpy target of the heater outlet and the boiler inlet";
  
  //Mass flow rate parameters
  parameter SI.MassFlowRate m_0 = 1e-8 "Minimum mass flow rate through any pipe";
  parameter SI.MassFlowRate m_min = 1e-8 "minimum mass flow rate to start"; //used to be 1e-7 for both
  parameter SI.MassFlowRate m_tol = 0.001*m_boiler_des "Hysteresis tolerance of 0.1% of design PB mass flow rate used for the controller";
  
  //Storage control parameters 
  parameter SI.Energy E_max = 12.0*3600.0*800.0e6 "Ideal storage capacity (J)"; 
  parameter SI.Time t_stor_cap = E_max/Q_boiler_des "Just the ideal storage capacity in terms of seconds (s)";
  parameter Real level_mid = 0.50 "Midpoint storage level determined via component-level analysis.";
  parameter Real util_storage_des = 0.50 "Storage utilisation determined via component-level analysis.";
  parameter SI.Time t_stor_start_dis = 0.25*3600.0 "Minimum effective storage seconds accumulated that is needed to start discharge (s)."; 
  parameter Real L_start_dis = t_stor_start_dis/t_stor_cap + level_mid-0.5*util_storage_des "Tank level that corresponds to t_stor_start_dis being available";
  
  //Boiler timer. It needs to wait t_wait number of seconds before being turned on again
  parameter SI.Time t_wait = 1.0*3600 "Waiting time between turning off PB and being able to turn on";
  SI.Time t_threshold(start=0.0) "if time passes this value, PB := true";
  
  

  //Relative flow magnitude State
  SI.MassFlowRate m_guess(start=0.0) "Required mass flow rate of the heater assuming it is on and not curtailed.";
  SI.MassFlowRate m_boiler_dem = m_boiler_des*Q_demand/Q_boiler_des "Instantaneous boiler demand mass flow rate at design enthalpy (kg/s)";
  Integer Flow_State(start = 0); //0 if heater < 0, 1 if heater between 0 and boiler, 2 if heater greater than boiler.
  Integer Control_State(start=6) "1-6 Determines which pumps are operational and whether curtailment is on";
  Boolean Chg(start=true) "Can the storage be charged?";
  Boolean Dis(start=false) "Can the storage be discharged based on temperature?";
  Boolean Boil(start=true) "Can the Boiler be turned on?";
  
  Modelica.Blocks.Interfaces.RealInput Level "Tank Storage Level 0-100"
    annotation (Placement(visible = true, transformation(extent = {{-124, 24}, {-84, 64}}, rotation = 0), iconTransformation(extent = {{-126, 6}, {-86, 46}}, rotation = 0)));
  
  Modelica.Blocks.Interfaces.RealInput T_top_tank "Temperature of the top of HTF in storage (K)"
    annotation (Placement(visible = true, transformation(extent = {{-124, -4}, {-84, 36}}, rotation = 0), iconTransformation(extent = {{-126, -36}, {-86, 4}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.RealInput T_bot_tank "Temperature of the bottom of HTF in storage (K)"
    annotation (Placement(visible = true, transformation(extent = {{-124, -30}, {-84, 10}}, rotation = 0), iconTransformation(extent = {{-126, -78}, {-86, -38}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.RealOutput m_boiler_signal(start=1.0e-8) "Signal output to the boiler loop pump (kg/s)." annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, -12}, {130, 28}}, rotation = 0))) ;
  
  Modelica.Blocks.Interfaces.RealOutput m_heater_signal(start=1.0e-8) "Signal output to the heater loop pump (kg/s)." annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, 38}, {130, 78}}, rotation = 0))) ;
  
  Modelica.Blocks.Interfaces.RealOutput Q_curtail(start=Q_boiler_des) "Signal output to the heater stating the required curtailed heating rate (W)" annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(origin = {-108, -90},extent = {{-20, -20}, {20, 20}}, rotation = 180))) ;
  
  Modelica.Blocks.Interfaces.BooleanOutput curtail(start=false) "Signal output to the heater stating the need for curtailment." annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, -60}, {130, -20}}, rotation = 0))) ;

  Modelica.Blocks.Interfaces.RealInput Q_heater_raw "Input signal from the heater stating the available heating rate (W)"
    annotation (Placement(visible = true, transformation(extent = {{-124, 52}, {-84, 92}}, rotation = 0), iconTransformation(extent = {{-126, 48}, {-86, 88}}, rotation = 0)));
  
  Modelica.Blocks.Interfaces.RealInput Q_demand "Input signal from the scheduler stating the heat-rate demanded by the boiler (W)"
    annotation (Placement(visible = true, transformation(extent = {{-124, -102}, {-84, -62}}, rotation = 0), iconTransformation(origin = {70, 112},extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  
  Modelica.Blocks.Interfaces.RealInput h_boiler_outlet "Cold outlet enthalpy from the boiler (J/kg)"
    annotation (Placement(visible = true, transformation(extent = {{-124, -54}, {-84, -14}}, rotation = 0), iconTransformation(origin = {14, 112},extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    
  Modelica.Blocks.Interfaces.RealInput h_tank_bot "Enthalpy of the HTF coming out of bottom of tank (J/kg)"
    annotation (Placement(visible = true, transformation(extent = {{-124, -78}, {-84, -38}}, rotation = 0), iconTransformation(origin = {-40, 112},extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    
  Modelica.Blocks.Interfaces.RealInput h_tank_top "Enthalpy of the HTF coming out of top of tank (J/kg)"
    annotation (Placement(visible = true, transformation(extent = {{-124, 80}, {-84, 120}}, rotation = 0), iconTransformation(origin = {-90, 112},extent = {{-20, -20}, {20, 20}}, rotation = -90)));

initial algorithm
  //Just determine if the storage can be discharged right at the start of simulation.
  if Level <= L_start_dis then
    Dis := false;
  else
    Dis := true;
  end if;

algorithm
  //Discharge Control
  when Level > L_start_dis then
    if T_top_tank > T_boiler_start then
      Dis := true;
    end if;
  end when;
  when T_top_tank > T_boiler_start then
    if Level > L_start_dis then
      Dis := true;
    end if;
  end when;
   when T_top_tank < T_boiler_min then 
    Dis := false;
  end when;
  
  //Charge Control
  when T_bot_tank > T_heater_max then 
    Chg := false;
  elsewhen T_bot_tank < T_heater_start then
    Chg := true;
  end when;
  
  //Boiler Timer Control
  when m_boiler_signal < 0.1 * m_boiler_des then //take this as shutdown
    Boil := false; //start the cooldown
    t_threshold := time + t_wait;
  end when;
  when time > t_threshold then
    Boil := true;
  end when;
  
  //System Flow State Control
  when m_guess > m_tol then
    if m_guess < m_boiler_dem + m_tol then //no overshoot
      Flow_State := 1;
    else //overshot
      Flow_State := 2;
    end if;
  end when;
  when m_guess < m_min then
    Flow_State := 0;
  end when;
  when m_guess > m_boiler_dem + m_tol then
    Flow_State := 2;
  end when;
  when m_guess < m_boiler_dem then
    if m_guess > m_tol then //no overshoot
      Flow_State := 1;
    else
      Flow_State := 0;
    end if;
  end when;

equation
  if Q_heater_raw >= Q_demand then
    m_guess = (Q_heater_raw + m_boiler_dem*(h_boiler_outlet-h_tank_bot))/(h_target-h_tank_bot);
  else
    m_guess = Q_heater_raw/(h_target-h_boiler_outlet);
  end if;

  if Flow_State == 0 then
    if Dis == true then
      if Boil == true then
        Control_State = 2;
      else
        Control_State = 6;
      end if;
    else
      Control_State = 6;
    end if;
  elseif Flow_State == 1 then
    if Dis == true then
      if Boil == true then
        Control_State = 4;
      else
        Control_State = 1;
      end if;
    else
      if Chg == true then
        Control_State = 1;
      else
        Control_State = 6;
      end if;
    end if;
  else //Flow_State == 2
    if Chg == true then
      if Boil == true then
        Control_State = 5;
      else
        Control_State = 1;
      end if;
    else
      if Boil == true then
        Control_State = 3;
      else
        Control_State = 6;
      end if;
    end if;
  end if;
   
  if Control_State == 1 then
    m_heater_signal = max(m_0,Q_heater_raw/(h_target-h_tank_bot));
    m_boiler_signal = m_0;
    curtail = false;
    Q_curtail = Q_boiler_des; //Not used anyway
    
  elseif Control_State == 2 then
    m_heater_signal = m_0;
    m_boiler_signal = max(m_0,m_boiler_dem*(h_target-h_boiler_outlet)/(h_tank_top-h_boiler_outlet));
    curtail = false;
    Q_curtail = Q_boiler_des; //Not used anyway

  elseif Control_State == 3 then
    m_heater_signal = max(m_0,m_boiler_dem);
    m_boiler_signal = max(m_0,m_boiler_dem);
    curtail = true;
    Q_curtail = m_boiler_dem*(h_target-h_boiler_outlet);

  elseif Control_State == 4 then
    m_heater_signal = max(m_0,Q_heater_raw/(h_target-h_boiler_outlet));
    m_boiler_signal = max(m_0, (m_heater_signal*(h_tank_top-h_target) + m_boiler_dem*(h_target-h_boiler_outlet))/(h_tank_top-h_boiler_outlet));
    curtail = false;
    Q_curtail = Q_boiler_des; //Not used anyway

  elseif Control_State == 5 then
    m_heater_signal = max(m_0,(Q_heater_raw + m_boiler_dem*(h_boiler_outlet-h_tank_bot))/(h_target-h_tank_bot));
    m_boiler_signal = max(m_0,m_boiler_dem);
    curtail = false;
    Q_curtail = Q_boiler_des; //Not used anyway

  else
    m_heater_signal = m_0;
    m_boiler_signal = m_0;
    curtail = false;
    Q_curtail = Q_boiler_des; //Not used anyway
  end if;
//Additional info about control states:
//1 = Recv is on, only charges the storage. PB is off.
//2 = PB is on, and only run by discharging the Storage. Recv is off.
//3 = PB is run only from receiver, storage is completely bypassed, receiver is defocused to exactly balance mass flows of recv and PB, obeying also the target receiver outlet temperature.
//4 = PB is on, and run by combining receiver outlet with a Storage discharge stream.
//5 = Receiver is on, and its outlet stream splits to charge the Storage and run the PB.
//6 = Everything is off.
  annotation(Documentation(revisions ="<html>
		<p>By Zebedee Kee on 03/12/2020</p>
		</html>",info="<html>
		<p>This component determines the mass flow rates of both the receiver and power block mass flow rates. The variable m_guess calculates the required receiver mass flow to achieve target outlet temperature T_target based on inlet enthalpy from either storage bottom outlet, PB outlet or a combination of both. Depending on whether the storage is allowed to charge, discharge and relative size of m_guess wrt minimum flowrate and PB design flowrate, one of the 6 operating states is chosen.</p>
		</html>"));
end WindPV_Thermocline_Control;