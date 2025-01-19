within SolarTherm.Models.Control;
//Power Block can only be turned on t-seconds after being shut down. Storage needs to accumulate x-amount of hours of discharge to start the power block.
model Thermocline_Controller_PBLimit_Level
  extends Icons.Control;
  replaceable package HTF = SolarTherm.Media.Sodium.Sodium_pT;
  parameter Real L_1 = 5.0 "Level to stop discharge, minimum PB temperature";
  parameter Real L_2 = 15.0 "Level to start discharge, slightly above L1";
  parameter Real L_3 = 85.0 "Level to start charge, slightly below L4";
  parameter Real L_4 = 95.0 "Level to stop charge, maximum recv temperature";
  parameter SI.Temperature T_target = 740.0 + 273.15 "Target receiver outlet temperature";
  
  parameter SI.MassFlowRate m_flow_PB_des = 100.0 "Reference power block mass flow rate";
  parameter SI.HeatFlowRate Q_des_blk = 200e6 "Power block design heat rate";
  parameter SI.SpecificEnthalpy h_target = HTF.specificEnthalpy(HTF.setState_pTX(101323.0, T_target)) "Target specific enthalpy of the receiver outlet";
  Integer Control_State(start=6) "1-6 Determines which pumps are flowing and whether defocus is on";
  Boolean Chg(start=true) "Can the storage be charged?";
  Boolean Disch(start=false) "Can the storage be discharged?";
  Boolean PB(start=false) "Can the PB be turned on?";

  //Timer to prevent PB from being turned on too many times
  //parameter SI.Time t_wait=1.0*3600 "Time you have to wait after shutdown before it can be turned on again";
  //SI.Time t_shutdown(start=0) "Time since it last shut down";
  //End timer
  parameter SI.HeatFlowRate Q_rcv_min = 0.10*Q_des_blk "Minimum receiver heat-rate to start mass flow to receiver";
  parameter SI.MassFlowRate m_0 = 1e-8 "Minimum mass flow rate through any pipe";
  parameter SI.MassFlowRate m_min = 1e-8 "minimum mass flow rate to start"; //used to be 1e-7 for both
  
  Modelica.Blocks.Interfaces.RealInput Level "Tank Storage Level 0-100"
    annotation (Placement(visible = true, transformation(extent = {{-126, -20}, {-86, 20}}, rotation = 0), iconTransformation(extent = {{-126, -8}, {-86, 32}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.RealOutput m_flow_PB(start=0.0) "Power block mass flow?" annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, -12}, {130, 28}}, rotation = 0))) ;
  
  Modelica.Blocks.Interfaces.RealOutput m_flow_recv(start=0.0) "Receiver mass flow?" annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, 38}, {130, 78}}, rotation = 0))) ;
  
  Modelica.Blocks.Interfaces.RealOutput Q_defocus(start=Q_des_blk) "Required defocus heat" annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(origin = {-108, -90},extent = {{-20, -20}, {20, 20}}, rotation = 180))) ;
  
  Modelica.Blocks.Interfaces.BooleanOutput defocus(start=false) "defocus receiver?" annotation (Placement(visible = true, transformation(extent = {{90, -20}, {130, 20}}, rotation = 0), iconTransformation(extent = {{90, -60}, {130, -20}}, rotation = 0))) ;

  Modelica.Blocks.Interfaces.RealInput Q_rcv_raw "The net receiver heat rate before curtailment"
    annotation (Placement(visible = true, transformation(extent = {{-124, 22}, {-84, 62}}, rotation = 0), iconTransformation(extent = {{-126, 48}, {-86, 88}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.RealInput h_PB_outlet "Enthalpy of the HTF coming out of PB"
    annotation (Placement(visible = true, transformation(extent = {{-128, -80}, {-88, -40}}, rotation = 0), iconTransformation(origin = {56, 112},extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    
  Modelica.Blocks.Interfaces.RealInput h_tank_outlet "Enthalpy of the HTF coming out of bottom of tank"
    annotation (Placement(visible = true, transformation(extent = {{-128, -100}, {-88, -60}}, rotation = 0), iconTransformation(origin = {-40, 112},extent = {{-20, -20}, {20, 20}}, rotation = -90)));

  SI.MassFlowRate m_guess(start=0.0) "Guess required flow rate of recv";
  parameter SI.Time t_wait = 2.0*3600 "Waiting time between turning off PB and being able to turn on";
  //SI.Time t_threshold(start=0.0) "if time passes this value, PB := true";
  
  //New model input parameters
  parameter SI.Energy E_max = 1200.0*3600.0e6 "theoretical storage capacity"; 
  parameter Real eff_storage_des = 0.60 "design storage utilisation";
  parameter SI.Time t_stor_startPB = 2.0*3600.0 "minimum hours of storage available to startup PB"; 
  
  //New calculated parameters
  parameter SI.Time t_stor_cap = E_max/Q_des_blk "design storage capacity in seconds of design PB operation";
  parameter Real L_startPB = t_stor_startPB/t_stor_cap + 0.5*(1.0-eff_storage_des) "Minimum tank level to start PB to ensure the required time of stored energy available";
  
  //Relative flow magnitude State
  Integer Flow_State(start = 0); //0 if receiver < 0, 1 if recv between 0 and PB, 2 if recv greater than PB
  parameter SI.MassFlowRate m_tol = 0.0001*m_flow_PB_des; //Tolerance of 0.1% of design PB mass flow rate

algorithm
  //Changing Storage State
  when Level < L_1 then 
    Disch := false;
  elsewhen Level > L_2 then 
    Disch := true;
  end when;
  
  when Level > L_4 then 
    Chg := false;
  elsewhen Level < L_3 then
    Chg := true;
  end when;
  
  /*
  //when m_flow_PB < 0.1*m_flow_PB_des then //take this as shutdown
  when m_flow_PB < 100.0*m_0 then
  //when pre(Control_State) > 5 or pre(Control_State) < 2 then
    if PB == true then
      PB := false; //start the cooldown
      t_threshold := time + t_wait;
    end if;
  end when;
  when time > t_threshold then
    PB := true;
  end when;
  */
  //New parts
  when m_guess > m_tol then
    if m_guess < m_flow_PB_des + m_tol then //no overshoot
      Flow_State := 1;
    else //overshot
      Flow_State := 2;
    end if;
  end when;
  when m_guess < m_min then
    Flow_State := 0;
  end when;
  
  when m_guess > m_flow_PB_des + m_tol then
    Flow_State := 2;
  end when;
  when m_guess < m_flow_PB_des then
    if m_guess > m_tol then //no overshoot
      Flow_State := 1;
    else
      Flow_State := 0;
    end if;
  end when;
  
  //New controls: PB triggered by shutdown and re-enabled after tank level increases past a certain value
  when m_flow_PB < 0.2 * m_flow_PB_des then
    if Level <= L_startPB then //try this condition
      PB := false;
    end if;
  end when;
  when Level > L_startPB then
    PB := true;
  end when;

equation
  //m_guess = Q_rcv_raw/(h_target-max(h_tank_outlet,h_PB_outlet));
  //m_guess =(Q_rcv_raw + m_flow_PB*(h_PB_outlet-h_tank_outlet))/(h_target-h_tank_outlet);
  m_guess =(Q_rcv_raw + m_flow_PB*(h_PB_outlet-h_tank_outlet))/(h_target-h_tank_outlet);

/* 
  if m_guess < m_min then
    if Disch == true then
      if PB == true then
        Control_State = 2;
      else
        Control_State = 6;
      end if;
    else
      Control_State = 6;
    end if;
  elseif m_guess >= m_min and m_guess <= m_flow_PB_des then
    if Disch == true then
      if PB == true then
        Control_State = 4;
      else
        Control_State = 6;
      end if;
    else
      if Chg == true then
        Control_State = 1;
      else
        Control_State = 6;
      end if;
    end if;
  else
    if Chg == true then
      if PB == true then
        Control_State = 5;
      else
        Control_State = 1;
      end if;
    else
      if PB == true then
        Control_State = 3;
      else
        Control_State = 6;
      end if;
    end if;
  end if;
*/
  if Flow_State == 0 then
    if Disch == true then
      if PB == true then
        Control_State = 2;
      else
        Control_State = 6;
      end if;
    else
      Control_State = 6;
    end if;
  elseif Flow_State == 1 then
    if Disch == true then
      if PB == true then
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
      if PB == true then
        Control_State = 5;
      else
        Control_State = 1;
      end if;
    else
      if PB == true then
        Control_State = 3;
      else
        Control_State = 6;
      end if;
    end if;
  end if;
   
  if Control_State == 1 then
    m_flow_recv = Q_rcv_raw/(h_target-h_tank_outlet);
    m_flow_PB = m_0;
    defocus = false;
    Q_defocus = m_flow_PB_des*(h_target-h_PB_outlet); //Not used anyway
    
  elseif Control_State == 2 then
    m_flow_recv = m_0;
    m_flow_PB = m_flow_PB_des;
    defocus = false;
    Q_defocus = m_flow_PB_des*(h_target-h_PB_outlet); //Not used anyway

  elseif Control_State == 3 then
    m_flow_recv = m_flow_PB_des;
    m_flow_PB = m_flow_PB_des;
    defocus = true;
    Q_defocus = m_flow_PB_des*(h_target-h_PB_outlet); //Not used anyway

  elseif Control_State == 4 then
    m_flow_recv = Q_rcv_raw/(h_target-h_PB_outlet);
    m_flow_PB = m_flow_PB_des; //whoops I switched these by mistake
    defocus = false;
    Q_defocus = m_flow_PB_des*(h_target-h_PB_outlet); //Not used anyway

  elseif Control_State == 5 then
    m_flow_recv = (Q_rcv_raw + m_flow_PB*(h_PB_outlet-h_tank_outlet))/(h_target-h_tank_outlet);
    m_flow_PB = m_flow_PB_des;
    defocus = false;
    Q_defocus = m_flow_PB_des*(h_target-h_PB_outlet); //Not used anyway

  else
    m_flow_recv = m_0;
    m_flow_PB = m_0;
    defocus = false;
    Q_defocus = m_flow_PB_des*(h_target-h_PB_outlet); //Not used anyway
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
end Thermocline_Controller_PBLimit_Level;