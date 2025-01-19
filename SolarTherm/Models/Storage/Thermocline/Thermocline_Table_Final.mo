within SolarTherm.Models.Storage.Thermocline;
  
model Thermocline_Table_Final
  extends SolarTherm.Interfaces.Models.StorageFluid_Thermocline;
  import SI = Modelica.SIunits;
  
  parameter String table_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Storage/Lookup_Table_SAR3.motab");
  Modelica.Blocks.Tables.CombiTable1Ds Table_Charging (tableOnFile=true, tableName="table_charging", columns=2:2, fileName=table_file);
  Modelica.Blocks.Tables.CombiTable1Ds Table_Discharging (tableOnFile=true, tableName="table_discharging", columns=2:2, fileName=table_file);
  
  //Metadata
  parameter Real[21] Meta = SolarTherm.Utilities.Metadata_Storage(table_file);
  parameter SI.Temperature T_min = Meta[1] + 273.15 "(K)";
  parameter SI.Temperature T_recv_start = Meta[2] + 273.15 "(K)";
  parameter SI.Temperature T_recv_max = Meta[3] + 273.15 "(K)";
  parameter SI.Temperature T_PB_min = Meta[4] + 273.15 "(K)";
  parameter SI.Temperature T_PB_start = Meta[5] + 273.15 "(K)";
  parameter SI.Temperature T_max = Meta[6] + 273.15 "(K)";
  parameter Real L_recv_start = Meta[7] "(Real)";
  parameter Real L_recv_max = Meta[8] "(Real)";
  parameter Real L_PB_min = Meta[9] "(Real)";
  parameter Real L_PB_start = Meta[10] "(Real)";
  parameter Real util = Meta[11] "Utilisation (Real)";
  parameter SI.Energy E_max = Meta[12]*3600.0e6 "Theoretical Energy Capacity (J)";
  parameter SI.Time t_charge = Meta[13]*3600.0 "Design Charging period (h)";
  parameter SI.Time t_discharge = Meta[14]*3600.0 "Design Discharging period (h)";
  parameter SI.CoefficientOfHeatTransfer U_loss = Meta[15] "Overall heatloss coefficient (W/m2K)";
  parameter SI.Area A_loss_total = Meta[16] "Total tank surface area for heat loss (m2)";
  parameter Real C_filler = Meta[17];
  parameter Real C_fluid = Meta[18];
  parameter Real C_tank = Meta[19];
  parameter Real C_insulation = Meta[20];
  parameter Real C_encapsulation = Meta[21];
  parameter Real C_total = C_filler + C_fluid + C_tank + C_insulation + C_encapsulation;

  //parameter Medium.AbsolutePressure p_amb = 101325;
  // Storage
  //replaceable package Storage = SolarTherm.Materials.Graphite "the storage medium";
  //parameter Modelica.SIunits.Temperature T_min = 500 + 273.15 "start temperature of the storage medium";
  //parameter Modelica.SIunits.Temperature T_max = 720 + 273.15 "the maximum temperature of the storage medium";
  parameter Modelica.SIunits.SpecificEnthalpy h_max = Fluid_Package.h_Tf(T_max, 1.0) "the maximum enthalpy of the storage medium";
  parameter Modelica.SIunits.SpecificEnthalpy h_min = Fluid_Package.h_Tf(T_min, 0.0) "the initial enthalpy of the storage medium";
  //parameter Modelica.SIunits.Mass m_st=E_max/(h_max - h_0) "total mass of the storage medium";
  //parameter Modelica.SIunits.Energy E_max;
  Modelica.Blocks.Interfaces.RealInput T_amb "Ambient Temperature" annotation (Placement(
        visible = true,transformation(
        
        origin={-50, 8.88178e-16},extent={{-10, -10}, {10, 10}},
        rotation=0), iconTransformation(
        
        origin={-46, 0},extent={{-6, -6}, {6, 6}},
        rotation=0)));
  Modelica.Blocks.Interfaces.RealInput p_amb "Ambient Pressure" annotation(
    Placement(visible = true, transformation(origin = {48, 8.88178e-16}, extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {46, 0}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
  //Modelica.SIunits.Temperature T_st"tempaerture of the storage medium";
  //Modelica.SIunits.SpecificEnthalpy h_st "enthlpy of the sotrage medium";
  // HTF
  redeclare replaceable package Medium = SolarTherm.Media.Sodium.Sodium_pT "heat transfer medium";
  replaceable package Fluid_Package = SolarTherm.Materials.Sodium_Table "fluid material package";
  Modelica.SIunits.Temperature T_in;
  Modelica.SIunits.Temperature T_out;
  Modelica.SIunits.SpecificEnthalpy h_in;
  Modelica.SIunits.SpecificEnthalpy h_out;
  Fluid_Package.State fluid_top "Property calculation of inlet";
  Fluid_Package.State fluid_bot "Property calculation of outlet";
  // energy balance
  Modelica.SIunits.MassFlowRate m_flow "mass flow with respect to height direction";
  Modelica.SIunits.Energy E_stored(start = L_start * E_max);
  Modelica.SIunits.EnergyFlowRate Q_loss "HeatLoss rate (W)";
  //Real e_out "outlet effectiveness";
  Modelica.Blocks.Interfaces.RealOutput h_bot_outlet "Enthaply at the bottom of the tank as an output signal (K)" annotation(
    Placement(visible = true, transformation(origin = {-40, -70}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-35, -65}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealOutput Level "Storage Level NonDimensional" annotation(
    Placement(visible = true, transformation(origin = {-41, 49}, extent = {{13, -13}, {-13, 13}}, rotation = 0), iconTransformation(origin = {-48, 38}, extent = {{8, -8}, {-8, 8}}, rotation = 0)));
  //Effectiveness
  parameter Real L_start = 0.0 "Starting storage level";
  Real L(start = L_start) "Tank Level";
  Real e_top "Nondimensional top temperature/effectiveness";
  Real e_bot "Nondimensional bot temperature/effectiveness if it were to flow out of bottom";
  //parameter Real C_total = 3.11232e7 "total cost";
  //Measured outlet temperature
  Real T_outlet_degC "Outlet temperature in degrees Celcius";
  //Fluid_Package.State fluid_bot;
  //Safety calculations
  //parameter Real e_top_min = (T_PB_min-T_min)/(T_max-T_min);
  //parameter Real e_bot_min = (T_max-T_recv_max)/(T_max-T_min);
equation
  //Table inputs
  Table_Charging.u = L;
  Table_Discharging.u = L;
  
  
//Theoretical bottom outlet effectiveness
  e_bot = min(1.0,Table_Charging.y[1]);
  e_top = min(1.0,Table_Discharging.y[1]);
  //e_bot = max(e_bot_min,min(1.0,Table_Charging.y[1]));
  //e_top = max(e_top_min,min(1.0,Table_Discharging.y[1]));
//Derived temperatures
  fluid_top.T = T_min + e_top*(T_max-T_min);
  fluid_bot.T = T_max - e_bot*(T_max-T_min);
  h_bot_outlet = fluid_bot.h;
  der(E_stored) = m_flow*(inStream(fluid_b.h_outflow)-fluid_a.h_outflow)-Q_loss;
  Q_loss = U_loss*A_loss_total*(T_min+Level*(T_max-T_min)-T_amb);


//effectiveness calculations
  m_flow = -1.0*fluid_a.m_flow; 
  if m_flow >= 0.0 then //flowing upwards so discharge
    h_in = inStream(fluid_b.h_outflow);
    fluid_b.h_outflow = fluid_bot.h;//h_in;
    fluid_a.h_outflow = fluid_top.h;
    h_out = fluid_top.h;
    T_out = fluid_top.T;//T_min + e_out*(T_max-T_min);
    T_in = fluid_bot.T;
    //h_bot_outlet = h_min;//h_min;
    //fluid_bot.T = T_min;
    //h_bot_outlet = Fluid_Package.h_Tf(T_max - e_bot*(T_max-T_min), 0.0);
  else //fluid flowing downwards so charging
    h_in = inStream(fluid_a.h_outflow);
    fluid_a.h_outflow = fluid_top.h;// h_in;
    fluid_b.h_outflow = fluid_bot.h;//h_out;
    h_out = fluid_bot.h;
    T_out = fluid_bot.T;//T_max - e_out*(T_max-T_min);
    T_in = fluid_top.T;
    //h_bot_outlet = h_out;
    //h_bot_outlet = fluid_bot.h;//Fluid_Package.h_Tf(T_out, 0.0);
    //fluid_bot.T = T_out;
  end if;
  //try constant e_out
  //e_out = 0.99;
  //fluid_bot.T = T_max - e_out*(T_max-T_min);
  
  
  L = E_stored/E_max;
  Level = L;
  
  //fluid_in.h = h_in;
  //T_in = fluid_in.T;
  
  //fluid_out.T = T_out;
  //h_out = fluid_out.h;
  
  fluid_a.m_flow + fluid_b.m_flow = 0;
  
  fluid_a.p = p_amb;
  fluid_b.p = p_amb;
  
  if m_flow > 1.0e-3 then //Discharging, outlet is the top
    T_outlet_degC = T_out - 273.15;
  elseif m_flow < -1.0e-3 then //Charging, outlet is the bottom
    T_outlet_degC = T_out - 273.15;
  else //No flow, output reference temperature
    T_outlet_degC = 25.0;
  end if;

  annotation(
    Icon(graphics = {Text(origin = {1, 4}, extent = {{-27, 32}, {27, -32}}, textString = "Reg")}));
end Thermocline_Table_Final;