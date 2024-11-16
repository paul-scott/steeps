within SolarTherm.Models.Storage.Thermocline.Parallel;

model Thermocline_HBS_LC_2P_MixedOutlet
  //Two-tank mixed outlet flow, charging B-BA-A discharging A-AB-B
  extends SolarTherm.Interfaces.Models.StorageFluid_Thermocline;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  //Initialise Material Packages
  replaceable package Medium = SolarTherm.Media.Air.Air_CoolProp_1bar;
  replaceable package Fluid_Package = SolarTherm.Materials.PartialMaterial;
  replaceable package Filler_Package_A = SolarTherm.Materials.PartialMaterial;
  replaceable package Filler_Package_B = SolarTherm.Materials.PartialMaterial;

  //Storage Parameter Settings
  parameter Integer Correlation = 1 "Internal Flow, Gas";
    //Storage Cpacity and Fractions
  parameter SI.Energy E_max = 144.0e9 "Maximum storage capacity of entire group (J)";
  parameter Real frac_1 = 1.0 / 2.0 "Fraction of storage capacity of Tank_A";
  
    //Aspect ratios (H/D) of tanks
  parameter Real ar_A = 2.0 "Aspect ratio of tank";
  parameter Real ar_B = ar_A "Aspect ratio of tank";
    //Porosity of tank filler materials
  parameter Real eta_A = 0.51 "Porosity";
  parameter Real eta_B = eta_A "Porosity";
    //Hole diameter of filler material
  parameter Real d_p_A = 0.03 "Filler hole hydraulic diameter (m)";
  parameter Real d_p_B = d_p_A "Filler hole hydraulic diameter (m)";
    
    //Discretization settings
  parameter Integer N_f_A = 100 "Number of fluid CVs in Tank_A";
  parameter Integer N_f_B = N_f_A "Number of fluid CVs in Tank_B";
  
    //Heat loss coefficient of tanks
  parameter SI.CoefficientOfHeatTransfer U_loss_top_A = 0.698 "Heat loss coefficient at the top of the tank (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_top_B = U_loss_top_A "W/m2K";
  
  parameter SI.CoefficientOfHeatTransfer U_loss_bot_A = 1.22 "Heat loss coefficient at the bottom of the tank (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_bot_B = U_loss_bot_A "W/m2K";
  
    //Temperature settings
  parameter SI.Temperature T_min = CV.from_deg(515) "Minimum temperature (design) also starting T";
  parameter SI.Temperature T_max = CV.from_deg(715) "Maximum design temperature (design)";
    //Internal control for temperature
  parameter SI.Temperature T_recv_set = 273.15 + 570.0 "Temperature of T_05 at which it switches to the next tank during charging";
  parameter SI.Temperature T_PB_set = 273.15 + 680.0 "Temperature of T_95 at which it switches to the previous tank durng discharging";
  //Internal Controls
  //Integer Active_Tank(start = 1) "Which tank is in use currently";
  Real f_chg_guess(start=0.0) "Fraction of mass flow of secondary tank during charging";
  Real f_disch_guess(start=0.0) "Fraction of mass flow of secondary tank during discharging";
  
  
  Real f_chg(start=0.0) "Fraction of mass flow of secondary tank during charging";
  Real f_disch(start=0.0) "Fraction of mass flow of secondary tank during discharging";
  
  Integer Charge_State(start=1) "1 = Only A, 2 = B assists A, 3 = Only B";
  Integer Discharge_State(start=1) "1 = Only A, 2 = B assists A, 3 = Only B";
  
  parameter SI.SpecificEnthalpy h_recv_set = Fluid_Package.h_Tf(T_recv_set,0.0);
  parameter SI.SpecificEnthalpy h_PB_set = Fluid_Package.h_Tf(T_PB_set,0.0);
  
  //Input and Output Ports
  Modelica.Blocks.Interfaces.RealOutput T_top_measured "Temperature at the top of the tank as an output signal (K)" annotation(
    Placement(visible = true, transformation(extent = {{40, 50}, {60, 70}}, rotation = 0), iconTransformation(origin = {45, 55}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput T_bot_measured "Temperature at the bottom of the tank as an output signal (K)" annotation(
    Placement(visible = true, transformation(extent = {{40, -70}, {60, -50}}, rotation = 0), iconTransformation(origin = {45, -55}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput h_bot_outlet "Enthaply at the bottom of the tank as an output signal (J/kg)" annotation(
    Placement(visible = true, transformation(extent = {{40, -72}, {60, -52}}, rotation = 0), iconTransformation(origin = {-27, -65}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealOutput h_top_outlet "Enthaply at the top of the tank as an output signal (J/kg)" annotation(
    Placement(visible = true, transformation(extent = {{40, -72}, {60, -52}}, rotation = 0), iconTransformation(origin = {-27, 65}, extent = {{5, -5}, {-5, 5}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealInput T_amb "Ambient Temperature" annotation(
    Placement(visible = true, transformation(origin = {-50, 8.88178e-16}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-46, 0}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput p_amb "Ambient Pressure" annotation(
    Placement(visible = true, transformation(origin = {48, 8.88178e-16}, extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {46, 0}, extent = {{6, -6}, {-6, 6}}, rotation = 0)));
  //Initialize Tank_A
  SolarTherm.Models.Storage.Thermocline.Thermocline_HBS_LC_Section_Final Tank_A(redeclare replaceable package Fluid_Package = Fluid_Package, redeclare replaceable package Filler_Package = Filler_Package_A, Correlation = Correlation, E_max = E_max * frac_1, ar = ar_A, eta = eta_A, d_p = d_p_A, T_min = T_min, T_max = T_max, N_f = N_f_A, U_loss_top = U_loss_top_A, U_loss_bot = U_loss_bot_A, z_offset = 0.0, eff_pump=eff_pump);
  //Initialize Tank_B
  SolarTherm.Models.Storage.Thermocline.Thermocline_HBS_LC_Section_Final Tank_B(redeclare replaceable package Fluid_Package = Fluid_Package, redeclare replaceable package Filler_Package = Filler_Package_B, Correlation = Correlation, E_max = E_max * (1.0 - frac_1), ar = ar_B, eta = eta_B, d_p = d_p_B, T_min = T_min, T_max = T_max, N_f = N_f_B, U_loss_top = U_loss_top_B, U_loss_bot = U_loss_bot_B, z_offset = 0.0, eff_pump=eff_pump);

  //Cost BreakDown
  parameter Real C_filler = Tank_A.C_filler + Tank_B.C_filler;
  parameter Real C_fluid = Tank_A.C_fluid + Tank_B.C_fluid;
  parameter Real C_total = Tank_A.C_section + Tank_B.C_section;
  parameter Real C_tank = Tank_A.C_tank + Tank_B.C_tank ;
  parameter Real C_insulation = Tank_A.C_insulation + Tank_B.C_insulation ;
  parameter Real C_encapsulation = Tank_A.C_encapsulation + Tank_B.C_encapsulation ;
  
  //Analytics
    //Tank Energy Levels
  Modelica.Blocks.Interfaces.RealOutput Level "Level of the entire storage" annotation(
    Placement(visible = true, transformation(extent = {{40, 50}, {60, 70}}, rotation = 0), iconTransformation(origin = {45, 21}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

    //Tank Non-dimensionalized vertical axis
  parameter Real ZDH_A[N_f_A] = Tank_A.ZDH;
  parameter Real ZDH_B[N_f_B] = Tank_B.ZDH;
    //Tank Temperature profiles in degC units
  Real T_f_A_degC[N_f_A](start = fill(T_min, N_f_A));
  Real T_f_B_degC[N_f_B](start = fill(T_min, N_f_B));
    //Tank temperature at 5% and 95% vertical positions

  parameter SI.MassFlowRate m_0 = 1.0e-7 "Minimum mass flow rate, just to avoid zero";
  
  //Analysis of fluid entering and exiting storage
  Fluid_Package.State fluid_top "Fluid entering/exiting top";
  Fluid_Package.State fluid_bot "Fluid entering/exiting bottom";
  
  //Enthaply tolerance to avoid chattering
  parameter SI.SpecificEnthalpy h_tol = 0.01*(Tank_A.h_f_max - Tank_A.h_f_min) "enthalpy tolerance, 1.0% of maximum difference";
  
  parameter SI.TemperatureDifference T_tol = 0.5;
  parameter Real f_tol = 0.05;
  
  //Total pumping losses
  SI.Power W_loss_pump = Tank_A.W_loss_pump + Tank_B.W_loss_pump;
  parameter Real eff_pump = 0.8 "Pumping efficiency, fed into physical model";

algorithm
  //Tank A assists Tank B during charging, therefore when Tank B needs help, modify the chargestate boolean
  when Tank_B.T_f[1] > T_recv_set then
    Charge_State := 2;
  elsewhen Tank_B.T_f[1] < T_recv_set - T_tol then
    Charge_State := 1;
  end when;
  
  when f_chg_guess > 1.0 then
    Charge_State := 3;
  elsewhen f_chg_guess < 1.0 - f_tol then
    Charge_State := 2;
  end when;
  
  when Tank_A.T_f[N_f_A] < T_PB_set then
    Discharge_State := 2;
  elsewhen Tank_A.T_f[N_f_A] > T_PB_set + T_tol then
    Discharge_State := 1;
  end when;
  
  when f_disch_guess > 1.0 then
    Discharge_State := 3;
  elsewhen f_disch_guess < 1.0 - f_tol then
    Discharge_State := 2;
  end when;
  
  //The assisting tank gets measured
  T_bot_measured := Tank_A.T_f[1];
  T_top_measured := Tank_B.T_f[N_f_B];


equation
  //Figure out assisted mass fraction
  //Tank B assists discharging
    if Tank_B.h_f[N_f_B] > Tank_A.h_f[N_f_A] + h_tol then
      f_disch_guess = (h_PB_set - Tank_A.h_f[N_f_A]) / (Tank_B.h_f[N_f_B] - Tank_A.h_f[N_f_A]);
    else //Cannot Assist
      f_disch_guess = 0.0;
    end if;
  //Tank A assists charging
    if Tank_A.h_f[1] < Tank_B.h_f[1] - h_tol then
      f_chg_guess = (h_recv_set - Tank_B.h_f[1]) / (Tank_A.h_f[1] - Tank_B.h_f[1] );
    else
      f_chg_guess = 0.0;
    end if;

  if Charge_State == 1 then
    f_chg = 0.0;
  elseif Charge_State == 2 then
    f_chg = f_chg_guess;
  else
    f_chg = 1.0;
  end if;
  
  if Discharge_State == 1 then
    f_disch = 0.0;
  elseif Discharge_State == 2 then
    f_disch = f_disch_guess;
  else
    f_disch = 1.0;
  end if;
  /*
  if Tank_B.h_f[N_f_B] > Tank_A.h_f[N_f_A] then
    f_disch = min(1.0, (h_PB_set - Tank_A.h_f[N_f_A]) / (Tank_B.h_f[N_f_B] - Tank_A.h_f[N_f_A] ));
  else
    f_disch = 0.0;
  end if;
  
  if Tank_B.h_f[1] < Tank_A.h_f[1] then
    f_chg = min(1.0, (h_recv_set - Tank_A.h_f[1]) / (Tank_B.h_f[1] - Tank_A.h_f[1] ));
  else
    f_chg = 0.0;
  end if;
  */
  if fluid_a.m_flow > 1e-6 then
    fluid_top.h = inStream(fluid_a.h_outflow);
    fluid_bot.h = fluid_b.h_outflow;
  elseif fluid_a.m_flow < -1e-6 then
    fluid_top.h = fluid_a.h_outflow;
    fluid_bot.h = inStream(fluid_b.h_outflow);
  else
    fluid_top.T = 298.15;
    fluid_bot.T = 298.15;
  end if;
  //Convert from Kelvin to degC for easier plotting
  T_f_A_degC = Tank_A.T_f .- 273.15;
  T_f_B_degC = Tank_B.T_f .- 273.15;
  //Calculate tank energy level
  Level = (1 - frac_1) * Tank_B.Level + frac_1 * Tank_A.Level;
  //Charging
  if fluid_a.m_flow > 0.0 then //Charging
      Tank_B.m_flow = -1.0 * fluid_a.m_flow * (1.0 - f_chg); //Swap Tank ID
      Tank_A.m_flow = -1.0 * fluid_a.m_flow * f_chg; //Swap Tank ID
      
      Tank_A.h_in = inStream(fluid_a.h_outflow);
      Tank_B.h_in = inStream(fluid_a.h_outflow);
      
      fluid_a.h_outflow = Tank_B.h_in * (1.0 - f_chg) + Tank_A.h_in * (f_chg); //Swap Tank ID
      fluid_b.h_outflow = Tank_B.h_out * (1.0 - f_chg) + Tank_A.h_out * (f_chg); //Swap Tank ID
      
      //h_top_outlet = (1.0 - f_chg)*Tank_B.h_f[N_f_B] + f_chg*Tank_B.h_f[N_f_B];
      
      
  else //Discharging
      Tank_A.m_flow = -1.0 * fluid_a.m_flow * (1.0 - f_disch);
      Tank_B.m_flow = -1.0 * fluid_a.m_flow * f_disch;
      
      Tank_A.h_in = inStream(fluid_b.h_outflow);
      Tank_B.h_in = inStream(fluid_b.h_outflow);
      
      fluid_a.h_outflow = Tank_A.h_out * (1.0 - f_disch) + Tank_B.h_out * (f_disch);
      fluid_b.h_outflow = Tank_A.h_in * (1.0 - f_disch) + Tank_B.h_in * (f_disch);
      
      
      
  end if;
  //h_bot_outlet = (1.0 - f_disch)*Tank_A.h_f[1] + f_disch*Tank_B.h_f[1];
  h_bot_outlet = (1.0 - f_chg)*Tank_B.h_f[1] + f_chg*Tank_A.h_f[1];
  h_top_outlet = (1.0 - f_disch)*Tank_A.h_f[N_f_A] + f_disch*Tank_B.h_f[N_f_B];
  //Connect pressure and ambient temp
  fluid_a.p = p_amb;
  fluid_a.p = fluid_b.p;
  T_amb = Tank_A.T_amb;
  T_amb = Tank_B.T_amb;
  
  
  //Steady state mass flows
  fluid_a.m_flow = -1.0 * fluid_b.m_flow;
  annotation(
    Icon(graphics = {Rectangle(origin = {9, 49}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-49, 11}, {31, -109}}), Text(origin = {-27, 33}, extent = {{-5, 5}, {5, -7}}, textString = "A"), Text(origin = {17, 33}, extent = {{-5, 5}, {5, -7}}, textString = "B"), Rectangle(origin = {7, 11}, fillColor = {144, 144, 144}, fillPattern = FillPattern.Solid, extent = {{-39, 7}, {-19, 1}}), Rectangle(origin = {7, -1}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, extent = {{-39, 7}, {-19, 1}}), Rectangle(origin = {-3, 17}, fillColor = {203, 203, 203}, fillPattern = FillPattern.Solid, extent = {{-29, 7}, {-9, 1}}), Rectangle(origin = {1, -17}, fillColor = {24, 24, 24}, fillPattern = FillPattern.Solid, extent = {{-33, -1}, {-13, -7}}), Rectangle(origin = {5, -13}, fillColor = {71, 71, 71}, fillPattern = FillPattern.Solid, extent = {{-37, 7}, {-17, 1}}), Rectangle(origin = {7, 5}, fillColor = {113, 113, 113}, fillPattern = FillPattern.Solid, extent = {{-39, 7}, {-19, 1}}), Rectangle(origin = {1, -11}, fillColor = {47, 47, 47}, fillPattern = FillPattern.Solid, extent = {{-33, -1}, {-13, -7}}), Line(origin = {0, 53}, points = {{0, 9}, {0, -7}, {0, -7}}), Line(origin = {0, -53}, points = {{0, -7}, {0, 7}, {0, 7}}), Ellipse(origin = {-5, -41}, extent = {{1, -1}, {9, -9}}, endAngle = 360), Ellipse(origin = {-5, 51}, extent = {{1, -1}, {9, -9}}, endAngle = 360), Rectangle(origin = {45, 17}, fillColor = {203, 203, 203}, fillPattern = FillPattern.Solid, extent = {{-33, 7}, {-13, 1}}), Rectangle(origin = {49, 11}, fillColor = {144, 144, 144}, fillPattern = FillPattern.Solid, extent = {{-37, 7}, {-17, 1}}), Rectangle(origin = {49, 5}, fillColor = {113, 113, 113}, fillPattern = FillPattern.Solid, extent = {{-37, 7}, {-17, 1}}), Rectangle(origin = {49, -1}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, extent = {{-37, 7}, {-17, 1}}), Rectangle(origin = {49, -13}, fillColor = {71, 71, 71}, fillPattern = FillPattern.Solid, extent = {{-37, 7}, {-17, 1}}), Rectangle(origin = {43, -11}, fillColor = {47, 47, 47}, fillPattern = FillPattern.Solid, extent = {{-31, -1}, {-11, -7}}), Rectangle(origin = {43, -17}, fillColor = {24, 24, 24}, fillPattern = FillPattern.Solid, extent = {{-31, -1}, {-11, -7}}), Rectangle(origin = {7, -7}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, extent = {{-39, 7}, {-19, 1}}), Rectangle(origin = {49, -7}, fillColor = {95, 95, 95}, fillPattern = FillPattern.Solid, extent = {{-37, 7}, {-17, 1}}), Polygon(origin = {-23, 6}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-3, 20}, {5, 20}, {11, 18}, {11, -30}, {5, -32}, {-3, -32}, {-9, -30}, {-9, 18}, {-3, 20}}), Line(origin = {-22, 23.1708}, points = {{-10, 0.82918}, {-4, -1.17082}, {4, -1.17082}, {10, 0.82918}}), Ellipse(origin = {-22, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {1.5, -2.5}}, endAngle = 360), Ellipse(origin = {-26, 24}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Ellipse(origin = {-26, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Ellipse(origin = {-20, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Ellipse(origin = {-20, 24}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Ellipse(origin = {-28, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {1.5, -2.5}}, endAngle = 360), Ellipse(origin = {-16, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {1.5, -2.5}}, endAngle = 360), Line(origin = {-26, -2}, points = {{0, 24}, {0, -24}}), Line(origin = {-18, -2}, points = {{0, 24}, {0, -24}}), Polygon(origin = {21, 6}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-3, 20}, {5, 20}, {11, 18}, {11, -30}, {5, -32}, {-3, -32}, {-9, -30}, {-9, 18}, {-3, 20}}), Line(origin = {22, 23.1325}, points = {{-10, 0.867544}, {-4, -1.13246}, {4, -1.13246}, {10, 0.867544}, {10, 0.867544}}), Line(origin = {18, -2}, points = {{0, 24}, {0, -24}}), Line(origin = {26, -2}, points = {{0, 24}, {0, -24}}), Ellipse(origin = {18, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Ellipse(origin = {24, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Ellipse(origin = {22, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {1.5, -2.5}}, endAngle = 360), Ellipse(origin = {16, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {1.5, -2.5}}, endAngle = 360), Ellipse(origin = {28, 26}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {1.5, -2.5}}, endAngle = 360), Ellipse(origin = {18, 24}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Ellipse(origin = {24, 24}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-0.5, -0.5}, {2.5, -1.5}}, endAngle = 360), Line(origin = {-22, 25}, points = {{0, 1}, {0, -1}}), Line(origin = {22, 25}, points = {{0, 1}, {0, -1}}), Polygon(origin = {-22, 10}, fillColor = {156, 156, 156}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {-22, 2}, fillColor = {138, 138, 138}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {-22, -6}, fillColor = {108, 108, 108}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {-22, -14}, fillColor = {80, 80, 80}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {-22, -22}, fillColor = {61, 61, 61}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {-22, 18}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Line(origin = {-26, -2}, points = {{0, 24}, {0, -24}}), Line(origin = {-18, -2}, points = {{0, 24}, {0, -24}}), Polygon(origin = {22, 10}, fillColor = {156, 156, 156}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {22, 18}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {22, 2}, fillColor = {138, 138, 138}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {22, -6}, fillColor = {108, 108, 108}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {22, -14}, fillColor = {80, 80, 80}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Polygon(origin = {22, -22}, fillColor = {61, 61, 61}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {-4, 4}, {4, 4}, {10, 6}, {10, -2}, {4, -4}, {-4, -4}, {-10, -2}, {-10, 6}}), Line(origin = {18, -2}, points = {{0, 24}, {0, -24}}), Line(origin = {26, -2}, points = {{0, 24}, {0, -24}}), Line(origin = {0, 36}, points = {{-22, -10}, {-22, 10}, {22, 10}, {22, -10}, {22, -10}}), Line(origin = {0, -36}, points = {{-22, 10}, {-22, -10}, {22, -10}, {22, 10}})}, coordinateSystem(initialScale = 0.1)));
end Thermocline_HBS_LC_2P_MixedOutlet;