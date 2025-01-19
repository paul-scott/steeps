within SolarTherm.Validation.Gan_HBS;

model HBS_6Layer_Tank

  extends SolarTherm.Interfaces.Models.StorageFluid_Thermocline;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  //Initialise Material Packages
  replaceable package Medium = SolarTherm.Validation.Gan_HBS.Media.Flue_Gas_107kPa;
  replaceable package Fluid_Package = SolarTherm.Materials.PartialMaterial;
  replaceable package Filler_Package_1 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_1;
  replaceable package Filler_Package_2 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_2;
  replaceable package Filler_Package_3 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_3;
  replaceable package Filler_Package_4 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_4;
  replaceable package Filler_Package_5 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_5;
  replaceable package Filler_Package_6 = SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_6;
  //replaceable package Encapsulation_Package = Filler_Package; //Defaults to filler material
  
  //Storage Parameter Settings
  parameter Integer Correlation = 1 "Interfacial convection correlation {1 = WakaoKaguei, 2 = MelissariArgyropoulos, 3 = Conservative}";
    //Storage CApacity
  parameter SI.Energy E_max = 144.0e9 "Maximum storage capacity";
  
    //Aspect ratios (H/D) of tank
  parameter Real ar = 2.0 "Aspect ratio of tank";
  
    //Porosity of tank filler materials
  parameter Real eta = 0.22 "Porosity";
  
    //Filler diameter of materials
  parameter SI.Length d_p = 0.02 "Diameter of the pores (m)";
  //parameter Real s_p = 0.04 "Separation of the pores (m)";

  
    //Encapsulation thickness
  //parameter SI.Length t_e = d_p/(2*N_p) "Encapsulation thickness"; //Defaults to equidistant radial discretization
  
    //Discretization Settings
  parameter Integer N_f = 100;
  parameter Integer N_p = 1;
  
  
  //Heat loss coefficient of tanks
  parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.0 "W/m2K";
  
  //Temperature Settings
  parameter SI.Temperature T_min = CV.from_deg(290) "Minimum temperature (design) also starting T";
  parameter SI.Temperature T_max = CV.from_deg(574) "Maximum design temperature (design)";


  //Input and Output Ports
  Modelica.Blocks.Interfaces.RealOutput T_top_measured "Temperature at the top of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          extent = {{40, 50}, {60, 70}}, rotation = 0), iconTransformation(origin = {45, 55}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput T_bot_measured "Temperature at the bottom of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          extent = {{40, -70}, {60, -50}}, rotation = 0), iconTransformation(origin = {45, -55}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
          
  Modelica.Blocks.Interfaces.RealOutput T_95_measured = SolarTherm.Models.Storage.Thermocline.Interpolate_Temperature(ZDH,T_f_degC,N_f,0.95) + 273.15 "Temperature at the 90% height of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          extent = {{40, 36}, {60, 56}}, rotation = 0), iconTransformation(origin = {45, 41}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput T_05_measured = SolarTherm.Models.Storage.Thermocline.Interpolate_Temperature(ZDH,T_f_degC,N_f,0.05) + 273.15"Temperature at the 10% height of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          extent = {{40, -54}, {60, -34}}, rotation = 0), iconTransformation(origin = {45, -39}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  
  Modelica.Blocks.Interfaces.RealOutput h_bot_outlet "Enthaply at the bottom of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          origin = {-40, -70},extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-35, -65}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
          
  Modelica.Blocks.Interfaces.RealOutput h_top_outlet "Enthaply at the top of the tank as an output signal (J/kg)" annotation (Placement(visible = true,transformation(
          origin = {-40, -74},extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-35, 65}, extent = {{5, -5}, {-5, 5}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealInput p_amb "Ambient Pressure" annotation (Placement(
        visible = true,transformation(
        
        origin={48, 8.88178e-16},extent={{10, -10}, {-10, 10}},
        rotation=0), iconTransformation(
        
        origin={-46, -20},extent={{-6, -6}, {6, 6}},
        rotation=0)));
  
  //Initialize Tank
  SolarTherm.Validation.Gan_HBS.HBS_6Layer Tank_A(redeclare replaceable package Fluid_Package = Fluid_Package, redeclare replaceable package Filler_Package_1 = Filler_Package_1, redeclare replaceable package Filler_Package_2 = Filler_Package_2, redeclare replaceable package Filler_Package_3 = Filler_Package_3, redeclare replaceable package Filler_Package_4 = Filler_Package_4, redeclare replaceable package Filler_Package_5 = Filler_Package_5, redeclare replaceable package Filler_Package_6 = Filler_Package_6, Correlation = Correlation, E_max = E_max, ar = ar, d_p = d_p, eta = eta, T_min = T_min, T_max = T_max, N_f = N_f, U_loss_tank = U_loss_tank);


  //Cost BreakDown
  parameter Real C_filler = Tank_A.C_filler;
  parameter Real C_fluid = Tank_A.C_fluid;
  parameter Real C_total = Tank_A.C_section;
  parameter Real C_tank = Tank_A.C_tank;
  parameter Real C_insulation = Tank_A.C_insulation;
  parameter Real C_encapsulation = Tank_A.C_encapsulation;
  //Theoretical Tank Level
  Modelica.Blocks.Interfaces.RealOutput Level "Theoretical Tank Level"
                                          annotation (Placement(visible = true,transformation(
          extent = {{40, 16}, {60, 36}}, rotation = 0), iconTransformation(origin = {45, 25}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));

  //Tank Non-dimensionalized vertical axis
  parameter Real ZDH[N_f] = Tank_A.ZDH;
  
  //Plotting Temperature degC
  Real T_f_degC[N_f](start=fill(T_min,N_f));
  
  //Analysis of fluid entering and exiting storage
  Fluid_Package.State fluid_top "Fluid entering/exiting top";
  Fluid_Package.State fluid_bot "Fluid entering/exiting bottom";
  Modelica.Blocks.Interfaces.RealInput T_amb annotation(
    Placement(visible = true, transformation(origin = {-50, 8.88178e-16}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-46, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  
equation
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
  T_f_degC = (Tank_A.T_f).-273.15;
  
  //Calculate tank energy level
  Level = Tank_A.Level;
  
  //Determine tank outlet enthalpy used by external control system
  h_bot_outlet = Tank_A.h_f[1];
  h_top_outlet = Tank_A.h_f[N_f];
  
  //Mass balance
  fluid_a.m_flow = -1.0*fluid_b.m_flow; //always true for a steady state component
  if fluid_a.m_flow > 0.0 then //mass is flowing into the top, direction is downwards so Tank_A.m_flow is (negative), charging
    Tank_A.m_flow = -1.0*fluid_a.m_flow;
    Tank_A.h_in = inStream(fluid_a.h_outflow);
    fluid_a.h_outflow = Tank_A.h_in;//Tank_A.h_f[N_f];//fluid_top.h;
    fluid_b.h_outflow = Tank_A.h_out;//Tank_A.h_f[1];//fluid_bot.h;
  else //discharging
    Tank_A.m_flow = -1.0*fluid_a.m_flow;
    Tank_A.h_in = inStream(fluid_b.h_outflow);
    fluid_a.h_outflow = Tank_A.h_out;//Tank_A.h_f[N_f];//fluid_top.h;
    fluid_b.h_outflow = Tank_A.h_in;//Tank_A.h_f[1];//fluid_bot.h;
  end if;

  fluid_a.p = p_amb;
  fluid_b.p = p_amb;
  T_amb = Tank_A.T_amb;
  T_top_measured = Tank_A.T_f[N_f];
  T_bot_measured = Tank_A.T_f[1];

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Text(origin = {-60, 31}, extent = {{-16, 7}, {6, -3}}, textString = "T_amb"), Text(origin = {-60, -10}, extent = {{-16, 10}, {8, -6}}, textString = "p_amb"), Text(origin = {-66, -67}, extent = {{-26, 5}, {28, -3}}, textString = "h_bot_outlet"), Text(origin = {66, -37}, extent = {{-18, 5}, {42, -25}}, textString = "T_bot_measured"), Text(origin = {69, 74}, extent = {{-15, 4}, {41, -20}}, textString = "T_top_measured"), Text(origin = {22, 84}, extent = {{-12, 4}, {24, -12}}, textString = "fluid_a"), Text(origin = {24, -76}, extent = {{-12, 4}, {26, -12}}, textString = "fluid_b"), Text(origin = {60, 51.5}, extent = {{-6, 2.5}, {14, -9.5}}, textString = "T_95%"), Text(origin = {54, -24.5}, extent = {{-6, 2.5}, {16, -11.5}}, textString = "T_05%"), Text(origin = {62, 33}, extent = {{-8, 3}, {12, -7}}, textString = "Level"), Text(origin = {-66, 65}, extent = {{-28, 5}, {28, -3}}, textString = "h_top_outlet"), Text(origin = {208, -136}, lineColor = {0, 0, 255}, extent = {{-150, 150}, {-20, 116}}, textString = "%name"), Rectangle(origin = {8, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-48, 60}, {32, -60}}), Polygon(origin = {0, 63}, fillColor = {190, 190, 190}, fillPattern = FillPattern.Solid, points = {{-16, 3}, {18, 3}, {40, -3}, {18, -9}, {-16, -9}, {-40, -3}, {-16, 3}, {-16, 3}}), Polygon(origin = {0, -51}, fillColor = {190, 190, 190}, fillPattern = FillPattern.Solid, points = {{-16, 3}, {16, 3}, {40, -3}, {16, -9}, {-16, -9}, {-40, -3}, {-16, 3}, {-16, 3}}), Polygon(origin = {-10, 42}, fillColor = {190, 190, 190}, fillPattern = FillPattern.Solid, points = {{-30, 6}, {-6, 0}, {28, 0}, {50, 6}, {50, -102}, {28, -108}, {-6, -108}, {-30, -102}, {-30, 6}}), Polygon(origin = {0, 40}, fillColor = {177, 177, 177}, fillPattern = FillPattern.Solid, points = {{-40, 20}, {-16, 14}, {18, 14}, {40, 20}, {40, 0}, {18, -6}, {-16, -6}, {-40, 0}, {-40, 20}}), Polygon(origin = {0, 22}, fillColor = {177, 177, 177}, fillPattern = FillPattern.Solid, points = {{-40, 20}, {-16, 14}, {18, 14}, {40, 20}, {40, 0}, {18, -6}, {-16, -6}, {-40, 0}, {-40, 20}}), Polygon(origin = {0, 2}, fillColor = {177, 177, 177}, fillPattern = FillPattern.Solid, points = {{-40, 20}, {-16, 14}, {18, 14}, {40, 20}, {40, 0}, {18, -6}, {-16, -6}, {-40, 0}, {-40, 20}}), Polygon(origin = {0, -18}, fillColor = {177, 177, 177}, fillPattern = FillPattern.Solid, points = {{-40, 20}, {-16, 14}, {18, 14}, {40, 20}, {40, 0}, {18, -6}, {-16, -6}, {-40, 0}, {-40, 20}}), Polygon(origin = {0, -40}, fillColor = {177, 177, 177}, fillPattern = FillPattern.Solid, points = {{-40, 22}, {-16, 16}, {18, 16}, {40, 22}, {40, 2}, {18, -4}, {-16, -4}, {-40, 2}, {-40, 22}}), Ellipse(origin = {0, 61}, fillColor = {208, 213, 255}, fillPattern = FillPattern.Solid, extent = {{-6, -2}, {6, 0}}, endAngle = 360), Ellipse(origin = {-24, 61}, fillColor = {208, 213, 255}, fillPattern = FillPattern.Solid, extent = {{-6, -2}, {6, 0}}, endAngle = 360), Ellipse(origin = {24, 61}, fillColor = {208, 213, 255}, fillPattern = FillPattern.Solid, extent = {{-6, -2}, {6, 0}}, endAngle = 360), Ellipse(origin = {10, 65}, fillColor = {208, 213, 255}, fillPattern = FillPattern.Solid, extent = {{-6, -2}, {6, 0}}, endAngle = 360), Ellipse(origin = {10, 57}, fillColor = {208, 213, 255}, fillPattern = FillPattern.Solid, extent = {{-6, -2}, {6, 0}}, endAngle = 360), Ellipse(origin = {-10, 65}, fillColor = {208, 213, 255}, fillPattern = FillPattern.Solid, extent = {{-6, -2}, {6, 0}}, endAngle = 360), Ellipse(origin = {-10, 57}, fillColor = {208, 213, 255}, fillPattern = FillPattern.Solid, extent = {{-6, -2}, {6, 0}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)),
    Documentation(revisions = "<html>
		<p>By Zebedee Kee on 03/12/2020</p>
		</html>", info = "<html>
		<p>This model contains the fluid_a (top) and fluid_b (bottom) ports, basically a complete CSP component. This model simply connects the Thermocline_Spheres_Section models to the correct ports.</p>
		</html>"));
end HBS_6Layer_Tank;