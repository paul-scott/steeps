within SolarTherm.Models.Storage.Thermocline;

model Thermocline_HBS_LC_SingleTank_NoVerticalConduction

  extends SolarTherm.Interfaces.Models.StorageFluid_Thermocline;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  //Initialise Material Packages
  replaceable package Medium = SolarTherm.Media.Sodium.Sodium_pT;
  replaceable package Fluid_Package = SolarTherm.Materials.PartialMaterial;
  replaceable package Filler_Package = SolarTherm.Materials.PartialMaterial;
  
  //Storage Parameter Settings
  parameter Integer Correlation = 1 "Interfacial convection correlation {1 = WakaoKaguei, 2 = MelissariArgyropoulos, 3 = Conservative}";
    //Storage CApacity
  parameter SI.Energy E_max = 144.0e9 "Maximum storage capacity";
  
    //Aspect ratios (H/D) of tank
  parameter Real ar = 2.0 "Aspect ratio of tank";
  
    //Porosity of tank filler materials
  parameter Real eta = 0.51 "Porosity";
  
    //Filler diameter of materials
  parameter SI.Length d_p = 0.03 "Channel hydraulic diameter";
  
    //Encapsulation thickness
  //parameter SI.Length t_e = d_p/(2*N_p) "Encapsulation thickness"; //Defaults to equidistant radial discretization
  
    //Discretization Settings
  parameter Integer N_f = 10;
  //parameter Integer N_p = 5;
  
  
  //Heat loss coefficient of tanks
  //parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.1 "W/m2K";
  parameter SI.CoefficientOfHeatTransfer U_loss_top = 0.698 "Heat loss coefficient at the top of the tank (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_bot = 1.22 "Heat loss coefficient at the bottom of the tank (W/m2K)";
  
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
          
  Modelica.Blocks.Interfaces.RealOutput T_95_measured = Interpolate_Temperature(ZDH,T_f_degC,N_f,0.95) + 273.15 "Temperature at the 90% height of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          extent = {{40, 36}, {60, 56}}, rotation = 0), iconTransformation(origin = {45, 41}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput T_05_measured = Interpolate_Temperature(ZDH,T_f_degC,N_f,0.05) + 273.15"Temperature at the 10% height of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          extent = {{40, -54}, {60, -34}}, rotation = 0), iconTransformation(origin = {45, -39}, extent = {{-5, -5}, {5, 5}}, rotation = 0)));
  
  Modelica.Blocks.Interfaces.RealOutput h_bot_outlet "Enthaply at the bottom of the tank as an output signal (K)"
                                          annotation (Placement(visible = true,transformation(
          origin = {-40, -70},extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-17, -73}, extent = {{-5, -5}, {5, 5}}, rotation = -90)));
          
  Modelica.Blocks.Interfaces.RealOutput h_top_outlet "Enthaply at the top of the tank as an output signal (J/kg)" annotation (Placement(visible = true,transformation(
          origin = {-40, -74},extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-17, 73}, extent = {{5, -5}, {-5, 5}}, rotation = -90)));

  Modelica.Blocks.Interfaces.RealInput T_amb "Ambient Temperature" annotation (Placement(
        visible = true,transformation(
        
        origin={-50, 8.88178e-16},extent={{-10, -10}, {10, 10}},
        rotation=0), iconTransformation(
        
        origin={-46, 20},extent={{-6, -6}, {6, 6}},
        rotation=0)));
  Modelica.Blocks.Interfaces.RealInput p_amb "Ambient Pressure" annotation (Placement(
        visible = true,transformation(
        
        origin={48, 8.88178e-16},extent={{10, -10}, {-10, 10}},
        rotation=0), iconTransformation(
        
        origin={-46, -20},extent={{-6, -6}, {6, 6}},
        rotation=0)));
  
  //Initialize Tank
  SolarTherm.Models.Storage.Thermocline.Thermocline_HBS_LC_Section_NoVerticalConduction Tank_A(redeclare replaceable package Fluid_Package = Fluid_Package, redeclare replaceable package Filler_Package = Filler_Package, Correlation = Correlation, E_max = E_max, ar = ar, eta = eta, d_p = d_p, T_min = T_min, T_max = T_max, N_f = N_f, U_loss_top = U_loss_top, U_loss_bot = U_loss_bot);


  //Cost BreakDown
  parameter Real C_filler = Tank_A.C_filler;
  parameter Real C_fluid = 0.0;
  parameter Real C_total = 0.0;
  parameter Real C_tank = Tank_A.C_tank;
  parameter Real C_insulation = Tank_A.C_insulation;
  parameter Real C_encapsulation = 0.0;
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
    Icon(graphics = {Text(origin = {-60, 31}, extent = {{-16, 7}, {6, -3}}, textString = "T_amb"), Text(origin = {-60, -10}, extent = {{-16, 10}, {8, -6}}, textString = "p_amb"), Text(origin = {-50, -77}, extent = {{-26, 5}, {28, -3}}, textString = "h_bot_outlet"), Text(origin = {66, -37}, extent = {{-18, 5}, {42, -25}}, textString = "T_bot_measured"), Text(origin = {69, 74}, extent = {{-15, 4}, {41, -20}}, textString = "T_top_measured"), Text(origin = {22, 84}, extent = {{-12, 4}, {24, -12}}, textString = "fluid_a"), Text(origin = {26, -76}, extent = {{-12, 4}, {26, -12}}, textString = "fluid_b"), Text(origin = {60, 51.5}, extent = {{-6, 2.5}, {14, -9.5}}, textString = "T_95%"), Text(origin = {54, -24.5}, extent = {{-6, 2.5}, {16, -11.5}}, textString = "T_05%"), Text(origin = {62, 33}, extent = {{-8, 3}, {12, -7}}, textString = "Level"), Text(origin = {-50, 77}, extent = {{-28, 5}, {28, -3}}, textString = "h_top_outlet"), Polygon(origin = {-31, 40}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{9, 28}, {53, 28}, {71, 20}, {71, -100}, {53, -108}, {9, -108}, {-9, -100}, {-9, 20}, {9, 28}}), Polygon(origin = {-30, 54}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {8, -2}, {52, -2}, {70, 6}, {70, -12}, {52, -20}, {8, -20}, {-10, -12}, {-10, 6}}), Polygon(origin = {-30, 36}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {8, -2}, {52, -2}, {70, 6}, {70, -12}, {54, -20}, {8, -20}, {-10, -12}, {-10, 6}}), Polygon(origin = {-30, 18}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {8, -2}, {52, -2}, {70, 6}, {70, -12}, {52, -20}, {8, -20}, {-10, -12}, {-10, 6}}), Polygon(origin = {-30, 0}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {8, -2}, {52, -2}, {70, 6}, {70, -12}, {52, -20}, {8, -20}, {-10, -12}, {-10, 6}}), Polygon(origin = {-30, -18}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {8, -2}, {52, -2}, {70, 6}, {70, -12}, {52, -20}, {8, -20}, {-10, -12}, {-10, 6}}), Polygon(origin = {-30, -36}, fillColor = {195, 195, 195}, fillPattern = FillPattern.Solid, points = {{-10, 6}, {8, -2}, {52, -2}, {70, 6}, {70, -10}, {52, -18}, {8, -18}, {-10, -10}, {-10, 6}}), Line(origin = {-22, -8}, points = {{0, 60}, {0, -60}}), Line(origin = {22, -8}, points = {{0, 60}, {0, -60}}), Ellipse(origin = {-4, 64}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-3.5, -1.5}, {11.5, -4.5}}, endAngle = 360), Line(origin = {0, 64}, points = {{0, 4}, {0, -4}}), Ellipse(origin = {-22, 68}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {13.5, -4.5}}, endAngle = 360), Ellipse(origin = {14, 68}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-5.5, -1.5}, {9.5, -4.5}}, endAngle = 360), Ellipse(origin = {-22, 60}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {13.5, -4.5}}, endAngle = 360), Ellipse(origin = {14, 60}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-5.5, -1.5}, {9.5, -4.5}}, endAngle = 360), Ellipse(origin = {-32, 64}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-1.5, -1.5}, {13.5, -4.5}}, endAngle = 360), Ellipse(origin = {24, 64}, fillColor = {211, 211, 255}, fillPattern = FillPattern.Solid, extent = {{-5.5, -1.5}, {9.5, -4.5}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)), Documentation(revisions ="<html>
		<p>By Zebedee Kee on 03/12/2020</p>
		</html>",info="<html>
		<p>This model contains the fluid_a (top) and fluid_b (bottom) ports, basically a complete CSP component. This model simply connects the Thermocline_Spheres_Section models to the correct ports.</p>
		</html>"));

end Thermocline_HBS_LC_SingleTank_NoVerticalConduction;