within SolarTherm.Models.CSP.CRS.Receivers;

model Basic_Heater_Energy
  //This model is an energy-based heater model used for direct charging of a radiatively heated TES component"
  parameter Real eff_heater = 0.99 "Electrical-to-heat conversion efficiency of the heater";
  parameter SI.HeatFlowRate Q_heater_des = 600.0e6 "Design maximum heater heat-rate output (W_th)";
  parameter SI.Power P_heater_des = Q_heater_des/eff_heater "Design maximum electrical power input to the heater (W_e)";
  SI.Power P_heater_out "Heater inlet power after P_supply is limited by P_heater_des";

  Modelica.Blocks.Interfaces.RealInput P_supply "Supplied electrical power" annotation(
    Placement(visible = true, transformation(origin = {108, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-111, -1}, extent = {{11, -11}, {-11, 11}}, rotation = 180)));
    
  Modelica.Blocks.Interfaces.RealOutput Q_out_raw "Heat-rate that the heater is currently able to output" annotation(
    Placement(visible = true, transformation(origin = {108, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {33, 43}, extent = {{11, -11}, {-11, 11}}, rotation = 180)));
    
  Modelica.Blocks.Interfaces.RealOutput Q_out "Actual heater heat output after considering on/off state" annotation(
    Placement(visible = true, transformation(origin = {108, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {111, -1}, extent = {{11, -11}, {-11, 11}}, rotation = 180)));

  Modelica.Blocks.Interfaces.BooleanInput on annotation (Placement(
        visible = true,transformation(extent = {{-126, -88}, {-86, -48}}, rotation = 0),iconTransformation(origin = {2, -112},extent = {{-12, -12}, {12, 12}}, rotation = 90)));

equation
  if on then
    Q_out = max(Q_out_raw,1.0e-6);
  else
    Q_out = 1.0e-6;
  end if;
  
  Q_out_raw = P_heater_out*eff_heater;
  P_heater_out = min(P_supply,P_heater_des);
 

  annotation (Documentation(info = "<html>
</html>", revisions = "<html>
<ul>
<li>by Zebedee Kee October 2023 </li>
</ul>
</html>"),
    Icon(graphics = {Line(origin = {-60, -2}, points = {{-40, 0}, {40, 0}}, thickness = 0.5), Line(origin = {0.62, 0.31}, points = {{6.00532, 33.005}, {-9.99468, 1.00502}, {10.0053, 1.00502}, {-7.99468, -32.995}}, thickness = 1.25), Line(origin = {61, -2}, points = {{-39, 0}, {39, 0}}, thickness = 0.5), Rectangle(origin = {2, 0}, lineThickness = 0.5, extent = {{-22, 54}, {20, -54}}), Line(origin = {2, -76}, points = {{0, -24}, {0, 22}}, color = {255, 0, 255}, pattern = LinePattern.Dash, thickness = 0.5)}, coordinateSystem(initialScale = 0.1)));
end Basic_Heater_Energy;