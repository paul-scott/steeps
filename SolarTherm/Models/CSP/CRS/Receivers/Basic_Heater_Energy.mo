within SolarTherm.Models.CSP.CRS.Receivers;
model Basic_Heater_Energy
  //This model is an energy-based heater model used for direct charging of a radiatively heated TES component"
  parameter Real eff_heater = 0.99 "Electrical-to-heat conversion efficiency of the heater";
  parameter SI.HeatFlowRate Q_heater_des = 600.0e6 "Design maximum heater heat-rate output (W_th)";
  parameter SI.Power P_heater_des = Q_heater_des/eff_heater "Design maximum electrical power input to the heater (W_e)";
  SI.Power P_heater_out "Heater inlet power after P_supply is limited by P_heater_des";
  //SI.HeatFlowRate Q_out_raw "Heat-rate before curtailment signal from the system controller";

    
  Modelica.Blocks.Interfaces.RealInput P_supply "Supplied electrical power" annotation(
    Placement(visible = true, transformation(origin = {108, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-111, -1}, extent = {{11, -11}, {-11, 11}}, rotation = 180)));
    
  Modelica.Blocks.Interfaces.RealOutput Q_out_raw "Heat-rate that the heater is currently able to output" annotation(
    Placement(visible = true, transformation(origin = {108, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {33, 43}, extent = {{11, -11}, {-11, 11}}, rotation = 180)));
    
  Modelica.Blocks.Interfaces.RealOutput Q_out "Actual heater heat output after considering on/off state" annotation(
    Placement(visible = true, transformation(origin = {108, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {111, -1}, extent = {{11, -11}, {-11, 11}}, rotation = 180)));

  Modelica.Blocks.Interfaces.BooleanInput on annotation (Placement(
        visible = true,transformation(extent = {{-126, -88}, {-86, -48}}, rotation = 0),iconTransformation(origin = {2, -112},extent = {{-12, -12}, {12, 12}}, rotation = 90)));
        
  //parameter SI.Area A_recv "Area of receiver";
  //parameter Medium.ThermodynamicState state_0=Medium.setState_pTX(1e5,T_0);
  //parameter SI.SpecificEnthalpy h_0=Medium.specificEnthalpy(state_0);
  //parameter SI.Temperature T_0=from_degC(290) "Start value of temperature";
  //parameter SI.Temperature T_max=from_degC(700) "Target outlet temperature";
  //parameter SI.CoefficientOfHeatTransfer h_conv=20.0 "Convective loss coefficient (W/m2K)";

  //SI.Temperature T_in;
  //SI.Temperature T_out;
  //parameter SI.SpecificEnthalpy h_out_ref = Medium.specificEnthalpy(Medium.setState_pTX(1e5,T_max));

//algorithm
  //T_4avg := (0.20*(T_in^4 + (T_in^3)*(T_out) + (T_in^2)*(T_out^2) + (T_in)*(T_out^3) + T_out^4))^0.25;
  //T_avg := 0.5*(T_in+T_out);
  //Medium.BaseProperties state_in;
equation
  if on then
    Q_out = max(Q_out_raw,1.0e-6);
  else
    Q_out = 1.0e-6;
  end if;
  
  Q_out_raw = P_heater_out*eff_heater;
  P_heater_out = min(P_supply,P_heater_des);
  
  //state_in.h = h_in;
  //state_in.p = 1e5;
  
  //fluid_a.m_flow + fluid_b.m_flow = 0.0;
  //fluid_a.p = fluid_b.p;
  //heat.T=0.5*(T_cold_set+T_hot_set);
  
  //h_in=inStream(fluid_a.h_outflow);

  //fluid_a.h_outflow = h_in;
  //fluid_b.h_outflow = h_out;
  
  //Q_heater_raw = Q_in_raw;
  //Q_in_raw = P_supply;
  //medium_in.h = h_in;
  //T_in = medium_in.T;
  //T_out = medium_out.T;
  //medium_in.p = fluid_a.p;
  
  //medium.h=0.5*(h_in+h_out);
  
  //New
  //T_4avg = (0.20*(T_in^4 + (T_in^3)*(T_out) + (T_in^2)*(T_out^2) + (T_in)*(T_out^3) + T_out^4))^0.25;
  //T_avg = 0.5*(T_in+T_out);
  //Q_rad = A_recv*sigma*em*((T_4avg^4)-(Tamb^4));
  //Q_conv = A_recv*h_conv*(T_avg-Tamb);
  //End New
  
  //Q_rad = A_recv*sigma*em*((medium.T^4)-(Tamb^4));
  //Q_conv = A_recv*h_conv*(medium.T-Tamb);
  

  
  //Q_rcv_raw = ab*heat.Q_flow-Q_rad-Q_conv; //Theoretical net receiver output before curtailment
  

  
  //medium_out.h = h_out;
  //medium_out.p = medium.p;
  annotation (Documentation(info = "<html>
</html>", revisions = "<html>
<ul>
<li>by Zebedee Kee, Receiver used in Packed Bed Storage Annual System Model. Only Radiative Losses. </li>
</ul>
</html>"),
    Icon(graphics = {Line(origin = {-60, -2}, points = {{-40, 0}, {40, 0}}, thickness = 0.5), Line(origin = {0.62, 0.31}, points = {{6.00532, 33.005}, {-9.99468, 1.00502}, {10.0053, 1.00502}, {-7.99468, -32.995}}, thickness = 1.25), Line(origin = {61, -2}, points = {{-39, 0}, {39, 0}}, thickness = 0.5), Rectangle(origin = {2, 0}, lineThickness = 0.5, extent = {{-22, 54}, {20, -54}}), Line(origin = {2, -76}, points = {{0, -24}, {0, 22}}, color = {255, 0, 255}, pattern = LinePattern.Dash, thickness = 0.5), Text(origin = {-114, 35}, extent = {{-36, 29}, {36, -29}}, textString = "P_supply"), Text(origin = {92, -15}, extent = {{-8, 5}, {36, -29}}, textString = "Q_out"), Text(origin = {-46, -99}, extent = {{-8, 5}, {30, -17}}, textString = "on"), Text(origin = {30, 89}, extent = {{-8, 5}, {66, -41}}, textString = "Q_out_raw")}, coordinateSystem(initialScale = 0.1)));
end Basic_Heater_Energy;