within SolarTherm.Models.CSP.CRS.Receivers;

model Basic_Heater
  extends Interfaces.Models.Electric_Heater;
  parameter Modelica.SIunits.Temperature T_cold_set = from_degC(290);
  parameter Modelica.SIunits.Temperature T_hot_set = from_degC(565);
  parameter Medium.ThermodynamicState state_cold_set = Medium.setState_pTX(Medium.p_default, T_cold_set) "Cold fluid thermodynamic state at design";
  parameter Medium.ThermodynamicState state_hot_set = Medium.setState_pTX(Medium.p_default, T_hot_set) "Hot fluid thermodynamic state at design";
  parameter Modelica.SIunits.SpecificEnthalpy h_cold_set = Medium.specificEnthalpy(state_cold_set) "Cold fluid specific enthalpy at design";  
  parameter Modelica.SIunits.SpecificEnthalpy h_hot_set = Medium.specificEnthalpy(state_hot_set) "Hot fluid specific enthalpy at design";
  //parameter SI.Power P_renewable_des "Design maximum ";
  parameter Real eff_heater = 0.99 "Electrical-to-heat conversion efficiency of the heater";
  parameter SI.HeatFlowRate Q_heater_des = 600.0e6 "Design maximum heater heat-rate output (W_th)";
  parameter SI.Power P_heater_des = Q_heater_des/eff_heater "Design maximum electrical power input to the heater (W_e)";
  
  SI.SpecificEnthalpy h_in "Specific enthalpy of fluid entering the heater";
  SI.SpecificEnthalpy h_out(start=h_hot_set) "Specific enthalpy of fluid exiting the heater";  
  SI.HeatFlowRate Q_out "Heat-rate going into the fluid after curtailment signal from the system controller";
  SI.HeatFlowRate Q_out_raw "Heat-rate before curtailment signal from the system controller";
  SI.Power P_heater_out "Heater inlet power after P_supply is limited by P_heater_des";
        
  Modelica.Blocks.Interfaces.RealInput Q_curtail "Required curtailment heat-rate signal from the system controller" annotation (Placement(
        visible = true,transformation(
        
        origin={44,84},extent={{-12,-12},{12,12}},
        rotation=-90), iconTransformation(
        
        origin={-112, 38},extent={{12, -12}, {-12, 12}},
        rotation=180)));
    
  Modelica.Blocks.Interfaces.RealOutput Q_heater_raw "Heat-rate that the heater is currently able to deliver to the fluid" annotation(
    Placement(visible = true, transformation(origin = {108, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {111, 75}, extent = {{11, -11}, {-11, 11}}, rotation = 180)));

  Modelica.Blocks.Interfaces.BooleanInput curtail "Does the heater output need to be curtailed based on system controller?" annotation (Placement(
        visible = true,transformation(extent = {{-126, -88}, {-86, -48}}, rotation = 0),iconTransformation(extent = {{-124, 84}, {-100, 108}}, rotation = 0)));
        
  Medium.BaseProperties state_in "Inlet fluid thermodynamic state";
equation
  state_in.h = h_in;
  state_in.p = 1e5;
  
  fluid_a.m_flow + fluid_b.m_flow = 0.0;
  fluid_a.p = fluid_b.p;
  
  h_in=inStream(fluid_a.h_outflow);

  fluid_a.h_outflow = h_in;
  fluid_b.h_outflow = h_out;
  
  Q_heater_raw = Q_out_raw;
  Q_out_raw = P_heater_out*eff_heater;
  P_heater_out = min(P_supply,P_heater_des);
  
  if fluid_a.m_flow > 1e-3 then
    Q_out = (if curtail == true then min(Q_out_raw,Q_curtail) else Q_out_raw);
    h_out = h_in + Q_out/max(1.0e-3,fluid_a.m_flow);

  else
    Q_out = 0.0;
    h_out = h_hot_set;

  end if;

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
    Icon(coordinateSystem(initialScale = 0.1), graphics = {Text(origin = {10, 10}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {129, -146}}, textString = "%name")}));
end Basic_Heater;