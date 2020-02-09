within SolarTherm.Models.Control;
model SB_Control
  extends Icons.Control;
  //replaceable package HTF = SolarTherm.Media.Sodium.Sodium_ph;
  parameter Real T_df_on=1123 "Level of start defocus"; //recv
  parameter Real T_df_off=1098 "Level of stop defocus"; //recv
  parameter Real T_on=1048 "Temp of start discharge"; //PB
  parameter Real T_off=1023 "Temp of stop discharge"; //PB
  
  parameter Real Q_flow_PB_des = 298.5e6 "Design heat flow rate into power block";
  Real h_fg;
  Real PB_load;
  
  Modelica.Blocks.Interfaces.RealInput T_stor(start=1023) "Temperature of the HTF in storage"
    annotation (Placement(visible = true, transformation(extent = {{-128, -80}, {-88, -40}}, rotation = 0), iconTransformation(extent = {{-128, -80}, {-88, -40}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.RealInput Q_recv_in "Heat input to receiver"
    annotation (Placement(visible = true, transformation(extent = {{-128, 40}, {-88, 80}}, rotation = 0), iconTransformation(extent = {{-128, 40}, {-88, 80}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.BooleanInput helio_on(start=false) "Is heliostat on?"
    annotation (Placement(visible = true, transformation(extent = {{-128, -20}, {-88, 20}}, rotation = 0), iconTransformation(extent = {{-128, -20}, {-88, 20}}, rotation = 0)));

  Modelica.Blocks.Interfaces.BooleanOutput defocus(start=false) annotation (Placement(visible = true, transformation(origin = {114, -2}, extent = {{-24, -20}, {24, 20}}, rotation = 0), iconTransformation(origin = {114, -2}, extent = {{-24, -20}, {24, 20}}, rotation = 0)));
    
  Modelica.Blocks.Interfaces.RealOutput m_flow_PB(start=0.0) "Power block mass flow?" annotation (Placement(visible = true, transformation(extent = {{90, -82}, {130, -42}}, rotation = 0), iconTransformation(extent = {{90, -82}, {130, -42}}, rotation = 0))) ;
  
  Modelica.Blocks.Interfaces.RealOutput m_flow_recv(start=0.0) "Receiver mass flow?" annotation (Placement(visible = true, transformation(extent = {{90, 40}, {130, 80}}, rotation = 0), iconTransformation(extent = {{90, 40}, {130, 80}}, rotation = 0))) ;


algorithm
  //Recv
  when T_stor > T_df_on then
    defocus := true;
  end when;
  when T_stor < T_df_off then
    defocus := false;
  end when;

  //PB
  when T_stor > T_on then
    PB_load := 1.0;
  end when;
  when T_stor < T_off then
    PB_load := 0.0;
  end when;

equation
  h_fg = SolarTherm.Media.Sodium.Sodium_utilities.h_fg_T(T_stor);
  m_flow_recv = (Q_recv_in/h_fg)*(1.0/0.2);
  m_flow_PB = PB_load*(Q_flow_PB_des/h_fg);

end SB_Control;