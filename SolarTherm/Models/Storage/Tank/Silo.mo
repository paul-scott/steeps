within SolarTherm.Models.Storage.Tank;

model Silo "A silo (from Ancient Greek σιρός (sirós) 'pit for holding grain') is a structure for storing bulk materials."
  extends Interfaces.Models.StorageFluid;
  replaceable package Filler_Package = SolarTherm.Materials.PartialMaterial "Filler Package";
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  Medium.BaseProperties medium;
  //Input Parameters
  parameter SI.Diameter D_tank=18.667 "Diameter";
  parameter SI.Height H_tank=20 "Height";
  parameter Real epsilon=0.26 "Filler void fraction"; //0.26 dense packing, 0.375-0.391 poured random packing, 0.359-0.375 close random packing.
  parameter SI.Temperature T_min = CV.from_degC(600) "Design cold Temperature of everything in the tank (K)";
  parameter SI.Temperature T_max = CV.from_degC(900) "Design hot Temperature of everything in the tank (K)";
  parameter SI.CoefficientOfHeatTransfer U_loss_tank=0.40 "Overall heat transfer coefficient to the ambient (W/m2K)";
  parameter SI.Temperature T_start = CV.from_degC(900) "Starting temperature of everything in the tank (K)";
  parameter SI.Temperature T_set=from_degC(800) "Tank Heater Temperature Set-Point";
  parameter Real L_start = 0.7 "Starting level of the packed bed, in decimal";
  parameter SI.Power P_max= 30e8 "Rating of auxilliary electrical heating (W)";
  parameter Real eff_heater=0.99 "Electrical-to-heat conversion efficieny of auxiliary heater";
  
  //Calculated Parameters
  parameter SI.Volume V_tank=0.25*CN.pi*D_tank*D_tank*H_tank "Volume of the cylindrical tank (m3)";
  parameter SI.Area A_cs=0.25*CN.pi*D_tank*D_tank "Cross sectional area of the tank (m2)";
  parameter SI.Area A_cyl=CN.pi*D_tank*H_tank "Cylindrical area of the tank (m2)";
  
  
  //Calculated filler properties
  parameter SI.SpecificEnthalpy h_p_max = Filler_Package.h_Tf(T_max, 1.0);
  parameter SI.SpecificEnthalpy h_p_min = Filler_Package.h_Tf(T_min, 0.0);
  parameter SI.Density rho_p_min = Filler_Package.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p_max = Filler_Package.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p = min(rho_p_min,rho_p_max) "Worst-case min. filler density (kg/m3)";
  
  //Variables
  Modelica.Blocks.Interfaces.RealOutput L=H_bed/H_tank "Level in decimals based on height (m/m)" annotation(
    Placement(visible = true, transformation(extent = {{96, 44}, {116, 64}}, rotation = 0), iconTransformation(extent = {{92, 48}, {112, 68}}, rotation = 0)));
  SI.Volume V_bed=V_p/(1.0-epsilon) "Volume of the packed bed, including the voids (m3)";
  SI.Volume V_p=m_p/rho_p "Volume of just the filler material (m3)";
  SI.Mass m_p(start=L_start*V_tank*rho_p*(1.0-epsilon)) "Mass of filler material (kg)";
  SI.SpecificEnthalpy h_p(start=Filler_Package.h_Tf(T_start, 1.0)) "Specific enthalpy of filler material (kg)";
  SI.Height H_bed=4.0*V_bed/(CN.pi*D_tank*D_tank) "Height of the packed bed (m)";
  
  SI.SpecificEnthalpy h_in "Inlet enthalpy";
  
  //Inlet and outlet enthalpies and temperatures
  //SI.SpecificEnthalpy h_in "Inlet Enthalpy depends on mass flow direction (J/kg)";
  //SI.SpecificEnthalpy h_out "Outlet Enthalpy depends on mass flow direction (J/kg)";
  //SI.Temperature T_in "Inlet Temperature depends on mass flow direction";
  //SI.Temperature T_out "Outlet Temperature depends on mass flow direction";
  SI.Power P_aux "Instantaneous power consumed by the auxiliary heater (W)";
  SI.HeatFlowRate Q_aux "Heat-rate delivered by the auxiliary heater (J/s)";
  SI.HeatFlowRate Q_loss "Instantaneous heat-loss-rate to the surroundings (J/s)";
  
  Real der_mh "Time derivative of the product of mass and enthalpy (J/s)";

  Modelica.Blocks.Interfaces.RealOutput h_fluid = h_p "Enthalpy of fluid (J/kg)" annotation (Placement(visible = true, transformation(origin = {107, 3}, extent = {{-11, -11}, {11, 11}}, rotation = 0), iconTransformation(origin = {102, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0))) ;
  
  Modelica.Blocks.Interfaces.RealOutput T_fluid = medium.T "Temperature of fluid (K)" annotation (Placement(visible = true, transformation(origin = {107, 3}, extent = {{-11, -11}, {11, 11}}, rotation = 0), iconTransformation(origin = {102, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0))) ;

  Modelica.Blocks.Interfaces.RealInput p "Pressure input signal" annotation (Placement(
        transformation(
        extent={{-16,-16},{16,16}},
        rotation=-90,
        origin={0,108}), iconTransformation(
        extent={{-11,-11},{11,11}},
        rotation=-90,
        origin={45,97})));

  Modelica.Blocks.Interfaces.RealInput T_amb "Ambient temperature input signal(K)" annotation (Placement(
        transformation(
        extent={{-16,-16},{16,16}},
        rotation=-90,
        origin={0,108}), iconTransformation(
        extent={{-11,-11},{11,11}},
        rotation=-90,
        origin={-41,97})));
  
  Filler_Package.State Filler_State;
algorithm
  when T_fluid < T_set then 
	P_aux:=min(Q_loss,P_max)*1.05;
  elsewhen T_fluid>T_set+1 then
	P_aux:=0.0;
  end when;
  der_mh := fluid_a.m_flow*h_in + fluid_b.m_flow*h_p + Q_aux - Q_loss;

equation
  h_in = inStream(fluid_a.h_outflow);

  medium.h = h_p;
  medium.p = p;
  
  Q_loss = (L*A_cyl+A_cs)*U_loss_tank*(T_fluid-T_amb);
  
  fluid_a.p = medium.p;
  fluid_b.p = medium.p;

  fluid_a.h_outflow=medium.h;
  fluid_b.h_outflow=medium.h;
  
  der(m_p)=fluid_a.m_flow+fluid_b.m_flow;
    
  der(h_p)=(der_mh - h_p*der(m_p))/m_p;
  
  Q_aux = eff_heater*P_aux;
  
  Filler_State.h = medium.h;

  annotation (Icon(coordinateSystem(preserveAspectRatio = false, initialScale = 0.1), graphics = {Text(origin = {3, -5}, lineColor = {255, 255, 255}, fillColor = {255, 255, 255}, extent = {{-55, 37}, {55, -37}}, textString = "Silo"), Text(origin = {14, 0}, lineColor = {0, 0, 255}, extent = {{-149, -114}, {129, -146}}, textString = "%name")}), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
<ul>
<li>Zebedee Kee:<br>Released first version. </li>
</ul>
</html>"));
end Silo;