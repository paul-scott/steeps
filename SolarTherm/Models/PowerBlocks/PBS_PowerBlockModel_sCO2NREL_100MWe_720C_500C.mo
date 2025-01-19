within SolarTherm.Models.PowerBlocks;
model PBS_PowerBlockModel_sCO2NREL_100MWe_720C_500C
  extends Interfaces.Models.PowerBlock;

  replaceable model Cycle =
     SolarTherm.Models.PowerBlocks.Correlation.sCO2NREL_100MWe_720C_500C constrainedby
    SolarTherm.Models.PowerBlocks.Correlation.Cycle
    annotation (Dialog(group="Regresion"),choicesAllMatching=true);
  Cycle cycle(T_in=T_in,load=load,T_amb=T_amb_internal,final T_in_ref=T_in_ref);
  
  parameter SI.HeatFlowRate W_des=100.0e6/nu_net "Design turbine gross output" annotation (Dialog(group="Design"));
  parameter Real nu_net=0.9 "Estimated gross to net conversion factor at the power block" annotation(Dialog(group="Parasities energy losses"));
  parameter SI.Temperature T_in_ref=from_degC(720) "HTF inlet temperature (design)" annotation (Dialog(group="Design"));
  parameter SI.Temperature T_out_ref=from_degC(500) "HTF outlet temperature (design)"
                                                                                     annotation (Dialog(group="Design"));
  parameter SI.AbsolutePressure p_bo=10e5 "Boiler operating pressure" annotation (Dialog(group="Design"));
  parameter SI.HeatFlowRate Q_flow_ref=232.109832e6 "Design thermal power" annotation (Dialog(group="Design"));
  parameter SI.MassFlowRate m_flow_ref=842.224 "Design HTF flow rate";
  
  parameter Real W_base=0.0 "Power consumed at all times" annotation(Dialog(group="Parasities energy losses"));
  parameter Real nu_min=0.60 "Minimum turbine operation" annotation (Dialog(group="Operating strategy"));
  SI.HeatFlowRate Q_flow( final start=0) "Cycle heat addition";
  SI.Temperature T_in;//=Medium.temperature(state_in);
  SI.Temperature T_out;//=Medium.temperature(state_out);
  
  parameter Boolean capped = false "Is the net output capped to the design rated output (not crediting excess power)";
  parameter SI.HeatFlowRate W_cap = W_des*nu_net "The design net rating of the PB";
  
  parameter Boolean enable_losses = true
    "= true enable thermal losses with environment"
      annotation (Dialog(group="Assumptions"), Evaluate=true, HideResult=true, choices(checkBox=true));
  parameter SI.Temperature T_des=from_degC(43) "Ambient temperature at design" annotation (Dialog(group="Assumptions",enable = enable_losses));

  parameter Boolean external_parasities = false
    "= true enable parasities as an input"
      annotation (Dialog(group="Parasities energy losses"), Evaluate=true, HideResult=true, choices(checkBox=true));



  Modelica.Blocks.Interfaces.RealOutput h_out_signal "outlet specific enthalpy" annotation (Placement(visible = true, transformation(origin = {-80, -72}, extent = {{20, -20}, {-20, 20}}, rotation = 0), iconTransformation(origin = { -61, -59}, extent = {{5, -5}, {-5, 5}}, rotation = 0))) ;

  Modelica.Blocks.Interfaces.RealInput T_amb if enable_losses annotation (Placement(
        transformation(extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={1.77636e-015,80}),                  iconTransformation(
        extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={-20,60})));


  replaceable model Cooling = SolarTherm.Models.PowerBlocks.Cooling.NoCooling constrainedby
    SolarTherm.Models.PowerBlocks.Cooling.Cooling
    annotation (Dialog(group="Cooling losses",enable = enable_losses),
      choicesAllMatching=true);
  Cooling cool(T_amb=T_amb_internal);
  Real load;
  SI.HeatFlowRate W_gross "Parasitic losses power";
  SI.HeatFlowRate W_loss "Parasitic losses power";
  //SI.HeatFlowRate W_net "Estimated net output at design";
  SI.Energy E_gross(final start=0, fixed=true, displayUnit="MW.h");
  SI.Energy E_net(final start=0, fixed=true, displayUnit="MW.h");

  Boolean logic;
  parameter Real C_PB_total = 143591149.0 "Total PB Cost (USD)";

   Modelica.Blocks.Interfaces.RealInput parasities if external_parasities annotation (Placement(
        transformation(extent={{-12,-12},{12,12}},
        rotation=-90,
        origin={1.77636e-015,80}),                  iconTransformation(
        extent={{-6,-6},{6,6}},
        rotation=-90,
        origin={20,60})));
  Real eff_PB_measured = k_q*k_w;
protected
  Modelica.Blocks.Interfaces.RealInput parasities_internal;
  Real k_q "In this case is the eta_Q of HX";
  Real k_w "In this case is the eta_gross of the power cycle";
  SI.SpecificEnthalpy h_in;
  SI.SpecificEnthalpy h_out;
  
  //Medium.ThermodynamicState state_in=Medium.setState_phX(fluid_a.p,inStream(fluid_a.h_outflow));
  //Medium.ThermodynamicState state_out=Medium.setState_phX(fluid_a.p,h_out);
  Medium.BaseProperties state_in;
  Medium.BaseProperties state_out;
  parameter Medium.ThermodynamicState state_in_ref=Medium.setState_pTX(1e5,T_in_ref);
  parameter Medium.ThermodynamicState state_out_ref=Medium.setState_pTX(1e5,T_out_ref);
  parameter SI.SpecificEnthalpy h_in_ref=Medium.specificEnthalpy(state_in_ref);
  parameter SI.SpecificEnthalpy h_out_ref=Medium.specificEnthalpy(state_out_ref);

  Modelica.Blocks.Interfaces.RealInput T_amb_internal;

  parameter Real nu_eps=0.1;

equation
  //States
  state_in.h = h_in;
  h_out = state_out.h;
  T_in = state_in.T;
  state_out.T = T_out;
  state_in.p = 1e5;
  state_out.p = 1e5;
  
  if enable_losses then
    connect(T_amb_internal,T_amb);
  else
    T_amb_internal=T_des;
  end if;
  if external_parasities then
    connect(parasities_internal,parasities);
  else
    parasities_internal=0;
  end if;

  logic=load>nu_min and fluid_a.m_flow > 1e-6;
  h_in=inStream(fluid_a.h_outflow);

 
  h_out_signal = h_out_ref;//h_out_ref;
  fluid_b.h_outflow = h_out_ref;//h_out_ref;//h_out;//
  
  fluid_a.h_outflow = h_in;//h_in_ref;//h_out_ref;
  fluid_a.m_flow+fluid_b.m_flow=0;
  fluid_a.p=fluid_b.p;

  load=max(nu_eps,fluid_a.m_flow/m_flow_ref); 

  if fluid_a.m_flow > 1e-6 then
    k_q=max(0.1,cycle.k_q);
    k_w=max(0.1,cycle.k_w);
    Q_flow = k_q*Q_flow_ref;
    h_out = h_in - Q_flow/fluid_a.m_flow;
    W_gross = k_w*Q_flow;
  else
    k_q=0.0;
    k_w=0.0;
    h_out=h_out_ref;
    Q_flow=fluid_a.m_flow*(h_in-h_out);
    W_gross=0.0;
  end if;

  //Q_flow/(cool.nu_q*Q_flow_ref*load)=k_q;
  //W_gross/(cool.nu_w*W_des*load)=k_w;

  der(E_gross)=W_gross;
  der(E_net)=W_net;
  W_loss=(1-nu_net)*W_gross+W_base+parasities_internal;
  if capped == true then
    W_net = min(W_cap,W_gross - W_loss);
  else
    W_net = W_gross - W_loss;
  end if;

end PBS_PowerBlockModel_sCO2NREL_100MWe_720C_500C;