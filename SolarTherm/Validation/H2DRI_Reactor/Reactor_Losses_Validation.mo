within SolarTherm.Validation.H2DRI_Reactor;

model Reactor_Losses_Validation
  //0D Model
  //In this example, the molar ratio of H2 is controlled such that 100% yield of Fe-DRI is guaranteed at least according to equilibrium calculations.
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  parameter SI.Pressure p_des = 100000.0 "Design pressure (Pa)";
  //Default 1 bar
  //parameter Real m_flow_Fe_MTperYear = 1.5 "Mega Metric-Tonnes per year of Iron DRI (Mt/yr)";
  parameter SI.MassFlowRate m_flow_Fe_des = 2.0 * SolarTherm.Models.Chemistry.ChemTable.Fe.M "Design outlet mass flow rate of Fe-DRI assuming 100% conversion";
  parameter SI.MassFlowRate m_flow_Fe2O3_stoi = 1.4297448 * m_flow_Fe_des "Stoichiometric mass flow rate of Fe2O3 required to ensure the design Fe-DRI output assuming 100% conversion rate";
  parameter SI.MassFlowRate m_flow_H2_stoi = 0.05414666 * m_flow_Fe_des "Stoichiometric mass flow rate of H2 required to ensure the design Fe-DRI output assuming 100% conversion rate";
  parameter SI.Temperature T_inlet = 925.0 + 273.15 "Inlet temperature (K)";
  //parameter SI.Temperature T_inlet_base = 600.0+273.15;
  parameter SI.MolarEnthalpy Q_loss_per_mol = 25000.0 "Heat loss per mole of Fe2O3 (J/mol)";
  //SI.Temperature T_products_des(start = T_products_des_start) "Design intended product temperature (K)";
  //SI.Temperature T_inlet "Reactant temperature before losses estimated using polynomial fit (K)";
  //SI.Temperature T_reactants "Reactant temperature after losses estimated using polynomial fit (K)";
  SI.Temperature T_Fe2O3_in = T_inlet "Inlet Fe2O3 reactant temperature (K)";
  SI.Temperature T_H2_in = T_inlet "Inlet H2 reactant temperature (K)";
  //Real r_min "Minimum inlet molar ratio of H2 to Fe2O3 required to ensure 100% yield of Fe-DRI at least according to equilibrium calculations, estimated using a polynomial fit (-)";
  parameter SI.MolarFlowRate n_flow_Fe_max = m_flow_Fe_des / SolarTherm.Models.Chemistry.ChemTable.Fe.M "Maximum number of moles of Fe that can be produced per second (mol/s)";
  //Yield of each iron species in decimal form
  Real Yield_Fe = Reactor.n_flow_Fe_out / n_flow_Fe_max;
  Real Yield_FeO = Reactor.n_flow_FeO_out / n_flow_Fe_max;
  //Up to the number of moles of Fe that can be produced
  Real Yield_Fe3O4 = Reactor.n_flow_Fe3O4_out / (1.0 / 3.0 * n_flow_Fe_max);
  //Up to one-third the moles of Fe that can be produced
  Real Yield_Fe2O3 = Reactor.n_flow_Fe2O3_out / (1.0 / 2.0 * n_flow_Fe_max);
  //Up to half the moles of Fe that can be produced
  //Enthalpy of reaction per mole of Fe2O3 in
  //Real Delta_H_per_mol_Fe2O3 = Reactor.Q_flow_absorbed/Reactor.n_flow_Fe2O3_in;
  //SI.MolarEnthalpy H_rxn_permolFe2O3 = 2.0*SolarTherm.Models.Chemistry.H2DRI.Isothermal.Overall_Rxn_Enthalpy(T_reactants,p_des);
  //In the event there are reactor losses
  //Reactor efficiency = (Reaction enthalpy - Q_loss)/(Reaction enthalpy) where reaction enthalpy is calculated based on reactant temperature
  //SI.MolarEnthalpy Delta_H_per_mol_Fe2O3_check = 2.0*SolarTherm.Models.Chemistry.H2DRI.Isothermal.Overall_Rxn_Enthalpy(T_reactants,p_des);
  //Loss Signal
  //parameter SI.MolarEnthalpy Q_loss_per_mol = 30000.0 "Heat loss per mole of Fe2O3 (J/mol)";
  SI.HeatFlowRate Q_flow_loss = Q_loss_per_mol * Reactor.n_flow_Fe2O3_in;
  //Real eff_reactor = (SolarTherm.Models.Chemistry.H2DRI.Isothermal.Overall_Rxn_Enthalpy(Reactor.T_reactants_afterloss,p_des)-0.5*Q_loss_per_mol)/SolarTherm.Models.Chemistry.H2DRI.Isothermal.Overall_Rxn_Enthalpy(Reactor.T_reactants_afterloss,p_des);
  //Real f_loss = Q_flow_loss/(2.0*Reactor.n_flow_Fe2O3_in*SolarTherm.Models.Chemistry.H2DRI.Isothermal.Overall_Rxn_Enthalpy(Reactor.T_reactants,p_des));
  Modelica.Fluid.Sources.Boundary_pT Fe2O3_source(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph, nPorts = 1, p = p_des, use_T_in = true) annotation(
    Placement(visible = true, transformation(origin = {-107, -23}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT H2_Source(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2, nPorts = 1, p = p_des, use_T_in = true) annotation(
    Placement(visible = true, transformation(origin = {-109, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  //Storage Model
  //Componenets and Connectors
  Modelica.Blocks.Sources.RealExpression m_flow_Fe2O3(y = m_flow_Fe2O3_signal) annotation(
    Placement(visible = true, transformation(origin = {-107, -65}, extent = {{-19, -17}, {19, 17}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression m_flow_H2(y = m_flow_H2_signal) annotation(
    Placement(visible = true, transformation(origin = {-110, 73}, extent = {{-20, -19}, {20, 19}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression p_amb(y = 100000) annotation(
    Placement(visible = true, transformation(origin = {77, 54}, extent = {{13, -16}, {-13, 16}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Fe2O3_pump(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph) annotation(
    Placement(visible = true, transformation(origin = {-54, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple H2_pump(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(origin = {-53, 25}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  //Mass flow Signals starts in charging state
  SI.MassFlowRate m_flow_Fe2O3_signal(start = m_flow_Fe2O3_stoi);
  SI.MassFlowRate m_flow_H2_signal;
  SolarTherm.Models.Reactors.Reactor_H2DRI_Overall_Losses Reactor annotation(
    Placement(visible = true, transformation(origin = {38, 0}, extent = {{-48, -48}, {48, 48}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink H2_Sink(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(origin = {139, 9}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink H2O_Sink(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) annotation(
    Placement(visible = true, transformation(origin = {129, -53}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression T_Fe2O3_signal(y = T_Fe2O3_in) annotation(
    Placement(visible = true, transformation(origin = {-166, -16}, extent = {{-18, -16}, {18, 16}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression T_H2_signal(y = T_H2_in) annotation(
    Placement(visible = true, transformation(origin = {-166, 32}, extent = {{-18, -16}, {18, 16}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Heat_Loss(y = Q_flow_loss) annotation(
    Placement(visible = true, transformation(origin = {-1, 54}, extent = {{-17, -16}, {17, 16}}, rotation = 0)));
initial equation

equation
  m_flow_Fe2O3_signal = m_flow_Fe2O3_stoi;
//Fixed based on plant size
//Sweep design product temperature
//T_products_des = T_products_des_start + (time/100.0)*(T_products_des_end-T_products_des_start);
//Estimate required reactant temperature
//T_reactants = 1.0862878*T_products_des+51.45030;
//T_reactants = 65.15226 - (14.0897*eff_reactor) + (1209.891*(T_products_des/1000.0)) - (123.165*eff_reactor*(T_products_des/1000.0));
//T_inlet = 51.02315 + (1086.761*T_products_des/1000.0) - (0.00838*Q_loss_per_mol/1000.0) + 2.16376*(T_products_des*Q_loss_per_mol*1.0e-6);
//Estimate required inlet molar ratio of H2 to Fe2O3
//r_min = 6.062457e1 + ((-7.646067e-2)*T_reactants) + ((2.853147e-5)*T_reactants*T_reactants);
//Vary molar ratios
  m_flow_H2_signal = m_flow_H2_stoi * ((3.1 + time) / 3.0);
//Connectors
  connect(H2_Source.ports[1], H2_pump.fluid_a) annotation(
    Line(points = {{-92, 25}, {-64, 25}}, color = {0, 127, 255}));
  connect(H2_pump.fluid_b, Reactor.fluid_H2_in) annotation(
    Line(points = {{-42, 25}, {-26, 25}, {-26, 9}, {-10, 9}}, color = {0, 127, 255}));
  connect(Reactor.fluid_H2_out, H2_Sink.port_a) annotation(
    Line(points = {{86, 9}, {118, 9}}, color = {0, 127, 255}));
  connect(T_H2_signal.y, H2_Source.T_in) annotation(
    Line(points = {{-146, 32}, {-129, 32}}, color = {0, 0, 127}));
  connect(Heat_Loss.y, Reactor.Q_flow_loss) annotation(
    Line(points = {{18, 54}, {24, 54}, {24, 31}}, color = {0, 0, 127}));
  connect(Reactor.fluid_H2O_out, H2O_Sink.port_a) annotation(
    Line(points = {{86, -10}, {108, -10}, {108, -53}}, color = {0, 127, 255}));
  connect(Fe2O3_pump.fluid_b, Reactor.fluid_Fe2O3_in) annotation(
    Line(points = {{-44, -22}, {-26, -22}, {-26, -10}, {-10, -10}}, color = {0, 127, 255}));
  connect(Fe2O3_source.ports[1], Fe2O3_pump.fluid_a) annotation(
    Line(points = {{-90, -22}, {-64, -22}, {-64, -22}, {-64, -22}}, color = {0, 127, 255}));
  connect(T_Fe2O3_signal.y, Fe2O3_source.T_in) annotation(
    Line(points = {{-146, -16}, {-128, -16}, {-128, -16}, {-128, -16}}, color = {0, 0, 127}));
  connect(m_flow_H2.y, H2_pump.m_flow) annotation(
    Line(points = {{-88, 74}, {-54, 74}, {-54, 34}, {-52, 34}}, color = {0, 0, 127}));
  connect(m_flow_Fe2O3.y, Fe2O3_pump.m_flow) annotation(
    Line(points = {{-86, -65}, {-74, -65}, {-74, -6}, {-54, -6}, {-54, -14}}, color = {0, 0, 127}));
  connect(p_amb.y, Reactor.p_reactor) annotation(
    Line(points = {{62, 54}, {38, 54}, {38, 40}}, color = {0, 0, 127}));
  annotation(
    experiment(StopTime = 16.8, StartTime = 0.0, Tolerance = 1e-3, Interval = 0.2),
    Diagram(coordinateSystem(extent = {{-150, -100}, {150, 100}}, preserveAspectRatio = false, initialScale = 0.1)),
    Icon(coordinateSystem(extent = {{-150, -100}, {150, 100}}, preserveAspectRatio = false)));
end Reactor_Losses_Validation;