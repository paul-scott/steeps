within SolarTherm.Models.Reactors;

model Reactor_H2DRI_Overall_v2 "No enthalpy correction"
  extends Icons.Reactor;
  //extends Icons.ChemEqn;
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  //Reactor Pressure Imposed
  Modelica.Blocks.Interfaces.RealInput p_reactor "Reactor Pressure (Pa)" annotation(
    Placement(visible = true, transformation(origin = {48, 8.88178e-16}, extent = {{10, -10}, {-10, 10}}, rotation = 0), iconTransformation(origin = {0, 86}, extent = {{6, -6}, {-6, 6}}, rotation = 90)));
  //Fluid Inlet Ports
  Modelica.Fluid.Interfaces.FluidPort_a fluid_Fe2O3_in(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-94, 18}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_a fluid_H2_in(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(extent = {{-110, 40}, {-90, 60}}, rotation = 0), iconTransformation(origin = {-94, -20}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  //Fluid Outlet Ports
  Modelica.Fluid.Interfaces.FluidPort_b fluid_H2_out(redeclare package Medium = Modelica.Media.IdealGases.SingleGases.H2) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {94, 68}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_H2O_out(redeclare package Medium = Modelica.Media.Water.WaterIF97_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {94, 44}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_Fe2O3_out(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe2O3_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {94, 20}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_Fe3O4_out(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe3O4_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {94, -4}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_FeO_out(redeclare package Medium = SolarTherm.Media.SolidParticles.FeO_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {94, -30}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  Modelica.Fluid.Interfaces.FluidPort_b fluid_Fe_out(redeclare package Medium = SolarTherm.Media.SolidParticles.Fe_ph) annotation(
    Placement(visible = true, transformation(extent = {{-110, -80}, {-90, -60}}, rotation = 0), iconTransformation(origin = {94, -54}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  //Outlet Pressure
  SI.Pressure p_H2_out "Partial pressure of hydrogen gas leaving the reactor";
  SI.Pressure p_H2O_out "Partial pressure of water vapour leaving the reactor (Assume superheated gas)";
  
  //Mass Flow
  SI.MassFlowRate m_flow_Fe2O3_in "Number of moles of Fe2O3 flowing into the reactor (kg/s)";
  SI.MassFlowRate m_flow_H2_in "Number of moles of H2 flowing into the reactor (kg/s)";
  SI.MassFlowRate m_flow_H2_out "Number of moles of H2 flowing out of the reactor (kg/s)";
  SI.MassFlowRate m_flow_H2O_out "Number of moles of H2O flowing out of the reactor (kg/s)";
  SI.MassFlowRate m_flow_Fe2O3_out "Number of moles of Fe2O3 flowing out of the reactor (kg/s)";
  SI.MassFlowRate m_flow_Fe3O4_out "Number of moles of Fe3O4 flowing out of the reactor (kg/s)";
  SI.MassFlowRate m_flow_FeO_out "Number of moles of FeO flowing out of the reactor (kg/s)";
  SI.MassFlowRate m_flow_Fe_out "Number of moles of Fe flowing out of the reactor (kg/s)";
  //Molar Flow
  SI.MolarFlowRate n_flow_Fe2O3_in "Number of moles of Fe2O3 flowing into the reactor (mol/s)";
  SI.MolarFlowRate n_flow_H2_in "Number of moles of H2 flowing into the reactor (mol/s)";
  SI.MolarFlowRate n_flow_H2_out "Number of moles of H2 flowing out of the reactor (mol/s)";
  SI.MolarFlowRate n_flow_H2O_out "Number of moles of H2O flowing out of the reactor (mol/s)";
  SI.MolarFlowRate n_flow_Fe2O3_out "Number of moles of Fe2O3 flowing out of the reactor (mol/s)";
  SI.MolarFlowRate n_flow_Fe3O4_out "Number of moles of Fe3O4 flowing out of the reactor (mol/s)";
  SI.MolarFlowRate n_flow_FeO_out "Number of moles of FeO flowing out of the reactor (mol/s)";
  SI.MolarFlowRate n_flow_Fe_out "Number of moles of Fe flowing out of the reactor (mol/s)";
  //Pre-reaction reactant mix
  SI.Temperature T_reactants "Temperature of reactant stream (K)";
  SI.SpecificEnthalpy h_H2_mix "Specific enthalpy of H2 after being mixed with Fe2O3 (J/kg)";
  SI.SpecificEnthalpy h_Fe2O3_mix "Specific enthalpy of Fe2O3 after being mixed with H2 (J/kg)";
  Real r_input "Ratio of H2 moles to Fe2O3 moles in the reactant stream";
  //Reaction scaling factor
  SI.Frequency f_reaction "Maximum possible moles of reaction per second assuming Fe2O3 is in excess = n_flow_Fe2O3_in/0.5 (1/s)";
  //Thermodynamic Outputs
  //Step 1
  SI.MolarEnthalpy H_1;
  //SI.MolarEnthalpy G_1;
  //SI.MolarEntropy S_1;
  //Real K_1;
  Real xi_1;
  //Step 2
  SI.MolarEnthalpy H_2;
  //SI.MolarEnthalpy G_2;
  //SI.MolarEntropy S_2;
  //Real K_2;
  Real xi_2;
  //Step 1
  SI.MolarEnthalpy H_3;
  //SI.MolarEnthalpy G_3;
  //SI.MolarEntropy S_3;
  //Real K_3;
  Real xi_3;
  //Reaction Temperature
  SI.Temperature T_rxn(start = 500.0 + 273.15);
  //Set a starting value
  //Reaction Extent Lookup Table
  SolarTherm.Models.Chemistry.H2DRI.Reaction_Extent_Table_Akima Extent(file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Reactions/H2DRI_3Step_Refined.motab"));
  //SolarTherm.Models.Chemistry.H2DRI.Reaction_Extent_Regression Extent;
  //Overall reaction enthalpy (based on the extents of each reaction
  SI.HeatFlowRate Q_flow_released "Heat released per second of reaction";
  //Products of reaction
  SI.Temperature T_products "Temperature of the products after considering reaction enthalpy";
  //Post-reaction mix
  SI.SpecificEnthalpy h_H2_out "Specific enthalpy of H2 coming out of the reactor (J/kg)";
  SI.SpecificEnthalpy h_H2O_out "Specific enthalpy of H2O coming out of the reactor (J/kg)";
  SI.SpecificEnthalpy h_Fe2O3_out "Specific enthalpy of Fe2O3 coming out of the reactor (J/kg)";
  SI.SpecificEnthalpy h_Fe3O4_out "Specific enthalpy of Fe3O4 coming out of the reactor (J/kg)";
  SI.SpecificEnthalpy h_FeO_out "Specific enthalpy of FeO coming out of the reactor (J/kg)";
  SI.SpecificEnthalpy h_Fe_out "Specific enthalpy of Fe coming out of the reactor (J/kg)";
  //Check Continuity
  SI.MassFlowRate m_flow_in_total;
  SI.MassFlowRate m_flow_out_total;
equation
//Mass-Molar flow rate conversions
  n_flow_Fe2O3_in = m_flow_Fe2O3_in / Chemistry.ChemTable.Fe2O3.M;
  n_flow_H2_in = m_flow_H2_in / Chemistry.ChemTable.H2.M;
  m_flow_H2_out = n_flow_H2_out * Chemistry.ChemTable.H2.M;
  m_flow_H2O_out = n_flow_H2O_out * Chemistry.ChemTable.H2O.M;
  m_flow_Fe2O3_out = n_flow_Fe2O3_out * Chemistry.ChemTable.Fe2O3.M;
  m_flow_Fe3O4_out = n_flow_Fe3O4_out * Chemistry.ChemTable.Fe3O4.M;
  m_flow_FeO_out = n_flow_FeO_out * Chemistry.ChemTable.FeO.M;
  m_flow_Fe_out = n_flow_Fe_out * Chemistry.ChemTable.Fe.M;
//Pre-reaction reactant mix
  m_flow_Fe2O3_in * inStream(fluid_Fe2O3_in.h_outflow) + m_flow_H2_in * inStream(fluid_H2_in.h_outflow) = m_flow_Fe2O3_in * h_Fe2O3_mix + m_flow_H2_in * h_H2_mix;
  h_Fe2O3_mix = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_reactants);
  h_H2_mix = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_reactor, T_reactants);
//Temperature of the reactant mixture should have been calculated by now.
  f_reaction = n_flow_Fe2O3_in / 0.5;
  r_input = n_flow_H2_in / n_flow_Fe2O3_in;
//"Ratio of gas to ore"
//Reaction thermodynamics
  //Try standard enthalpy of reaction
  H_1 = SolarTherm.Models.Chemistry.H2DRI.Isothermal.Step_1_Enthalpy(T_reactants, p_reactor);
  H_2 = SolarTherm.Models.Chemistry.H2DRI.Isothermal.Step_2_Enthalpy(T_reactants, p_reactor);
  H_3 = SolarTherm.Models.Chemistry.H2DRI.Isothermal.Step_3_Enthalpy(T_reactants, p_reactor);
  //H_1 = SolarTherm.Models.Chemistry.H2DRI.Adiabatic.Step_1_Enthalpy(T_reactants,p_reactor,T_products,p_reactor);
  //H_2 = SolarTherm.Models.Chemistry.H2DRI.Adiabatic.Step_2_Enthalpy(T_reactants,p_reactor,T_products,p_reactor);
  //H_3 = SolarTherm.Models.Chemistry.H2DRI.Adiabatic.Step_3_Enthalpy(T_reactants,p_reactor,T_products,p_reactor);
//Reaction extent table
  Extent.r_input = r_input;
  Extent.T_rxn = T_reactants;
  xi_1 = Extent.xi_1;
  xi_2 = Extent.xi_2;
  xi_3 = Extent.xi_3;
//Rate of heat released by reaction
  Q_flow_released = -1.0 * (xi_1 * H_1 + xi_2 * H_2 + xi_3 * H_3) * f_reaction;
//Check this
//Number of moles output for each product in the stream
  n_flow_H2_out = n_flow_H2_in - f_reaction * (1.0 / 6.0 * xi_1 + 1.0 / 3.0 * xi_2 + xi_3);
  n_flow_H2O_out = f_reaction * (1.0 / 6.0 * xi_1 + 1.0 / 3.0 * xi_2 + xi_3);
  n_flow_Fe2O3_out = n_flow_Fe2O3_in - f_reaction * (0.5 * xi_1);
  n_flow_Fe3O4_out = f_reaction * (1.0 / 3.0 * xi_1 - 1.0 / 3.0 * xi_2);
  n_flow_FeO_out = f_reaction * (xi_2 - xi_3);
  n_flow_Fe_out = f_reaction * xi_3;
//Figure out the temperature of reactants before enthalpy of reaction, then add reaction enthalpy
  //m_flow_Fe2O3_in*h_Fe2O3_mix + m_flow_H2_in*h_H2_mix + Q_flow_released = m_flow_H2_out * h_H2_out + m_flow_H2O_out * h_H2O_out + m_flow_Fe2O3_out * h_Fe2O3_out + m_flow_Fe3O4_out * h_Fe3O4_out + m_flow_FeO_out * h_FeO_out + m_flow_Fe_out * h_Fe_out;
  
  m_flow_H2_out * Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_reactor, T_reactants) + m_flow_H2O_out * Modelica.Media.Water.IF97_Utilities.h_pT(p_reactor, T_reactants) + m_flow_Fe2O3_out * SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_reactants) + m_flow_Fe3O4_out * SolarTherm.Media.SolidParticles.Fe3O4_utilities.h_T(T_reactants) + m_flow_FeO_out * SolarTherm.Media.SolidParticles.FeO_utilities.h_T(T_reactants) + m_flow_Fe_out * SolarTherm.Media.SolidParticles.Fe_utilities.h_T(T_reactants) + Q_flow_released = m_flow_H2_out * h_H2_out + m_flow_H2O_out * h_H2O_out + m_flow_Fe2O3_out * h_Fe2O3_out + m_flow_Fe3O4_out * h_Fe3O4_out + m_flow_FeO_out * h_FeO_out + m_flow_Fe_out * h_Fe_out;
  
  
  h_H2_out = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_reactor, T_products);
  h_H2O_out = Modelica.Media.Water.IF97_Utilities.h_pT(p_reactor, T_products);
  h_Fe2O3_out = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_products);
  h_Fe3O4_out = SolarTherm.Media.SolidParticles.Fe3O4_utilities.h_T(T_products);
  h_FeO_out = SolarTherm.Media.SolidParticles.FeO_utilities.h_T(T_products);
  h_Fe_out = SolarTherm.Media.SolidParticles.Fe_utilities.h_T(T_products);
//Partial pressure of gases
  p_H2_out = p_reactor*(n_flow_H2_out)/(n_flow_H2_out+n_flow_H2O_out);
  p_H2O_out = p_reactor*(n_flow_H2O_out)/(n_flow_H2_out+n_flow_H2O_out);
  
  
//Fluid Ports
  m_flow_Fe2O3_in = fluid_Fe2O3_in.m_flow;
  m_flow_H2_in = fluid_H2_in.m_flow;
  fluid_H2_out.m_flow = -1.0 * m_flow_H2_out;
  fluid_H2O_out.m_flow = -1.0 * m_flow_H2O_out;
  fluid_Fe2O3_out.m_flow = -1.0 * m_flow_Fe2O3_out;
  fluid_Fe3O4_out.m_flow = -1.0 * m_flow_Fe3O4_out;
  fluid_FeO_out.m_flow = -1.0 * m_flow_FeO_out;
  fluid_Fe_out.m_flow = -1.0 * m_flow_Fe_out;
  fluid_Fe2O3_in.h_outflow = inStream(fluid_Fe2O3_in.h_outflow);
  fluid_H2_in.h_outflow = inStream(fluid_H2_in.h_outflow);
  fluid_H2_out.h_outflow = h_H2_out;
  fluid_H2O_out.h_outflow = h_H2O_out;
  fluid_Fe2O3_out.h_outflow = h_Fe2O3_out;
  fluid_Fe3O4_out.h_outflow = h_Fe3O4_out;
  fluid_FeO_out.h_outflow = h_FeO_out;
  fluid_Fe_out.h_outflow = h_Fe_out;
  fluid_Fe2O3_in.p = p_reactor;
  fluid_H2_in.p = p_reactor;
  fluid_H2_out.p = p_H2_out;
  fluid_H2O_out.p = p_H2O_out;
  fluid_Fe2O3_out.p = p_reactor;
  fluid_Fe3O4_out.p = p_reactor;
  fluid_FeO_out.p = p_reactor;
  fluid_Fe_out.p = p_reactor;
//More equations
//T_rxn = T_reactants;
  T_rxn = (T_products + T_reactants) / 2.0;
//T_rxn = T_products;
//Check continuity
  m_flow_in_total = m_flow_H2_in + m_flow_Fe2O3_in;
  m_flow_out_total = m_flow_H2_out + m_flow_H2O_out + m_flow_Fe_out + m_flow_FeO_out + m_flow_Fe3O4_out + m_flow_Fe2O3_out;
  annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),
    Icon(graphics = {Line(origin = {-67, -20}, points = {{-21, 0}, {21, 0}}), Text(origin = {73, 81}, extent = {{1, -1}, {13, -13}}, textString = "H2"), Line(origin = {-67, 18}, points = {{-21, 0}, {21, 0}}), Text(origin = {75, 49}, extent = {{-5, 7}, {13, -7}}, textString = "H2O"), Text(origin = {-75, -7}, extent = {{-1, -3}, {13, -13}}, textString = "H2"), Text(origin = {-71, 23}, extent = {{-11, 7}, {13, -7}}, textString = "Fe2O3"), Line(origin = {65.2285, 43.861}, points = {{27.139, 24.139}, {3.13904, 24.139}, {-16.861, -27.861}}), Line(origin = {67.8201, 29.8201}, points = {{20.1799, 14.1799}, {0.1799, 14.1799}, {-19.8201, -17.8201}}), Line(origin = {66.8015, 13.8015}, points = {{21.1985, 8.19845}, {1.19845, 8.19845}, {-20.8015, -7.80155}}), Line(origin = {67.814, 0.185994}, points = {{20.186, -6.18599}, {0.185994, -6.18599}, {-19.814, 5.81401}}), Line(origin = {66.8251, -16.8251}, points = {{23.1749, -13.1749}, {1.17492, -13.1749}, {-20.8251, 20.8251}}), Line(origin = {67.2259, -29.8584}, points = {{21.1416, -24.1416}, {1.14157, -24.1416}, {-20.8584, 29.8584}}), Text(origin = {77, 27}, extent = {{-11, 7}, {13, -7}}, textString = "Fe2O3"), Text(origin = {75, -1}, extent = {{-11, 7}, {13, -7}}, textString = "Fe3O4"), Text(origin = {77, -23}, extent = {{-9, 3}, {13, -7}}, textString = "FeO"), Text(origin = {77, -47}, extent = {{-9, 3}, {13, -7}}, textString = "Fe")}, coordinateSystem(initialScale = 0.1)));
end Reactor_H2DRI_Overall_v2;