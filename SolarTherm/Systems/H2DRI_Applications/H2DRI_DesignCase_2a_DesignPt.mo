within SolarTherm.Systems.H2DRI_Applications;

model H2DRI_DesignCase_2a_DesignPt
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  
  //Medium Packages
  replaceable package Medium_H2 = Modelica.Media.IdealGases.SingleGases.H2;
  replaceable package Medium_H2O = Modelica.Media.Water.WaterIF97_ph;
  replaceable package Medium_Fe2O3 = SolarTherm.Media.SolidParticles.Fe2O3_ph;
  
  //Material Packages
  replaceable package Material_Fe2O3 = SolarTherm.Materials.Fe2O3;
  
  //Chemical Constants
  parameter SI.MolarMass M_Fe2O3 = SolarTherm.Models.Chemistry.ChemTable.Fe2O3.M;
  parameter SI.MolarMass M_H2O = SolarTherm.Models.Chemistry.ChemTable.H2O.M;
  parameter SI.MolarMass M_H2 = SolarTherm.Models.Chemistry.ChemTable.H2.M;
//Plant Constants e.g. reference/feedstock conditions
  parameter SI.Pressure p_des = 100000.0 "Design pressure of the plant (Pa)";
  parameter SI.Temperature T_H2_feedstock_des = 25.0+273.15 "Design feedstock H2 temperature (K)";
  parameter SI.Temperature T_Fe2O3_feedstock_des = 25.0+273.15 "Design feedstock Fe2O3 temperature (K)";
  parameter SI.SpecificEnthalpy h_H2_feedstock_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_H2_feedstock_des) "Specific Enthalpy of feedstock H2 (J/kg)";
  
  //Iterative solve for this:
  parameter SI.Temperature T_Fe2O3_hot_des = 842.04 + 273.15 "Actual temperature of H2 entering the GGXH cold stream after being mixed with topup stream (J/kg)";
  
  
  //Design Inputs
  parameter Real m_flow_Fe_Mtperyr = 1.0 "Plant Size in MegaTons per year (Mt/yr)";
  parameter SI.MassFlowRate m_flow_Fe_des = m_flow_Fe_Mtperyr*(1.0e9)/(31536000.0) "Design maximum mass flow rate output of Iron (kg/s)"; 
  //We assume a year is 365.0 days
  //47.6
  //parameter SI.Temperature T_reactants_des = 757.0 + 273.15 "Design reactant temperature (K)";
  parameter SI.Temperature T_products_des = 600.0 + 273.15 "Design product temperature (K)";
  parameter SI.Temperature T_reactants_des = 1.0862878*T_products_des+51.45030 "Design estimated reactant temperature (K)";
  parameter SI.Temperature T_inlet_des = 51.02315 + (1086.761*T_products_des/1000.0) - (0.00838*Q_loss_per_mol_des/1000.0) + 2.16376*(T_products_des*Q_loss_per_mol_des*1.0e-6);
  
  //Variable
  SI.Temperature T_mix_actual;
  
  parameter Real r_min_des = 6.062457e1 + ((-7.646067e-2)*T_reactants_des) + ((2.853147e-5)*T_reactants_des*T_reactants_des);
  
  parameter SI.SpecificEnthalpy h_Fe2O3_hot_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_Fe2O3_hot_des);
  
  parameter SI.MolarEnthalpy Q_loss_per_mol_des = 34.57e3 "Design Heat loss per mole of Fe2O3 (J/mol) Note factor of 2.0 due to stoichiometry";
  
  parameter Real eff_reactor = 0.8 "Thermal efficiency of the reactor at design point (-).";
  
  parameter SI.MolarEnthalpy H_rxn_per_mol_des = 2.0*SolarTherm.Models.Chemistry.H2DRI.Isothermal.Overall_Rxn_Enthalpy(T_reactants_des,p_des);
 
  parameter SI.MolarEnthalpy Q_loss_per_mol_des2 = Q_loss_des/n_flow_Fe2O3_des "Design Heat loss per mole of Fe2O3 (J/mol) Note factor of 2.0 due to stoichiometry";
  parameter SI.HeatFlowRate Q_loss_des = n_flow_Fe2O3_des*Q_loss_per_mol_des;
  
  
  //parameter SI.SpecificEnthalpy h_Fe2O3_hot_des = ( m_flow_H2_des*(h_H2_reactant_des-h_H2_pre2_des) + m_flow_Fe2O3_des*h_Fe2O3_reactant_des )/m_flow_Fe2O3_des;
  
  parameter Real eff_allHX_des = 0.80 "Effectiveness of all HX (-)";
  parameter Real eff_GGHX_des = eff_allHX_des "Effectiveness of the GGHX (-)";
  parameter Real eff_PGHX1_des = eff_allHX_des "Effectiveness of PGHX1 (-)";
  parameter Real eff_PGHX2_des = eff_allHX_des "Effectiveness of PGHX2 (-)";
  
  //Design Reactor Inlet Mixed Enthalpies
  parameter SI.SpecificEnthalpy h_H2_inlet_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_inlet_des);
  parameter SI.SpecificEnthalpy h_Fe2O3_inlet_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_inlet_des);
  
  parameter SI.Temperature T_condenser_out_des = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p_H2O_offgas_des) "Design saturation temperature of the condenser (K)";
  parameter SI.SpecificEnthalpy h_H2_condenser_out_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_condenser_out_des) "Design specific enthalpy of H2 leaving the condenser (J/kg)";
  
  

  parameter SI.SpecificEnthalpy h_H2_mix_des = (m_flow_H2_stoi*h_H2_feedstock_des + m_flow_H2_excess*h_H2_condenser_out_des)/(m_flow_H2_stoi+m_flow_H2_excess) "Design specific enthalpy of H2 entering the GGXH cold stream after being mixed with topup stream (J/kg)";
  parameter SI.Temperature T_H2_mix_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_mix_des) "Design temperature of H2 entering the GGXH cold stream after being mixed with topup stream (J/kg)";
  
  
  //parameter SI.Temperature T_Fe2O3_hot_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_hot_des) "Design temperature of Fe2O3 in the hot TES (K)";
  
  
  //parameter SI.SpecificEnthalpy h_Fe2O3_hot_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_Fe2O3_hot_des);

  //GGHX Outlet design temperature
  parameter SI.SpecificEnthalpy h_H2_GGHX_max=Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_products_des) "Design specific enthalpy of H2 leaving the GGHX cold stream assuming it had been heated to the maximum possible temperature i.e. reactor product temperature (J/kg)";
  parameter SI.SpecificEnthalpy h_H2_PGHX2_max=Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_Fe2O3_hot_des);
  
  
  
  //GGHX Calculations where f_split_GGHX worth of offgas gets diverted to
  parameter SI.SpecificEnthalpy h_H2_pre1_des=h_H2_mix_des + eff_GGHX_des*(h_H2_GGHX_max-h_H2_mix_des) "Design of H2 leaving the GGHX cold stream after considering HX effectiveness and assuming C_C is C_min (J/kg)";
  
  parameter SI.SpecificHeatCapacity cp_Hot_H2O_GGHX = 0.5*(Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, T_products_des) + Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, max(100.0, T_H2_mix_des))) "cp of the hot H2O entering the GGHX (J/kgK)";
  
  parameter SI.SpecificHeatCapacity cp_Hot_H2_GGHX = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_products_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_H2_mix_des)) / (T_products_des - T_H2_mix_des) "cp of the hot H2 entering the GGHX (J/kgK)";
  
  parameter SI.SpecificHeatCapacity cp_Cold_H2_GGHX = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_products_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_H2_mix_des)) / (T_products_des - T_H2_mix_des) "cp of the cold H2 entering the GGHX (J/kgK)";
  
  //parameter SI.ThermalConductance C_GGHX_H = m_flow_H2_excess*cp_GGHX_Hot_H2 + m_flow_H2O_des*cp_GGHX_Hot_H2O "Heat capacity rate of the hot reactor offgas stream exiting the reactor and entering the GGHX (W/K)";
  
  parameter SI.ThermalConductance C_C_GGHX = m_flow_H2_des*cp_Cold_H2_GGHX "Heat capacity rate of the cold H2 stream exiting the H2 mixer and entering the GGHX (W/K)"; //this is the one that C_GGHX_H must match.
  
  parameter Real f_split_GGHX = C_C_GGHX/(m_flow_H2_excess*cp_Hot_H2_GGHX + m_flow_H2O_des*cp_Hot_H2O_GGHX) "Fraction of off-gas flow that is diverted to the GGHX in order to match C_C and C_H of that HX component";  //results in C_r_GGHX
  
  
  //PGHX1 Calculations, where (1-f_split_GGHX) worth of gas gets diverted to.
  //The mass flow rates of hot H2 and hot H2O are (1-f_split_GGHX)*m_flow_H2_excess and (1-f_split_GGHX)*m_flow_H2O_des respectively
  //The hot and cold temperatures in this case are T_products_des and T_Fe2O3_feedstock_des
  //T_Fe2O3_feedstock_des
  parameter SI.SpecificHeatCapacity cp_Hot_H2O_PGHX1 = 0.5*(Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, T_products_des) + Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, max(100.0, T_Fe2O3_feedstock_des))) "cp of the hot H2O entering the PGHX1 (J/kgK)";
  parameter SI.SpecificHeatCapacity cp_Hot_H2_PGHX1 = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_products_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_Fe2O3_feedstock_des)) / (T_products_des - T_H2_mix_des) "cp of the hot H2 entering the PGHX1 (J/kgK)";
  
  parameter SI.SpecificHeatCapacity cp_Cold_Fe2O3_PGHX1 = (Material_Fe2O3.h_Tf(T_products_des,0.0) - Material_Fe2O3.h_Tf(T_Fe2O3_feedstock_des,0.0))/(T_products_des - T_Fe2O3_feedstock_des) "cp of the cold Fe2O3 entering the PGHX1 (J/kgK)";
  
  parameter SI.ThermalConductance C_H_PGHX1 = (1-f_split_GGHX)*m_flow_H2_excess*cp_Hot_H2_PGHX1 + (1-f_split_GGHX)*m_flow_H2O_des*cp_Hot_H2O_PGHX1 "Heat capacity rate of the hot off-gas entering PGHX1 (W/K)";
  
  //The mass flow rate of Fe2O3 feedstock entering the PGHX1 is adjusted to match C_C_PGHX1 to C_H_PGHX1
  parameter SI.MassFlowRate m_flow_Fe2O3_feedstock = C_H_PGHX1/cp_Cold_Fe2O3_PGHX1 "Design mass flow rate of cold Fe2O3 feedstock entering PGHX1 (kg/s)";
  
  //What is the enthalpy of Fe2O3 exiting the PGHX1 and entering Med Tank?
  parameter SI.SpecificEnthalpy h_Fe2O3_feedstock_des = Material_Fe2O3.h_Tf(T_Fe2O3_feedstock_des,0.0) "Specific enthalpy of Fe2O3 feedstock (J/kg)";
  parameter SI.SpecificEnthalpy h_Fe2O3_PGHX1_max = Material_Fe2O3.h_Tf(T_products_des,0.0) "Maximum possible specific enthalpy of Fe2O3 leaving PGHX1 assuming 100% effectiveness (J/kg)";
  parameter SI.SpecificEnthalpy h_Fe2O3_med1_des = h_Fe2O3_feedstock_des + eff_PGHX1_des*(h_Fe2O3_PGHX1_max - h_Fe2O3_feedstock_des) "Specific enthalpy of Fe2O3 leaving PGHX1 and entering the med tank (J/kg)";
  parameter SI.Temperature T_Fe2O3_med1_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_med1_des) "Temperature of Fe2O3 leaving PGHX1 and entering the med tank (K)";
  
  
  
  parameter SI.SpecificEnthalpy h_H2_pre2_des=h_H2_pre1_des + eff_PGHX2_des*(h_H2_PGHX2_max-h_H2_pre1_des) "Design of H2 leaving the PGHX2 cold stream after considering HX effectiveness and assuming C_H2 is C_min (J/kg)";

  //Storage Temperatures
  parameter SI.Temperature T_H2_pre1_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_pre1_des) "Design temperature of the H2 leaving the GGHX";
  
  parameter SI.Temperature T_H2_pre2_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_pre2_des) "Design temperature of the H2 leaving PGHX2";
  
  //Flow Rates
  parameter SI.MassFlowRate m_flow_Fe2O3_des = 1.4297448*m_flow_Fe_des "Design required Fe2O3 flow rate (kg/s)";
  
  parameter SI.MolarFlowRate n_flow_Fe2O3_des = m_flow_Fe2O3_des/M_Fe2O3 "Design required molar Fe2O3 flow rate (mol/s)";
  parameter SI.MolarFlowRate n_flow_H2_stoi = 3.0*n_flow_Fe2O3_des "Design stoichiometric molar flow rate of H2 consumed assuming 100% yield of iron, also required top-up moles (mol/s)";
  parameter SI.MolarFlowRate n_flow_H2O_des = n_flow_H2_stoi "Design molar flow rate of H2O produced in the reactor assuming 100% yield of iron (mol/s)";
  parameter SI.MassFlowRate m_flow_H2O_des = n_flow_H2O_des*M_H2O "Design mass flow rate of H2O produced in the reactor assuming 100% yield of iron (kg/s)";
  parameter SI.MassFlowRate m_flow_H2_stoi = n_flow_H2_stoi*M_H2 "Design stoichiometric mass flow rate of H2 consumed assuming 100% yield of iron, also required top-up mass flow (kg/s)";
  
  parameter SI.MolarFlowRate n_flow_H2_des = n_flow_Fe2O3_des*r_min_des "Design molar flow rate of H2 required to ensure 100% yield of iron (mol/s)";
  parameter SI.MassFlowRate m_flow_H2_des = n_flow_H2_des*M_H2 "Design mass flow rate of H2 required to ensure 100% yield of iron (kg/s)";
  parameter SI.MassFlowRate m_flow_H2_excess = m_flow_H2_des - m_flow_H2_stoi "Excess mass flow rate of H2 exiting the reactor (kg/s)";
  
  //Reactor Outlet Pressures
  parameter SI.Pressure p_H2_offgas_des = p_des*(n_flow_H2_des - n_flow_H2_stoi)/(n_flow_H2_des - n_flow_H2_stoi + n_flow_H2O_des) "Design partial pressure of H2 in the offgas (Pa)";
  parameter SI.Pressure p_H2O_offgas_des = p_des*(n_flow_H2O_des)/(n_flow_H2_des - n_flow_H2_stoi + n_flow_H2O_des) "Design partial pressure of H2O in the offgas (Pa)";
  
  //Reference enthalpies
  parameter SI.SpecificEnthalpy h_H2_ref = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(1.0e5,298.15);
  parameter SI.SpecificEnthalpy h_H2O_ref = Modelica.Media.Water.IF97_Utilities.h_pT(1.0e5,298.15);
  
  
  //PGHX2 Calculations
  parameter SI.MassFlowRate m_flow_Fe2O3_PGHX2 = m_flow_H2_des*cp_H2_PGHX2/cp_Fe2O3_PGHX2;
  parameter SI.SpecificHeatCapacity cp_Fe2O3_PGHX2 = (SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_Fe2O3_hot_des)-SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_H2_pre1_des))/(T_Fe2O3_hot_des-T_H2_pre1_des);
  parameter SI.SpecificHeatCapacity cp_H2_PGHX2 = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_Fe2O3_hot_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_H2_pre1_des)) / (T_Fe2O3_hot_des - T_H2_pre1_des);
  
  parameter SI.SpecificEnthalpy h_Fe2O3_med2_des = h_Fe2O3_hot_des - eff_PGHX2_des*(h_Fe2O3_hot_des - SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_H2_pre1_des)) "Specific enthalpy of Fe2O3 leaving PGHX2 and entering the Medium temperature tank (J/kg)";
  parameter SI.Temperature T_Fe2O3_med2_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_med2_des) "Temperature of Fe2O3 leaving PGHX2 and entering the Medium temperature tank (K)";
  
  //What is the expected temperature of materials entering the med tank assuming they get mixed well?
  parameter SI.SpecificEnthalpy h_Fe2O3_med_des = (m_flow_Fe2O3_feedstock*h_Fe2O3_med1_des + m_flow_Fe2O3_PGHX2*h_Fe2O3_med2_des)/(m_flow_Fe2O3_feedstock + m_flow_Fe2O3_PGHX2) "Expected specific enthalpy of the med tank (J/kg)";
  parameter SI.Temperature T_Fe2O3_med_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_med_des) "Expected specific enthalpy of the med tank (K)";
  
  
  //Linear Guess
  parameter SI.HeatFlowRate H_flow_inlet_des = m_flow_Fe2O3_des*h_Fe2O3_inlet_des + m_flow_H2_des*h_H2_inlet_des;
  
  parameter SI.Temperature T_Fe2O3_hot_guess = (H_flow_inlet_des+m_flow_H2_des*(h_H2_pre1_des*eff_PGHX2_des - h_H2_pre1_des) +261521.6*m_flow_Fe2O3_des+762976.9*m_flow_H2_des*eff_PGHX2_des)/(895.7669*m_flow_Fe2O3_des+15214.199*m_flow_H2_des*eff_PGHX2_des);
  
  //Reactor Physical Design
  //Input parameters and assumptions
  //parameter SI.Diameter d_p = 250.0e-6 "Particle diameter of iron ore (m)";
  //parameter Real f_fluidisation = 3.0 "The gas velocity is f times the minimum fluidisation velocity";
 //parameter Real epsilon_pb = 0.32 "Packed Bed porosity, assuming static conditions  (-)";
  //parameter Real phi_s = 0.86 "Sphericity of the Fe2O3 ore particles (-)";
  //parameter SI.Time t_residence_min = 1.0*3600.0 "Min Residence time of the reactor (s)";
  //parameter SI.Length t_shell = 0.0127 "Minimum shell thickness of the reactor (m)";
  //parameter SI.CoefficientOfHeatTransfer U_loss_reactor = 29.47;

  
  
  
  
  //Properties
  //parameter SI.Density rho_s_reactor = SolarTherm.Media.SolidParticles.Fe2O3_utilities.rho_T(T_products_des) "Density of Fe2O3 at the reactor product temperature";
  //parameter SI.Density rho_g_reactor = p_des*SolarTherm.Models.Chemistry.ChemTable.H2.M/(CN.R*T_products_des) "Density of H2 at the reactor product temperature";
  //parameter SI.DynamicViscosity mu_g_reactor = Modelica.Media.IdealGases.SingleGases.H2.dynamicViscosity(Modelica.Media.IdealGases.SingleGases.H2.setState_pT(p_des,T_products_des)) "Dynamic visocsity of H2 gas at the reactor product temperature";
  //parameter Real V_flow_Fe2O3_reactor = m_flow_Fe2O3_des/(rho_s_reactor*(1.0-epsilon_mf)*(1.0-epsilon_reactor)) "Mass of Fe2O3 in the reactor (kg)";
  //parameter Real epsilon_material = 0.20 "Material porosity of the Fe2O3 ore (-)";
  //parameter Real epsilon_reactor = 0.25 "Void space above the reacting bed and half of the tube that is not occupied by the ore (-)";
  //parameter SI.Length D_reactor = 6.40 "Diameter of reactor (m)"; //Maximum diameter for vertical vessel is 21 ft.
  //parameter SI.Length H_reactor = 12.192 "Height of reactor (m)";
  //parameter SI.
  //parameter SI.Density rho_shell = 7861.09;// SolarTherm.Materials.Inconel625.rho_Tf(298.15,0.0) "Density of Inconel625 at 25degC";
  
  //Calculated Parameters
  //parameter SI.Velocity u_mf = (d_p*d_p*(rho_s_reactor-rho_g_reactor)*9.81)/(1650.0*mu_g_reactor) "Minimum fluidisation velocity (superficial) (m/s)";
  //parameter SI.Velocity u_gas = f_fluidisation*u_mf "Superficial velocity of the H2 gas (m/s)";
  //parameter Real epsilon_mf = (1.0/(14.0*phi_s))^(1.0/3.0) "Bed porosity at minimum fluidisation (-)";
  //parameter SI.Area A_plate = m_flow_H2_des/(rho_g_reactor*u_gas) "Rectangular area of the perforated plate where H2 is flowed through (m2)";
  //parameter SI.Mass m_shell = CN.pi*(D_reactor+t_shell)*(H_reactor+0.8*D_reactor)*t_shell*rho_shell "Mass of empty reactor shell must be less than 417305 kg";
  //parameter SI.Length L_reactor_H2 = A_plate/D_reactor;
  //parameter Real N_reactors = A_plate/(CN.pi*D_reactor*D_reactor*0.25);
  
  
  //parameter SI.Volume V_reactor_total = N_reactors*(CN.pi*D_reactor*D_reactor*0.25*H_reactor) "Total volume of the reactors (m3)";
  //parameter SI.Length L_reactor = 4.0*V_reactor/(CN.pi*D_reactor*D_reactor) "Minimum length of reactor to satisfy the residence time of Fe2O3 ore (m)";
  //parameter SI.Area A_surface_reactor_total = N_reactors*(CN.pi*D_reactor*H_reactor + 0.5*CN.pi*D_reactor*D_reactor);
  //parameter SI.HeatFlowRate Q_loss_des = U_loss_reactor*A_surface_reactor_total*(0.5*(T_inlet_des+T_products_des) - 298.15);
  //parameter Real C_vessel = (816.0/500.0)*exp(7.0132 + 0.18255*(log(2.20462*m_shell)) + 0.02297*(log(2.20462*m_shell))^2.0);
  //parameter Real C_pl = (816.0/500.0)*361.8*((3.28084*D_reactor)^0.7396)*((3.28084*H_reactor)^0.70684);
  //parameter Real C_reactor = N_reactors*(F_M*C_vessel+C_pl);
  //parameter Real F_M = 3.9 "Inconel-600 Table 22 Seider";
  
  //parameter Real BFM = 4.16;
  //parameter Real IDC = 2.05;
  //parameter Real C_reactor_final = BFM*IDC*C_reactor;
  
  //parameter SI.Time t_residence = V_reactor_total/V_flow_Fe2O3_reactor;
equation

  m_flow_Fe2O3_des*h_Fe2O3_hot_des + m_flow_H2_des*h_H2_pre2_des = m_flow_Fe2O3_des*SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_mix_actual) + m_flow_H2_des*Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(1.0e5,T_mix_actual);
  

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-150, -100}, {150, 100}})),
    Icon(coordinateSystem(extent = {{-150, -100}, {150, 100}}, preserveAspectRatio = false)));
end H2DRI_DesignCase_2a_DesignPt;