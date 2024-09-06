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
  parameter SI.MolarMass M_Fe3O4 = SolarTherm.Models.Chemistry.ChemTable.Fe3O4.M;
  parameter SI.MolarMass M_FeO = SolarTherm.Models.Chemistry.ChemTable.FeO.M;
  parameter SI.MolarMass M_Fe = SolarTherm.Models.Chemistry.ChemTable.Fe.M;
  parameter SI.MolarMass M_H2O = SolarTherm.Models.Chemistry.ChemTable.H2O.M;
  parameter SI.MolarMass M_H2 = SolarTherm.Models.Chemistry.ChemTable.H2.M;
//Plant Constants e.g. reference/feedstock conditions
  parameter SI.Pressure p_des = 100000.0 "Design pressure of the plant (Pa)";
  parameter SI.Temperature T_amb_des = 25.0 + 273.15 "Design ambient temperature (K)";
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
  
  //parameter Real eff_reactor = 0.8 "Thermal efficiency of the reactor at design point (-).";
  
  parameter SI.MolarEnthalpy H_rxn_per_mol_des = 2.0*SolarTherm.Models.Chemistry.H2DRI.Isothermal.Overall_Rxn_Enthalpy(T_reactants_des,p_des);
 
  //parameter SI.MolarEnthalpy Q_loss_per_mol_des2 = Q_loss_des/n_flow_Fe2O3_des "Design Heat loss per mole of Fe2O3 (J/mol) Note factor of 2.0 due to stoichiometry";
  parameter SI.HeatFlowRate Q_loss_des = n_flow_Fe2O3_des*Q_loss_per_mol_des;
  
  
  //parameter SI.SpecificEnthalpy h_Fe2O3_hot_des = ( m_flow_H2_des*(h_H2_reactant_des-h_H2_pre2_des) + m_flow_Fe2O3_des*h_Fe2O3_reactant_des )/m_flow_Fe2O3_des;
  
  parameter Real eff_allHX_des = 0.80 "Effectiveness of all HX (-)";
  parameter Real eff_GGHX_des = eff_allHX_des "Effectiveness of the GGHX (-)";
  parameter Real eff_PGHX1_des = eff_allHX_des "Effectiveness of PGHX1 (-)";
  parameter Real eff_PGHX2_des = eff_allHX_des "Effectiveness of PGHX2 (-)";
  
  //Design Reactor Inlet Mixed Enthalpies
  parameter SI.SpecificEnthalpy h_H2_inlet_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_inlet_des);
  parameter SI.SpecificEnthalpy h_Fe2O3_inlet_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_inlet_des);
  
  parameter SI.Temperature T_condenser_out_des = Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p_H2O_offgas_des) - 1.0 "Design saturation temperature of the condenser (K)";
  parameter SI.SpecificEnthalpy h_H2_condenser_out_des = Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_condenser_out_des) "Design specific enthalpy of H2 leaving the condenser (J/kg)";
  
  

  parameter SI.SpecificEnthalpy h_H2_mix_des = (m_flow_H2_stoi*h_H2_feedstock_des + m_flow_H2_excess*h_H2_condenser_out_des)/(m_flow_H2_stoi+m_flow_H2_excess) "Design specific enthalpy of H2 entering the GGXH cold stream after being mixed with topup stream (J/kg)";
  parameter SI.Temperature T_H2_mix_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_mix_des) "Design temperature of H2 entering the GGXH cold stream after being mixed with topup stream (J/kg)";
  
  
  //parameter SI.Temperature T_Fe2O3_hot_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_hot_des) "Design temperature of Fe2O3 in the hot TES (K)";
  
  
  //parameter SI.SpecificEnthalpy h_Fe2O3_hot_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_Fe2O3_hot_des);

  //GGHX Outlet design temperature
  parameter SI.SpecificEnthalpy h_H2_GGHX_max=Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_products_des) "Design specific enthalpy of H2 leaving the GGHX cold stream assuming it had been heated to the maximum possible temperature i.e. reactor product temperature (J/kg)";
  
  
  
  
  //GGHX Calculations where f_split_GGHX worth of offgas gets diverted to
  parameter SI.SpecificEnthalpy h_H2_pre1_des=h_H2_mix_des + eff_GGHX_des*(h_H2_GGHX_max-h_H2_mix_des) "Design of H2 leaving the GGHX cold stream after considering HX effectiveness and assuming C_C is C_min (J/kg)";
  
  parameter SI.SpecificHeatCapacity cp_Hot_H2O_GGHX = 0.5*(Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, T_products_des) + Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p_H2O_offgas_des)+1.0)) "cp of the hot H2O entering the GGHX (J/kgK)";
  
  parameter SI.SpecificHeatCapacity cp_Hot_H2_GGHX = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_products_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_H2_mix_des)) / (T_products_des - T_H2_mix_des) "cp of the hot H2 entering the GGHX (J/kgK)";
  
  parameter SI.SpecificHeatCapacity cp_Cold_H2_GGHX = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_products_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_H2_mix_des)) / (T_products_des - T_H2_mix_des) "cp of the cold H2 entering the GGHX (J/kgK)";
  
  //parameter SI.ThermalConductance C_GGHX_H = m_flow_H2_excess*cp_GGHX_Hot_H2 + m_flow_H2O_des*cp_GGHX_Hot_H2O "Heat capacity rate of the hot reactor offgas stream exiting the reactor and entering the GGHX (W/K)";
  
  parameter SI.ThermalConductance C_C_GGHX = m_flow_H2_des*cp_Cold_H2_GGHX "Heat capacity rate of the cold H2 stream exiting the H2 mixer and entering the GGHX (W/K)"; //this is the one that C_GGHX_H must match.
  
  parameter Real f_split_GGHX = C_C_GGHX/(m_flow_H2_excess*cp_Hot_H2_GGHX + m_flow_H2O_des*cp_Hot_H2O_GGHX) "Fraction of off-gas flow that is diverted to the GGHX in order to match C_C and C_H of that HX component";  //results in C_r_GGHX
  
  parameter SI.HeatFlowRate Q_flow_GGHX_des = m_flow_H2_des*eff_GGHX_des*(h_H2_GGHX_max-h_H2_mix_des) "Design heat transfer rate inside GGHX (J/s)";
  
  parameter SI.ThermalConductance C_min_GGHX_des = C_C_GGHX;
  
  //Calculate temperature of cooled off-gas exiting GGHX
  SI.Temperature T_condenser2_in_des(start = 500.0) "Temperature of cooled offgas exiting GGHX and entering condenser2 (K)";
  
  
  //PGHX1 Calculations, where (1-f_split_GGHX) worth of gas gets diverted to.
  //The mass flow rates of hot H2 and hot H2O are (1-f_split_GGHX)*m_flow_H2_excess and (1-f_split_GGHX)*m_flow_H2O_des respectively
  //The hot and cold temperatures in this case are T_products_des and T_Fe2O3_feedstock_des
  //T_Fe2O3_feedstock_des
  parameter SI.SpecificHeatCapacity cp_Hot_H2O_PGHX1 = 0.5*(Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, T_products_des) + Modelica.Media.Water.IF97_Utilities.cp_pT(p_H2O_offgas_des, max(Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p_H2O_offgas_des)+1.0,T_H2_mix_des))) "cp of the hot H2O entering the PGHX1 (J/kgK)";
  parameter SI.SpecificHeatCapacity cp_Hot_H2_PGHX1 = (Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_products_des) - Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_Fe2O3_feedstock_des)) / (T_products_des - T_H2_mix_des) "cp of the hot H2 entering the PGHX1 (J/kgK)";
  
  parameter SI.SpecificHeatCapacity cp_Cold_Fe2O3_PGHX1 = (Material_Fe2O3.h_Tf(T_products_des,0.0) - Material_Fe2O3.h_Tf(T_Fe2O3_feedstock_des,0.0))/(T_products_des - T_Fe2O3_feedstock_des) "cp of the cold Fe2O3 entering the PGHX1 (J/kgK)";
  
  parameter SI.ThermalConductance C_H_PGHX1 = (1-f_split_GGHX)*m_flow_H2_excess*cp_Hot_H2_PGHX1 + (1-f_split_GGHX)*m_flow_H2O_des*cp_Hot_H2O_PGHX1 "Heat capacity rate of the hot off-gas entering PGHX1 (W/K)";
  
  //The mass flow rate of Fe2O3 feedstock entering the PGHX1 is set to m_flow_Fe2O3_des
  parameter SI.MassFlowRate m_flow_Fe2O3_feedstock = m_flow_Fe2O3_des "Design mass flow rate of cold Fe2O3 feedstock entering PGHX1 (kg/s)";//C_H_PGHX1/cp_Cold_Fe2O3_PGHX1 ;
  
  //It is likely that C_C_PGHX1 is the C_min_PGHX, check that anyway
  parameter SI.ThermalConductance C_C_PGHX1 = m_flow_Fe2O3_feedstock*cp_Cold_Fe2O3_PGHX1;
  
  parameter SI.ThermalConductance C_min_PGHX1 = min(C_C_PGHX1,C_H_PGHX1);
  
  parameter SI.ThermalConductance C_max_PGHX1 = max(C_C_PGHX1,C_H_PGHX1);
  
  parameter Real C_r_PGHX1 = C_min_PGHX1/C_max_PGHX1;
  
  parameter SI.HeatFlowRate Q_flow_PGHX1_des = eff_PGHX1_des*C_min_PGHX1*(T_products_des-T_Fe2O3_feedstock_des);
  
  parameter Real NTU_PGHX1_des = (1.0/(C_r_PGHX1-1.0))*log((eff_PGHX1_des-1.0)/(eff_PGHX1_des*C_r_PGHX1-1.0)) "For PGHX1, C_r is unlikely to be one";
  
  //What is the enthalpy of Fe2O3 exiting the PGHX1 and entering Med Tank?
  parameter SI.SpecificEnthalpy h_Fe2O3_feedstock_des = Material_Fe2O3.h_Tf(T_Fe2O3_feedstock_des,0.0) "Specific enthalpy of Fe2O3 feedstock (J/kg)";
  //parameter SI.SpecificEnthalpy h_Fe2O3_PGHX1_max = Material_Fe2O3.h_Tf(T_products_des,0.0) "Maximum possible specific enthalpy of Fe2O3 leaving PGHX1 assuming 100% effectiveness (J/kg)";
  parameter SI.SpecificEnthalpy h_Fe2O3_med1_des = h_Fe2O3_feedstock_des + Q_flow_PGHX1_des/m_flow_Fe2O3_feedstock "Specific enthalpy of Fe2O3 leaving PGHX1 and entering the med tank (J/kg)";
  
  //h_Fe2O3_feedstock_des + eff_PGHX1_des*(h_Fe2O3_PGHX1_max - h_Fe2O3_feedstock_des) ;
  parameter SI.Temperature T_Fe2O3_med1_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_med1_des) "Temperature of Fe2O3 leaving PGHX1 and entering the med tank (K)";
  
  //What is the temperature of off-gas entering Condenser 1?
  //First calculate Q_flow_PGHX1_des
  
  //parameter SI.HeatFlowRate Q_flow_PGHX1_des = ((1-f_split_GGHX)*m_flow_H2_excess*cp_Hot_H2_PGHX1 + (1-f_split_GGHX)*m_flow_H2O_des*cp_Hot_H2O_PGHX1)*(T_products_des-T_Fe2O3_feedstock_des)*eff_PGHX1_des;
  //parameter SI.HeatFlowRate Q_flow_PGHX1_des = m_flow_Fe2O3_feedstock*eff_PGHX1_des*(h_Fe2O3_PGHX1_max - h_Fe2O3_feedstock_des) "Design heat transfer rate in PGHX1 (J/s)";
  //Now the model will try to find the temperature of cooled offgas leaving PGHX1
  SI.Temperature T_condenser1_in_des(start = 500.0) "Temperature of cooled offgas leaving PGHX1 and entering Condenser1 (K)";
   
  
  
  
  
  parameter SI.SpecificEnthalpy h_H2_pre2_des=h_H2_pre1_des + eff_PGHX2_des*(h_H2_PGHX2_max-h_H2_pre1_des) "Design of H2 leaving the PGHX2 cold stream after considering HX effectiveness and assuming C_H2 is C_min (J/kg)";

  //Storage Temperatures
  parameter SI.Temperature T_H2_pre1_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_pre1_des) "Design temperature of the H2 leaving the GGHX";
  
  
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
  
  parameter SI.SpecificEnthalpy h_H2_PGHX2_max=Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des,T_Fe2O3_hot_des) "Design specific enthalpy of H2 leaving the PGHX2 cold stream assuming it had been heated to the maximum possible temperature i.e. reactor product temperature (J/kg)";
  
  parameter SI.Temperature T_H2_pre2_des = Modelica.Media.IdealGases.SingleGases.H2.temperature_ph(p_des,h_H2_pre2_des) "Design temperature of the H2 leaving PGHX2 after considering PGHX2 effectiveness";
  
  //What is the expected temperature of materials entering the med tank assuming they get mixed well?
  parameter SI.SpecificEnthalpy h_Fe2O3_med_des = (m_flow_Fe2O3_feedstock*h_Fe2O3_med1_des + m_flow_Fe2O3_PGHX2*h_Fe2O3_med2_des)/(m_flow_Fe2O3_feedstock + m_flow_Fe2O3_PGHX2) "Expected specific enthalpy of the med tank (J/kg)";
  parameter SI.Temperature T_Fe2O3_med_des = SolarTherm.Media.SolidParticles.Fe2O3_utilities.T_h(h_Fe2O3_med_des) "Expected specific enthalpy of the med tank (K)";
  
  parameter SI.HeatFlowRate Q_flow_PGHX2_des = m_flow_Fe2O3_PGHX2*eff_PGHX2_des*(h_Fe2O3_hot_des - SolarTherm.Media.SolidParticles.Fe2O3_utilities.h_T(T_H2_pre1_des)) "Design heat transfer rate across PGHX2 (J/s)";
  
  parameter Real NTU_PGHX2_des = eff_PGHX2_des/(1.0-eff_PGHX2_des) "For PGHX2, C_r is one";
  
  parameter SI.ThermalConductance C_H_PGHX2 = m_flow_Fe2O3_PGHX2*cp_Fe2O3_PGHX2;
  
  
  //Linear Guess
  parameter SI.HeatFlowRate H_flow_inlet_des = m_flow_Fe2O3_des*h_Fe2O3_inlet_des + m_flow_H2_des*h_H2_inlet_des;
  
  parameter SI.Temperature T_Fe2O3_hot_guess = (H_flow_inlet_des+m_flow_H2_des*(h_H2_pre1_des*eff_PGHX2_des - h_H2_pre1_des) +261521.6*m_flow_Fe2O3_des+762976.9*m_flow_H2_des*eff_PGHX2_des)/(895.7669*m_flow_Fe2O3_des+15214.199*m_flow_H2_des*eff_PGHX2_des);
  
  //Condenser1
  SI.HeatFlowRate Q_flow_cooling_condenser1 = (1.0-f_split_GGHX)*m_flow_H2_excess*(Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_condenser1_in_des)-Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_condenser_out_des)) + (1.0-f_split_GGHX)*m_flow_H2O_des*(Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_condenser1_in_des) - Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_condenser_out_des));
  
  SI.HeatFlowRate Q_flow_cooling_condenser2 = f_split_GGHX*m_flow_H2_excess*(Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_condenser2_in_des)-Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_des, T_condenser_out_des)) + f_split_GGHX*m_flow_H2O_des*(Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_condenser2_in_des) - Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_condenser_out_des));
  
  //Electrical Heater
  //The design mass flow rate of the electrical heater is (m_flow_Fe2O3_des+m_flow_Fe2O3_PGHX2)
  //There is a heater multiple, assume 2.0 for now.
  //Heater efficiency, let's say 0.95 for now.
  //Initial temperature should be T_Fe2O3_med_des and it is heated up to T_Fe2O3_hot_des.\
  parameter Real HM = 2.0 "Heater Multiple";
  parameter Real eff_heater = 0.95 "Heater electrical to thermal efficiency";
  
  parameter SI.HeatFlowRate Q_flow_heater_curtail = (m_flow_Fe2O3_des+m_flow_Fe2O3_PGHX2)*(h_Fe2O3_hot_des - h_Fe2O3_med_des) "Heating rate required during steady operation where the hot tank is full (J/s)";
  parameter SI.HeatFlowRate Q_flow_heater_max = HM*Q_flow_heater_curtail "Maximum heat output rate of the heater (J/s)";
  
  parameter SI.Power P_heater_curtail = Q_flow_heater_curtail/eff_heater "Heater power input when curtailed (W)";
  
  parameter SI.Power P_heater_max = Q_flow_heater_max/eff_heater "Maximum heater power input (W)";
  
    
    
  //Cost and Sizing Parameters
  parameter Real CEPCI = 816.0 "CEPCI index of the year used in the study e.g. 816.0 for year 2022";
  //PGHX1 assume C_r is 1.0
  parameter SI.CoefficientOfHeatTransfer U_GGHX_des = 30.0 "Overall heat transfer coefficient in the GGHX (W/m2K)";
   
  parameter Real NTU_GGHX_des = eff_GGHX_des/(1.0 - eff_GGHX_des) "Design NTU of the GGHX (-)";
  parameter SI.ThermalConductance UA_GGHX_des = C_min_GGHX_des*NTU_GGHX_des "Design UA of the GGHX (W/K)";
  parameter SI.Area A_GGHX_des = UA_GGHX_des/U_GGHX_des "Design heat transfer area of the GGHX (m2)"; 
  parameter Real FOB_GGHX = div(A_GGHX_des, 185.8)*(CEPCI/500.0)*150945.38 + (CEPCI/500.0)*6200.0*((10.764*rem(A_GGHX_des,185.8))^0.42) "Free-on-board cost of the GGHX component, based on a Spiral-Plate HX (USD_year)";
  parameter Real FCI_GGHX = FOB_GGHX*1.05*2.4541;
   
  //H2 Blower cost, based on cast iron, 3psig (1.22 bar absolute)
  parameter SI.Density rho_H2_mix_des = Modelica.Media.IdealGases.SingleGases.H2.density_pT(p_des,T_H2_mix_des) "Density of H2 after the mixer (kg/m3)";
  parameter SI.VolumeFlowRate V_flow_H2_mix_des = m_flow_H2_des/rho_H2_mix_des "Volumetric flow rate of H2 in the blower fan (m3/s)";
  parameter Real N_blower_quo = div(V_flow_H2_mix_des, 23.6) "Number of maximum-sized blower fans (max 23.6 m3/s vol flow rate)";
  parameter SI.VolumeFlowRate V_blower_rem = rem(V_flow_H2_mix_des, 23.6) "Remaining volumetric flow rate of the blowers (m3/s)";
  
  parameter Real FOB_blower = if V_blower_rem > 2.36 then N_blower_quo*(CEPCI/460.0)*23.0*(2119.0*23.6)^0.82 + (CEPCI/460.0)*23.0*(2119.0*V_blower_rem)^0.82 else N_blower_quo*(CEPCI/460.0)*23.0*(2119.0*23.6)^0.82 + (CEPCI/460.0)*516.0*(2119.0*V_blower_rem)^0.46 "Free-on-board cost of the blowers, cast iron, (USD_year)";
  
  parameter Real FCI_blower = FOB_blower*1.05*2.5;
  
  //Hot Tank Cost
  parameter SI.Time t_stor_Fe2O3_hot = 10.0*3600.0 "Number of seconds of storage of hot-tank Fe2O3 (s)";
  //parameter Real ar_tank_hot = 2.0 "Aspect ratio of hot tank (-)";
  
  parameter Real porosity_Fe2O3_material = 1.0 - (4300.0/Material_Fe2O3.rho_Tf(900.0+273.15,0.0)) "Material porosity of Fe2O3 based on measurements at 900degC obtained from Mahdiar (-)";
  parameter Real porosity_Fe2O3_packing = 0.20 "Packing porosity of Fe2O3 packed bed (-)";
  parameter SI.Density rho_Fe2O3_hot = Material_Fe2O3.rho_Tf(T_Fe2O3_hot_des,0.0) "Density of pure, dense Fe2O3 at the hot tank temperature (kg/m3)";
  
  parameter SI.Volume V_tank_hot = (m_flow_Fe2O3_des+m_flow_Fe2O3_PGHX2)*t_stor_Fe2O3_hot/(rho_Fe2O3_hot*(1.0-porosity_Fe2O3_material)*(1.0-porosity_Fe2O3_packing));
  
  parameter Real FOB_tank_hot = (CEPCI/500)*2.1*570.0*(35.315*V_tank_hot)^0.46 "FOB cost of the hot Fe2O3 storage tank (USD_year)";
  parameter Real FCI_tank_hot = FOB_tank_hot*1.05*2.744;
  
  //Med Tank Cost
  parameter SI.Time t_stor_Fe2O3_med = 10.0*3600.0 "Number of seconds of storage of med-tank Fe2O3 (s)";
  parameter SI.Density rho_Fe2O3_med = Material_Fe2O3.rho_Tf(T_Fe2O3_med_des,0.0) "Density of pure, dense Fe2O3 at the hot tank temperature (kg/m3)";
  parameter SI.Volume V_tank_med = (m_flow_Fe2O3_des+m_flow_Fe2O3_PGHX2)*t_stor_Fe2O3_med/(rho_Fe2O3_med*(1.0-porosity_Fe2O3_material)*(1.0-porosity_Fe2O3_packing));
    
  parameter Real FOB_tank_med = (CEPCI/500)*1.0*570.0*(35.315*V_tank_hot)^0.46 "FOB cost of the medium temp Fe2O3 storage tank (USD_year)";
  parameter Real FCI_tank_med = FOB_tank_med*1.05*4.0;
  //Reactor Cost
  //Formula is (816/1000)*14000000*(m_flow_Fe2O3_des^0.53) where m_flow is up to 1000kg/s (for the first stage)
  //Formula is (816/1000)*14000000*(m_flow_Fe3O4_des^0.53) where m_flow is up to 1000kg/s (for the second stage)
  //Formula is (816/1000)*14000000*(m_flow_FeO_des^0.53) where m_flow is up to 1000kg/s (for the third stage)
  parameter SI.MassFlowRate m_flow_Fe3O4_des = (2.0*m_flow_Fe2O3_des*M_Fe3O4)/(3.0*M_Fe2O3);
  parameter SI.MassFlowRate m_flow_FeO_des = (2.0*m_flow_Fe2O3_des*M_FeO)/(M_Fe2O3);
  
  parameter Real FOB_reactor = (CEPCI/1000.0)*(14.0e6*(m_flow_Fe2O3_des^0.53 + m_flow_Fe3O4_des^0.53 + m_flow_FeO_des^0.53)) "FOB Cost of the reactors (USD_year)";
  parameter Real FCI_reactor = FOB_reactor*1.05*2.0;
  
  //Heater Cost
  parameter Real pri_heater = 0.150 "Cost per W of heater. Assumed to be 200 USD_2022/kW";
  parameter Real FOB_heater = (CEPCI/816.0)*pri_heater*P_heater_max "FOB cost of the electrical heater (USD_year)";
  parameter Real FCI_heater = FOB_heater*1.05*4.0;
  
  //Condenser Q_flow_cooling = U_condenser_des*A_condenser1*(T_condenser1_in_des-T_amb_des)
  parameter SI.CoefficientOfHeatTransfer U_condenser1_des = 700.0 "W/m2K";
  parameter SI.CoefficientOfHeatTransfer U_condenser2_des = 700.0 "W/m2K";
  SI.Area A_condenser1 = Q_flow_cooling_condenser1/(U_condenser1_des*(T_condenser1_in_des-T_amb_des)) "Required heat transfer area of air-cooled condenser (m2)";
  SI.Area A_condenser2 = Q_flow_cooling_condenser2/(U_condenser2_des*(T_condenser2_in_des-T_amb_des)) "Required heat transfer area of air-cooled condenser (m2)";
  
  Real FOB_condenser1 = (CEPCI/500.0)*10000.0*((10.764*A_condenser1)^0.40);
  Real FOB_condenser2 = (CEPCI/500.0)*10000.0*((10.764*A_condenser2)^0.40);
  
  Real FCI_condenser1 = FOB_condenser1*1.05*2.4541;
  Real FCI_condenser2 = FOB_condenser2*1.05*2.4541;
  
  //Cost of PGHX1
  parameter SI.ThermalConductance U_PGHX1_des = 373.96 "Overall heat transfer coefficient of PGHX1 (W/K)";
  parameter SI.Area A_PGHX1_des = C_min_PGHX1*NTU_PGHX1_des/U_PGHX1_des;
  parameter Real FOB_PGHX1 = (CEPCI/708.8)*(1909.21*A_PGHX1_des);
  parameter Real FCI_PGHX1 = FOB_PGHX1*1.05*1.8688;
  
  //Cost of PGHX2
  parameter SI.ThermalConductance U_PGHX2_des = 533.96 "Overall heat transfer coefficient of PGHX2 (W/K)";
  parameter SI.Area A_PGHX2_des = C_H_PGHX2*NTU_PGHX2_des/U_PGHX2_des;
  parameter Real FOB_PGHX2 = (CEPCI/708.8)*(1856.0*A_PGHX2_des);
  parameter Real FCI_PGHX2 = FOB_PGHX2*1.05*1.8688;
  
  
  
  
  
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
  
  //Calculate Temperature of H2 and H2O leaving PGHX1 and entering condenser 1:
  Q_flow_PGHX1_des = (1.0 - f_split_GGHX)*m_flow_H2_excess*(Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_products_des)-Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_condenser1_in_des)) + (1.0 - f_split_GGHX)*m_flow_H2O_des*(Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_products_des)-Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_condenser1_in_des));
  
  //Calculate Temperature of H2 and H2O leaving GGHX and entering condenser 2
  Q_flow_GGHX_des = f_split_GGHX*m_flow_H2_excess*(Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_products_des)-Modelica.Media.IdealGases.SingleGases.H2.specificEnthalpy_pT(p_H2_offgas_des, T_condenser2_in_des)) +  f_split_GGHX*m_flow_H2O_des*(Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_products_des)-Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p_H2O_offgas_des, T_condenser2_in_des));

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-150, -100}, {150, 100}})),
    Icon(coordinateSystem(extent = {{-150, -100}, {150, 100}}, preserveAspectRatio = false)));
end H2DRI_DesignCase_2a_DesignPt;