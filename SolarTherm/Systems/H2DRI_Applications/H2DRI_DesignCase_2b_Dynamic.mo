within SolarTherm.Systems.H2DRI_Applications;

model H2DRI_DesignCase_2b_Dynamic
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  
  //Material Packages
  replaceable package Material_IOE_OreD = SolarTherm.Materials.IOE_Dehydroxylated;
  replaceable package Medium_Ore_Dehydroxylated = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_ph;
  
  parameter Real Plant_Scale = 1.0 "Plant DRI design production rate in (Mt_DRI/yr)";
  parameter String PV_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/PV_Pilbara_1MW.motab");
  parameter String Wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/Wind_Pilbara_320MW.motab");
  
  parameter Real CEPCI = 816.0 "CEPCI index of the year used in the study e.g. 816.0 for year 2022";
  
  //Free Parameter Inputs
  parameter Real RM = 1.2 "Renewable Multiple (pre-transmission oversizing)";
  parameter Real HM = 2.0 "Heater Multiple";
  parameter SI.Time t_storage = 70.0*3600.0 "Seconds of storage (h)";
  parameter Real PV_fraction = 0.4 "PV_fraction";
  
  //DRI Properties
  parameter Real DRI_f_Fe = 0.8530 "Mass fraction of Fe in the DRI (-)";
  
  parameter SI.HeatFlowRate Q_process_des = 6.105e7;// 6.14926e7;
  parameter SI.MassFlowRate m_flow_ore_des = 131.4894*(Plant_Scale/1.0) "Mass flow of dehydroxylated ore out of the hot tank and into the med tank if running at design point (kg/s)";
  
  parameter SI.MassFlowRate m_flow_ore_reactor_des = 43.34*(Plant_Scale/1.0) "Mass flow of dehydroxylated ore into the reactor if running at design point (kg/s)";
  parameter SI.MassFlowRate m_flow_DRI_des = 31.71*(Plant_Scale/1.0) "Mass flow rate of DRI produced by the reactor at design point, this mass includes Al2O3 and SiO2 (kg/s)";
  parameter SI.MassFlowRate m_flow_Fe_des = DRI_f_Fe*m_flow_DRI_des "Mass flow rate of metallic iron produced by the reactor at design point (kg/s)";
  parameter SI.MassFlowRate m_flow_H2_consumed_des = 1.5*(m_flow_Fe_des/55.845e-3)*(2.01588e-3) "Mass flow rate of stoichiometric H2 consumed by the reactor at design point (kg/s)"; //1.5 moles of H2 are needed per mol of Fe produced".
 //1.5 moles of H2 are needed per mol of Fe produced".
  
  parameter Real eff_heater = 0.95 "Electrical-to-heat conversion efficiency of the heater";
  //Renewable Parameters
  parameter SI.Power P_renewable_des = RM * P_heater_des;
  parameter SI.HeatFlowRate Q_heater_des = HM * Q_process_des;
  parameter SI.Power P_heater_des = Q_heater_des / eff_heater;
  parameter SI.Power PV_ref_size = 1.0e6; //1MW reference size
  parameter SI.Power Wind_ref_size = 320.0e6; //320MW reference size
  parameter Real LOF_PV = 1.300; //Loss oversize factor, so we actually get the required power from PV field
  parameter Real LOF_Wind = 1.131; //Loss oversize factor, so we actually get the required power from Wind field
  parameter SI.Power P_wind_net = (1.0 - PV_fraction)*P_renewable_des;
  parameter SI.Power P_PV_net = PV_fraction*P_renewable_des;
  parameter SI.Power P_wind_gross = P_wind_net*LOF_Wind;
  parameter SI.Power P_PV_gross = P_PV_net*LOF_PV;
  
  //Tank Parameters
  parameter SI.Temperature T_hot_set = 757.6 + 273.15;
  parameter SI.Temperature T_med_set = 292.7 + 273.15;//289.0 + 273.15;
  
  parameter Real ar = 2.0;
  
  //Sizing Parameters
    //Storage Bins
  parameter Real eps_packing = 0.20;
  parameter Real eps_material = 0.00;
  parameter Real epsilon = eps_packing + eps_material - eps_packing*eps_material "Effective porosity";
  parameter SI.Density rho_ore_med = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_utilities.rho_T(T_med_set) "Density of pure ore (kg/m3)";
  parameter SI.Density rho_ore_hot = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_utilities.rho_T(T_hot_set) "Density of pure ore (kg/m3)";
  
  parameter SI.Mass m_ore_tank_hot = m_flow_ore_des*t_storage "Maximum mass capacity of the hot tank (kg)";
  parameter SI.Mass m_ore_tank_med = m_flow_ore_des*t_storage "Maximum mass capacity of the hot tank (kg)";
  
  parameter SI.Volume V_tank_hot = m_ore_tank_hot/(rho_ore_hot*(1.0-epsilon));
  parameter SI.Volume V_tank_med = m_ore_tank_med/(rho_ore_med*(1.0-epsilon));
  
  parameter SI.Length D_tank_hot = (4.0*V_tank_hot/(CN.pi*ar))^(1.0/3.0);
  
  parameter SI.Length D_tank_med = (4.0*V_tank_med/(CN.pi*ar))^(1.0/3.0);
  
  parameter SI.Length H_tank_hot = D_tank_hot*ar;
  parameter SI.Length H_tank_med = D_tank_med*ar;
  
  parameter SI.Length D_tank_max = ((4.0*2832.0)/(CN.pi*ar))^(1.0/3.0) "Maximum tank diameter (m3)";
  
  //Cold Tank
  parameter Real N_quo_med = div(V_tank_med,2832.0);
  parameter SI.Volume V_rem_med = rem(V_tank_med,2832.0);
  parameter SI.Length D_rem_med = ((4.0*V_rem_med)/(CN.pi*ar))^(1.0/3.0) "Remainder cold tank diameter (m3)";
  parameter SI.Area A_loss_med_total = N_quo_med*(CN.pi*D_tank_max*D_tank_max*ar + 0.5*CN.pi*D_tank_max*D_tank_max) + CN.pi*D_rem_med*D_rem_med*ar + 0.5*CN.pi*D_rem_med*D_rem_med;
  
  //Hot Tank
  parameter Real N_quo_hot = div(V_tank_hot,2832.0);
  parameter SI.Volume V_rem_hot = rem(V_tank_hot,2832.0);
  parameter SI.Length D_rem_hot = ((4.0*V_rem_hot)/(CN.pi*ar))^(1.0/3.0) "Remainder hot tank diameter (m3)";
  parameter SI.Area A_loss_hot_total = N_quo_hot*(CN.pi*D_tank_max*D_tank_max*ar + 0.5*CN.pi*D_tank_max*D_tank_max) + CN.pi*D_rem_hot*D_rem_hot*ar + 0.5*CN.pi*D_rem_hot*D_rem_hot;
  
  parameter Real FOB_tank_med = (CEPCI/500)*1.0*(N_quo_med*113736.27 + 570.0*(35.315*V_rem_med)^0.46) "FOB cost of the medium temp Fe2O3 storage tank (USD_year)";
  
  parameter Real FOB_tank_hot = (CEPCI/500)*2.1*(N_quo_hot*113736.27 + 570.0*(35.315*V_rem_med)^0.46) "FOB cost of the hot temp Fe2O3 storage tank (USD_year)";
  
  //parameter Real FOB_tank_hot = (CEPCI/500)*2.1*570.0*(35.315*V_tank_hot)^0.46 "FOB cost of the hot Fe2O3 storage tank (USD_year)";
  
  //parameter Real FOB_tank_med = (CEPCI/500)*1.0*570.0*(35.315*V_tank_med)^0.46 "FOB cost of the medium temp Fe2O3 storage tank (USD_year)";
  
    //Fluidised Bed Heater
  parameter SI.HeatFlux q_flow_heater_max = 60000.0 "Maximum radiant heat flux of the fluidised bed heater (W/m2)"; //Placeholder
  parameter SI.Area A_cs_FB = Q_heater_des/q_flow_heater_max "Minimum cross sectional area of the fluidised bed (m2)";
  parameter SI.Diameter d_p = 2.5e-4 "Assumed ore particle diameter (m)"; //250 micrometres
  parameter SI.Density rho_p_FB = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_utilities.rho_T(0.5*(T_hot_set+T_med_set));
  
  parameter SI.Density rho_g_FB = SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.rho_T(0.5*(T_hot_set+T_med_set));
  parameter SI.DynamicViscosity mu_g_FB = SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.mu_T(0.5*(T_hot_set+T_med_set));
  
  parameter SI.Velocity u_mf_FB = (d_p*d_p*(rho_p_FB-rho_g_FB)*9.81)/(1650.0*mu_g_FB);
  parameter SI.Velocity u_g_FB = 3.0*u_mf_FB;
  
  parameter SI.MassFlowRate m_flow_air_FB = rho_g_FB*A_cs_FB*u_g_FB "Mass flow rate through FB blower (kg/s)";
  
  parameter SI.VolumeFlowRate V_flow_air_FB = m_flow_air_FB/SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.rho_T(298.15) "Volumetric flow rate of ambient temp air into FB blower (m3/s)";
  
  parameter SI.Power P_C_FB = (0.9855*(1.41/0.41)*(V_flow_air_FB*1.0e5/0.75)*(((1.1/1.0)^(0.41/1.41))-1.0))/0.9 "Sizing power of blower of FB heater (W)";
  parameter Real FOB_blower_FB = (CEPCI/500.0)*1.0*exp(6.8929+0.79*log(P_C_FB/745.7)) "FOB cost of the blower in FB heater (USD_year)";
  //PGHX2's max temp range is T_amb_des to T_OreD_hot_des effectiveness is assumed to be 0.80
  parameter SI.HeatFlowRate Q_flow_recup_FB = 0.8*m_flow_air_FB*(SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.h_T(T_hot_set)-SolarTherm.Media.Air.Air_CoolProp_1bar_utilities.h_T(298.15));
  parameter SI.ThermalConductance U_recup_FB = 30.0 "Gas-gas heat transfer coefficient (W/m2K)";
  parameter SI.Area A_recup_FB = Q_flow_recup_FB/(0.80*U_recup_FB*(T_hot_set-298.15)) "Required gas-gas heat recuperator heat exchanger area of FB heater (m2), assumed 0.80 effectiveness";
  parameter Real FOB_recup_FB = div(A_recup_FB, 185.8)*(CEPCI/500.0)*150945.38 + (CEPCI/500.0)*6200.0*((10.764*rem(A_recup_FB,185.8))^0.42);
  
  
  
  
  //Cost parameters
  parameter Real pri_PV = 1.10725 "Cost per W_gross of PV plant, already including 3% contingency (USD_2022/W)";
  parameter Real pri_Wind = 1.50607 "Cost per W_gross of PV plant already including 3% contingency (USD_2022/W)";
  parameter Real pri_FB_Heating = 0.1539 "Cost per W of heater (USD_2022/W)";
  
  //Fixed-Size Capital Costs
  parameter Real FCI_Reactor = 506546275.0 "Reactor FCI cost (USD_2022)";
  parameter Real FCI_GGHX = 41675600.0 "GGHX FCI cost (USD_2022)";
  parameter Real FCI_Blower_H2 = 1404500.0 "H2 Blower FCI cost (USD_2022)";
  parameter Real FCI_Condenser_1 = 702567.0 "Condenser 1 cost (USD_2022)";
  parameter Real FCI_Condenser_2 = 2178690.0 "Condenser 2 cost (USD_2022)";
  parameter Real FCI_PGHX1 = 1128020.0 + 595218.0 + 147284.0 "PGHX1 cost (USD_2022)";
  parameter Real FCI_PGHX2 = 23195900.0 + 6255300.0 + 977190.0 "PGHX2 cost (USD_2022)";

  //Variable-Sized Capital Costs
  parameter Real FCI_heating_FB = pri_FB_Heating*P_heater_des;
  parameter Real FCI_blower_FB = FOB_blower_FB*1.05*3.5*0.7012 "FCI cost of the blower in FB heater (USD_year)";
  parameter Real FCI_recup_FB = FOB_recup_FB*1.05*3.5*0.7012;
  parameter Real FCI_tank_hot = FOB_tank_hot*1.05*2.744;
  parameter Real FCI_tank_med = FOB_tank_med*1.05*4.0;
  
  //PV and Wind Capital Costs
  parameter Real FCI_PV = pri_PV*P_PV_gross;
  parameter Real FCI_Wind = pri_Wind*P_wind_gross;
  
  //Crushing Capital Cost
  parameter Real FCI_Crushing = 50548009.0*Plant_Scale/1.0;
  
  //DRI Plant Capital Costs
  parameter Real FCI_Plant = FCI_Reactor + FCI_GGHX + FCI_Blower_H2 + FCI_Condenser_1 + FCI_Condenser_2 + FCI_PGHX1 + FCI_PGHX2 + FCI_heating_FB + FCI_blower_FB + FCI_recup_FB + FCI_tank_hot + FCI_tank_med;
  
  //FCI_Plant, FCI_PV and FCI_Wind make up total CAPEX.
  
  //Fixed Annual Costs
  parameter Real AC_Plant = 0.04*FCI_Plant "Annual cost due to plant O&M, 4% of Plant CAPEX (USD/year)";
  parameter Real AC_Labour = 6613253.0*Plant_Scale/1.0 "Annual cost due to labour needed to run plant (USD/year)";
  parameter Real AC_PV = 0.01268*P_PV_gross "Annual O&M costs for PV plant (USD/year)";
  parameter Real AC_Wind = 0.01868*P_wind_gross "Annual O&M costs for Wind plant (USD/year)";
  
  //Variable Annual Costs
  Real AC_H2 = 183801405.0*CapF_Process*(Plant_Scale/1.0) "Variable annual costs due to stoichiometric consumption of H2 (USD/yr)";
  
  //m_flow_H2_consumed_des*(86400.0*365.0)*CapF_Process*3.50*(816.0/708.8);
  
  //218178221.0*0.8530*CapF_Process*(Plant_Scale/1.0) "Variable annual costs due to stoichiometric consumption of H2 (USD/yr)"; //0.8530kg of Fe per 1.0kg of DRI
  Real AC_Mining = 23254357.0*CapF_Process*(Plant_Scale/1.0) "Variable annual costs due to stoichiometric mining of iron ore (USD/yr)";
  Real AC_Electric = 4183700.0*CapF_Process*(Plant_Scale/1.0) "Variable annual costs due to electricity cost of processing iron ore (USD/yr)";
  
  
  //Sum Up Everything
  parameter Real C_capital = FCI_Plant + FCI_PV + FCI_Wind + FCI_Crushing;
  Real C_annual = AC_Plant + AC_Labour + AC_PV + AC_Wind + AC_H2 + AC_Mining + AC_Electric;
  
  //Some Financial Parameters
  parameter Real n = 30.0 "Plant lifetime (years)";
  parameter Real r_nom = 0.07 "Nominal discount rate (-)";
  parameter Real r_inf = 0.025 "Inflation rate (-)";
  
  parameter Real r = ((1.0+r_nom)/(1.0+r_inf)) - 1 "Real discount rate (-)";
  parameter Real f = (r*((1.0+r)^n))/(((1.0+r)^n)-1.0) "Annuity factor on LCOD calculation (-)";
  
  
  
  
  
  
  Modelica.Blocks.Sources.CombiTimeTable PV_input(fileName = PV_file, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "Power", tableOnFile = true) annotation(
    Placement(visible = true, transformation(origin = {-88, 54}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.CombiTimeTable Wind_input(fileName = Wind_file, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative, tableName = "Power", tableOnFile = true) annotation(
    Placement(visible = true, transformation(origin = {-88, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Add Grid_Sum(k1 = P_PV_gross / PV_ref_size, k2 = P_wind_gross / Wind_ref_size) annotation(
    Placement(visible = true, transformation(origin = {-51, 39}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
  SolarTherm.Models.CSP.CRS.Receivers.Basic_Heater Heater(redeclare package Medium = Medium_Ore_Dehydroxylated, P_heater_des = P_heater_des, Q_heater_des = Q_heater_des, T_cold_set = T_med_set, T_hot_set = T_hot_set, eff_heater = eff_heater) annotation(
    Placement(visible = true, transformation(origin = {14, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Storage.Tank.Silo Cold_Tank(redeclare package Medium = Medium_Ore_Dehydroxylated, redeclare package Filler_Package = Material_IOE_OreD, T_min = T_med_set - 10.0, T_max = T_med_set + 10.0, U_loss_tank = 0.0, T_start = T_med_set, T_set = T_med_set - 10.0, L_start = 0.50, epsilon = epsilon, H_tank = H_tank_med, D_tank = D_tank_med) annotation(
    Placement(visible = true, transformation(origin = {-32, -8}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  SolarTherm.Models.Storage.Tank.Silo Hot_Tank(redeclare package Medium = Medium_Ore_Dehydroxylated, redeclare package Filler_Package = Material_IOE_OreD, T_min = T_hot_set - 10.0, T_max = T_hot_set + 10.0, U_loss_tank = 0.0, T_start = T_hot_set, T_set = T_hot_set - 10.0, L_start = 0.50, epsilon = epsilon, H_tank = H_tank_hot, D_tank = D_tank_hot) annotation(
    Placement(visible = true, transformation(origin = {40, -6}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Heater_Lift(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {-10, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  SolarTherm.Models.Fluid.Sources.FluidSink2 Sink(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {90, -12}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT OreD_source(redeclare package Medium = Medium_Ore_Dehydroxylated, T = T_med_set, nPorts = 1, p = 100000, use_T_in = false) annotation(
    Placement(visible = true, transformation(origin = {-88, -12}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple Cold_Lift(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {-56, -14}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure Reactor_Lift(redeclare package Medium = Medium_Ore_Dehydroxylated) annotation(
    Placement(visible = true, transformation(origin = {60, -12}, extent = {{-6, -6}, {6, 6}}, rotation = 0)));

//Controller
  parameter SI.Time t_wait = 1.0*3600 "Waiting time between turning off process and being able to turn on";
  SI.Time t_threshold(start=0.0) "if time passes this value, Process := true";

  Boolean Chg(start = true);
  Boolean Dis(start = true);
  
  Boolean Process(start = true);
  
  Integer State(start = 1);
  
  SI.Mass m_DRI_produced(start=0);
  SI.Mass m_DRI_target(start=0);
  Real CapF_Process(start=0);
  Real CapF_Heater(start=0);
  //Boolean Defocus(start = false);
  
  
  //if Heater.Q_heater_raw <= 0.0 then
    //if Dis == true then
      //State = 4;
    //end if;
  //end if;
  
  //Energy Accounting
  SI.Energy E_PV_out(start=0);
  SI.Energy E_Wind_out(start=0);
  SI.Energy E_renewable_raw(start=0);
  SI.Energy E_heater_raw(start=0);
  SI.Energy Q_heater_raw(start=0);
  SI.Energy Q_heater_out(start=0);
  
  SI.Energy Q_heater_target(start=0); //maximum possible heater output at 100% operation
  
  //LCOD Calculation
  Real LCOD_2022 "Levelised cost per ton of DRI (USD_2022/tDRI)";
  Real LCOD_2021 = LCOD_2022*(708.8/816.0) "Levelised cost per ton of DRI (USD_2021/tDRI)"; 
  Real LCOD_2021_AUD = LCOD_2021*1.330849 "Levelised cost per ton of DRI (AUD_2021/tDRI)";
  Real LCOD_2022_numerator = f*C_capital + C_annual "Numerator of the LCOD formula (USD/yr)";
  
  //These are annual costs
  Real LCOD_frac_Hydrogen = AC_H2/LCOD_2022_numerator;
  Real LCOD_frac_Mining = AC_Mining/LCOD_2022_numerator;
  Real LCOD_frac_Electricity = AC_Electric/LCOD_2022_numerator;
  
  Real LCOD_frac_PlantOM = AC_Plant/LCOD_2022_numerator;
  Real LCOD_frac_Labour = AC_Labour/LCOD_2022_numerator;
  
  //These are capital costs
  Real LCOD_frac_Plant = f*FCI_Plant/LCOD_2022_numerator;
  Real LCOD_frac_Crushing = f*FCI_Crushing/LCOD_2022_numerator;
  
  //These are a combination
  Real LCOD_frac_PV = (f*FCI_PV + AC_PV)/LCOD_2022_numerator;
  Real LCOD_frac_Wind = (f*FCI_Wind + AC_Wind)/LCOD_2022_numerator;
  
  //Check if they sum to 1.0
  Real Sum_frac = LCOD_frac_Hydrogen + LCOD_frac_Mining + LCOD_frac_Electricity + LCOD_frac_PlantOM + LCOD_frac_Labour + LCOD_frac_Plant + LCOD_frac_PV + LCOD_frac_Wind + LCOD_frac_Crushing;
  
algorithm


  //Boiler Timer Control
  when Reactor_Lift.m_flow < 0.01 * m_flow_ore_des then //take this as shutdown
    Process := false; //start the cooldown
    t_threshold := time + t_wait;
  end when;
  when time > t_threshold then
    Process := true;
  end when;


  when Hot_Tank.L > 0.98 then
    Chg := false;
  end when;
  when Hot_Tank.L < 0.95 then
    Chg := true;
  end when;
  
  when Hot_Tank.L < 0.02 then
    Dis := false;
  end when;
  when Hot_Tank.L > 0.05 then
    Dis := true;
  end when;

equation
  if time > 86400.0*10.0 then //10 days initialisation has finished
    der(E_PV_out) = (Grid_Sum.u1/PV_ref_size)*P_PV_gross;
    der(E_Wind_out) = (Grid_Sum.u2/Wind_ref_size)*P_wind_gross;

    der(E_renewable_raw) = Grid_Sum.y;
    der(E_heater_raw) = Heater.P_heater_out;
    der(Q_heater_raw) = Heater.Q_heater_raw;
    der(Q_heater_out) = Heater.Q_out;

    der(Q_heater_target) = Q_heater_des;
    der(m_DRI_target) = m_flow_ore_des*(50.8/154.1)*(0.8924*2.0*55.845/159.6882)*(1.0/DRI_f_Fe);
    der(m_DRI_produced) = Sink.port_a.m_flow*(50.8/154.1)*(0.8924*2.0*55.845/159.6882)*(1.0/DRI_f_Fe);
  else
    der(E_PV_out) = 0.0;
    der(E_Wind_out) = 0.0;
    der(E_renewable_raw) = 0.0;
    der(E_heater_raw) = 0.0;
    der(Q_heater_raw) = 0.0;
    der(Q_heater_out) = 0.0;

    der(Q_heater_target) = 0.0;
    der(m_DRI_target) = 0.0;
    der(m_DRI_produced) = 0.0;
  end if;
  
  if time < 10.0 + 86400.0*10.0 then
    CapF_Process = 0.0;
    CapF_Heater = 0.0;
    LCOD_2022 = 0.0;
  else
    CapF_Process = m_DRI_produced/m_DRI_target;
    CapF_Heater = Q_heater_out/Q_heater_target;
    LCOD_2022 = (f*C_capital + C_annual)/(1.0e-3*m_DRI_produced);
  end if;


  if Heater.Q_heater_raw <= 1.0e-6 then
    if Dis == true then
      if Process == true then
        State = 4;
      else
        State = 6;
      end if;
    else
      State = 6;
    end if;
  elseif Heater.Q_heater_raw >= Q_process_des then
    if Process == true then
      if Chg == true then
        State = 1;
      else
        State = 3;
      end if;
    else
      if Chg == true then
        State = 5;
      else
        State = 6;
      end if;
    end if;
  else
    if Dis == true then
      if Process == true then
        State = 2;
      else
        State = 6;
      end if;
    else
      if Chg == true then
        State = 5;
      else
        State = 6;
      end if;
    end if;
  end if;
  if State == 1 then
    Heater_Lift.m_flow = Heater.Q_heater_raw / Q_process_des * m_flow_ore_des;
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des;
//Not used anyway
    Cold_Lift.m_flow = m_flow_ore_des;
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 2 then
    Heater_Lift.m_flow = Heater.Q_heater_raw / Q_process_des * m_flow_ore_des;
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des;
//Not used anyway
    Cold_Lift.m_flow = m_flow_ore_des;
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 3 then
    Heater_Lift.m_flow = m_flow_ore_des;
    Heater.curtail = true;
    Heater.Q_curtail = Q_process_des;
//Curtailed
    Cold_Lift.m_flow = m_flow_ore_des;
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 4 then
    Heater_Lift.m_flow = 1.0e-9;
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des;
//Not used anyway
    Cold_Lift.m_flow = m_flow_ore_des;
    Reactor_Lift.m_flow = m_flow_ore_des;
  elseif State == 5 then
    Heater_Lift.m_flow = Heater.Q_heater_raw / Q_process_des * m_flow_ore_des;
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des;
//Not used anyway
    Cold_Lift.m_flow = 1.0e-9;
    Reactor_Lift.m_flow = 1.0e-9;
  else
    Heater_Lift.m_flow = 1.0e-9;
    Heater.curtail = false;
    Heater.Q_curtail = Q_process_des;
//Not used anyway
    Cold_Lift.m_flow = 1.0e-9;
    Reactor_Lift.m_flow = 1.0e-9;
  end if;
  Cold_Tank.T_amb = 298.15;
  Hot_Tank.T_amb = 298.15;
  Cold_Tank.p = 100000.0;
  Hot_Tank.p = 100000.0;
  connect(Wind_input.y[1], Grid_Sum.u2) annotation(
    Line(points = {{-76, 26}, {-70, 26}, {-70, 34}, {-62, 34}}, color = {0, 0, 127}));
  connect(PV_input.y[1], Grid_Sum.u1) annotation(
    Line(points = {{-76, 54}, {-69.5, 54}, {-69.5, 44}, {-62, 44}}, color = {0, 0, 127}));
  connect(Cold_Lift.fluid_b, Cold_Tank.fluid_a) annotation(
    Line(points = {{-50, -14}, {-40, -14}, {-40, -4}}, color = {0, 127, 255}));
  connect(OreD_source.ports[1], Cold_Lift.fluid_a) annotation(
    Line(points = {{-80, -12}, {-62, -12}, {-62, -14}}, color = {0, 127, 255}));
  connect(Reactor_Lift.fluid_b, Sink.port_a) annotation(
    Line(points = {{66, -12}, {76, -12}}, color = {0, 127, 255}));
  connect(Hot_Tank.fluid_b, Reactor_Lift.fluid_a) annotation(
    Line(points = {{48, -12}, {54, -12}}, color = {0, 127, 255}));
  connect(Grid_Sum.y, Heater.P_supply) annotation(
    Line(points = {{-41, 39}, {-6, 39}, {-6, 1}, {3, 1}}, color = {0, 0, 127}));
  connect(Heater.fluid_b, Hot_Tank.fluid_a) annotation(
    Line(points = {{24, -6}, {32, -6}, {32, -2}}, color = {0, 127, 255}));
  connect(Cold_Tank.fluid_b, Heater_Lift.fluid_a) annotation(
    Line(points = {{-24, -14}, {-16, -14}}, color = {0, 127, 255}));
  connect(Heater_Lift.fluid_b, Heater.fluid_a) annotation(
    Line(points = {{-4, -14}, {4, -14}, {4, -6}}, color = {0, 127, 255}));

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)), experiment(StopTime = 3.24e+07, StartTime = 0, Tolerance = 1.0e-5, Interval = 300, maxStepSize = 60, initialStepSize = 60));
end H2DRI_DesignCase_2b_Dynamic;