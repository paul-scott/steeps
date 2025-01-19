within examples;
//SM=3.0, 25 ideal storage hours, Power block 720-500 degC. 759153 MWhe per year
model PBS_Surround_SCO2NREL_Table_Series
    //Note here we used the optimised isp value
  function opt_file_naming
    input String prefix;
    //"modelica://SolarTherm/Data/Optics/SodiumBoiler/surround/Ref/"
    input String phi_pct_string;
    input Real SM_guess;
    input Real HT_pct_guess;
    input Real f_recv_guess;
    output String opt_file;
  protected
    Integer phi;
    Integer SM;
    Integer HT_pct;
    Integer f_recv;
    String SM_string;
    String HT_pct_string;
    String f_recv_string;
  algorithm
    SM := max(14, min(38, 1 * round(SM_guess * 10)));
//Actually SM*10"
    HT_pct := max(70, min(130, 5 * round(HT_pct_guess * 0.2)));
    f_recv := max(70, min(130, 5 * round(f_recv_guess * 20.0)));
    SM_string := String(SM);
    HT_pct_string := String(HT_pct);
    f_recv_string := String(f_recv);
    opt_file := Modelica.Utilities.Files.loadResource(prefix + SM_string + "dSM/isp_optimised/" + phi_pct_string + "%phi_" + HT_pct_string + "%HT_" + f_recv_string + "%Arecv_optics.motab");
  end opt_file_naming;

  function stor_file_naming
    input String prefix;
    //"modelica://SolarTherm/Data/Storage/Storage_Table_Series42_252_"
    input Real t_storage;
    output String stor_file;
  protected
    String t_stor_string;
  algorithm
    t_stor_string := String(round(t_storage*100.0));
    stor_file := Modelica.Utilities.Files.loadResource(prefix + t_stor_string + "ch.motab");
  end stor_file_naming;

  function round
    input Real number;
    output Integer int;
  protected
    Integer quotient;
    Real remainder;
  algorithm
    quotient := integer(number);
    remainder := number - floor(number);
    if remainder >= 0.5 then
      int := 1 + quotient;
    else
      int := quotient;
    end if;
  end round;

  import SolarTherm.{Models,Media};
  import Modelica.SIunits.Conversions.from_degC;
  import SI = Modelica.SIunits;
  import nSI = Modelica.SIunits.Conversions.NonSIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import FI = SolarTherm.Models.Analysis.Finances;
  import SolarTherm.Types.Solar_angles;
  import SolarTherm.Types.Currency;
  extends Modelica.Icons.Example;
  // Input Parameters
  // *********************
  //To be optimised
  //--Tank
  //Concept Parameters
  //Tolerance of PB and receiver inlet
  parameter SI.TemperatureDifference T_tol_recv = 40.0 "Temperature tolerance above design receiver input temperature before receiver is shut off";
  parameter SI.TemperatureDifference T_tol_PB = 20.0 "Temperature tolerance below design PB input temperature before PB is shut off";
  //Temperature controls
  parameter SI.Temperature T_max = 720.0 + 273.15 "Ideal high temperature of the storage";
  parameter SI.Temperature T_PB_start = T_max - T_tol_PB + 10.0 "Temperature at top of tank where PB can start";
  parameter SI.Temperature T_PB_min = T_max - T_tol_PB "Temperature at top of tank where PB must stop";
  parameter SI.Temperature T_recv_max = T_min + T_tol_recv "Maximum temperature at bottom of tank when it can no longer be pumped into the receiver";
  parameter SI.Temperature T_recv_start = T_min + T_tol_recv - 10.0 "Temperature at bottom of tank when it can start being pumped into the receiver again";
  parameter SI.Temperature T_min = 500.0 + 273.15 "Ideal low temperature of the storage";
  //Additional Power Block design temperature
  parameter SI.Temperature T_PB_in_des = 720.0 + 273.15 "Power Block design inlet temperature";
  parameter SI.Temperature T_PB_out_des = 500.0 + 273.15 "Power Block design outlet temperature";
  //Controller Parameters
  parameter SI.Time t_stor_startPB = 1.0 * 3600.0 "minimum hours of storage available to startup PB";
  //Change this to sweep
  parameter Real eff_storage_des = Tank.util;// 0.8317; //Override and set it to 0.55 Tank.util "design storage utilisation";
  parameter String storage_file = stor_file_naming("modelica://SolarTherm/Data/Storage/Storage_Table_Series42_252_", t_storage);
  parameter Real t_storage(unit = "h") = 25.0 "Hours of storage";
  //parameter Real L_recv_max = 0.94599;
  //L_4
  //parameter Real L_recv_start = 0.28678;
  //L_3
  //parameter Real L_PB_start = 0.067738;
  //L_2
  //parameter Real L_PB_min = 0.062044;
  //L_1
  //Material and Media Packages
  replaceable package Medium = SolarTherm.Media.Sodium.Sodium_pT "Medium props for molten salt";
  replaceable package Fluid = SolarTherm.Materials.Sodium_Table "Material model for Sodium Chloride PCM";
  replaceable package Filler = SolarTherm.Materials.MgO_Constant "Tank filler";
  //Storage Design
  parameter SI.Energy E_max = Q_flow_ref_blk * t_storage * 3600.0 "Theoretical max capacity of storage";
  //parameter Integer N_f = 50 "Number of discretizations in vertical fluid phase";
  //parameter Integer N_p = 10 "Number of discretizations in radial filler phase";
  //parameter SI.Length d_p = 0.10 "Tank filler diameter";
  //parameter Real eta = 0.26 "Packed bed void fraction (porosity)";
  //parameter Real ar = 2.0 "Aspect ratio (H/D) of tank";
  
  //parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.00 "Heat loss coefficient of all tanks";
  //Weather Data (Alice Springs)
  parameter String wea_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Weather/example_TMY3.motab");
  parameter Real wdelay[8] = {0, 0, 0, 0, 0, 0, 0, 0} "Weather file delays";
  parameter nSI.Angle_deg lon = 133.889 "Longitude (+ve East)";
  parameter nSI.Angle_deg lat = -23.795 "Latitude (+ve North)";
  parameter nSI.Time_hour t_zone = 9.5 "Local time zone (UCT=0)";
  parameter Integer year = 1996 "Meteorological year";
  //Optical Model Parameters
  parameter Solar_angles angles = Solar_angles.dec_hra "Angles used in the lookup table file";
  parameter String field_type = "surround";
  parameter String opt_file_prefix = "modelica://SolarTherm/Data/Optics/surround/50MWe/4555c%/883K/1000kWpm2/";
  parameter String phi_pct_string = "100";
  parameter Real SM_guess = 3.0;
  parameter Real HT_pct_guess = 100;
  parameter Real f_recv = 1.0;
  parameter String opt_file = opt_file_naming(opt_file_prefix, phi_pct_string, SM_guess, HT_pct_guess, f_recv);
  parameter String SM_string = String(2 * SolarTherm.Utilities.Round(SM_guess * 5)) "Solar Multiple rounded to the nearest 0.2, multiplied by 10 and converted to string";
  parameter Real he_av_design = 1.00 "Helisotats availability";
  parameter SI.Energy helio_E_start = 90e3 * A_heliostat / 144.375 "Heliostat startup energy consumption";
  parameter SI.Power helio_W_track = 0.0553 * A_heliostat / 144.375 "Heliostat tracking power";
  parameter SI.Angle ele_min = 0.13962634015955 "Heliostat stow deploy angle";
  parameter Boolean use_wind = true "true if using wind stopping strategy in the solar field";
  parameter SI.Velocity Wspd_max = 15 if use_wind "Wind stow speed";
  parameter Real nu_start = 0.4 "Minimum energy start-up fraction to start the receiver";
  parameter Real nu_min_sf = 0.3 "Minimum turn-down energy fraction to stop the receiver";
  parameter Real nu_defocus = 1 "Energy fraction to the receiver at defocus state";
  parameter Real[8] MetaA = SolarTherm.Utilities.Metadata_Optics(opt_file);
  parameter Integer n_heliostat = SolarTherm.Utilities.Round(MetaA[1]) * 2 "Number of heliostats";
  //Two-tower system
  parameter SI.Area A_heliostat = MetaA[2] "Area of one heliostat";
  parameter Real eff_opt_des = MetaA[3];
  parameter SI.Length H_recv = MetaA[4];
  parameter SI.Length D_recv = MetaA[5];
  parameter SI.Length H_tower = MetaA[6] "Height of the tower";
  parameter SI.Area A_field = A_heliostat * n_heliostat "Area of the entire field (reflective area)";
  parameter SI.Area A_land = land_mult * A_field "Land area occupied by the plant";
  //Receiver Parameters
  parameter SI.Area A_recv = if field_type == "polar" then 2.0 * H_recv * D_recv else 2.0 * H_recv * D_recv * CN.pi "Total Receiver area";
  //Two-tower system
  //parameter Integer N_pa_recv = 20 "Number of panels in receiver";
  //parameter SI.Thickness t_tb_recv = 1.25e-3 "Receiver tube wall thickness";
  //parameter SI.Diameter D_tb_recv = 40e-3 "Receiver tube outer diameter";
  parameter SI.Efficiency ab_recv = 0.9872 "Receiver coating absorptance";
  parameter SI.Efficiency em_recv = 0.9408 "Receiver coating emissivity";
  parameter SI.CoefficientOfHeatTransfer h_conv_recv = 20.0 "W/m2K";
  parameter SI.Temperature T_4avg_recv = (0.20 * (T_min ^ 4 + T_min ^ 3 * T_max + T_min ^ 2 * T_max ^ 2 + T_min * T_max ^ 3 + T_max ^ 4)) ^ 0.25;
  parameter SI.Temperature T_avg_recv = 0.50 * (T_min + T_max);
  parameter SI.Temperature T_amb_recv_des = 25.0 + 273.15 "Design ambient temperature (K)";
  // Power block
  parameter SI.Power P_gross_des = P_name / eff_net_des "Power block gross rating at design point";
  parameter SI.Power P_name = 100e6;
  parameter SI.Temperature T_pb_cool_des = 35.0 + 273.15 "Design cooling temperature ambient of PB, is 6 degrees below 41 degC";
  parameter SI.Efficiency eff_net_des = 0.90 "Power block net efficiency rating";
  parameter SI.Efficiency eff_blk_des = 0.4555 "Power block efficiency at design point";
  parameter SI.Time t_PB_wait = 1.0 * 3600.0 "Wait time between shutdown and turning back on";
  //Pumping and Parasitics
  parameter SI.SpecificEnergy k_loss_cold = 0.15e3 "Cold pump parasitic power coefficient";
  parameter SI.SpecificEnergy k_loss_hot = 0.55e3 "Hot pump parasitic power coefficient";
  //Design Calculations
  parameter SI.Irradiance dni_des = SolarTherm.Utilities.DNI_Models.Meinel(abs(lat)) "Design point DNI value based on meinel calculation";
  parameter Real SM = Q_flow_rec_des / Q_flow_ref_blk "Real solar multiple";
  //Enthalpies
  parameter SI.SpecificEnthalpy h_in_ref_blk = Medium.specificEnthalpy(Medium.setState_pTX(101323.0, T_PB_in_des)) "Specific enthalpy of sodium entering PB at design pt";
  parameter SI.SpecificEnthalpy h_out_ref_blk = Medium.specificEnthalpy(Medium.setState_pTX(101323.0, T_PB_out_des)) "Specific enthalpy of sodium leaving PB at design pt";
  parameter SI.SpecificEnthalpy h_in_ref_recv = Medium.specificEnthalpy(Medium.setState_pTX(101323.0, T_min)) "Specific enthalpy of sodium entering receiver at design pt";
  parameter SI.SpecificEnthalpy h_out_ref_recv = Medium.specificEnthalpy(Medium.setState_pTX(101323.0, T_max)) "Specific enthalpy of sodium leaving receiver at design pt";
  //Heat Flow Rates
  parameter SI.HeatFlowRate Q_flow_ref_blk = P_gross_des / eff_blk_des "design heat input rate into the PB";
  parameter SI.HeatFlowRate Q_flow_rec_loss_des = CN.sigma * em_recv * A_recv * (T_4avg_recv ^ 4 - T_amb_recv_des ^ 4) + h_conv_recv * A_recv * (T_avg_recv - T_amb_recv_des) "Receiver design heat loss rate";
  parameter SI.HeatFlowRate Q_flow_rec_des = dni_des * he_av_design * eff_opt_des * A_field * ab_recv - Q_flow_rec_loss_des "Receiver Thermal power output at design";
  parameter SI.HeatFlowRate Q_flow_defocus = (Q_flow_ref_blk + Q_flow_rec_loss_des) / ab_recv "Solar field thermal power at defocused state (not used)";
  //Mass flow rates
  parameter SI.MassFlowRate m_flow_blk_des = Q_flow_ref_blk / (h_in_ref_blk - h_out_ref_blk) "Design point mass flow rate of sodium vapor condensing into the power block";
  parameter SI.MassFlowRate m_flow_recv_des = Q_flow_rec_des / (h_out_ref_recv - h_in_ref_recv) "Design mass flow rate into recv";
  //Cost & Financial Parameters
  parameter String pri_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Prices/aemo_vic_2014.motab") "Electricity price file";
  parameter Currency currency = Currency.USD "Currency used for cost analysis";
  parameter Real land_mult = 5.0 "Land area multiplier";
  // Cost data in USD (default) or AUD
  parameter Real r_disc = 0.044 "Real discount rate";
  parameter Real r_i = 0.00 "Inflation rate";
  parameter Integer t_life = 27 "Lifetime of plant";
  parameter Integer t_cons = 3 "Years of construction";
  parameter Real r_cur = 0.71 "The currency rate from AUD to USD";
  // Valid for 2019. See https://www.rba.gov.au/
  parameter Real f_Subs = 0 "Subsidies on initial investment costs";
  parameter FI.AreaPrice pri_field = if currency == Currency.USD then 75.00 else 75.00 / r_cur "Field cost per design aperture area";
  // SAM 2018 cost data: 177*(603.1/525.4) in USD. Note that (603.1/525.4) is CEPCI index from 2007 to 2018
  parameter FI.AreaPrice pri_site = if currency == Currency.USD then 16.00 else 16.00 / r_cur "Site improvements cost per area";
  // SAM 2018 cost data: 16
  //parameter FI.EnergyPrice pri_storage = if currency == Currency.USD then 37 / (1e3 * 3600) else 37 / (1e3 * 3600) / r_cur "Storage cost per energy capacity";
  // SAM 2018 cost data: 22 / (1e3 * 3600)
  //parameter FI.PowerPrice pri_block = if currency == Currency.USD then 1040.00 / 1e3 else 1040.00 / 1e3 / r_cur "Power block cost per gross rated power";
  parameter FI.PowerPrice pri_block = powerBlock.C_PB_total / P_gross_des;
  // SAM 2018 cost data: 1040
  parameter FI.PowerPrice pri_bop = if currency == Currency.USD then 0.29 else 0.29 "Balance of plant cost per gross rated power";
  //not used
  //SAM 2018 cost data: 290
  parameter FI.AreaPrice pri_land = if currency == Currency.USD then 2.471 else 2.471 "Land cost per area";
  parameter Real pri_om_name(unit = "$/W/year") = if currency == Currency.USD then 40.0 / 1e3 else 40.0 / 1e3 / r_cur "Fixed O&M cost per nameplate per year";
  //SAM 2018 cost data: 66
  parameter Real pri_om_prod(unit = "$/J/year") = if currency == Currency.USD then 3.0 / (1e6 * 3600) else 3.0 / (1e6 * 3600) / r_cur "Variable O&M cost per production per year";
  //SAM 2018 cost data: 3.5
  parameter FI.Money C_field = pri_field * A_field "Field cost";
  parameter FI.Money C_site = pri_site * A_field "Site improvements cost";
  parameter FI.Money C_tower = 3000000.0 * exp(0.0113 * H_tower) * 2.0 "Tower cost";
  //two-tower system
  parameter FI.Money C_receiver = (4780420.0 + 21480230.0 * (D_recv / 14.0) * (H_recv / 14.5) ^ 0.6) * 2.0 "Receiver cost";
  //two-tower system
  // SAM 2018 cost data: 103e6 * (A_receiver / 1571) ^ 0.7
  parameter FI.Money C_storage = Tank.C_total;
  //tankHot.C_Storage "Storage cost";
  //parameter FI.Money C_block = pri_block * P_gross_des "Power block cost";
  parameter FI.Money C_block = powerBlock.C_PB_total "Power block cost";
  parameter FI.Money C_bop = 16216995.0 * (P_name / 100e6) + 11919765.0 * (2.0 * H_tower / 300.0) "Balance of plant cost";
  parameter FI.Money C_cap_dir_sub = (1 - f_Subs) * (C_field + C_site + C_tower + C_receiver + C_storage + C_block + C_bop) "Direct capital cost subtotal";
  // i.e. purchased equipment costs
  parameter FI.Money C_contingency = 0.10 * C_cap_dir_sub "Contingency costs";
  parameter FI.Money C_cap_dir_tot = C_cap_dir_sub + C_contingency "Direct capital cost total";
  parameter FI.Money C_EPC = 0.09 * C_cap_dir_tot "Engineering, procurement and construction(EPC) and owner costs";
  // SAM 2018 cost data: 0.13
  parameter FI.Money C_land = pri_land * A_land "Land cost";
  parameter FI.Money C_cap = C_cap_dir_tot + C_EPC + C_land "Total capital (installed) cost";
  parameter FI.MoneyPerYear C_year = pri_om_name * P_name "Fixed O&M cost per year";
  parameter Real C_prod(unit = "$/J/year") = pri_om_prod "Variable O&M cost per production per year";
  // System components
  // *********************
  //Weather data
  SolarTherm.Models.Sources.DataTable.DataTable data(lon = lon, lat = lat, t_zone = t_zone, year = year, file = wea_file) annotation(
    Placement(visible = true, transformation(extent = {{-120, 82}, {-90, 110}}, rotation = 0)));
  //DNI_input
  Modelica.Blocks.Sources.RealExpression DNI_input(y = data.DNI) annotation(
    Placement(visible = true, transformation(extent = {{-114, 60}, {-94, 80}}, rotation = 0)));
  //Tamb_input
  Modelica.Blocks.Sources.RealExpression Tamb_input(y = data.Tdry) annotation(
    Placement(visible = true, transformation(extent = {{122, 70}, {102, 90}}, rotation = 0)));
  //WindSpeed_input
  Modelica.Blocks.Sources.RealExpression Wspd_input(y = data.Wspd) annotation(
    Placement(visible = true, transformation(extent = {{-118, 38}, {-92, 58}}, rotation = 0)));
  //pressure_input
  Modelica.Blocks.Sources.RealExpression Pres_input(y = data.Pres) annotation(
    Placement(visible = true, transformation(extent = {{120, 86}, {100, 106}}, rotation = 0)));
  //parasitic inputs
  Modelica.Blocks.Sources.RealExpression parasities_input(y = heliostatsField.W_loss + pumpHot.W_loss + pumpCold.W_loss) annotation(
    Placement(visible = true, transformation(origin = {121, 64}, extent = {{-13, -10}, {13, 10}}, rotation = 180)));
  // Or block for defocusing
  //Sun
  SolarTherm.Models.Sources.SolarModel.Sun sun(lon = data.lon, lat = data.lat, t_zone = data.t_zone, year = data.year, redeclare function solarPosition = Models.Sources.SolarFunctions.PSA_Algorithm) annotation(
    Placement(transformation(extent = {{-82, 60}, {-62, 80}})));
  // Solar heliostat field
  SolarTherm.Models.CSP.CRS.HeliostatsField.HeliostatsField heliostatsField(redeclare model Optical = Models.CSP.CRS.HeliostatsField.Optical.Table_Full(angles = angles, file = opt_file), A_h = A_heliostat, Q_design = Q_flow_defocus, Wspd_max = Wspd_max, ele_min(displayUnit = "deg") = ele_min, he_av = he_av_design, lat = data.lat, lon = data.lon, n_h = n_heliostat, nu_defocus = nu_defocus, nu_min = nu_min_sf, nu_start = nu_start, use_defocus = false, use_on = true, use_wind = true) annotation(
    Placement(transformation(extent = {{-88, 2}, {-56, 36}})));
  // Hot Pump (power block)
  // Cold pump (receiver)
  //Market
  SolarTherm.Models.Analysis.Market market(redeclare model Price = Models.Analysis.EnergyPrice.Constant) annotation(
    Placement(visible = true, transformation(origin = {144, 20}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  //Receiver
  SolarTherm.Models.CSP.CRS.Receivers.PBS_Receiver receiver(redeclare package Medium = Medium, H_rcv = H_recv, D_rcv = D_recv, A_recv = A_recv, ab = ab_recv, em = em_recv, T_0 = T_min, Q_des_blk = Q_flow_ref_blk, T_max = T_max) annotation(
    Placement(visible = true, transformation(origin = {-30, 24}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  //Storage
  //Loop Breakers
  //Cold Controller (Receiver)
  //Hot Controller (Power Block)
  //Power Block
  SolarTherm.Models.PowerBlocks.PBS_PowerBlockModel_sCO2NREL_100MWe_720C_500C powerBlock(redeclare package Medium = Medium, redeclare model Cooling = SolarTherm.Models.PowerBlocks.Cooling.NoCooling, capped = false, Q_flow_ref = Q_flow_ref_blk, T_out_ref = T_PB_out_des, W_base = 0.0055 * P_gross_des, m_flow_ref = m_flow_blk_des, nu_net = eff_net_des) annotation(
    Placement(visible = true, transformation(origin = {101, 21}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
  //Annual Simulation variables
  SI.Power P_elec "Output power of power block";
  SI.Energy E_elec(start = 0, fixed = true, displayUnit = "MW.h") "Generate electricity";
  FI.Money R_spot(start = 0, fixed = true) "Spot market revenue";
  //Boolean constrained(start = false);
  //Real distance(start = 0);
  /*
                                    //Analytics
                                    //Accumulated energy
                                    SI.Energy E_resource(start = 0) "Integral of DNI with time if greater than zero";
                                    SI.Energy E_helio_incident(start = 0) "Cumulative heat energy incident on heliostats after curtailment (low-DNI/high-wind)";
                                    SI.Energy E_helio_raw(start = 0) "Cumulative heat energy delivered by field to receiver after he_av losses + optical losses";
                                    SI.Energy E_helio_net(start = 0) "Cumulative heat energy delivered by field to receiver after defocusing losses";
                                    SI.Energy E_recv_absorbed(start = 0) "Cumulative heat energy absorbed by the receiver before re-emission and convection";
                                    SI.Energy E_recv_output(start = 0) "Cumulative heat energy outputted by the receiver after thermal losses";
                                    SI.Energy E_PB_input(start = 0) "Cumulative heat energy inputted into the power block";
                                    SI.Energy E_PB_gross(start = 0) "Cumulative gross electrical energy produced by the power block";
                                    SI.Energy E_PB_net(start = 0) "Cumulative electrical output of the power block after parasitics and generator losses";
                                    Real sum_shading(start = 0) "Shading efficiency multiplied by time when heliostats are on";
                                    Real sum_cosine(start = 0) "Shading efficiency multiplied by time when heliostats are on";
                                    Real sum_reflection(start = 0) "Shading efficiency multiplied by time when heliostats are on";
                                    Real sum_blocking(start = 0) "Shading efficiency multiplied by time when heliostats are on";
                                    Real sum_attenuation(start = 0) "Shading efficiency multiplied by time when heliostats are on";
                                    Real sum_intercept(start = 0) "Shading efficiency multiplied by time when heliostats are on";
                                    Real sum_timehelio(start = 0) "Sum of time when heliostat is on";
                                    Real eta_shading;
                                    Real eta_cosine;
                                    Real eta_reflection;
                                    Real eta_blocking;
                                    Real eta_attenuation;
                                    Real eta_intercept;
                                    //Annual efficiencies
                                    Real eta_curtail_off "Curtailment: Heliostat off";
                                    Real eta_he_av "Heliostat Availability";
                                    Real eta_optical "Field optical efficiency including spillage";
                                    Real eta_curtail_defocus "Curtailment: Full Storage";
                                    Real eta_recv_abs "Receiver Absorptivity";
                                    Real eta_recv_thermal "Receiver thermal efficiency";
                                    Real eta_storage "Storage thermal efficiency";
                                    Real eta_pb_gross "Power block gross efficiency";
                                    Real eta_pb_net "Power block net efficiency";
                                    Real eta_solartoelec "Solar to electric";
                                  */
  //Storage Uitlization
  //Real E_max_today(start=0.0) "Today's max energy stored";
  //Real E_min_today(start=0.0) "Today's min energy stored";
  //Real eta_util_ytd(start=0.0) "yesterday's total utilization pct of storage";
  //End Analytics
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction_LoopBreaker Splitter_bot(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {26, 13}, extent = {{-10, -9}, {10, 9}}, rotation = 180)));
  SolarTherm.Models.Fluid.Valves.PBS_TeeJunction Splitter_top(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {26, 59}, extent = {{-10, -9}, {10, 9}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpCold(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {-2, 6}, extent = {{4, -4}, {-4, 4}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.PumpSimple_EqualPressure pumpHot(redeclare package Medium = Medium) annotation(
    Placement(visible = true, transformation(origin = {66, 68}, extent = {{-4, -4}, {4, 4}}, rotation = 0)));
  SolarTherm.Models.Control.Thermocline_Controller_PBLimit_Level pBS_Controller_PBLimit_Level(redeclare package HTF = Medium, T_target = T_max, m_flow_PB_des = m_flow_blk_des, Q_des_blk = Q_flow_ref_blk, L_1 = Tank.L_PB_min, L_2 = Tank.L_PB_start, L_3 = Tank.L_recv_start, L_4 = Tank.L_recv_max, t_wait = t_PB_wait, eff_storage_des = Tank.util, t_stor_startPB = t_stor_startPB, E_max = Tank.E_max) annotation(
    Placement(visible = true, transformation(origin = {38, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Storage.Thermocline.Thermocline_Table_Final Tank(redeclare package Medium = Medium, redeclare package Fluid_Package = Fluid, table_file = storage_file, L_start = 0.0) annotation(
    Placement(visible = true, transformation(origin = {26, 36}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
algorithm
/*when rem(time,86400) > 86399 then //reset the storage utilization
    if time > 432000 then //after 5 days
      //eta_util_ytd := (E_max_today - E_min_today)/E_max;
      E_max_today := Tank.Tank_A.E_stored;
      E_min_today := Tank.Tank_A.E_stored;
    end if;
  end when;
  when rem(time,600) > 599 then //update the Max and Min energy stored for the day if appropriate;
  if time > 432000 then //after 5 days
    if Tank.Tank_A.E_stored > E_max_today then
      E_max_today := Tank.Tank_A.E_stored;
    end if;
    if Tank.Tank_A.E_stored < E_min_today then
      E_min_today := Tank.Tank_A.E_stored;
    end if;
  end if;
  end when;
  */
/*if time > 60.0 then
      if tankHot.m_avail < 1.0e-6 then
      constrained := true;
      distance := 1000 + 1000 * ((31536000.0 - time) / 31536000.0);
      terminate("Receiver Dries Out");
    elseif tankHot.L_total > z_tray then
      constrained := true;
      distance := 1000 + 1000 * ((31536000.0 - time) / 31536000.0);
      terminate("Liquid PCM spills out the tray");
    elseif tankHot.T_HTF < 600.0 then
      constrained := true;
      distance := 1000 + 1000 * ((31536000.0 - time) / 31536000.0);
      terminate("Tank gets too cold to continue");
    elseif tankHot.x_HTF < 0.0 or tankHot.x_HTF > 1.0 then
      constrained := true;
      distance := 1000 + 1000 * ((31536000.0 - time) / 31536000.0);
      terminate("Invalid Quality");
    end if;
  end if;
*/
/*
//Annual Efficieny
  if time > 31449600 then
    eta_curtail_off := E_helio_incident / E_resource;
    eta_optical := E_helio_raw / (E_helio_incident * he_av_design);
    eta_he_av := he_av_design;
    eta_curtail_defocus := E_helio_net / E_helio_raw;
    eta_recv_abs := E_recv_absorbed / E_helio_net;
    eta_recv_thermal := E_recv_output / E_recv_absorbed;
    eta_storage := E_PB_input / E_recv_output;
    eta_pb_gross := E_PB_gross / E_PB_input;
    eta_pb_net := E_PB_net / E_PB_gross;
    eta_solartoelec := E_PB_net / E_resource;
    eta_shading := sum_shading / sum_timehelio;
    eta_cosine := sum_cosine / sum_timehelio;
    eta_reflection := sum_reflection / sum_timehelio;
    eta_blocking := sum_blocking / sum_timehelio;
    eta_attenuation := sum_attenuation / sum_timehelio;
    eta_intercept := sum_intercept / sum_timehelio;
  end if;
*/
//Optics
equation
/*
//Analytics
//Cumulative heat
  der(E_resource) = max(sun.dni * A_field, 0.0);
  der(E_helio_incident) = if heliostatsField.on_internal then heliostatsField.n_h * heliostatsField.A_h * max(0.0, heliostatsField.solar.dni) else 0.0;
  der(E_helio_raw) = heliostatsField.Q_raw;
  der(E_helio_net) = heliostatsField.Q_net;
  der(E_recv_absorbed) = receiver.ab * receiver.heat.Q_flow;
  der(E_recv_output) = receiver.Q_rcv;
  der(E_PB_input) = powerBlock.Q_flow;
  der(E_PB_gross) = powerBlock.W_gross;
  der(E_PB_net) = powerBlock.W_net;
  if heliostatsField.on_internal then
//cumulative optics efficiency
    der(sum_shading) = max(sun.dni * A_field, 0.0) * heliostatsField.optical.nu_shading;
    der(sum_cosine) = max(sun.dni * A_field, 0.0) * heliostatsField.optical.nu_cosine;
    der(sum_reflection) = max(sun.dni * A_field, 0.0) * heliostatsField.optical.nu_reflection;
    der(sum_blocking) = max(sun.dni * A_field, 0.0) * heliostatsField.optical.nu_blocking;
    der(sum_attenuation) = max(sun.dni * A_field, 0.0) * heliostatsField.optical.nu_attenuation;
    der(sum_intercept) = max(sun.dni * A_field, 0.0) * heliostatsField.optical.nu_intercept;
    der(sum_timehelio) = max(sun.dni * A_field, 0.0);
  else
    der(sum_shading) = 0.0;
    der(sum_cosine) = 0.0;
    der(sum_reflection) = 0.0;
    der(sum_blocking) = 0.0;
    der(sum_attenuation) = 0.0;
    der(sum_intercept) = 0.0;
    der(sum_timehelio) = 0.0;
  end if;
  */
//Connections from data
  connect(DNI_input.y, sun.dni) annotation(
    Line(points = {{-93, 70}, {-83, 70}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Wspd_input.y, heliostatsField.Wspd) annotation(
    Line(points = {{-91, 48}, {-100.35, 48}, {-100.35, 30}, {-88, 30}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
// Fluid connections
// controlCold connections
// controlHot connections
//Solar field connections i.e. solar.heat port and control
  connect(sun.solar, heliostatsField.solar) annotation(
    Line(points = {{-72, 60}, {-72, 36}}, color = {255, 128, 0}));
//PowerBlock connections
  P_elec = powerBlock.W_net;
  E_elec = powerBlock.E_net;
  R_spot = market.profit;
  connect(heliostatsField.heat, receiver.heat) annotation(
    Line(points = {{-56, 27.5}, {-46, 27.5}, {-46, 29}}, color = {191, 0, 0}));
  connect(parasities_input.y, powerBlock.parasities) annotation(
    Line(points = {{107, 64}, {107, 38}}, color = {0, 0, 127}));
  connect(powerBlock.T_amb, Tamb_input.y) annotation(
    Line(points = {{95, 38}, {95, 80}, {101, 80}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(powerBlock.W_net, market.W_net) annotation(
    Line(points = {{116, 20}, {132, 20}}, color = {0, 0, 127}));
  connect(Tamb_input.y, receiver.Tamb) annotation(
    Line(points = {{101, 80}, {-30, 80}, {-30, 36}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
//connect(Wspd_input.y, receiver.Wspd) annotation(
//  Line(points = {{-113, 48}, {-32, 48}, {-32, 36}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Splitter_bot.fluid_b, pumpCold.fluid_a) annotation(
    Line(points = {{18, 6}, {2, 6}}, color = {0, 127, 255}));
  connect(pumpCold.fluid_b, receiver.fluid_a) annotation(
    Line(points = {{-6, 6}, {-14, 6}, {-14, 10}, {-27, 10}}, color = {0, 127, 255}));
  connect(receiver.fluid_b, Splitter_top.fluid_a) annotation(
    Line(points = {{-25, 32}, {-18, 32}, {-18, 68}, {18, 68}}, color = {0, 127, 255}));
  connect(pumpHot.fluid_b, powerBlock.fluid_a) annotation(
    Line(points = {{70, 68}, {82, 68}, {82, 30}, {88, 30}}, color = {0, 127, 255}));
  connect(Splitter_top.fluid_b, pumpHot.fluid_a) annotation(
    Line(points = {{34, 68}, {62, 68}}, color = {0, 127, 255}));
  connect(Splitter_top.fluid_c, Tank.fluid_a) annotation(
    Line(points = {{26, 54}, {26, 46}}, color = {0, 127, 255}));
  connect(Tank.fluid_b, Splitter_bot.fluid_c) annotation(
    Line(points = {{26, 26}, {26, 14}}, color = {0, 127, 255}));
  connect(Pres_input.y, Tank.p_amb) annotation(
    Line(points = {{100, 96}, {48, 96}, {48, 36}, {32, 36}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Tank.h_bot_outlet, pBS_Controller_PBLimit_Level.h_tank_outlet) annotation(
    Line(points = {{22, 28}, {22, 20}, {34, 20}, {34, -11}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(receiver.Q_rcv_raw, pBS_Controller_PBLimit_Level.Q_rcv_raw) annotation(
    Line(points = {{-26, 20}, {14, 20}, {14, -10}, {27, -10}, {27, -15}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Tank.Level, pBS_Controller_PBLimit_Level.Level) annotation(
    Line(points = {{20, 41}, {10, 41}, {10, -21}, {27, -21}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(pBS_Controller_PBLimit_Level.Q_defocus, receiver.Q_defocus) annotation(
    Line(points = {{27, -31}, {-12, -31}, {-12, 16}, {-26, 16}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(pBS_Controller_PBLimit_Level.defocus, receiver.defocus) annotation(
    Line(points = {{49, -26}, {86, -26}, {86, -36}, {-44, -36}, {-44, 20}, {-36, 20}}, color = {255, 0, 255}, pattern = LinePattern.Dash));
  connect(pBS_Controller_PBLimit_Level.m_flow_recv, pumpCold.m_flow) annotation(
    Line(points = {{49, -16}, {92, -16}, {92, -44}, {-8, -44}, {-8, 16}, {-2, 16}, {-2, 10}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(pBS_Controller_PBLimit_Level.m_flow_PB, pumpHot.m_flow) annotation(
    Line(points = {{49, -21}, {76, -21}, {76, 78}, {66, 78}, {66, 72}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(powerBlock.fluid_b, Splitter_bot.fluid_a) annotation(
    Line(points = {{84, 8}, {62, 8}, {62, 6}, {34, 6}, {34, 6}}, color = {0, 127, 255}));
  connect(pBS_Controller_PBLimit_Level.h_PB_outlet, powerBlock.h_out_signal) annotation(
    Line(points = {{44, -10}, {66, -10}, {66, 4}, {84, 4}, {84, 4}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
  connect(Tamb_input.y, Tank.T_amb) annotation(
    Line(points = {{101, 80}, {14, 80}, {14, 36}, {20, 36}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
protected
  annotation(
    Diagram(coordinateSystem(extent = {{-140, -120}, {160, 140}}, initialScale = 0.1), graphics = {Text(origin = {0, 6}, extent = {{-52, 8}, {-4, -12}}, textString = "Receiver", fontSize = 12, fontName = "CMU Serif"), Text(origin = {12, 2}, extent = {{-110, 4}, {-62, -16}}, textString = "Heliostats Field", fontSize = 12, fontName = "CMU Serif"), Text(origin = {-16, 10}, extent = {{-80, 86}, {-32, 66}}, textString = "Sun", fontSize = 12, fontName = "CMU Serif"), Text(origin = {0, -12}, extent = {{80, 12}, {128, -8}}, textString = "Power Block", fontSize = 12, fontName = "CMU Serif"), Text(origin = {8, -2}, extent = {{112, 16}, {160, -4}}, textString = "Market", fontSize = 12, fontName = "CMU Serif"), Text(origin = {-66, -36}, extent = {{80, 12}, {128, -8}}, textString = "Controller", fontSize = 12, fontName = "CMU Serif")}),
    Icon(coordinateSystem(extent = {{-140, -120}, {160, 140}})),
    experiment(StopTime = 3.1536e+07, StartTime = 0, Tolerance = 0.001, Interval = 300, maxStepSize = 60, initialStepSize = 60),
    __Dymola_experimentSetupOutput,
    Documentation(revisions = "<html>
	<ul>
	<li>Zebedee Kee (June 2020) :<br>Released first version. </li>
	</ul>

	</html>"));
end PBS_Surround_SCO2NREL_Table_Series;
