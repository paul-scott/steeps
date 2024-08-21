within examples;

model Hybrid_CSP_PV_Particle_HotH2_Burner_OptDispatch2
  //*********************** As per December 7 2020, the tower cost function is changed to the Latest Tower Cost Function
  //*********************** based on the email by J.Sment (Sandia) Sat 05/12/2020 05:48
  import SolarTherm.{Models,Media};
  import Modelica.SIunits.Conversions.from_degC;
  import SI = Modelica.SIunits;
  import nSI = Modelica.SIunits.Conversions.NonSIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import FI = SolarTherm.Models.Analysis.Finances;
  import Utils = SolarTherm.Media.SolidParticles.CarboHSP_utilities;
  import SolarTherm.Types.Solar_angles;
  import SolarTherm.Types.Currency;
  import Modelica.Math;
  import SolarTherm.Utilities.Kriging_OnTheFly.*;
  import SolarTherm.Utilities.ANN_OnTheFly.*;
  import metadata = SolarTherm.Utilities.Metadata_Optics_windy;
  import MedGas = Modelica.Media.IdealGases.Common.Functions;
  import MedSteam = Modelica.Media.Water.StandardWater;
  extends SolarTherm.Media.CO2.PropCO2;
  extends Modelica.Icons.Example;
  parameter Real PV_fraction = 1.49937772041 "Fraction of the hybrid system that is PV nameplate";
  parameter SI.Power P_hybrid_system = 1968485792.43 "Hybrid system nameplate [W] --> allow to be oversized";
  parameter Real t_storage_H2(unit = "h") = 157.1263624096 "[ST] H2 Storage capacity";
  parameter Real t_storage(unit = "h") = 582.182264329 "[ST] Storage capacity";
  parameter Real PB_fraction = 0;
  parameter Real upstreammultiplier_change = 1 "The multiplier controlled by the optimisation, and passed";
  parameter Boolean set_no_thermal = false "If true, then turn off PB, electrical heater, indirect particle gas HX";
  parameter Boolean set_dispatch_optimiser_dual_tank = false;
  parameter Real CSP_fraction = 1 "Fraction of the hybrid system that is CSP nameplate";
  parameter Boolean use_g3p3_field = true "If true, then always use the pre-LCOE-optimised 100MWe based G3P3 OELT";
  parameter Boolean use_TES = true "If true, use_TES";
  parameter SI.MassFlowRate H2_mdot_HX = H2_mdot_numberup + H2_mdot_recycle * 6 "Mass flow rate of the H2 going to the reactor";
  parameter Real factor_multiplier = 1 "Multiplier to multiply heat ratio to LHV of cold H2";
  parameter SI.MassFlowRate H2_mdot_recycle = if T_out_gas_target_H2 <= 1073.15 then (1.2718496962994789e+002 - 2.2858457240116964e-001 * T_out_gas_target_H2 + 1.0357499995171332e-004 * T_out_gas_target_H2 ^ 2) * factor_multiplier else 0.116168 * factor_multiplier "H2 needed to push the reaction and make the heat duty in reactor zero as a function of H2 out temperature. Based on 0.25 mT Fe/a";
  //Calculation of the H2 for burner
  //****************************** Simulation Set-up
  parameter Boolean set_pri_field_wspd = false "[H&T] true = using wspd_max dependent heliostat cost based on Emes 2017 https://is.gd/xSgpMV ";
  parameter Boolean set_detail_field_om = false "[H&T] true = use detail washing and field O&M cost";
  parameter Boolean set_const_dispatch = true "[CTRL] Constant dispatch of energy";
  parameter Boolean set_feedforward = true "[CTRL] if true then defocus heat is calculated using expensive PB model. False = defocus heat is a function of Q_in_rcv";
  parameter Boolean set_dispatch_optimiser = false;
  parameter Boolean set_new_storage_calc = false "[ST] if true = use tuffcrete and microporous insulation (Zebedee Kee's idea)";
  parameter Boolean set_use_neural_network = true;
  parameter Boolean set_dome_storage = true "[ST] If true then use dome storage calculation provided by Jeremy Sment from Sandia";
  parameter Boolean set_simple_PB_cost = true "[PB] sub system (excluding the primary heat exchanger) are evaluated using gen3_cost";
  parameter Boolean set_tnk_use_p_top = true "true if tank pressure is to connect to weather file";
  parameter Boolean set_tnk_enable_losses = true "true if the tank heat loss calculation is enabled";
  parameter Boolean set_external_storage = true "[ST] true if storage bins are not integrated with tower";
  parameter Boolean set_SAM_tower_cost = true "[H&T] true tower cost is evaluated to match SAM";
  parameter Boolean set_single_field = true "[H&T] True for single field, false for multi tower";
  parameter Boolean set_external_parasities = true "[PB] True = net power calculation in the PB model will consider parasitic losses";
  parameter Boolean set_use_wind = true "True if using wind stopping strategy in the solar field";
  parameter Boolean set_swaying_optical_eff = false "[H&T] True if optical efficiency depends on the wind speed due to swaying effect";
  parameter Boolean set_absolute_tower_cost = false "[H&T] False if tower cost is an absolute value, false means using SBP/SAM tower cost";
  parameter Boolean get_optics_breakdown = false "if true, the breakdown of the optical performance will be processed";
  parameter Boolean set_optics_verbose = false "[H&T] true if to save all the optical simulation details";
  parameter Boolean set_optics_view_scene = false "[H&T] true if to visualise the optical simulation scene (generate vtk files)";
  parameter Boolean set_scheduler = false "if true dispatched is scheduled (not optimised) based on the time alone. Must be OFF when dispatch
                                                                             optimiser is on";
  parameter Boolean HX_always_off = true "Whether to use industrial HX to draw heat from TES";
  parameter Boolean set_SMR_always_off = true "Whether SMR is on or off";
  //****************************** Importing medium and external files
  replaceable package Particle_Package = SolarTherm.Media.SolidParticles.CarboHSP_utilities;
  replaceable package Medium = SolarTherm.Media.SolidParticles.CarboHSP_ph "Medium props for Carbo HSP 40/70";
  replaceable package MedPB = SolarTherm.Media.CarbonDioxide_ph "Medium props for sCO2";
  //parameter String pri_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Prices/g3p3_TOD.motab") "[FN] Electricity or Time-of-day factor file";
  parameter String pri_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Prices/aemo_vic_2014.motab") "[FN] Electricity or Time-of-day factor file";
  parameter String sch_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Schedules/daily_sch_0.motab") if not set_const_dispatch "Discharging schedule from a file";
  parameter String wea_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Weather/pilbara3.motab") "[SYS] Weather file";
  parameter String DNI_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Weather/pilbara3.motab") "[CTRL] Weather file for dispatch optimisation - there was a bug in Modelica s.t. when I parse wea_file to dispatch optimisation, the 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            file path changes to .home/philgun..... which makes C program threw back segfault";
  parameter String price_file = pri_file;
  //************************ Weather Data Properties -  based on Dagget TMY2 sent by Luis G from UP Madrid to Philipe Gunawan at 7 April 2020
  parameter nSI.Angle_deg lon = 119.443359;
  //-116.783 "[SYS] Longitude (+ve East) TMY2 Dagget 1967 Location ID 23161";
  parameter nSI.Angle_deg lat = -24.166802;
  //34.86667 "[SYS] Lati1tude (+ve North) TMY2 Dagget 1967 Location ID 23161";
  parameter nSI.Time_hour t_zone = 8;
  //-8 "[SYS] Local time zone (UTC=0) TMY2 Dagget 1967 Location ID 23161";
  parameter Integer year = 1967 "[SYS] Meteorological year TMY2 Dagget 1967 Location ID 23161";
  //***************************** The fixed parameters -- H2 mdot feedstock --> dictattes how much Fe produced
  parameter Real upstreammultiplier(fixed = false) "The multiplier to be divided to H2_mdot_numberup and used to size the upstream sub-system";
  parameter SI.MassFlowRate H2_mdot_numberup = 2.5755 "H2 mass flow rate that will be consumed by the reactor ---> equivalent to 1.5 million tonne of Fe per year";
  parameter SI.MassFlowRate H2_mdot_feedstock = 0.42924999999999996 "H2 mass flow rate taht will be consumed by the reactor regardless the operating temperature of the reactor is. Equivalent to 0.25 million ton of Fe per year";
  parameter SI.Pressure p_H2O_excess_partial = H2O_product * 6 / 18 / (H2O_product * 6 / 18 + H2_mdot_recycle * 6 / 2) * 7e5 "Partial pressure of the H2O";
  parameter SI.Temperature T_H2O_excess_condensed = Modelica.Media.Water.WaterIF97_base.saturationTemperature(p_H2O_excess_partial);
  //****************************** Field simulation parameters
  parameter String opt_file(fixed = false);
  parameter String casefolder = "." "[H&T] Folder to which the OELT_Solstice look-up table will be stored";
  parameter Real metadata_list[9] = metadata(opt_file);
  parameter Solar_angles angles = Solar_angles.dec_hra "[SYS] Angles used in the lookup table file";
  parameter String sunshape = "buie" "Buie sunshape (buie) or pillbox sunshape (pillbox)";
  parameter Real buie_csr = 0.02 "circum solar ratio for Buie sunshape";
  parameter String field_type = "polar" "[H&T] Other options are : surround";
  parameter SI.Area A_helio = 144.375 "[H&T] Heliostat mirror area (m^2)";
  parameter Real AR_helio = 1 "[H&T] ratio between H_helio/W_helio";
  parameter SI.Length W_helio = sqrt(A_helio / AR_helio) "width of heliostat in m";
  parameter SI.Length H_helio = AR_helio * W_helio "height of heliostat in m";
  parameter SI.Angle slope_error = 1.53e-3 "[H&T] slope error of heliostats, in radiance";
  parameter Real windy_slope_error_factor = 2 / 1.53 "[H&T] Multiplier to the slope_error such that the product is the windy_slope_error (2 mrad at nominal)";
  parameter SI.Angle slope_error_windy = windy_slope_error_factor * slope_error "A larger optical error of heliostats under windy conditions, in radiance";
  parameter SI.Length H_tower = 200.676394 "[H&T] Tower height";
  parameter SI.Length H_tower_final(fixed = false);
  parameter SI.Length R_tower(fixed = false) "Inner tower radius";
  parameter SI.Length R1 = 102.599577 "[H&T] distance between the first row heliostat and the tower";
  parameter Real fb = 0.655736 "[H&T] factor to grow the field layout";
  parameter Real he_av_design = 0.99 "[H&T] Helisotats availability";
  //parameter Integer n_rays = 10000 "[H&T] number of rays for solstice";
  parameter Real n_row_oelt = 5 "[H&T] number of rows of the look up table (simulated days in a year)";
  parameter Real n_col_oelt = 22 "[H&T] number of columns of the lookup table (simulated hours per day)";
  parameter Real n_procs = 0 "[H&T] number of processors to run the MCRT, 0 is using maximum available num cpu, 1 is 1 CPU,i.e run in series mode";
  parameter SI.Efficiency helio_rho = 0.95 "[H&T] Reflectivity of heliostat. 0.95 is the default value in SolarPILOT";
  parameter SI.Efficiency helio_sf_ratio = 0.97 "[H&T] Reflective surface ratio. 0.97 is the default value in SolarPILOT";
  parameter SI.Efficiency helio_soil = 0.95 "[H&T] Heliostat soiling factor. 0.95 is the default value in SolarPILOT";
  parameter SI.Efficiency helio_uncertain_factor = 1 "[H&T] Uncertainty multiplier to the effective heliostat reflectance. The uncertain range is made by making the effective reflectance in the range of 0.8 to 0.95";
  parameter SI.Efficiency helio_refl = helio_rho * helio_sf_ratio * helio_soil * helio_uncertain_factor "The effective heliostat reflectance (product of helio_soil, helio_sf_ratio and helio_rho and the helio_uncertain_factor)";
  //****************************** Design condition of the Hybrid Plant
  parameter SI.Energy LHV_H2 = 120.21e6 "LHV of H2 source: https://h2tools.org/hyarc/hydrogen-data/lower-and-higher-heating-values-hydrogen-and-other-fuels";
  parameter SI.Temperature T_out_gas_target_H2 = 700 + 273.15 "Target inlet temperature of the gas to the reactor [K]";
  parameter SI.MassFlowRate H2O_product = 76.6575 / 20 "Steam as a product of the reaction in the reactor [kg/s] ---> equivalent to 250k tonne of Fe / annum";
  parameter Real f_TES = 1 "Fraction of the heat needed to heat up the H2_mdot_total to T_out_gas_DP_H2 provided by TES";
  parameter SI.TemperatureDifference delta_T_H2 = T_out_gas_target_H2 - T_in_gas_DP_H2 "Delta T of the gas, from the AEL to before the reactor [K]";
  parameter SI.Temperature T_out_gas_DP_H2_HX = T_in_gas_DP_H2 + f_TES * delta_T_H2 "Since the heat can be provided by either TES and/or burner, there are scenarios where the outlet of the gas from the HX is not at T_out_gas_target_H2, the rest of the heat is provided by the burner";
  parameter SI.Temperature T_out_intercooler_target = 20 + 273.14;
  parameter Modelica.Media.IdealGases.Common.DataRecord gas_data_H2 = Modelica.Media.IdealGases.Common.SingleGasesData.H2 "Used together with MedGas to calculate thermodynamic properties of H2";
  parameter SI.Efficiency eta_burner = 0.76 "Eta burner https://doi.org/10.3390/en13205526";
  parameter SI.HeatFlowRate Q_burner = H2_mdot_HX * (MedGas.h_T(gas_data_H2, T_out_gas_target_H2) - MedGas.h_T(gas_data_H2, T_out_gas_DP_H2_HX)) / eta_burner "Heat flow rate that has to be supplied by the burner";
  parameter SI.MassFlowRate H2_mdot_burner = Q_burner / LHV_H2 "H2 mass flow rate that will be burned to as a suplementary burner to reach T_out_gas_target_H2 if f_TES < 1";
  parameter SI.HeatFlowRate Q_burner_emergency = H2_mdot_HX * (MedGas.h_T(gas_data_H2, T_out_gas_target_H2) - MedGas.h_T(gas_data_H2, T_in_gas_DP_H2)) / eta_burner "Heat flow rate that has to be supplied by the burner";
  parameter SI.MassFlowRate H2_mdot_burner_emergency = Q_burner_emergency / LHV_H2 "H2 mass flow rate that will be burned to as a suplementary burner to reach T_out_gas_target_H2 if f_TES < 1";
  parameter SI.MassFlowRate H2_mdot_tank = H2_mdot_numberup + H2_mdot_burner "Mass flow rate that has to be drawn from the hydrogen storage tank -- a sum of H2 feedstock and H2 for burner (if burner is used)";
  parameter SI.MassFlowRate H2_mdot_tank_emergency = H2_mdot_numberup + H2_mdot_burner_emergency "Mass flow rate that has to be drawn from the hydrogen storage tank -- a sum of H2 feedstock and H2 for burner (if burner is used)";
  parameter SI.Area A_recup(fixed = false);
  parameter SI.CoefficientOfHeatTransfer U_recup = 65 "gas to gas overall U value. source: Cao, Eduardo. 2010. “TYPICAL HEAT TRANSFER COEFFICIENTS.” Chap. G in Heat Transfer in Process Engineering. 1st ed. New York: . https://www.accessengineeringlibrary.com/content/book/9780071624084/back-matter/appendix7";
  parameter SI.Power P_hybrid_system_final(fixed = false) "Hybrid system nameplate [W]";
  parameter SI.MassFlowRate H2_mdot_target = H2_mdot_numberup "Hydrogen production per second in that has to be supplied by the tank [kg/s]";
  parameter SI.MassFlowRate O2_mdot_target(fixed = false) "Oxygen production per second in AEL calculated using initial equation in AEL [kg/s]";
  parameter Real CSP_fraction_final(fixed = false);
  parameter Real PV_fraction_final(fixed = false);
  parameter SI.Power P_CSP = if use_g3p3_field then 906924733.561042 else CSP_fraction_final * P_hybrid_system_final / upstreammultiplier "[PB] Power block net rating at design point [W]";
  parameter Boolean on_CSP = if P_CSP > 0 then true else false "Boolean to control CSP block";
  parameter SI.Power PV_Target = PV_fraction_final * P_hybrid_system_final "PV array nameplate in W";
  parameter SI.Power P_heater = PV_Target - P_hybrid_system_final "Rating of the electrical heater [W]";
  parameter SI.Efficiency eta_heater = 0.99 "Heater electric to thermal efficiency https://doi.org/10.3390/en14123437";
  //****************************** Design condition of the CSP plant
  parameter Real SMR_reaction_conversion = 1 "Steam methane reforming conversion extent";
  parameter Real W_consumption_SMR = 1492000 * (H2_mdot_target / (8994 / 3600)) ^ scaler_n "Electricity needs to be fed to the SMR [W]";
  parameter Real H2_mol_target = H2_mdot_target * 1000 / 2 "Mol target of H2 per second [mol/s]";
  parameter SI.HeatFlowRate Q_in_rcv = if use_g3p3_field then 906924733.561042 else if P_CSP > 0 then (P_CSP / eff_blk + Q_HX_industrial) / eta_rcv_assumption * SM else (P_gross / eff_blk + Q_HX_industrial) / eta_rcv_assumption * SM "Incident heat flow rate to the receiver at design point [Wth]";
  parameter String rcv_type = "particle" "[RCV] other options are : flat, cylindrical, stl";
  parameter SI.Area A_rcv(fixed = false) "Receiver aperture area is calculated during the initialisation";
  parameter nSI.Angle_deg tilt_rcv = 0 "[RCV] tilt of receiver in degree relative to tower axis";
  parameter Real SM = 2.5 "[SYS] Solar multiple";
  parameter SI.Power P_net_default_value = 123456789 "Default value to handle P_net = 0 [W]";
  parameter Real PB_fraction_final(fixed = false);
  parameter SI.Power P_net = if PB_fraction_final < 1e-3 or set_no_thermal then P_net_default_value else PB_fraction_final * P_hybrid_system_final / upstreammultiplier;
  /*
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      if CSP_fraction > 0 then 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            P_CSP 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      else 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            if P_heater > 0 then 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  if P_heater*0.5 < P_hybrid_system then 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        P_heater*0.5  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  else P_hybrid_system 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            else 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  P_net_default_value
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "Power of the PB to size the components [W] ---> if the PB size is the result of oversizing PV, then the size of the PB has to be divided by two since the eff. of the PB is around 0.4";*/
  //********************* PB size == CSP size since CSP size >0
  //********************* If CSP size == 0 (No CSP)
  //********************* There is heater, therefore must have PB to harness the power
  //********************* If P_heater < electrolyser name plate, size the PB according to heater
  //********************* If P_heater > electrolyser name plate, size the PB according to electrolyser rating
  //********************* If CSP size == 0 (No CSP) and no heater, put the default value to size the PB
  parameter SI.Power P_gross = P_net / (1 - par_fr) "The mechanical power of the PB (Turbine power - all of compressors) before cooling and parasities losses";
  parameter SI.Efficiency eff_blk(fixed = false) "Power block efficiency at design point";
  parameter SI.Temperature T_in_rec = T_cold_set "Particle inlet temperature to particle receiver at design";
  parameter SI.Irradiance dni_des = 950 "[SYS] DNI at design point - common values that everyone use?";
  parameter SI.Efficiency eta_rcv_assumption = 0.88 "[RCV] Receiver net thermal efficiency guess value to size the receiver";
  parameter Real n_helios = metadata_list[1] "Number of heliostats";
  parameter SI.Temperature T_amb_des_rcv = from_degC(10) "[RCV] Design point ambient temp of the receiver (K)";
  parameter SI.Area A_field = metadata_list[1] * metadata_list[2] "Heliostat field reflective area (m^2)";
  parameter Real A_land = metadata_list[9] "Land area consumed by the plant - calculated value by Solstice (m^2)";
  parameter SI.Velocity Wspd_max = 15.65 if set_use_wind "[CTRL] Wind stow speed - based on DOE suggestion (m/s)";
  parameter SI.Efficiency packing_factor = 0.6 "[RCV] Based on EES model by Sandia / Luis";
  //****************************** Design condiction of the PV array
  parameter nSI.Angle_deg azi_s = 180 "Surface azimuth angle";
  parameter nSI.Angle_deg ele_s = 30 "Surface elevation angle";
  parameter Integer N_parallel_final_PV(fixed = false) "Number of PV - Inverter unit";
  parameter Integer N_series_final_PV(fixed = false) "Number of array in series";
  //****************************** Industrial HX sizing
  replaceable package Med_PCL = SolarTherm.Media.SolidParticles.CarboHSP_ph;
  replaceable package Med_AIR = Modelica.Media.Air.ReferenceAir.Air_ph;
  replaceable package AIR = Modelica.Media.Air.ReferenceAir.Air_Utilities;
  parameter Integer num_seg = 6 "Segmentation of heat exchanger";
  parameter SI.Length W_HX = 8 "Heat exchanger air manifold width in m";
  parameter SI.Length H_HX = 8 "Heat exchanger air manifold height in m";
  parameter SI.Temperature T_in_AIR_DP = 25 + 273.15 "Air inlet temperature at design point (K)";
  parameter SI.Temperature T_out_AIR_DP = 500 + 273.15 "Air outlet temperature at design point (K)";
  parameter SI.MassFlowRate m_dot_AIR_DP = 1e-5 "Desired mass flow rate of air at design point (kg/s)";
  parameter SI.MassFlowRate m_dot_PCL_industrial(fixed = false);
  parameter SI.Area A_HX_industrial(fixed = false);
  parameter SI.Power Q_HX_industrial = 0;
  /*Design of the indirect particle gas hx*/
  replaceable package MedGas2_O2 = Modelica.Media.IdealGases.Common.SingleGasNasa(data = Modelica.Media.IdealGases.Common.SingleGasesData.O2, mediumName = "Hydrogen", fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.O2});
  replaceable package MedGas2_H2 = Modelica.Media.IdealGases.Common.SingleGasNasa(data = Modelica.Media.IdealGases.Common.SingleGasesData.H2, mediumName = "Hydrogen", fluidConstants = {Modelica.Media.IdealGases.Common.FluidData.H2});
  parameter Modelica.Media.IdealGases.Common.DataRecord gas_data_O2 = Modelica.Media.IdealGases.Common.SingleGasesData.O2;
  parameter SI.Temperature T_out_gas_target_O2 = 600 + 273.15 "Design inlet temperature of the gas [K]";
  parameter SI.Temperature T_out_pcl_off = T_cold_set "Outlet temperature of the particle [K]";
  /*On-design parameters*/
  parameter SI.Temperature T_in_gas_DP_O2 = 80 + 273.15 "Design inlet temperature of the gas [K]";
  parameter SI.Temperature T_out_gas_DP_O2 = 600 + 273.15 "Design outlet temperature of the gas [K]";
  parameter SI.Temperature T_in_pcl_DP_O2 = 800 + 273.15 "Design inlet temperature of the particle [K]";
  parameter SI.Temperature T_out_pcl_DP_O2 = 550 + 273.15 "Design outlet temperature of the pcl [K]";
  parameter SI.Length d_gas_pipe_O2 = 0.5e-3;
  parameter SI.Length w_ch_O2 = 6e-3;
  parameter SI.Length th_wall_O2 = 1e-3;
  parameter SI.Area A_plate_O2 = 0.1;
  parameter SI.Temperature T_in_gas_DP_H2 = 590 + 273.15 "Design inlet temperature of the gas [K]";
  parameter SI.Temperature T_in_pcl_DP_H2 = T_hot_set "Design inlet temperature of the particle [K]";
  parameter SI.Temperature T_out_pcl_DP_H2 = T_cold_set "Design outlet temperature of the pcl [K]";
  parameter SI.Length d_gas_pipe_H2 = 0.5e-3;
  parameter SI.Length w_ch_H2 = 6e-3;
  parameter SI.Length th_wall_H2 = 1e-3;
  parameter SI.Area A_plate_H2 = 0.1;
  parameter SI.Area A_HX_O2(fixed = false) "A HX [m2]";
  parameter SI.Area A_HX_H2(fixed = false) "A HX [m2]";
  parameter SI.HeatFlowRate Q_HX_H2(fixed = false) "Q HX [Wt]";
  parameter String saved_model_dir_H2 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/SurrogateModels/IndirectParticleGasHX/H2/3/surrogate_model");
  parameter String saved_model_dir_O2 = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/SurrogateModels/IndirectParticleGasHX/H2/surrogate_model");
  parameter Real[14] X_max_H2 = {873.147608736259, 1261.78929155903, 4.99874743580914, 1273.13856914357, 973.140145246266, 0.004999810880909, 0.059989212153262, 0.009999872519837, 0.999594690837943, 4.99848048118983, 1273.12769197071, 873.148855294655, 9528.33599417378, 3344.86768214966};
  parameter Real[14] X_min_H2 = {863.152208048605, 883.216016729446, 0.300902007162029, 1023.35331623183, 773.22731757903, 0.000100141721352, 0.002038008911878, 0.000501216437113, 0.050087863490287, 0.301800982674091, 1023.19701709183, 863.151268150423, 34.633804858446, 176.628395023518};
  parameter Real out_max_H2 = 49.98120625;
  parameter Real out_min_H2 = 1.50082196;
  parameter Real[14] X_max_O2 = {8.63151005e+02, 8.83170753e+02, 3.00902007e-01, 1.02335332e+03, 7.73227318e+02, 1.00141721e-04, 2.03800891e-03, 5.01216437e-04, 5.00878635e-02, 3.01800983e-01, 1.02319702e+03, 8.63151268e+02, -3.77393901e+02, 1.76628395e+02};
  parameter Real[14] X_min_O2 = {343.15195730620184, 823.2152898282988, 3.000271703042015, 1023.1544618085376, 773.1647588827092, 0.0001000259962433, 0.0020098650004692, 0.0005002794822282, 0.0500448367343232, 3.000284480139946, 973.1577692357974, 353.05000781902413, 353.09829689205435, 477.3614596432316};
  parameter Real out_max_O2 = 19.99971179818061;
  parameter Real out_min_O2 = 1.5008219564222385;
  /*Off-design conditions*/
  parameter SI.Temperature T_in_gas_off_H2(fixed = false);
  // = T_in_gas_DP_H2 "Design inlet temperature of the gas [K]";
  parameter SI.Temperature T_in_gas_off_O2 = 80 + 273.15 "Design inlet temperature of the gas [K]";
  //****************************** Design condiction of the Electrolyser
  parameter SI.Temperature T_electrolyser = 80 + 273.15 "Working temperature of the electrolyser";
  //****************************** Design condition of the Particle Receiver
  parameter Real ar_rec = 1.067451 "[RCV] Height to diameter aspect ratio of receiver aperture";
  parameter Real rec_fr(fixed = false) "Receiver loss fraction at design point";
  inner parameter SI.Efficiency eta_rec_th_des = 1 - rec_fr "Receiver thermal efficiency (Q_pcl / Q_sol)";
  parameter Real f_loss = 0.000001 "[RCV] Fraction of particles flow lost in receiver";
  parameter SI.Length th_w = 0.05 "[RCV] Backwall thickness of the receiver";
  parameter SI.ThermalConductance k_w = 0.2 "[RCV] Thermal conductance of the back wall of the receiver";
  parameter SI.CoefficientOfHeatTransfer h_conv_curtain = 32. "[RCV] Constant convective heat transfer coefficient (curtain) [W/m^2-K]";
  parameter SI.CoefficientOfHeatTransfer h_conv_backwall = 10. "[RCV] Constant Convective heat transfer coefficient (backwall) [W/m^2-K]";
  parameter SI.Length d_p = 0.00035 "[RCV] Particle diameter [m]";
  parameter SI.SpecificHeatCapacity cp_s = 1200 "[RCV] particle specific heat capacity [J/kgK]";
  parameter SI.Density rho_s = 3300 "[RCV] Particle intrinsic density [kg/m3]";
  parameter Real eps_s = 0.9 "[RCV] Particle emmisivity - based on EES code by Sandia";
  parameter Real abs_s = 0.9 "[RCV] Particle absorptivity - based on EES code by Sandia";
  parameter Real F = 0.54 "[RCV] View Factor of the particle curtain - based on EES code by Sandia";
  parameter Real eps_w = 0.8 "[RCV] Receiver wall emmisivity - based on EES code by Sandia";
  parameter Integer discretisations = 30 "Discretisation of the particle receiver during sizing";
  //****************************** NN Based Receiver Parameter
  parameter Integer inputsize_rcv = 7;
  parameter SI.Length H_drop_max_rcv = 45.0 "[RCV] maximum drop height where the receiver surrogate model is still valid [m]";
  parameter Real ar_rec_max_rcv = 3.49968553847121 "[RCV] maximum receiver aspect ratioto the receiver surrogate model [-]";
  parameter SI.Temperature T_out_max_rcv = 1273.11962227051 "[RCV] maximum particle outlet temperature [K]";
  parameter SI.HeatFlowRate Q_in_max_rcv = 3016745046.48895 "[RCV] maximum incident heat to the receiver [W]";
  parameter SI.Temperature T_in_max_rcv = 1048.00233180374 "[RCV] maximum particle inlet temperature to the surrogate model [K]";
  parameter SI.Temperature T_amb_max_rcv = 313.139295255538 "[RCV] maximum ambient temperature to the surrogate model [K]";
  parameter Real F_wind_max_rcv = 6.54742745445228 "[RCV] maximum wind factor to the surrogate model [-]";
  parameter SI.Length H_drop_min_rcv = 5.0 "[RCV] minimum drop height where the receiver surrogate model is still valid [m]";
  parameter Real ar_rec_min_rcv = 0.253011055900887 "[RCV] minimum receiver aspect ratioto the receiver surrogate model [-]";
  parameter SI.Temperature T_out_min_rcv = 1073.33624339311 "[RCV] minimum particle outlet temperature [K]";
  parameter SI.HeatFlowRate Q_in_min_rcv = 17816023.817307 "[RCV] minimum incident heat to the receiver [W]";
  parameter SI.Temperature T_in_min_rcv = 823.276338690234 "[RCV] minimum particle inlet temperature to the surrogate model [K]";
  parameter SI.Temperature T_amb_min_rcv = 253.173406522872 "[RCV] minimum ambient temperature to the surrogate model [K]";
  parameter Real F_wind_min_rcv = 1.00000038230021 "[RCV] minimum wind factor to the surrogate model [-]";
  parameter Real y_max_rcv = 0.999960147834589 "[RCV] maximum thermal efficiency";
  parameter Real y_min_rcv = 0.025096078157653 "[RCV] minimum thermal efficiency";
  parameter Real[inputsize_rcv] X_max_rcv = {H_drop_max_rcv, ar_rec_max_rcv, T_out_max_rcv, Q_in_max_rcv, T_in_max_rcv, T_amb_max_rcv, F_wind_max_rcv};
  parameter Real[inputsize_rcv] X_min_rcv = {H_drop_min_rcv, ar_rec_min_rcv, T_out_min_rcv, Q_in_min_rcv, T_in_min_rcv, T_amb_min_rcv, F_wind_min_rcv};
  parameter String saved_model_dir_rcv = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/SurrogateModels/ParticleReceiver/single_aperture_hot_hydrogen") "[RCV] path to which the static particle receiver surrogate model is stored";
  //****************************** OnTheFlySurrogate Power Block Parameters
  /************************************************************************************************************** /
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /       NREL PB and CEA are sizing the power block based on cycle power. In this code it is called P_gross   /
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /      -cycle power: W_turb_des - W_comp_des - W_recomp_des                                                  / 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /      -net power : (cycle_power - W_cooling_fan) * eta_motor * (1-f_fixed_load)                             /
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /      All power above are before parasities_input => heliostat field, pump/lift power consumption           /
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /      The output of the on the fly surrogates are: eta_gross and eta Q                                      /
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /      eta_gross: (W_cycle-W_cooling) / Q_HX                                                                 /
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /      eta_Q: (Q_HX / Q_HX_des)                                                                              /
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        /      The power block initalisation will produce Q_HX_des, regardless which PB model is used                /   
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ***************************************************************************************************************/
  //******************************** OnTheFlySurrogate PB Simulation Set-up
  parameter String base_path = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/SurrogateModels/PowerBlock") "[PB] Base path that points to which folder the C program located";
  parameter String SolarTherm_path = Modelica.Utilities.Files.loadResource("modelica://SolarTherm") "[PB] Base path that points to which folder SolarTherm libs are located";
  parameter Integer inputsize_PB = 3 "[PB] the input dimenstion of the surrogate model PB";
  parameter Integer outputsize_PB = 2 "[PB] the output dimension of the surrogate model PB";
  parameter Real tolerance_kriging = 6e-3 "[PB] Error tolerance of the on-the-fly surrogate model Kriging (Mean Absolute Error)";
  parameter Real tolerance_ANN = 6e-3 "[PB] Error tolerance of the on-the-fly surrogate model ANN (Mean Absolute Error)";
  parameter Integer which_surrogate = 2 "1 for Kriging, 2 for ANN";
  parameter Real eta_gross_base(fixed = false) "By product of PB initialisation, regardless which PB model is chosen e.g CEA or SAM";
  parameter Real eta_Q_base(fixed = false) "By product of PB initialisation, regardless which PB model is chosen e.g CEA or SAM";
  //******************************** General PB Parameters
  parameter SI.ThermodynamicTemperature T_low = 41 + 273.15 "[PB] Compressor inlet temperature at design point (K)";
  parameter Integer N_HTR = 15 "[PB] Discretisation of high temperature recuperator - CEA PB Parameter";
  parameter Integer N_LTR_parameter = 15 "Discretisation of low temperature recuperator";
  parameter Real f_fixed_load = 0.0055 "[PB] Fixed load constantly consumed by PB when it operates";
  parameter SI.ThermodynamicTemperature T_high = 700 + 273.15 "[PB] Turbine inlet temperature [K] - CEA PB Parameters";
  parameter Real gamma = 0.28 "[PB] Part of the mass flow going to the recompression directly - CEA PB parameters";
  parameter Real par_fr = 0.02 "[PB] Guess value of the parasitics fraction of power block rating at design point";
  //******************************** CEA PB Parameters
  parameter SI.Temperature T_in_ref_blk = T_hot_set "Particle inlet temperature to particle heat exchanger at design";
  parameter SI.AbsolutePressure p_high = 25e6 "[PB] Compressor outlet pressure (Pa)";
  parameter Real PR = 2.78 "[PB] Pressure ratio - CEA PB parameters";
  parameter SI.TemperatureDifference pinch_exchanger = 25 "[PB] pinch point of the heat exchanger - CEA PB parameters";
  parameter SI.TemperatureDifference dTemp_HTF_PHX = T_hot_set - T_cold_set;
  parameter SI.Temperature blk_T_amb_des = T_low - dT_mc_approach "[PB] Ambient temperature at design for power block";
  //******************************** NREL SAM PB Parameters - some parameters are defined in the CEA's section
  parameter Integer htf_choice = 50 "[PB] 50 is user defined fluid properties";
  parameter SI.TemperatureDifference dT_PHX_hot_approach(fixed = false) "[PB] difference between hot HTF and TIT (C/K). Product of iteration s.t. the T_HTF_cold_design == T_cold_set";
  parameter SI.TemperatureDifference dT_PHX_cold_approach = 15 "[PB] Temp. difference between cold HTF and cold CO2 PHX inlet (C/K)";
  parameter SI.Efficiency eta_comp_main = 0.89 "[PB] main compressor isentropic efficiency";
  parameter SI.Efficiency eta_comp_re = 0.89 "[PB] re-compressor isentropic efficiency";
  parameter SI.Efficiency eta_turb = if P_gross > 3e7 then 0.93 else 0.85 "turbine isentropic efficiency";
  parameter SI.TemperatureDifference dT_mc_approach = 6.0 "[PB] Temp. difference between main compressor CO2 inlet and ambient (C/K)";
  parameter Integer which_PB_model = 1 "0 is for CEA power block, 1 is for NREL-SAM power block";
  parameter String HTF_name = "CarboHSP" "[PB]";
  parameter Real NREL_PB_configurations[13](each fixed = false) if which_PB_model == 1 "array to store design point sizing result";
  //****************************** Storage Parameters
  parameter SI.ThermalInsulance U_value_hot_tank = 0.25 "[ST] Desired U_value for the tanks";
  parameter SI.ThermalInsulance U_value_cold_tank = 0.25 "[ST] Desired U value for the tanks";
  parameter Real NS_particle = 0.05 "[ST] Fraction of additional non-storage particles";
  parameter SI.Temperature T_cold_set = T_in_gas_DP_H2 + 10 "[ST] Cold tank target temperature ==  HTF outlet temperature from PB at design point (K)";
  parameter SI.Temperature T_hot_set = 1200 "[ST] Hot tank target temperature == HTF inlet temperature to the PB at design point (K)";
  parameter SI.Temperature T_cold_start = T_cold_set "Cold tank starting temperature";
  parameter SI.Temperature T_hot_start = T_hot_set "Hot tank starting temperature";
  /*Thermophysical of the particle*/
  parameter Real split_cold = (100 - hot_tnk_empty_ub + 1) / 100 "Starting medium fraction in cold tank, must be the function of the upper bound trigger level of the hot tank 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           so the simulation wont crash at t=0, since the control logic use t_on - t_start etc";
  parameter SI.Density rho_cold_set = Particle_Package.rho_T(T_cold_set) "Cold particles density at design";
  parameter SI.Density rho_hot_set = Particle_Package.rho_T(T_hot_set) "Hot particles density at design";
  parameter SI.MassFlowRate m_dot_pcl_DP_H2_HX(fixed = false) "mass flow rate of the pcl to heat up H2 in the indirect particle gas HX";
  //parameter SI.Mass m_max = E_max / (h_hot_set - h_cold_set) "Max particles mass in tanks [kg]";
  parameter SI.Mass m_max = m_dot_pcl_DP_H2_HX * t_storage * 3600 "Max particles mass in tanks [kg]";
  parameter SI.Energy E_max = t_storage * 3600 * (m_max * (h_hot_set - h_cold_set)) "Maximum tank stored energy [J]";
  //m_flow_blk * t_storage * 3600 "Max particles mass in tanks [kg]";
  parameter SI.Volume V_max = m_max / (rho_s) / packing_factor "Volume needed to host particles in the tank with certain packing factor value";
  parameter Real tank_ar = 1.17 "[ST] storage aspect ratio";
  parameter SI.Diameter D_storage = H_storage / tank_ar "Storage tank diameter";
  parameter Real X_dome = 1 / (4 * tank_ar ^ 2) - 1 / (8 * tank_ar ^ 3) + 1 / (12 * tank_ar ^ 3);
  parameter SI.Length H_storage = if set_dome_storage then ceil((V_max / (CN.pi * X_dome)) ^ (1 / 3)) else ceil((4 * V_max * tank_ar ^ 2 / CN.pi) ^ (1 / 3)) "Storage tank total height [m]";
  parameter SI.Area SA_storage = if set_dome_storage then CN.pi * D_storage * (H_storage - 0.5 * D_storage) + 2 * CN.pi * (D_storage / 2) ^ 2 else CN.pi * D_storage * H_storage "Storage tank surface area";
  //************************ Surface area of the hemisphere and the cylinder
  //************************ Assuming cylindrical storage
  parameter SI.Length Th_refractory_hot_tank = 1.003 "[ST] Thickness of the refractory of the hot tank [m] - Dome Storage only";
  parameter SI.Length Th_refractory_cold_tank = 0.895 "[ST] Thickness of the refractory of the cold tank [m] - Dome Storage only";
  //****************************** Power Block Technical Parameters - CEA Power Block
  /*Heat Exchanger Parameters*/
  parameter Real pinch_recuperator = 15 "Pinch point of the recuperators - CEA PB Parameters";
  parameter SI.Temperature T_out_ref_blk = T_cold_set "Particle outlet temperature from particle heat exchanger at design";
  parameter SI.Temperature T_in_ref_co2 = T_cold_set - pinch_exchanger "CO2 inlet temperature to particle heat exchanger at design";
  parameter SI.Temperature T_out_ref_co2 = T_high "CO2 outlet temperature from particle heat exchanger at design";
  parameter Integer N_exch_parameter = 15 "Heat exchanger discretisation - CEA PB Parameters";
  parameter Real nu_min_blk = 0.01 "[PB] Minimum allowed part-load mass flow fraction to power block";
  parameter SI.CoefficientOfHeatTransfer h_conv_CO2 = 2000 "[PB] According to Sandia EES Code";
  parameter Real A_HX(fixed = false) "Heat exchanger total surface area --> UA_HX / U_HX";
  parameter Real UA_HX(fixed = false) "By product of initialisation of the PB model regardless which PB model is chosen";
  parameter Real U_HX(fixed = false) "Calculated parameter in the initial equation";
  //****************************** Particle Lift Parameters
  parameter SI.Height dh_liftRC = H_tower_final "Vertical displacement in receiver lift";
  parameter SI.Height dh_liftHX = 10 "[H&T] Vertical displacement in heat exchanger lift";
  parameter SI.Height dh_LiftCold = H_storage "Vertical displacement in cold storage lift";
  parameter SI.Efficiency eff_lift = 0.8 "[H&T] Lift total efficiency";
  //****************************** Power Plant Control
  parameter SI.Time t_start = 3600 "[CTRL] Start-up traking delay";
  parameter SI.Angle ele_min = 0.0872665 "[CTRL] Heliostat stow deploy angle = 5 degree";
  parameter SI.HeatFlowRate Q_flow_defocus_calculated(fixed = false);
  parameter SI.HeatFlowRate Q_flow_defocus = if set_feedforward == true then Q_flow_defocus_calculated else Q_flow_des / (1 - rec_fr) "Solar field thermal power at defocused state";
  parameter Real nu_start = 0.29 "[CTRL] Fraction of Q_in_rcv above which the receiver stop working - SAM 2020.2.29 value minus hysterisis";
  parameter Real nu_min_sf = 0.3 "[CTRL] Fraction of Q_in_rcv below which the receiver stop working - SAM 2020.2.29 value";
  parameter Real nu_defocus = 1 "[CTRL] Energy fraction to the receiver at defocus state -> nu_defocus * Q_flow_defocus_calculated";
  parameter Real hot_tnk_empty_lb = 5 "[CTRL] Hot tank empty trigger lower bound";
  parameter Real hot_tnk_empty_ub = 10 "[CTRL] Hot tank trigger to start dispatching";
  parameter Real hot_tnk_full_lb = 93 "[CTRL] Hot tank full trigger lower bound";
  parameter Real hot_tnk_full_ub = 95 "[CTRL] Hot tank full trigger upper bound";
  parameter Real cold_tnk_defocus_lb = 5 "[CTRL] Cold tank defocus trigger lower bound";
  parameter Real cold_tnk_defocus_ub = 7 "[CTRL] Cold tank defocus trigger upper bound";
  parameter Real cold_tnk_crit_lb = 10 "[CTRL] Cold tank critically empty trigger lower bound";
  parameter Real cold_tnk_crit_ub = 15 "[CTRL] Cold tank critically empty trigger upper bound";
  //****************************** Dispatch Optimiser Parameters
  parameter Integer horizon = 24 "[CTRL] Forecast horizon of the dispatch optimisation algoriuthm in GLPK";
  parameter Real dt = 0.5 "[CTRL] delta t of the optimisation. For each how many hour the dispatch optimiser is gonna be called";
  parameter Real etaG(fixed = false) "By product of power block intialisation";
  parameter Real etaC = metadata_list[3];
  parameter Real SLminrel = hot_tnk_empty_lb / 100 "hot tank empt trigger point";
  parameter Real Ahelio = A_field;
  parameter Real const_t = -dt * 3600;
  //****************************** Dual Tank Dispatch Optimiser Parameters
  parameter Real c0 = 20 "Rewards for producing hot hydrogen";
  parameter Real c1 = 1 "Penalty for burning hydrogen in the hydrogen tank";
  parameter Real c2 = 0 "Penalty for dumping heat from solar field";
  parameter Real c3 = 1 "Penalty for dumping electricity from PV";
  parameter Integer n_horizon = 24 "Horizon , how many look-ahead time steps";
  parameter Real dt_dual_tank = 1 "Time step in [h]";
  parameter Real c_ratio = H2_mdot_HX / H2_mdot_numberup "Ratio between H_recycled and H2_reactor [-]";
  parameter Real etaRCV(fixed = false) "Receiver efficiency at design point [-]";
  parameter Real etaPV = 0.2 "PV efficiency at design point [-]";
  parameter Real etaBurner = eta_burner "Hydrogen burner efficiency at design point [-]";
  parameter Real etaResistive = eta_heater "Resistive heater efficiency [-]";
  //----------------- upstream subsystem
  parameter Real P_AEL_nameplate = P_hybrid_system_final / 1e6 "Electrolyser nameplate [MWe]";
  parameter Real A_sf = A_field * upstreammultiplier "Total area of the solar field [m2] -- has to be multiplied by upstreammultiplier";
  parameter Real PB_size = if PB_fraction < 1e-3 then 0 else P_name / 1e6 * upstreammultiplier "PB name plate [MWe] -- has to be multiplied by upstreammultiplier";
  parameter Real TES_capacity = if use_TES then E_max * 2.77778e-10 * upstreammultiplier else 0 "TES capacity [MWhth] -- has to be multiplied by upstreammultiplier";
  parameter Real Q_TES_HX_max = if use_TES then Q_burner_emergency * eta_burner / 1e6 else 0 "Maximum dispatchable thermal power from the TES to the TES-heater to heat-up H2_reactor + H2_recycle [MWth] (Q_burner_emergency is already multiplied by upstream multiplier)";
  parameter Real E_TES_min = if use_TES then hot_tnk_empty_lb / 100 * TES_capacity else 0 "Minimum TES level [MWhth] (TES_capacity is already multiplied by upstream multiplier)";
  parameter Real const_t_dual_tank = -dt_dual_tank * 3600;
  //----------------- downstream subsystem
  parameter Real Q_H2_reactor_max = H2_mdot_numberup * LHV_H2 / 1e6 "Maximum heat from the H₂ storage to the reactor (MWth) -- caried by the H2_reactor";
  parameter Real LHV = LHV_H2 / 1e6 "Lower heating value of the hydrogen [MJ]";
  parameter Real delta_h_H2 = Q_burner_emergency * eta_burner / H2_mdot_HX / 1e6 "Enthalpy change of the H2 before and after the pre-heating [MJ/kg]";
  parameter Real H2stg_capacity = E_H2_tank * 2.77778e-10 "H2 storage capacity [MWhth] -- already multiplied by upstreammultiplier";
  parameter Real H2stg_min = hot_tnk_empty_lb / 100 * H2stg_capacity "Minimum level of H2 storage [MWhth]";
  //----------------- calculated parameters
  parameter Real etaPB(fixed = false) "Power block efficiency at design point [-]";
  parameter Real etaAEL(fixed = false) "Electrolyser efficiency at design point [-]";
  //----------------- variables
  Real time_simul_dual_tank "Starting time [s]";
  Real H2stg_init "Initial level of H2 storage [MWhth]";
  Real E_TES_init "Initial level of TES [MWhth] -- has to be multiplied by upstreammultiplier";
  Real[9] optimalSolution "Array to store the the optimal solution of dual tank dispatch optimisation";
  //****************************** H2 and O2 tank parameters
  parameter Real t_storage_H2_final(fixed = false, unit = "h") "[ST] H2 Storage capacity";
  parameter Real t_storage_H2_threshold(unit = "h") = 3 "[ST] H2 Storage capacity";
  parameter Real t_storage_O2(unit = "h") = 10.0 "[ST] O2 Storage capacity";
  parameter SI.Energy E_H2_tank = H2_mdot_tank * t_storage_H2_final * 3600 * LHV_H2 "Capacity of the energy-based H2 tank [J]";
  parameter SI.Mass M_H2_tank = H2_mdot_tank * t_storage_H2_final * 3600 "Total mass of the hydrogen contained in the H2 tank [kg]";
  parameter SI.Mass M_O2_tank = O2_mdot_target * t_storage_O2 * 3600 "Capacity of the mass-based O2 tank [kg]";
  //******************************* Other Parameters
  parameter SI.HeatFlowRate Q_flow_des(fixed = false) "Heat to power block at design--> result of PB initialisation (W)";
  parameter SI.Length H_rcv = if use_g3p3_field then 30.72272 else sqrt(A_rcv * ar_rec) "Receiver aperture height";
  parameter SI.Length W_rcv = if use_g3p3_field then 28.781387 else A_rcv / H_rcv "Receiver aperture width";
  parameter SI.Length L_rcv = 1 "[RCV] Receiver length (depth) (m)";
  parameter SI.SpecificEnthalpy h_cold_set = Particle_Package.h_T(T_cold_set) "Cold particles specific enthalpy at design";
  parameter SI.SpecificEnthalpy h_hot_set = Particle_Package.h_T(T_hot_set) "Hot particles specific enthalpy at design";
  parameter SI.MassFlowRate m_flow_fac(fixed = false);
  parameter SI.MassFlowRate m_flow_rec_max = 1.5 * m_flow_fac "Maximum mass flow rate to receiver";
  parameter SI.MassFlowRate m_flow_rec_start = 0.8 * m_flow_fac "Initial https://pubs.acs.org/doi/pdf/10.1021/jp206115por guess value of mass flow rate to receiver in the feedback controller";
  parameter SI.MassFlowRate m_flow_blk(fixed = false);
  parameter SI.Power P_name = if abs(P_net - P_net_default_value) < 1 then 0 else P_net "Nameplate rating of power block";
  //******************************* Financial parameters
  parameter Currency currency = Currency.USD "[FN] Currency used for cost analysis 2 is USD 1 is AUD";
  parameter Real r_i = 0.025 "[FN] Inflation rate";
  parameter Real r_disc_nom = 0.0701 "[FN] Nominal discount rate";
  parameter Real r_disc = (1 + r_disc_nom) / (1 + r_i) - 1 "[FN] Real discount rate : Calculated using Fisher equation";
  parameter Integer t_life = 30 "[FN] Lifetime of plant";
  parameter Integer t_cons = 0 "[FN] Years of construction";
  parameter Real r_cur = 0.71 "[FN] The currency rate from AUD to USD valid for 2019. See https://www.rba.gov.au/";
  parameter Real r_contg = 0.1 "[FN] Contingency rate";
  parameter Real r_cons = 0.09 "[FN] Construction cost rate";
  //******************************* Specific cost for field, particle receiver, tower and lift
  parameter FI.AreaPrice pri_field = if set_pri_field_wspd == true then FI.heliostat_specific_cost_w_spd(Wspd_max = Wspd_max, A_helio = A_helio) * 0.3716 else 75 "Heliostat field specific price (USD/m^2)";
  //********** Field price depends on the wind speed and mirror size Emes 2017 https://is.gd/xSgpMV
  //******************************* Field price is based on Gen3 Roadmap Report
  parameter FI.AreaPrice pri_site = 10 "Site improvements cost per area [USD/m2 of heliostat]";
  parameter FI.AreaPrice pri_land = 10000 / 4046.86 "Land cost per area [USD/m2]";
  parameter Real pri_tower_fix_SAM = 3e6 "[H&T] Fix tower cost according to SAM cost function v2020.2.29";
  parameter Real pri_tower_scalar_exp_SAM = 0.0113 "[H&T] Scaler for tower cost according to SAM v2020.2.29";
  parameter Real pri_lift = 58.37 "Lift cost per rated mass flow per height Kevin Albrecht, 2019 https://is.gd/3VN0O7";
  parameter FI.AreaPrice pri_receiver = 37400 "Falling particle receiver cost per design aperture area Kevin Albrect, 2019 https://is.gd/3VN0O7";
  //******************************** Specific cost for the storage system
  parameter FI.MassPrice pri_particle = 1.0 "[ST] Unit cost of particles per kg";
  parameter Real pri_bin = 1230 "[ST] bin specific cost to tank surface area based on Kevin Albrect, 2019 https://is.gd/3VN0O7";
  parameter Real pri_bin_linear = 0.369 "[ST] bin linear scaler to tank temperature [USD/m^2] based on Kevin Albrect, 2019 https://is.gd/3VN0O7";
  parameter Real pri_bin_multiplier = 1 "[ST] number of each type of tank - use 1 if only have 1 hot tank and 1 cold tank";
  /*Starting from here, specific cost is for Jeremy Sment storage system cost function*/
  parameter Real pri_HRC = 850 "[ST] Specific cost of high resistance concrete [USD/m^3]";
  parameter Real pri_portland = 229 "[ST] Specific cost of portland concrete [USD/m^3]";
  parameter Real pri_RF = 2700 "[ST] Specific cost of refractory [USD/m^3]";
  parameter Real pri_filler_floor = 150 "[ST] Specific cost of floor filler material [USD/m^3]";
  parameter Real pri_roof = 457.4662 "[ST] Specific cost for roof (integrated storage)";
  parameter Real[6] c_storage = {pri_HRC, pri_portland, pri_RF, pri_filler_floor, pri_particle, pri_roof} "Specific cost of the storage system in an array";
  parameter Real D_outlet = if E_max * 2.777778e-10 < 600 then 0.21 else 0.5 "If stg capacity >= 600 MWh then D_outlet is 0.5 else 0.21 m";
  //******************************* Specific Cost parameters for Power Block Components
  parameter Real pri_recuperator = 5.2 "[PB] Specific cost of turbine (USD-K^0.8933/W^0.8933) based on Albrecht 2019 https://is.gd/3VN0O7";
  parameter Real pri_turbine = 9923.7 "[PB] Specific cost of turbine (USD/kW^0.5886) based on Albrecht 2019 https://is.gd/3VN0O7";
  parameter Real pri_compressor = 643.15 "[PB] Specific cost of compressor (USD/kW^0.9142) based on Albrecht 2019 https://is.gd/3VN0O7";
  parameter Real pri_cooler = 2.3 "[PB] Main cooler specific cost:
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         >  Based on Albrecht 2019 https://is.gd/3VN0O7 the specific cost is 76.25 (USD-K^0.8919/W^0.8919)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         >  Based on NREL sCO2 PB model used in SAM --> 2.3 USD-K/W
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          If we use UA_cooler from SAM Simulation Core sCO2 model and use Albrecht cost function, the cooler cost can reach up to 100 M.USD.
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          I believe we can just use NREL sCO2 PB model cooler cost function --> 2.3 x UA_cooler [W/K]";
  parameter Real pri_generator = 108900 "[PB] Generator cost (USD/MWe^0.5463) based on Weiland 2019 https://is.gd/uTaFkD";
  parameter Real pri_PHX_BOP_CO2 = 4753 "[PB] Primary Heat Exchanger sCO2 Line Cost (USD-s/kg) - G3P3 conversation email by Cliff 11 Nov 2020";
  parameter Real pri_PHX_BOP_s = 9153 "[PB] Primary Heat Exchanger Particle Cost [USD-s/kg] - G3P3 conversation email by Cliff 11 Nov 2020";
  parameter Real pri_PHX_per_area = 6594.5 "[PB] Primary Heat Exchanger Material+Manufacture Cost (USD/m2) - G3P3 conversation email by Cliff 11 Nov 2020";
  parameter FI.Money pri_exchanger = 150 "[PB] price of the primary exchanger in (USD/(kW_th). Value from v.9 EES sandia result c_hx";
  parameter FI.PowerPrice pri_bop = 290 / 1040 * 600 / 1000 "USD/We Balance of plant cost per gross rated power. 290--> Maximum BOP cost per MWe from SAM. 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             1040 is the maximum power block cost per MWe at SAM. 600 is the specific cost of the power block in USD/kWe according to DOE guidline";
  parameter FI.PowerPrice pri_block = 600 "sCO2 PB cost USD per kWe net based on the G3P3 Roadmap Report";
  //******************************* O&M & Washing Heliostat Specific Cost
  // Source : Heliostat Cost Reduction Study Gregory J. Kolb, page 138 Table 1
  parameter Real pri_om_name(unit = "$/W/year") = 40 / 1e3 * upstreammultiplier "Fixed O&M cost per nameplate per year";
  parameter Real pri_om_prod(unit = "$/J/year") = 0.003 / (1e6 * 3.6) * upstreammultiplier "Variable O&M cost per production per year";
  parameter Real C1_washing = 98 "[H&T] cleanliness of the mirror after 1 cleaning pass for method 1 (%)";
  parameter Real C2_washing = 96.5 "[H&T] cleanliness of the mirror after 1 cleaning pass for method 2 (%)";
  parameter Real C_target = helio_rho * 100 "annual reflectivity target";
  parameter Real R_soil = 0.49 "[H&T] reflectivity loss from the design reflectivity due to soiling";
  parameter Real P_w = 12 "[H&T] Month(s) in a year that the cleaning will be conducted";
  parameter Real omega_n = 1.5 "[H&T] Natural washing frequency (rain, snow, etc)";
  parameter Real omega_twister(fixed = false) "Supplementary washing (with machine) frequency per year";
  parameter Real pri_washing_deluge_method = 0.0027 * 1.3 "[H&T] USD/m.sq field annually. 1.3 is a factor of conversion from USD 2007 to 2020";
  parameter Real pri_washing_twister_method = 0.0076 * 1.3 "[H&T] USD/m.sq field annually. 1.3 is a factor of conversion from USD 2007 to 2020";
  parameter Real omega_deluge = 2 * omega_twister "this approach uses KJC cleaning method (1 Twister and 2 Deluge truck in between)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Source : Heliostat Cost Reduction Study Gregory J. Kolb, page 121 Table A-8";
  parameter Real pri_om_field = 52.8815449319 * A_helio ^ (-1.0359277351) "O&M field based on number of heliostat in USD / unit. The price is multiplied by 1.5 to converT it to USD 2020 from USD 2000";
  //*************************** Price of H2 storage
  parameter Real pri_H2_storage_ps = 430 "Specific cost of the H2 pipe storage USD/kg-H2";
  parameter Real pri_H2_storage_sc = 10 ^ (0.212669 * log10(M_H2_tank / 1e3) ^ 2 - 1.638654 * log10(M_H2_tank / 1e3) + 4.403100) "Specific cost of the H2 salt cavern storage USD/kg-H2";
  //***************************** Price of indirect particle gas HX (according to https://doi.org/10.1016/j.solener.2018.11.065 Fig 10)
  parameter Real pri_indirect_HX = 150 "150 USD/kWt --> according to https://doi.org/10.1016/j.solener.2018.11.065 which is based on DOE 2020 target";
  //***************************** Price of gas to gas H2 recuperator (according to https://doi.org/10.1016/j.solener.2019.01.078 )
  parameter Real pri_recup_gas_to_gas = 0;
  //1250 "Specific price of gas to gas recuperator as a function of the conductance";
  //***************************** Price of the gas burner
  parameter Real pri_gas_burner = 750.00 "USD/kWt based on the boiler price of Biomass single owner model in SAM 2020.11.29";
  //******************************* Cost of the PB components ---> product of PB initalisation
  parameter FI.Money C_HTR(fixed = false) "cost of the high temperature heat recuperator";
  parameter FI.Money C_LTR(fixed = false) "cost of the low temperature heat recuperator";
  parameter FI.Money C_turbine(fixed = false) "cost of the turbine";
  parameter FI.Money C_mainCompressor(fixed = false) "cost of the main compressor";
  parameter FI.Money C_reCompressor(fixed = false) "cost of the re compressor";
  parameter FI.Money C_exchanger(fixed = false) "cost of the exchanger";
  parameter FI.Money C_generator(fixed = false) "cost of the generator";
  parameter FI.Money C_cooler(fixed = false) "cost of the cooler";
  parameter FI.Money C_block(fixed = false) "Power block cost";
  //******************************* Cost of the indirect particle gas HX
  parameter FI.Money C_exchanger_indirect(fixed = false);
  // = if f_TES < 1e-3 then 0 else pri_indirect_HX * Q_HX_H2 / 1000 "Price of indirect particle gas HX [USD]";
  //******************************* Cost Washing calculation
  parameter FI.Money C_washing = (omega_twister * pri_washing_twister_method + omega_deluge * pri_washing_deluge_method) * A_field "Washing cost [USD/year]";
  parameter FI.Money C_om_field = pri_om_field * A_field "OnM field exclude washing cost [USD/year]";
  //******************************* Cost of solar field calculation
  parameter FI.Money C_field = if CSP_fraction < 1e-3 then 0 else A_field * pri_field * upstreammultiplier "Field cost";
  parameter FI.Money C_site = if CSP_fraction < 1e-3 then 0 else A_field * pri_site * upstreammultiplier "Site improvements cost";
  parameter FI.Money C_land = if CSP_fraction < 1e-3 then 0 else A_land * pri_land * upstreammultiplier "Land cost";
  parameter FI.Money C_field_total = C_field + C_site "Heliostat field plus site preparation costs";
  //******************************* Cost of tower sub-system (receiver + tower + receiver lift)
  //******************************* As per December 7 2020, the tower cost function is changed to the Latest Tower Cost Function
  parameter Real Euro_to_USD_exchange_rate = 1.09 "[USD/Euro] Google, accessed 8 April 2022";
  parameter Real AUD_to_USD_exchange_rate = 0.75 "[USD/AUD] Google, accessed 8 April 2022";
  parameter FI.Money C_extra_structure(fixed = false);
  parameter FI.Money C_tower_absolute = 83060926 "Absolute tower cost [USD]";
  /*Latest Tower Cost Function Based on the email by J.Sment (Sandia) Sat 05/12/2020 05:48 */
  parameter FI.Money C_tower = if CSP_fraction < 1e-3 then 0 else if set_SAM_tower_cost then C_extra_structure - 1.992 * H_tower_final ^ 2.747 + 523100 + pri_tower_fix_SAM * Modelica.Math.exp(pri_tower_scalar_exp_SAM * (H_tower_final + 0.5 * H_helio - H_rcv / 2)) - 28000 * Euro_to_USD_exchange_rate * H_tower_final + 1573 * H_tower_final else C_extra_structure - 1.992 * H_tower_final ^ 2.747 + 523100 + (0.7452 * H_tower_final ^ 3 - 148.25 * H_tower_final ^ 2 + 37204 * H_tower_final - 731236) * Euro_to_USD_exchange_rate + 1573 * H_tower_final "Cost of tower based on J.Sment (Sandia) email to G3P3 Team at Sat 05/12/2020 05:48
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              > Tim Harvey structure only cost model is a function of tower height [H_tower] and maximum particle mass in one storage tank [m_max]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    - Regression model for Tim Harvey cost : 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      ----> online tool https://stats.blue/Stats_Suite/multiple_linear_regression_calculator.html:
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      C_harvey = 2293496.5853409-45954.7293032756*H_tower+
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0.1048843661*m_max+256.311306896*H_tower^2+0.0015436937*m_max*H_tower-0.0000000021*m_max^2 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              > The Upper Boundary cost is the one with SAM cost function:
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  C_tower = Tim Harvey cost [USD] - SBP Material Cost [USD] + SAM Tower Cost [USD] - Piping Cost [Euro] * USD_to_Euro + Ducting cost [USD]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              > The Lower Boundary cost:
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  C_tower = Tim Harvey cost [USD] - SBP Material Cost [USD] + SBP Tower Cost (no pipe) [USD] - Ducting cost [USD]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              As per December 7 2020, the tower cost function is changed to the Latest Tower Cost Function";
  //*********************************** Evaluating tower cost using SAM tower correlation - piping cost + ducting cost + extra structure cost
  //*********************************** Based on the email by J.Sment (Sandia) Wed 09/12/2020 19:35
  /*C_tower = Tim Harvey structure only cost [USD]- SBP Materials [USD]+ Sam Tower Cost [USD]- Piping Cost [in Euro] + Ducting cost [USD]*/
  //==================> SBP Material Cost
  //==========> SAM Cost v2020.2.29
  //================> Piping cost
  //==================> Ducting Cost
  //*********************************** Evaluating tower cost using SAM tower correlation - piping cost + ducting cost + extra structure cost
  //*********************************** Based on the email by J.Sment (Sandia) Wed 09/12/2020 19:35
  /*C_tower = Tim Harvey structure only cost [USD]- SBP Materials [USD]+ SBP Tower Cost without Piping [Euro] + Ducting cost [USD]*/
  //==================> SBP Material Cost
  //======> SBP no pipe cost
  //======> Ducting Cost
  parameter FI.Money C_fpr = if CSP_fraction < 1e-3 then 0 else pri_receiver * A_rcv "Falling particle receiver cost";
  parameter FI.Money C_lift_rec = if CSP_fraction < 1e-3 then 0 else pri_lift * dh_liftRC * m_flow_fac "Receiver lift cost";
  parameter FI.Money C_receiver = if CSP_fraction < 1e-3 then 0 else if set_absolute_tower_cost == true then (C_fpr + C_tower_absolute + C_lift_rec) * upstreammultiplier else (C_fpr + C_tower + C_lift_rec) * upstreammultiplier "Total receiver sub-system cost";
  //******************************* Cost of storage sub-system (bins + cold tank lift + particles + PHX lift + insulation)
  parameter FI.Money C_lift_cold = if use_TES == false then 0 else if set_external_storage then pri_lift * dh_LiftCold * m_flow_blk else 0 "Cold storage tank lift cost";
  //******************************* Storage bin cost calculation based on Kevin Albrect, 2019 https://is.gd/3VN0O7
  parameter FI.Money C_bins = if use_TES == false then 0 else if set_new_storage_calc then 750 * CN.pi * (D_storage - 1 - 1) * H_storage else pri_bin_multiplier * (pri_bin + pri_bin_linear * (T_hot_set - 600) / 400) * SA_storage + pri_bin_multiplier * (pri_bin + pri_bin_linear * (T_cold_set - 600) / 400) * SA_storage;
  //******************************* Storage bin cost calculation based on Jeremy Sment of Sandia study (g3p3 project)
  parameter FI.Money C_bins_dome = if use_TES == false then 0 else if set_external_storage then SolarTherm.Utilities.G3P3StorageCostFunction(Th_refractory_hot_tank, H_storage, D_storage, D_outlet, t_storage, m_max, c_storage) + SolarTherm.Utilities.G3P3StorageCostFunction(Th_refractory_cold_tank, H_storage, D_storage, D_outlet, t_storage, m_max, c_storage) else SolarTherm.Utilities.G3P3StorageCostFunction_Integrated(Th_refractory_cold_tank, H_tower_final, R_tower * 2, m_max, D_outlet, packing_factor, c_storage) + SolarTherm.Utilities.G3P3StorageCostFunction_Integrated(Th_refractory_hot_tank, H_tower_final, R_tower * 2, m_max, D_outlet, packing_factor, c_storage) "Storage bin for dome storage --> based on the type of storage";
  parameter FI.Money C_insulation = if use_TES == false then 0 else if U_value_hot_tank == 0 and U_value_cold_tank == 0 then 0 else SA_storage * (131.0426 / U_value_hot_tank + 23.18) + SA_storage * (131.0426 / U_value_cold_tank + 23.18) "Insulation cost based on Kevin Albrecht 2019 https://is.gd/3VN0O7";
  parameter FI.Money C_particles = if use_TES == false then 0 else (1 + NS_particle) * pri_particle * m_max "Cost of particles";
  parameter FI.Money C_lift_hx = if abs(P_net - P_net_default_value) < 1 then 0 else if set_external_storage then pri_lift * dh_liftHX * m_flow_blk else 0 "Heat exchanger lift cost";
  /******************************************************************************************************
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        FIXME: There are 2 u_values now, implement it in the tuffcrete x microporous analysis
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        (131.0426 / U_value + 23.18) ======> cost function insulation of Tuffcrete, Microporous and Concrete
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        (873.11/U_value) - 322.202 ======> cost function insulation of Tuffcrete, Pumplite60 and Concrete
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0.03293006 / U_value + 0.01518 =====> thickness function of Pumplite60;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0.32368 / U_value - 0.146096   =====> thickness function of Microporous;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        parameter SI.Length t_mp = 0.32368 / (U_value_hot_tank + U_value_cold_tank) - 0.146096;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        parameter SI.Length t_tuffcrete47 = 0.01;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ******************************************************************************************************/
  //parameter FI.Money C_storage = if abs(P_net - P_net_default_value) < 1 or set_no_thermal then 0 else if set_dome_storage then (C_bins_dome + C_particles + C_lift_hx + C_lift_cold + 0 + f_loss * t_life * pri_particle * 0) * upstreammultiplier else (C_bins + C_particles + C_lift_hx + C_lift_cold + C_insulation + f_loss * t_life * pri_particle * 1.753e10) * upstreammultiplier "Total storage cost. Dome storage bin cost calculation already considers insulation (refractory) s.t. C_insulation = 0";
  parameter FI.Money C_storage = if use_TES == false then 0 else if set_dome_storage then (C_bins_dome + C_particles + C_lift_hx + C_lift_cold + 0 + f_loss * t_life * pri_particle * 0) * upstreammultiplier else (C_bins + C_particles + C_lift_hx + C_lift_cold + C_insulation + f_loss * t_life * pri_particle * 1.753e10) * upstreammultiplier "Total storage cost. Dome storage bin cost calculation already considers insulation (refractory) s.t. C_insulation = 0";
  //******************************* Cost of H2 storage
  parameter FI.Money C_storage_H2 = if t_storage_H2_final < t_storage_H2_threshold then 0 else min(pri_H2_storage_sc * M_H2_tank, pri_H2_storage_ps * M_H2_tank) "Cost of the hydrogen storage";
  //******************************* Cost of H2 recuperators (gas to gas HX)
  parameter Real A_recup_ft = A_recup * 10.7639 "Recuperator area in ft^2";
  parameter Real p_recup_psi = 7e5 * 0.000145038 "Pressure of the recuperator in psi";
  parameter Real Fp = 0.9083 + 0.018 * (p_recup_psi / 100) + 0.0017 * (p_recup_psi / 100) ^ 2 "Seider et al p.571";
  parameter Real Fl = 1 "Assuming tube length > 20 ft Seider et al p.571";
  parameter Real Fm = 2.70 + (A_recup_ft / 100) ^ 0.07 "Assuming stainless steel material Seider et al p.571";
  parameter FI.Money C_recuperator_H2_base = exp(11.0545 - 0.9228 * log(A_recup_ft) + 0.09005 * log(A_recup_ft) ^ 2) "Price of the H2 recuperators assuming Fixed-Head S&T HX Seider et al p.571";
  parameter FI.Money C_recuperator_H2 = C_recuperator_H2_base * Fp * Fm * Fl;
  //******************************* Price of the H2 burner
  parameter FI.Money C_burner = pri_gas_burner * Q_burner_emergency / 1e3 "Price of gas burner";
  //******************************* Cost of BOP
  parameter FI.Money C_bop = if abs(P_net - P_net_default_value) < 1 then 0 else P_gross * pri_bop * upstreammultiplier "Balance of plant cost";
  parameter FI.Money C_prod = if abs(P_net - P_net_default_value) < 1 then 0 else pri_om_prod "Variable O&M cost per production per year";
  //******************************* PV capital and OM cost
  parameter Real pri_PV = 340 "TIC cost of PV per $/kWe 2030 [340 - 834] https://www.irena.org/-/media/Files/IRENA/Agency/Publication/2019/Nov/IRENA_Future_of_Solar_PV_2019.pdf";
  parameter Real pri_om_PV = 12.682 "OnM cost for PV in USD/kWe installed capacity in 2030 ---> depends on the capacity factor https://ec.europa.eu/energy/sites/ener/files/documents/2018_06_27_technology_pathways_-_finalreportmain2.pdf page 45";
  parameter FI.Money C_PV = PV_Target / 1e3 * pri_PV "PV cost in $";
  parameter FI.MoneyPerYear C_year_PV = 0 "Added later in the post processing since the PV OM depends on the PV CF";
  //pri_om_PV * PV_Target / 1e3 "Fixed OM cost for PV";
  //******************************* Industrial electric heater cost
  parameter Real pri_heater = 140 "Euro /KWth 2030 [UB] https://doi.org/10.1080/15567249.2020.1843565 table 1. Other values are 140 USD/kWe TIC https://doi.org/10.1016/j.energy.2020.118472; 220 USD/kWe https://doi.org/10.1016/j.enconman.2020.113779";
  parameter FI.Money C_heater = if P_heater <= 0 or set_no_thermal then 0 else pri_heater * P_heater / 1e3 "Price of electric heater using scaling formula [USD]";
  //******************************* Capital cost of SMR
  parameter Real pri_SMR = 117232000 * Euro_to_USD_exchange_rate "cost of SMR component for 'Standalone (Merchant) H2 plant' as per https://ieaghg.org/exco_docs/2017-02.pdf.         
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      It is scalled using scaler with exchange rate 1.1 USD/euro (Google, accessed on 29 March 2022)";
  parameter Real scaler_n = 0.7;
  parameter FI.Money C_SMR = if set_SMR_always_off then 0 else pri_SMR * (H2_mdot_target / (8994 / 3600)) ^ scaler_n;
  //******************************* Captial and OM cost of Electrolyser
  parameter Real pri_electrolyser = 317 "Electrolyser price 2030 in (317 - 817) USD/kWe https://doi.org/10.3390/en14123437";
  parameter Real pri_om_electrolyser = 0.035 "Fraction of electrolyser OnM based on C_electrolyser https://doi.org/10.3390/en14123437";
  parameter FI.Money C_electrolyser = if CSP_fraction < 1e-3 and PV_fraction < 1e-3 then 0 else P_hybrid_system_final / 1000 * pri_electrolyser "cost of electrolyser";
  parameter FI.MoneyPerYear C_year_electrolyser = C_electrolyser * pri_om_electrolyser "Fixed OM cost for electrolyser in USD/year";
  //******************************* Cost per kg of Natural gas consumed by SMR
  parameter Real pri_natural_gas = 0.279 "Cost of natural gas USD per kg. Cost is 8 AUD/GJ (https://doi.org/10.1016/j.ijhydene.2021.04.104) NSW price, LHV is 46.5 MJ/kg (https://ieaghg.org/exco_docs/2017-02.pdf)";
  //******************************* Cost per kg of water consumed by SMR
  parameter Real pri_water_SMR = 0.00375 "Cost of processed water per kg consumed by Electrolyser https://doi.org/10.1016/j.ijhydene.2021.04.104";
  //******************************* Cost per kg of water consumed by Electrolyser
  parameter Real pri_water_ele = 0.00375 "Cost of processed water per kg consumed by Electrolyser https://doi.org/10.1016/j.ijhydene.2021.04.104";
  //******************************* Cost per kg of H2 production
  parameter Real pri_h2 = 0.0 "Variable cost to produce H2 $/kg";
  //******************************* Cost of Carbon
  parameter Real pri_carbon_tax = 78.5 * Euro_to_USD_exchange_rate "Euro per kg Carbon";
  //******************************* Cost of H2 conditioning -- compressing and expansion (https://doi.org/10.1115/ES2019-3893)
  parameter FI.Money C_compressor_H2 = 643.15 * (P_compressor_H2 / 1e3) ^ 0.9142;
  parameter FI.Money C_expander_H2 = 0;
  //0 * P_expander_H2 ^ 0.5886 "Screw -- not significant and can be neglected";
  parameter SI.Power P_compressor_H2(fixed = false);
  parameter SI.Power P_expander_H2(fixed = false);
  //******************************* Cost of O&M (fixed + varied)
  //parameter FI.MoneyPerYear C_year = if set_detail_field_om then P_name * pri_om_name + C_om_field + C_washing + C_year_PV + C_year_electrolyser else P_name * pri_om_name + C_year_PV + C_year_electrolyser "Fixed O&M cost per year (PV, CSP, Electrolyser)";
  /*****************************************************  	
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            	Specific cost comes from https://www.nrel.gov/docs/fy04osti/34440.pdf
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            *****************************************************/
  parameter Real labour_PB = 15 "Number of labour for operating PB (constant at 25 according to Sargent and Lundy)";
  parameter Real labour_mirror = if A_field >= 231000 then ceil(1.75965439123334 * (A_field / 1e6) ^ 2 + 9.32364996575522 * A_field / 1e6 + 5.7523399399) else max(A_field / 231000 * 8, 1);
  parameter SI.Volume heliostat_cleaning_water = 0.022 * A_field "Water consumed to clean the heliostat per year";
  parameter Real pri_om_labour = if on_CSP == false then 71000 else if A_field < 71100 then 71000 * upstreammultiplier else max((2.11549413053375e-9 * A_field ^ 2 - 0.014387055541508 * A_field + 65188.593) * upstreammultiplier, 42000 * upstreammultiplier);
  parameter Real pri_om_helio_wear_and_tear = 0.243 * upstreammultiplier "USD/m2, wear and tear of the heliostat";
  parameter Real pri_om_helio_utilities = if A_field >= 231000 then (1.0535e-7 * A_field + 0.09432704889) * upstreammultiplier else A_field / 231000 * 0.12 * upstreammultiplier "Utilities and eequipment cost to operate heliostat, spread over 5 years of equipment lifetime";
  parameter Real pri_om_cleaning_water = 0.32 * upstreammultiplier "Cleaning water cost per m3";
  parameter FI.MoneyPerYear C_year_heliostat_labour = if CSP_fraction < 1e-3 then 0 else labour_mirror * pri_om_labour * 1.58 * upstreammultiplier;
  parameter FI.MoneyPerYear C_year_heliostat_wear_and_tear = if CSP_fraction < 1e-3 then 0 else pri_om_helio_wear_and_tear * A_field * 289.4 / 132;
  parameter FI.MoneyPerYear C_year_heliostat_equipment = if CSP_fraction < 1e-3 then 0 else pri_om_helio_utilities * t_life / 5 * A_field * 289.4 / 132;
  parameter FI.MoneyPerYear C_year_heliostat_cleaning_water = if CSP_fraction < 1e-3 then 0 else pri_om_cleaning_water * heliostat_cleaning_water;
  parameter FI.MoneyPerYear C_year_pb_labour = labour_PB * pri_om_labour * 1.58;
  parameter FI.MoneyPerYear C_year = C_year_PV + C_year_electrolyser + C_year_heliostat_labour + C_year_heliostat_wear_and_tear + C_year_heliostat_equipment + C_year_heliostat_cleaning_water + C_year_pb_labour;
  //******************************* Total cost calculation
  parameter FI.Money C_cap_CSP(fixed = false) "Equipment cost of the CSP plant";
  parameter FI.Money C_cap_PV(fixed = false) "Equipment cost of the PV plant";
  parameter FI.Money C_cap_SMR(fixed = false) "Equipment cost of the SMR";
  parameter FI.Money C_cap_electrolyser(fixed = false) "Equipment cost of the electrolyser";
  parameter FI.Money C_cap_total(fixed = false) "equipment cost";
  parameter FI.Money C_direct(fixed = false) "Direct capital costs";
  parameter FI.Money C_indirect(fixed = false) "Indirect capital costs";
  parameter FI.Money C_cap(fixed = false) "Capital costs";
  //********************* System Components Instantiation
  //********************* Weather data
  SolarTherm.Models.Sources.DataTable.DataTable data(lon = lon, lat = lat, t_zone = t_zone, year = year, file = wea_file) annotation(
    Placement(visible = true, transformation(origin = {-129, -40}, extent = {{-7, -8}, {7, 8}}, rotation = 0)));
  //********************* DNI_input
  Modelica.Blocks.Sources.RealExpression DNI_input(y = data.DNI) annotation(
    Placement(visible = true, transformation(origin = {-127, 70}, extent = {{-13, -10}, {13, 10}}, rotation = 0)));
  //********************* Weather inputs
  Modelica.Blocks.Sources.RealExpression Tamb_input(y = data.Tdry) annotation(
    Placement(visible = true, transformation(origin = {139, 80}, extent = {{19, -10}, {-19, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Wspd_input(y = data.Wspd) annotation(
    Placement(transformation(extent = {{-140, 20}, {-114, 40}})));
  Modelica.Blocks.Sources.BooleanExpression always_on(y = on_CSP) annotation(
    Placement(visible = true, transformation(origin = {-128, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Wind_dir(y = data.Wdir) annotation(
    Placement(visible = true, transformation(origin = {-129, 51}, extent = {{-11, -13}, {11, 13}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Pres_input(y = data.Pres) annotation(
    Placement(transformation(extent = {{76, 18}, {56, 38}})));
  SolarTherm.Utilities.WspdScaler wspdScaler(H_tower = H_tower_final) annotation(
    Placement(visible = true, transformation(origin = {-51, 45}, extent = {{-7, -7}, {7, 7}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression DHI_input(y = data.DHI) annotation(
    Placement(visible = true, transformation(origin = {-125, -68}, extent = {{-13, -10}, {13, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression Albedo_input(y = data.Albedo) annotation(
    Placement(visible = true, transformation(origin = {-125, -88}, extent = {{-13, -10}, {13, 10}}, rotation = 0)));
  //********************* Parasitic inputs
  Modelica.Blocks.Sources.RealExpression parasities_input(y = if set_external_storage then heliostatsField.W_loss + liftHX.W_loss + liftRC.W_loss + tankHot.W_loss + tankCold.W_loss else heliostatsField.W_loss + 0 + liftRC.W_loss + tankHot.W_loss + tankCold.W_loss) annotation(
    Placement(visible = true, transformation(origin = {147, 58}, extent = {{-13, -10}, {13, 10}}, rotation = -90)));
  //*********************  Or block for defocusing
  Modelica.Blocks.Logical.Or or1 annotation(
    Placement(transformation(extent = {{-102, 4}, {-94, 12}})));
  //********************* Sun Model
  SolarTherm.Models.Sources.SolarModel.Sun sun(lon = lon, lat = lat, t_zone = t_zone, year = year, redeclare function solarPosition = Models.Sources.SolarFunctions.PSA_Algorithm) annotation(
    Placement(transformation(extent = {{-82, 60}, {-62, 80}})));
  //********************* Solar field
  SolarTherm.Models.CSP.CRS.HeliostatsField.HeliostatsFieldSolstice_1stApproach heliostatsField(lon = lon, lat = lat, ele_min(displayUnit = "deg") = ele_min, use_wind = set_use_wind, t_start = t_start, Wspd_max = Wspd_max, he_av = he_av_design, use_on = true, use_defocus = true, A_h = A_helio, nu_defocus = nu_defocus, nu_min = nu_min_sf, Q_design = Q_flow_defocus, nu_start = nu_start, Q_in_rcv = Q_in_rcv, H_rcv = H_rcv, W_rcv = W_rcv, tilt_rcv = tilt_rcv, field_type = field_type, W_helio = W_helio, H_helio = H_helio, H_tower = H_tower_final, R_tower = R_tower, R1 = R1, fb = fb, helio_refl = helio_refl, slope_error = slope_error, slope_error_windy = slope_error_windy, n_row_oelt = n_row_oelt, n_col_oelt = n_col_oelt, n_procs = n_procs, psave = casefolder, wea_file = wea_file, sunshape = sunshape, buie_csr = buie_csr, optics_verbose = set_optics_verbose, optics_view_scene = set_optics_view_scene, get_optics_breakdown = get_optics_breakdown, set_swaying_optical_eff = set_swaying_optical_eff, on_CSP = on_CSP, use_g3p3_field = use_g3p3_field) annotation(
    Placement(transformation(extent = {{-88, 2}, {-56, 36}})));
  //********************* Washing calculator
  SolarTherm.Models.CSP.CRS.HeliostatsField.WashingFrequencyCalculator washingFrequencyCalculator(C_tw = C1_washing, C_dl = C2_washing, R_soil = R_soil, P_w = P_w, omega_n = omega_n, C_target = C_target) annotation(
    Placement(visible = true, transformation(origin = {-130, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //********************* Surrogate Receiver
  SolarTherm.Models.CSP.CRS.Receivers.ParticleReceiver particleReceiver(H_drop_design = H_rcv, use_neural_network = set_use_neural_network, redeclare package Medium = Medium, X_max = X_max_rcv, X_min = X_min_rcv, inputsize = inputsize_rcv, out_max = y_max_rcv, out_min = y_min_rcv, saved_model_dir = saved_model_dir_rcv, T_0 = T_cold_set, T_out = T_hot_set, m_flow_rec_max = m_flow_rec_max, ar_rec = ar_rec) annotation(
    Placement(visible = true, transformation(origin = {-26, 30}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  //********************* Receiver Design Condition Sizing Calculator - Using expensive model
  SolarTherm.Models.CSP.CRS.Receivers.ParticleReceiver1DCalculator particleReceiver1DCalculator(Q_in = Q_in_rcv, T_out_design = T_in_ref_blk, T_in_design = T_in_rec, T_amb_design = T_amb_des_rcv, Wspd_design = 0, Wspd_dir = 0, dni_des = dni_des, h_conv_backwall = h_conv_backwall, h_conv_curtain = h_conv_curtain, phi_max = packing_factor, eps_w = eps_w, th_w = th_w, k_w = k_w, F = F, d_p = d_p, cp_s = cp_s, rho_s = rho_s, eps_s = eps_s, abs_s = abs_s, ar_rec = ar_rec, with_detail_h_ambient = true, with_wind_effect = true, test_mode = false, fixed_geometry = false, iterate_Q_flow = false, with_iterate_mdot = false, with_pre_determined_eta = false, N = discretisations, with_catch_and_release_mechanism = false) annotation(
    Placement(visible = true, transformation(origin = {150, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //********************* Defocus calculator
  SolarTherm.Models.CSP.CRS.Receivers.ParticleReceiver1DCalculator defocuscalculator(T_out_design = T_in_ref_blk, T_in_design = T_in_rec, T_amb_design = T_amb_des_rcv, Wspd_design = 0, Wspd_dir = 0, dni_des = dni_des, h_conv_backwall = h_conv_backwall, h_conv_curtain = h_conv_curtain, phi_max = packing_factor, eps_w = eps_w, th_w = th_w, k_w = k_w, F = F, d_p = d_p, cp_s = cp_s, rho_s = rho_s, eps_s = eps_s, abs_s = abs_s, ar_rec = ar_rec, with_detail_h_ambient = true, with_wind_effect = true, fixed_geometry = true, iterate_Q_flow = true, with_iterate_mdot = false, m_design = m_flow_blk, H_drop_design = H_rcv, with_pre_determined_eta = false, N = discretisations, with_catch_and_release_mechanism = false) annotation(
    Placement(visible = true, transformation(origin = {12, 130}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //********************* Hot tank
  //SolarTherm.Models.Storage.Tank.Tank_Particle tankHot(redeclare package Medium = Medium, D = D_storage, H = H_storage, dome_storage = set_dome_storage, Th_refractory = Th_refractory_hot_tank, H_tower = H_tower, D_tower = 2 * R_tower, m_tot = m_max, D_outlet = D_outlet, t_storage = t_storage, T_start = T_hot_start, L_start = (1 - split_cold) * 100, use_p_top = set_tnk_use_p_top, enable_losses = set_tnk_enable_losses, use_L = true, W_max = 0, T_set = 0, U_value = U_value_hot_tank, packing_factor = packing_factor) annotation(
  //Placement(transformation(extent = {{16, 54}, {36, 74}})));
  SolarTherm.Models.Storage.Tank.Tank_Particle_3_outlet tankHot(redeclare package Medium = Medium, D = D_storage, H = H_storage, dome_storage = set_dome_storage, Th_refractory = Th_refractory_hot_tank, H_tower = H_tower_final, D_tower = 2 * R_tower, m_tot = m_max, D_outlet = D_outlet, t_storage = t_storage, T_start = T_hot_start, L_start = (1 - split_cold) * 100, use_p_top = set_tnk_use_p_top, enable_losses = set_tnk_enable_losses, use_L = true, W_max = 0, T_set = 0, U_value = U_value_hot_tank, packing_factor = packing_factor) annotation(
    Placement(visible = true, transformation(origin = {27, 64}, extent = {{-11, -10}, {11, 10}}, rotation = 0)));
  //********************* Cold tank
  SolarTherm.Models.Storage.Tank.Tank_Particle_2_outlet tankCold(redeclare package Medium = Medium, D = D_storage, H = H_storage, dome_storage = set_dome_storage, Th_refractory = Th_refractory_cold_tank, H_tower = H_tower_final, D_tower = 2 * R_tower, m_tot = m_max, D_outlet = D_outlet, T_start = T_cold_start, L_start = split_cold * 100, use_p_top = set_tnk_use_p_top, enable_losses = set_tnk_enable_losses, use_L = true, W_max = 0, T_set = 0, U_value = U_value_cold_tank, packing_factor = packing_factor, t_storage = t_storage) annotation(
    Placement(visible = true, transformation(origin = {54, -18}, extent = {{8, -10}, {-8, 10}}, rotation = 0)));
  //********************* Tower Inner Diameter Calculator - for integrated storage concept only
  SolarTherm.Utilities.TowerInnerDiameterCalculator_IntegratedG3P3Storage towerInnerDiameterCalculator(rho_particle = (rho_cold_set + rho_hot_set) / 2, H_tower = H_tower_final, m_max = m_max, Th_refractory = if Th_refractory_cold_tank == Th_refractory_hot_tank then Th_refractory_cold_tank else max(Th_refractory_cold_tank, Th_refractory_hot_tank), D_outlet = D_outlet, phi = packing_factor);
  //********************* Receiver lift
  SolarTherm.Models.Fluid.Pumps.LiftSimple liftRC(redeclare package Medium = Medium, cont_m_flow = true, use_input = true, dh = dh_liftRC, CF = 0, eff = eff_lift) annotation(
    Placement(visible = true, transformation(origin = {-1, -27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  //********************* Simple Electrical Heater
  SolarTherm.Models.UtilitiesComponent.SimpleElectricalHeater simpleElectricalHeater(W_dumped = if P_heater > 1 then P_heater else 0, cold_tnk_empty_lb = cold_tnk_crit_lb, cold_tnk_empty_ub = cold_tnk_crit_ub, eta = eta_heater) annotation(
    Placement(visible = true, transformation(origin = {16, -190}, extent = {{-44, -44}, {44, 44}}, rotation = 0)));
  //********************* Electric heater lift
  //********************* Heat exchanger lift
  SolarTherm.Models.Fluid.Pumps.LiftSimple liftHX(redeclare package Medium = Medium, cont_m_flow = true, use_input = true, dh = dh_liftHX, CF = 0, eff = eff_lift) annotation(
    Placement(visible = true, transformation(origin = {76, 42}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  //********************* Cold storage tank lift
  SolarTherm.Models.Fluid.Pumps.LiftSimple LiftCold(redeclare package Medium = Medium, cont_m_flow = false, use_input = false, dh = dh_LiftCold, CF = 0, eff = eff_lift) annotation(
    Placement(visible = true, transformation(origin = {106, -24}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  //********************* Temperature sensor
  //********************* PowerBlockControl
  SolarTherm.Models.Control.PowerBlockControl_PVCSP_Particle_H2O2Storage controlHot(m_flow_on = m_flow_blk, L_on = hot_tnk_empty_ub, L_off = hot_tnk_empty_lb, L_df_on = hot_tnk_full_ub, L_df_off = hot_tnk_full_lb, logic.dispatch_optimiser = set_dispatch_optimiser, logic.set_scheduler = set_scheduler, P_net = P_hybrid_system_final, logic.upstreammultiplier = upstreammultiplier, CSP_name_plate = P_net, Q_HX_des = Q_flow_des, eta_gross_base = eff_blk, t_storage = t_storage, logic.P_net = P_hybrid_system_final, optimal_dispatch_dual_tank = set_dispatch_optimiser_dual_tank) annotation(
    Placement(visible = true, transformation(extent = {{46, 76}, {58, 62}}, rotation = 0)));
  //********************* Power block
  SolarTherm.Models.PowerBlocks.sCO2PB_ConstantEfficiency powerBlock(P_gross = P_gross, T_in_ref_blk = T_hot_set, p_high = p_high, PR = PR, pinch_PHX = pinch_exchanger, dTemp_HTF_PHX = dTemp_HTF_PHX, T_amb_base = blk_T_amb_des, htf_choice = htf_choice, dT_PHX_hot_approach = dT_PHX_hot_approach, dT_PHX_cold_approach = dT_PHX_cold_approach, eta_isen_mc = eta_comp_main, eta_isen_rc = eta_comp_re, eta_isen_t = eta_turb, dT_mc_approach = dT_mc_approach, which_PB_model = which_PB_model, load_base = 1, eta_gross_base = eta_gross_base, eta_Q_base = eta_Q_base, Q_HX_des = Q_flow_des, m_HTF_des = m_flow_blk, base_path = base_path, SolarTherm_path = SolarTherm_path, inputsize = inputsize_PB, outputsize = outputsize_PB, tolerance_kriging = tolerance_kriging, tolerance_ANN = tolerance_ANN, which_surrogate = which_surrogate, test_mode = false, eta_motor = 1, f_fixed_load = f_fixed_load, external_parasities = set_external_parasities, P_net = P_net, P_net_default_value = P_net_default_value, nu_min = nu_min_blk, T_cold_set = T_cold_set) annotation(
    Placement(transformation(extent = {{88, 4}, {124, 42}})));
  //*********************Power Block Calculator
  //SolarTherm.Models.PowerBlocks.sCO2PBCalculator_Using_JPidea sCO2PBDesignPointCalculator(redeclare package Medium = Medium, P_net = P_net, T_in_ref_blk = T_in_ref_blk, p_high = p_high, PR = PR, pinch_PHX = pinch_exchanger, dTemp_HTF_PHX = dTemp_HTF_PHX, T_HTF_in = T_in_ref_blk, T_amb_input = blk_T_amb_des, load = 1, f_fixed_load = f_fixed_load, blk_T_amb_des = blk_T_amb_des, T_low = T_low, nu_min_blk = nu_min_blk, N_exch_parameter = N_exch_parameter, N_LTR_parameter = N_LTR_parameter, pri_recuperator = pri_recuperator, pri_turbine = pri_turbine, pri_compressor = pri_compressor, pri_cooler = pri_cooler, pri_generator = pri_generator, pri_exchanger = pri_exchanger, eta_motor = 1, pinch_recuperator = pinch_recuperator, par_fr = par_fr, test_mode = true, external_parasities = set_external_parasities) annotation(
  //  Placement(visible = true, transformation(origin = {-176, 110}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  //********************* PV Array
  //********************* Electrolyser
  SolarTherm.Models.Electrochemical.Simple_Electrolyser electrolyser(P_electro_requested = P_hybrid_system_final, with_storage = true, upstreammultiplier = upstreammultiplier, set_no_thermal = set_no_thermal) annotation(
    Placement(visible = true, transformation(origin = {218, -56}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  //Steam methane reformer
  SolarTherm.Models.ChemicalComponent.Simple_SMR SMR(CH4_reaction_extent = SMR_reaction_conversion, H2_mol_target = H2_mol_target, H2_mdot_target = H2_mdot_target, W_consumption = W_consumption_SMR, set_SMR_always_off = set_SMR_always_off) annotation(
    Placement(visible = true, transformation(origin = {131, -143}, extent = {{-33, -33}, {33, 33}}, rotation = 0)));
  //********************* Price
  SolarTherm.Models.Analysis.Market market(redeclare model Price = Models.Analysis.EnergyPrice.Table(file = pri_file)) annotation(
    Placement(visible = true, transformation(extent = {{196, 34}, {216, 54}}, rotation = 0)));
  SolarTherm.Models.Sources.Schedule.Scheduler sch if not set_const_dispatch;
  //********************* Tim Harvey Extra Cost for the tower structure
  SolarTherm.Utilities.TowerExtraCostTimHarvey structureExtraCost(H_tower = H_tower_final, D_inner_tower = R_tower * 2);
  //********************* HX and industrial heat
  Modelica.Blocks.Sources.BooleanExpression booleanExpression(y = false) annotation(
    Placement(visible = true, transformation(origin = {-128, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Sum sum1(nin = 2) annotation(
    Placement(visible = true, transformation(origin = {152, -56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.LiftSimple lift(cont_m_flow = true, use_input = true, dh = 10, CF = 0, eff = 1, redeclare package Medium = Med_PCL) annotation(
    Placement(visible = true, transformation(origin = {190, 100}, extent = {{16, -16}, {-16, 16}}, rotation = 0)));
  SolarTherm.Models.Fluid.HeatExchangers.HeatExchanger_ParticleGas heatExchanger_ParticleGas(W_HX = W_HX, H_HX = H_HX, level_on = hot_tnk_empty_ub, level_off = hot_tnk_empty_lb, T_in_PCL_DP = T_hot_set, T_out_PCL_DP = T_cold_set, T_in_AIR_DP = T_in_AIR_DP, T_out_AIR_DP = T_out_AIR_DP, m_dot_AIR_DP = m_dot_AIR_DP, HX_always_off = HX_always_off) annotation(
    Placement(visible = true, transformation(origin = {240, 16}, extent = {{-26, -26}, {26, 26}}, rotation = -90)));
  //********************* Instantiate the optical, price and weather table, for looping for the next forecast horizon
  Modelica.Blocks.Types.ExternalCombiTable2D opt_eff_calm = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = "optics", fileName = opt_file, table = fill(0.0, 0, 2), smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
  Modelica.Blocks.Types.ExternalCombiTable2D opt_eff_windy = Modelica.Blocks.Types.ExternalCombiTable2D(tableName = if set_swaying_optical_eff == true then "optics_windy" else "optics", fileName = opt_file, table = fill(0.0, 0, 2), smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
  Modelica.Blocks.Types.ExternalCombiTable1D wea_table = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "data", fileName = wea_file, table = fill(0.0, 0, 2), columns = 1:11, smoothness = Modelica.Blocks.Types.Smoothness.ContinuousDerivative);
  Modelica.Blocks.Types.ExternalCombiTable1D pri_table = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "prices", fileName = pri_file, table = fill(0.0, 0, 2), columns = 1:2, smoothness = Modelica.Blocks.Types.Smoothness.ConstantSegments);
  Modelica.Blocks.Math.Sum sumH2(nin = 2) annotation(
    Placement(visible = true, transformation(origin = {298, -116}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //********************* Variables
  SI.Power P_elec "Net output power of power block";
  SI.Energy E_elec(start = 0, fixed = true, displayUnit = "MW.h") "Generate electricity";
  FI.Money R_spot(start = 0, fixed = true) "Spot market revenue";
  SI.Energy E_resource(start = 0, fixed = true);
  SI.Energy E_resource_after_optical_eff(start = 0, fixed = true);
  SI.Energy E_helio_incident(start = 0, fixed = true);
  SI.Energy E_helio_raw(start = 0, fixed = true);
  SI.Energy E_helio_net(start = 0, fixed = true);
  SI.Energy E_helio_spillage(start = 0, fixed = true);
  SI.Energy E_helio_cosine(start = 0, fixed = true);
  SI.Energy E_recv_incident(start = 0, fixed = true);
  SI.Energy E_recv_net(start = 0, fixed = true);
  SI.Energy E_pb_input(start = 0, fixed = true);
  SI.Energy E_pb_gross(start = 0, fixed = true);
  SI.Energy E_pb_net(start = 0, fixed = true);
  SI.Energy E_pv_net_to_electrolyser(start = 0, fixed = true);
  SI.Energy E_pv(start = 0, fixed = true);
  SI.Energy E_pv_net_to_heater(start = 0, fixed = true);
  SI.Energy E_elec_net_to_electrolyser(start = 0, fixed = true);
  SI.Energy E_losses_curtailment(start = 0, fixed = true);
  SI.Energy E_losses_availability(start = 0, fixed = true);
  SI.Energy E_losses_optical(start = 0, fixed = true);
  SI.Energy E_losses_defocus(start = 0, fixed = true);
  SI.Energy E_check;
  //******************** New terms for Sankey diagram
  SI.Energy E_solars(start = 0) "The energy of the sun received by the heliostat, whether it turns on or not";
  SI.Energy E_loss_curtailments(start = 0) "Loss energy due to wind speed, low solar angles, low DNI";
  SI.Energy E_resources_after_curtailments(start = 0) "The energy that goes to the heliostat field after curtailment losses";
  SI.Energy E_loss_availabilitys(start = 0) "Loss energy due to heliostat availability";
  SI.Energy E_helio_raws(start = 0) "The energy that goes into the heliostat field";
  SI.Energy E_loss_opticals(start = 0) "Loss energy due to optical efficiency";
  SI.Energy E_helio_nets(start = 0) "Energy after optical efficiency";
  SI.Energy E_loss_defocuss(start = 0) "Loss energy due to defocus";
  SI.Energy E_incident_rcv(start = 0) "Energy that comes into the receiver aperture";
  SI.Energy E_loss_rcvs(start = 0) "Energy loss ion the rcv";
  SI.Energy E_absorbed_rcv(start = 0) "Energy that is absorbed by the rcv";
  SI.Energy E_loss_tank(start = 0) "Energy loss in the tank";
  SI.Energy E_PB_input(start = 0) "PB input";
  SI.Energy E_PB_loss_thermal_to_electric(start = 0) "PB loss (efficiency thermal to electric)";
  SI.Energy E_PB_gross(start = 0) "PB gross energy before parasitics";
  SI.Energy E_PB_loss_parasitics(start = 0) "PB loss (parasitics)";
  SI.Energy E_PB_net(start = 0) "PB net energy after parasitics";
  SI.Energy E_checks(start = 0);
  Real eta_curtail_off(start = 0);
  Real eta_optical(start = 0);
  Real eta_he_av(start = 0);
  Real eta_curtail_defocus(start = 0);
  Real eta_recv_abs(start = 0);
  Real eta_recv_thermal(start = 0);
  Real eta_storage(start = 0);
  Real eta_pb_gross(start = 0);
  Real eta_pb_net(start = 0);
  Real eta_solartoelec(start = 0);
  Integer counter_shut_down(start = 0, fixed = true);
  //*********************Dispatch optimiser variables
  Real SLinit;
  Real counter(start = const_t);
  Real counter_dual_tank(start = const_t_dual_tank);
  Real time_simul;
  Real optimalDispatch;
  Real DEmax(start = Q_flow_des / 1e6) "Maximum dispatchable heat from the storage in MWth";
  Real SLmax(start = E_max * 2.77778e-10) "Storage capacity in MWh th";
  Real dummyRatio;
  Real accumulated_m;
  Real TOD_W(start = 0, fixed = true) "Product of Time-of-day factor and instant of electric power";
  //**********************
  SI.Mass M_pcl_rcv "mass pcl from the rcv";
  SI.Energy E_PV_AEL_direct "PV electricity directly converted by the AEL";
  SI.Energy E_dumped_PV;
  SI.Energy E_dumped_PB;
  SI.Energy E_PB_AEL "PB electricity used in the conversion in AEL";
  //********************** Instantiate arrays (double *) to store foresights
  SI.Angle dec_horizon[integer(horizon / dt)] "Solar declination angle for the next forecast horizon";
  SI.Angle hra_horizon[integer(horizon / dt)] "Solar hour angle for the next forecast horizon";
  SI.Efficiency eta_opt_horizon_calm[integer(horizon / dt)] "Optical efficiency for calm condition for the next horizon";
  SI.Efficiency eta_opt_horizon_windy[integer(horizon / dt)] "Optical efficiency for windy  condition for the next horizon";
  SI.Efficiency eta_opt_horizon[integer(horizon / dt)] "Optical efficiency for the next horizon";
  SI.Velocity Wspd_horizon[integer(horizon / dt)] "Wind speed for the next horizon";
  SI.Angle slope_error_runtime[integer(horizon / dt)] "Slope error as a function of wind speed during runtime";
  Real pri_horizon[integer(horizon / dt)] "Price for the next horizon";
  Real DNI_horizon[integer(horizon / dt)] "DNI for the next horizon";
  SI.Angle dec_horizon_dual_tank[integer(n_horizon / dt_dual_tank)] "Solar declination angle for the next forecast horizon";
  SI.Angle hra_horizon_dual_tank[integer(n_horizon / dt_dual_tank)] "Solar hour angle for the next forecast horizon";
  SI.Efficiency eta_opt_horizon_calm_dual_tank[integer(n_horizon / dt_dual_tank)] "Optical efficiency for calm condition for the next horizon";
  SI.Efficiency eta_opt_horizon_windy_dual_tank[integer(n_horizon / dt_dual_tank)] "Optical efficiency for windy  condition for the next horizon";
  SI.Efficiency eta_opt_horizon_dual_tank[integer(n_horizon / dt_dual_tank)] "Optical efficiency for the next horizon";
  SI.Velocity Wspd_horizon_dual_tank[integer(n_horizon / dt_dual_tank)] "Wind speed for the next horizon";
  SI.Angle slope_error_runtime_dual_tank[integer(n_horizon / dt_dual_tank)] "Slope error as a function of wind speed during runtime";
  SI.Efficiency eta_rcv_horizon_dual_tank[integer(n_horizon / dt_dual_tank)] "Rcv efficiency for the next horizon";
  Real dni_dual_tank[integer(n_horizon / dt_dual_tank)] "DNI for the next forecast horizon [W/m2]";
  Real dhi_dual_tank[integer(n_horizon / dt_dual_tank)] "GHI for the next forecast horizon [W/m2]";
  Real t_forecast_dual_tank[integer(n_horizon / dt_dual_tank)];
  SI.Power P_PV_in_z[integer(n_horizon / dt_dual_tank)] "PV net output for the next horizon";
  Real albedo_horizon_dual_tank[integer(n_horizon / dt_dual_tank)];
  SI.Temperature Tdry_horizon_dual_tank[integer(n_horizon / dt_dual_tank)];
  SI.Angle zen_horizon_dual_tank[integer(n_horizon / dt_dual_tank)];
  SI.Angle azi_horizon_dual_tank[integer(n_horizon / dt_dual_tank)];
  SI.Angle ele_horizon_dual_tank[integer(n_horizon / dt_dual_tank)];
  //*************** Variables to capture the opt dispatch dual tank
  Real fraction_P_PV_AEL "Fraction of PV power that goes to the AEL";
  Real fraction_P_PV_TES "Fraction of PV power that goes to the TES";
  Real fraction_P_PV_dumped "Fraction of PV power that is dumped";
  Real fraction_Q_TES_PB "Fraction of the Q_flow_des that is sent to the PB by TES in form of hot particle";
  Real fraction_Q_TES_HX "Fraction of the Q_TES_HX_max that is sent to the heat exchanger by TES";
  Real fraction_Q_H2_reactor "Fraction of the H2 LHV to the reactor that is sent by the H2 tank";
  Real fraction_Q_H2_burner "Fraction of the heat H2 burner sent by the H2 tank";
  Real fraction_Q_SF_TES "Fraction of SF heat goes to TES";
  Real fraction_Q_SF_dumped "Fraction of SF heat being dumped";
  SI.Mass H2_mass "Accummulate mass of hydrogen produed by the plant";
  SI.Volume H2_volume "Accummulated H2 volumen in m3 produced by the plant";
  SI.Mass NaturalGas "Accummulated backup natural gas used annually consumed by SMR";
  SI.Mass H2_mass_SMR "Accummulated H2 mass produced by SMR";
  SI.Mass H2O_mass_SMR "Acummulated H2O mass consumed by the SMR";
  SI.Mass H2_mass_electrolyser "Accummulated H2 mass produced by Electrolyser";
  SI.Mass H2O_mass_electrolyser "Acummulated H2O mass consumed by the electrolyser";
  replaceable package MedH2 = Modelica.Media.IdealGases.SingleGases.H2;
  /*Energy balance of the system*/
  // PV energy balance
  SI.Energy E_PV_in(start = 0) "Energy delivered by the PV to the system before PV power splitter [J]";
  SI.Energy E_PV_dumped(start = 0) "Dumped PV electricity [J]";
  // TES
  SI.Energy E_sf_after_opt(start = 0) "Solar field delivered energy after optical efficiency and environment curtailment [J]";
  SI.Energy E_loss_rcv(start = 0) "Energy dumped from the system due to losses in receiver [J]";
  SI.Energy E_dumped_sf(start = 0) "Energy dumped from the system due to defocussing solar field (Tank is full) [J]";
  SI.Energy E_loss_TES(start = 0) "Energy loss in the system due to losses in the tank [J]";
  SI.Energy E_TES_HX(start = 0) "Energy provided by the TES to heat up H2 using TES HX [J]";
  // PB
  SI.Energy E_loss_PB(start = 0) "Energy loss in the PB [J]";
  SI.Energy E_loss_PB_ramping(start = 0) "Energy loss in the PB during ramping [J]";
  // AEL
  SI.Energy E_loss_AEL(start = 0) "AEL loss [J]";
  //H2 stg
  SI.Energy E_H2_reactor(start = 0) "Hot hydrogen produced by the system [J]";
  SI.Energy E_H2_burner(start = 0) "Wasted H2 for burning [J]";
  // HTR
  SI.Energy E_dumped_HTR(start = 0) "Dumped heater energy [J]";
  SI.Energy E_loss_HTR(start = 0) "Loss heater energy [J]";
  // Delta E of the system
  SI.Energy E_delta_system(start = 0) "Difference in the TES and H2stg level relative to the initial condition [J]";
  SI.Energy E_balance(start = 0);
  MedH2.ThermodynamicState state_H2;
  SI.Mass M_hot_air(start = 0) "Accummulated hot air mass for the industry";
  Modelica.Blocks.Math.Sum sumMdot(nin = 2) annotation(
    Placement(visible = true, transformation(origin = {24, 38}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  SolarTherm.Models.Storage.Tank.Tank_H2 tank_H2(E_tank_capacity = E_H2_tank, H2_demand = H2_mdot_tank, H2_demand_emergency = H2_mdot_tank_emergency, LHV_H2 = LHV_H2, t_storage = t_storage_H2_final, on_dispatch_optimiser_dual_tank = set_dispatch_optimiser_dual_tank) annotation(
    Placement(visible = true, transformation(origin = {409, -71}, extent = {{-41, -41}, {41, 41}}, rotation = 0)));
  SolarTherm.Models.Fluid.HeatExchangers.HeatExchanger_ParticleGasIndirect heatExchanger_ParticleGasIndirectH2(T_in_gas_DP = T_in_gas_DP_H2, T_out_gas_DP = T_out_gas_DP_H2_HX, m_dot_gas_DP = H2_mdot_numberup / 6 + H2_mdot_recycle, T_in_pcl_DP = T_in_pcl_DP_H2, T_out_pcl_DP = T_out_pcl_DP_H2, d_gas_pipe = d_gas_pipe_H2, w_ch = w_ch_H2, th_wall = th_wall_H2, A_plate = A_plate_H2, T_in_gas_off = T_in_gas_off_H2, T_out_pcl_off = T_cold_set, saved_model_dir = saved_model_dir_H2, gas_data = gas_data_H2, X_max = X_max_H2, X_min = X_min_H2, out_max = out_max_H2, out_min = out_min_H2, redeclare replaceable package MedGas2 = MedGas2_H2, m_dot_gas_recycle = H2_mdot_recycle, on_dispatch_optimiser_dual_tank = set_dispatch_optimiser_dual_tank) annotation(
    Placement(visible = true, transformation(origin = {186, -226}, extent = {{-14, 14}, {14, -14}}, rotation = 0)));
  SolarTherm.Models.Fluid.Pumps.LiftSimple liftH2HX annotation(
    Placement(visible = true, transformation(origin = {78, -136}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.UtilitiesComponent.H2_splitter h2_splitter(mdot_split = H2_mdot_burner, mdot_split_emergency = H2_mdot_burner_emergency, m_dot_HX = H2_mdot_numberup) annotation(
    Placement(visible = true, transformation(origin = {369, -259}, extent = {{29, -29}, {-29, 29}}, rotation = 0)));
  SolarTherm.Models.UtilitiesComponent.SimpleGasBurner simpleGasBurner(eta_burner = eta_burner, LHV_fuel = LHV_H2) annotation(
    Placement(visible = true, transformation(origin = {112, -318}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  SolarTherm.Models.UtilitiesComponent.H2_recycle h2_recycle(H2_mdot_recycle = H2_mdot_recycle * 6, H2_feedstock_des = H2_mdot_numberup) annotation(
    Placement(visible = true, transformation(origin = {342, -198}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  Modelica.Blocks.Math.Sum sum(nin = 2) annotation(
    Placement(visible = true, transformation(origin = {266, -192}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression upstream_multiplier(y = upstreammultiplier) annotation(
    Placement(visible = true, transformation(origin = {88, 120}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  SolarTherm.Models.UtilitiesComponent.DecompressingUnit decompressingUnit(T_tank = T_electrolyser, redeclare replaceable package MedGas2 = MedGas2_H2, m_dot_H2_AEL_DP = H2_mdot_tank) annotation(
    Placement(visible = true, transformation(origin = {524, -246}, extent = {{56, -56}, {-56, 56}}, rotation = 0)));
  SolarTherm.Models.UtilitiesComponent.ReactorRecuperator reactorRecuperator(redeclare replaceable package MedGas2 = MedGas2_H2, redeclare replaceable package MedWater = Modelica.Media.Water.StandardWater, m_dot_excess_H2 = H2_mdot_recycle * 6, m_dot_excess_H2O = H2O_product * 6, T5 = T_out_intercooler_target, T4 = T_H2O_excess_condensed) annotation(
    Placement(visible = true, transformation(origin = {514, 154}, extent = {{-48, -48}, {48, 48}}, rotation = 0)));
  SolarTherm.Models.PV.PVarray PVArray(PV_Target = PV_Target, azi_s = azi_s, ele_min = ele_min, ele_s = ele_s, lat = lat) annotation(
    Placement(visible = true, transformation(origin = {-68, -86}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  SolarTherm.Models.PV.PV_Power_Splitter pV_Power_Splitter annotation(
    Placement(visible = true, transformation(origin = {-12, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Sum sumTESHeater(nin = 2) annotation(
    Placement(visible = true, transformation(origin = {-12, -126}, extent = {{-6, -6}, {6, 6}}, rotation = -90)));
  SolarTherm.Models.UtilitiesComponent.SimpleExchanger simpleExchanger(T_in_HTF = T_out_gas_DP_H2_HX, T_in_HTF_emergency = T_in_gas_DP_H2, T_out_HTF_target = T_out_gas_target_H2, m_dot_HTF_recycle = H2_mdot_recycle * 6) annotation(
    Placement(visible = true, transformation(origin = {25, -253}, extent = {{31, -31}, {-31, 31}}, rotation = 0)));
  Modelica.Blocks.Math.Product numberup annotation(
    Placement(visible = true, transformation(origin = {142, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.RealExpression upstream_numberdown(y = 1 / upstreammultiplier) annotation(
    Placement(visible = true, transformation(origin = {-134, -118}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Math.Product pv_divider annotation(
    Placement(visible = true, transformation(origin = {-76, -136}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  SolarTherm.Models.UtilitiesComponent.CompressorTrainH2 compressorTrainH2(T_in = T_electrolyser, m_dot_H2_AEL_DP = H2_mdot_tank, redeclare replaceable package MedGas2 = MedGas2_H2) annotation(
    Placement(visible = true, transformation(origin = {344, -2}, extent = {{-16, -16}, {16, 16}}, rotation = -90)));
algorithm
/*
//********************************** Energy balance dual tank
// PV energy balance
  der(E_PV_in) = pV_Power_Splitter.P_PV_in * upstreammultiplier;
  der(E_PV_dumped) = pV_Power_Splitter.P_PV_dumped * upstreammultiplier;
// TES
  der(E_sf_after_opt) = particleReceiver.heat.Q_flow * upstreammultiplier;
  der(E_loss_rcv) = particleReceiver.Q_loss * upstreammultiplier;
  der(E_dumped_sf) = particleReceiver.Q_dumped * upstreammultiplier;
  der(E_loss_TES) = (der(tankHot.E_losses_total) + der(tankCold.E_losses_total)) * upstreammultiplier;
  der(E_TES_HX) = heatExchanger_ParticleGasIndirectH2.particle_port_in.m_flow * upstreammultiplier * (inStream(heatExchanger_ParticleGasIndirectH2.particle_port_in.h_outflow) - inStream(heatExchanger_ParticleGasIndirectH2.particle_port_out.h_outflow));
// PB
  der(E_loss_PB) = if powerBlock.ramping then 0 else (powerBlock.Q_HX - powerBlock.W_gross) * upstreammultiplier;
  der(E_loss_PB_ramping) = if powerBlock.ramping then powerBlock.Q_HX * powerBlock.PB_output_scaling_factor * upstreammultiplier else 0;
// AEL
  der(E_loss_AEL) = if electrolyser.H2_mdot_out > 1e-3 then (electrolyser.W_electrolyser_final - electrolyser.H2_mdot_out * electrolyser.LHV) * upstreammultiplier else 0;
//H2 stg
  der(E_H2_reactor) = (tank_H2.H2_out - tank_H2.mdot_h2_burner) * tank_H2.LHV_H2;
  der(E_H2_burner) = tank_H2.mdot_h2_burner * tank_H2.LHV_H2;
// HTR
  der(E_dumped_HTR) = simpleElectricalHeater.W_dumped_htr * upstreammultiplier "Dumped heater energy [J]";
  der(E_loss_HTR) = simpleElectricalHeater.W_loss * upstreammultiplier;
//********************************** End energy balance dual tank
*/
  if time > 3.15359e7 then
    E_delta_system := (tankHot.L - tankHot.L_start) / 100 * E_max * upstreammultiplier + (tank_H2.L - tank_H2.L_start) / 100 * E_H2_tank;
    E_balance := E_sf_after_opt + E_PV_in + E_delta_system - E_PV_dumped - E_loss_rcv - E_dumped_sf - E_dumped_HTR - E_loss_HTR - E_loss_TES - E_loss_PB - E_loss_PB_ramping - E_loss_AEL - E_TES_HX - E_H2_burner - E_H2_reactor;
    if on_CSP then
      eta_curtail_off := E_helio_incident / E_resource;
      eta_optical := E_resource_after_optical_eff / E_resource;
      eta_he_av := he_av_design;
      eta_curtail_defocus := E_helio_net / E_helio_raw;
      eta_recv_abs := E_recv_incident / E_helio_net;
      eta_recv_thermal := E_recv_net / E_recv_incident;
      eta_storage := E_pb_input / E_recv_net;
      eta_pb_gross := E_pb_gross / E_pb_input;
      eta_pb_net := E_pb_net / E_pb_input;
      eta_solartoelec := E_pb_net / E_resource;
      E_check := E_resource - E_losses_availability - E_losses_curtailment - E_losses_defocus - E_losses_optical - E_helio_net;
    else
      eta_curtail_off := 0;
      eta_optical := 0;
      eta_he_av := 0;
      eta_curtail_defocus := 0;
      eta_recv_abs := 0;
      eta_recv_thermal := 0;
      eta_storage := 0;
      eta_pb_gross := 0;
      if E_pb_input > 1000 then
        eta_pb_net := E_pb_net / E_pb_input;
      else
        eta_pb_net := 0;
      end if;
      eta_solartoelec := 0;
      E_check := 0;
    end if;
  end if;
initial equation
//T_in_gas_DP_H2 = reactorRecuperator.T_H2_after_recup;
  H_tower_final = if use_g3p3_field then 234.023630 else H_tower "[H&T] Tower height";
  etaRCV = particleReceiver1DCalculator.particleReceiver1D.eta_rec;
  upstreammultiplier = upstreammultiplier_change;
  P_compressor_H2 = compressorTrainH2.W_comp_des1 + compressorTrainH2.W_comp_des2 + compressorTrainH2.W_comp_des3;
  P_expander_H2 = decompressingUnit.W_exp_total;
  etaAEL = electrolyser.H2_mdot_design_point * LHV_H2 / P_hybrid_system_final;
  A_recup = reactorRecuperator.A_recup_total;
  T_in_gas_off_H2 = reactorRecuperator.T8;
  t_storage_H2_final = t_storage_H2;
  m_dot_pcl_DP_H2_HX = heatExchanger_ParticleGasIndirectH2.m_dot_pcl_DP;
  A_HX_H2 = heatExchanger_ParticleGasIndirectH2.A_HX;
  Q_HX_H2 = heatExchanger_ParticleGasIndirectH2.designPoint.Q_rating_HX;
  P_hybrid_system_final = P_hybrid_system;
  PB_fraction_final = PB_fraction;
  CSP_fraction_final = CSP_fraction;
  PV_fraction_final = PV_fraction;
//H2_mdot_target = electrolyser.H2_mdot_design_point;
  O2_mdot_target = electrolyser.O2_mdot_design_point;
  m_dot_PCL_industrial = heatExchanger_ParticleGas.m_dot_PCL_DP;
  A_HX_industrial = heatExchanger_ParticleGas.A_HX;
  N_parallel_final_PV = PVArray.N_parallel_final;
  N_series_final_PV = PVArray.N_series;
  C_extra_structure = structureExtraCost.C_extra_structure_cost;
  if set_external_storage then
    R_tower = W_rcv / 2;
//**********Assuming tower diameter = receiver width
  else
    R_tower = max(towerInnerDiameterCalculator.D_inner_tower / 2 + 1.83, 25 / 2);
  end if;
  omega_twister = ceil(washingFrequencyCalculator.omega);
//opt_file = heliostatsField.optical.tablefile;
  opt_file = heliostatsField.optical.tablefile;
  A_rcv = if use_g3p3_field then 28.757315 * 30.980541 else particleReceiver1DCalculator.particleReceiver1D.A_ap;
  rec_fr = 1 - particleReceiver1DCalculator.particleReceiver1D.eta_rec;
  m_flow_fac = particleReceiver1DCalculator.particleReceiver1D.mdot;
  Q_flow_defocus_calculated = defocuscalculator.particleReceiver1D.q_solar * A_rcv;
//********************* Retrieve the PB design points result from the PB model intialisation
  if which_PB_model == 0 then
//********************* CEA - Power Block
    eta_gross_base = -1;
    eta_Q_base = -1;
    Q_flow_des = -1;
    etaG = -1;
    etaPB = -1;
    eff_blk = -1;
    m_flow_blk = -1;
    UA_HX = -1;
  elseif which_PB_model == 1 then
//********************* NREL Power Block - Initialisation
    NREL_PB_configurations = SolarTherm.Utilities.designNRELPB(P_gross, T_in_ref_blk, p_high, dT_PHX_cold_approach, eta_comp_main, eta_comp_re, eta_turb, dT_mc_approach, blk_T_amb_des, HTF_name, htf_choice, SolarTherm_path, T_cold_set);
//********************* NREL Power Block - Retrieve initialisation result
    Q_flow_des = NREL_PB_configurations[10] "Heat transfer of the PHX at the design point";
    eta_gross_base = NREL_PB_configurations[11] "After cooling power, before other parasities and fixed self-power consumption";
    eta_Q_base = NREL_PB_configurations[12] "";
    etaG = eta_gross_base "after cooling power, motor eff and fixed load parasities";
    etaPB = eta_gross_base;
    eff_blk = etaG;
    m_flow_blk = NREL_PB_configurations[9] "HTF mass flow rate at the design point";
    dT_PHX_hot_approach = NREL_PB_configurations[13];
    UA_HX = NREL_PB_configurations[7];
  else
//********************* Wrong which_PB_model choice
    eta_gross_base = -1;
    eta_Q_base = -1;
    Q_flow_des = -1;
    etaG = -1;
    etaPB = -1;
    eff_blk = -1;
    m_flow_blk = -1;
    dT_PHX_hot_approach = -1;
    UA_HX = -1;
  end if;
//************************************ Heat Exchanger Area Calculation
  U_HX = 1 / (1 / h_conv_CO2 + 1 / SolarTherm.Utilities.h_particle(SolarTherm.Media.SolidParticles.CarboHSP_utilities.lamda_T(0.5 * (T_cold_set + T_hot_set)), 0.008, 1, 1));
  A_HX = UA_HX / U_HX "Heat exchanger total area";
//************************************ Power Block Component Cost Calculation
  if which_PB_model == 1 then
//************************************ Based on NREL Power Block - Cost function is based on Kevin Albrect, 2019 https://is.gd/3VN0O7
    C_HTR = pri_recuperator * NREL_PB_configurations[1] ^ 0.8933;
    C_LTR = pri_recuperator * NREL_PB_configurations[2] ^ 0.8933;
    C_turbine = pri_turbine * (NREL_PB_configurations[3] / 10 ^ 3) ^ 0.5886;
    C_mainCompressor = pri_compressor * (NREL_PB_configurations[4] / 10 ^ 3) ^ 0.9142;
    C_reCompressor = pri_compressor * (NREL_PB_configurations[5] / 10 ^ 3) ^ 0.9142;
    C_cooler = pri_cooler * NREL_PB_configurations[6];
    C_exchanger = pri_PHX_per_area * A_HX + pri_PHX_BOP_s * m_flow_blk + pri_PHX_BOP_CO2 * NREL_PB_configurations[8];
    C_generator = pri_generator * (P_gross / 10 ^ 6) ^ 0.5463;
  elseif which_PB_model == 0 then
//************************************ Based on NREL Power Blockbased on CEA Power Block*/
    C_HTR = sCO2PBDesignPointCalculator.powerBlock.C_HTR;
    C_LTR = sCO2PBDesignPointCalculator.powerBlock.C_LTR;
    C_turbine = sCO2PBDesignPointCalculator.powerBlock.C_turbine;
    C_mainCompressor = sCO2PBDesignPointCalculator.powerBlock.C_mainCompressor;
    C_reCompressor = sCO2PBDesignPointCalculator.powerBlock.C_reCompressor;
    C_cooler = sCO2PBDesignPointCalculator.powerBlock.C_cooler;
    C_exchanger = sCO2PBDesignPointCalculator.powerBlock.C_exchanger;
    C_generator = sCO2PBDesignPointCalculator.powerBlock.C_generator;
  else
    C_HTR = -1;
    C_LTR = -1;
    C_turbine = -1;
    C_mainCompressor = -1;
    C_reCompressor = -1;
    C_cooler = -1;
    C_exchanger = -1;
    C_generator = -1;
  end if;
//********************* Power Block Cost Calculation
  if which_PB_model == 0 then
    if set_simple_PB_cost then
      C_block = if abs(P_net - P_net_default_value) < 1 then 0 else (pri_block * P_gross / 1000 + C_exchanger) * upstreammultiplier;
    else
      C_block = if abs(P_net - P_net_default_value) < 1 then 0 else sCO2PBDesignPointCalculator.powerBlock.C_PB * upstreammultiplier;
    end if;
  elseif which_PB_model == 1 then
    if set_simple_PB_cost then
      C_block = if abs(P_net - P_net_default_value) < 1 then 0 else (pri_block * P_gross / 1000 + C_exchanger) * upstreammultiplier;
    else
      C_block = if abs(P_net - P_net_default_value) < 1 then 0 else (C_HTR + C_LTR + C_turbine + C_mainCompressor + C_reCompressor + C_cooler + C_exchanger + C_generator) * upstreammultiplier;
    end if;
  else
    C_block = if abs(P_net - P_net_default_value) < 1 then 0 else -1;
  end if;
//********************************* Cost of the indirect particle-gas-HX between particle and H2
  C_exchanger_indirect = if f_TES < 1e-3 or set_no_thermal then 0 else (pri_PHX_per_area * A_HX_H2 + pri_PHX_BOP_s * heatExchanger_ParticleGasIndirectH2.designPoint.m_dot_pcl_DP + pri_PHX_BOP_CO2 * heatExchanger_ParticleGasIndirectH2.designPoint.m_dot_gas_DP) * upstreammultiplier;
//************************************ CapitalCost Calculation
  C_cap_CSP = C_field + C_site + C_receiver + C_storage + C_block + C_bop + C_heater / (1 + r_contg + r_cons) "Converting heater from TCI to bare module";
  C_cap_PV = C_PV / (1 + r_contg + r_cons) "Converting TCI to bare module";
  C_cap_SMR = C_SMR;
  C_cap_electrolyser = C_electrolyser;
  C_cap_total = C_field + C_site + C_receiver + C_storage + C_cap_PV + C_block + C_bop + C_SMR + C_electrolyser + C_heater / (1 + r_contg + r_cons) + C_storage_H2 + C_exchanger_indirect + C_recuperator_H2 + C_burner + C_compressor_H2 + C_expander_H2 "Total equipment cost";
  C_direct = (1 + r_contg) * C_cap_total;
  C_indirect = r_cons * C_direct + C_land;
  C_cap = C_direct + C_indirect;
equation
//********************************** Energy balance dual tank
// PV energy balance
  der(E_PV_in) = pV_Power_Splitter.P_PV_in;
  der(E_PV_dumped) = pV_Power_Splitter.P_PV_dumped;
// TES
  der(E_sf_after_opt) = particleReceiver.heat.Q_flow * upstreammultiplier;
  der(E_loss_rcv) = if particleReceiver.logic then particleReceiver.Q_loss * upstreammultiplier else 0;
  der(E_dumped_sf) = if particleReceiver.logic then 0 else particleReceiver.Q_dumped * upstreammultiplier;
  der(E_loss_TES) = (der(tankHot.E_losses_total) + der(tankCold.E_losses_total)) * upstreammultiplier;
  der(E_TES_HX) = heatExchanger_ParticleGasIndirectH2.particle_port_in.m_flow * upstreammultiplier * (inStream(heatExchanger_ParticleGasIndirectH2.particle_port_in.h_outflow) - inStream(heatExchanger_ParticleGasIndirectH2.particle_port_out.h_outflow));
// PB
  der(E_loss_PB) = (powerBlock.Q_HX - powerBlock.W_gross) * upstreammultiplier;
  der(E_loss_PB_ramping) = 0;
// AEL
  der(E_loss_AEL) = if electrolyser.H2_mdot_out > 1e-3 then electrolyser.W_electrolyser_final - electrolyser.H2_mdot_out * electrolyser.LHV else 0;
//H2 stg
  der(E_H2_reactor) = (tank_H2.H2_out - tank_H2.mdot_h2_burner) * tank_H2.LHV_H2;
  der(E_H2_burner) = tank_H2.mdot_h2_burner * tank_H2.LHV_H2;
// HTR
  der(E_dumped_HTR) = simpleElectricalHeater.W_dumped_htr;
  der(E_loss_HTR) = simpleElectricalHeater.W_loss;
//********************************** End energy balance dual tank
  when simpleExchanger.m_dot_hot_HTF < 0 then
    counter_shut_down = pre(counter_shut_down) + 1;
  end when;
//******************************************************************
//************************* TO CLOSE THE EQUATION (EMERGENCY BURNER ON AND OFF) FIXME: DRAW BOOLEAN OUTPUT FROM THE TANK TO THE SPLITTER
  h2_splitter.emergency_burner = tank_H2.emergency_burner;
  simpleExchanger.emergency_burner = tank_H2.emergency_burner;
//************************* ENDING
  der(E_dumped_PV) = electrolyser.W_dumped_PV;
  der(E_dumped_PB) = electrolyser.W_dumped_PB;
  der(M_pcl_rcv) = particleReceiver.m_flow_out "mass going to the hot tank that is accummulated";
  der(E_PV_AEL_direct) = min(PVArray.W_net, P_hybrid_system);
  der(E_PB_AEL) = powerBlock.W_net;
  der(E_pv) = PVArray.P_net;
  der(E_pv_net_to_electrolyser) = electrolyser.W_electrolyser_PV - electrolyser.W_dumped_PV;
  der(E_pv_net_to_heater) = electrolyser.W_dumped_PV;
  der(E_elec_net_to_electrolyser) = electrolyser.W_electrolyser_final;
  state_H2.p = electrolyser.p_electrolyser;
  state_H2.T = electrolyser.T_electrolyser;
  der(NaturalGas) = der(SMR.CH4_SMR);
  der(H2_volume) = der(H2_mass) / MedH2.density(state_H2);
  if HX_always_off then
    der(M_hot_air) = 0;
  else
    der(M_hot_air) = if heatExchanger_ParticleGas.on then heatExchanger_ParticleGas.m_dot_air else 0;
  end if;
  der(H2_mass) = max(0, simpleExchanger.m_dot_hot_HTF) "H2 yield of the plant [kg]";
  der(H2_mass_SMR) = SMR.H2_SMR "Accummulated H2 mass produced by SMR";
  der(H2O_mass_SMR) = der(SMR.H2O_SMR) "Acummulated H2O mass consumed by the SMR";
  der(H2_mass_electrolyser) = der(H2_mass) "Accummulated H2 mass produced by Electrolyser";
  der(H2O_mass_electrolyser) = der(electrolyser.H2O_mass) "Acummulated H2O mass consumed by the electrolyser";
//************************************ Sankey diagram calculation
  der(E_solars) = max(sun.dni * n_helios * A_helio, 0.0);
  der(E_loss_curtailments) = if heliostatsField.on_internal == false then max(sun.dni * n_helios * A_helio, 0.0) else 0;
  der(E_resources_after_curtailments) = if heliostatsField.on_internal == true then max(sun.dni * n_helios * A_helio, 0.0) else 0;
  der(E_loss_availabilitys) = if heliostatsField.on_internal == true then max(sun.dni * n_helios * A_helio * (1 - he_av_design), 0.0) else 0;
  der(E_helio_raws) = if heliostatsField.on_internal == true then max(sun.dni * n_helios * A_helio * he_av_design, 0.0) else 0;
  der(E_loss_opticals) = if heliostatsField.on_internal == true then max(sun.dni * n_helios * A_helio * he_av_design * (1 - heliostatsField.nu), 0.0) else 0;
  der(E_helio_nets) = if heliostatsField.on_internal == true then max(sun.dni * n_helios * A_helio * he_av_design * heliostatsField.nu, 0.0) else 0;
  der(E_loss_defocuss) = if heliostatsField.on_internal == true then if heliostatsField.defocus_internal == true then max(heliostatsField.Q_raw - heliostatsField.Q_defocus, 0) else 0 else 0;
  der(E_incident_rcv) = particleReceiver.Q_incident;
  der(E_loss_rcvs) = particleReceiver.Q_loss;
  der(E_absorbed_rcv) = particleReceiver.Q_rcv;
  E_loss_tank = abs(tankCold.E_losses_total) + abs(tankHot.E_losses_total);
  der(E_PB_input) = powerBlock.Q_HX;
  der(E_PB_loss_thermal_to_electric) = powerBlock.Q_HX - powerBlock.W_gross;
  der(E_PB_gross) = powerBlock.W_gross;
  der(E_PB_loss_parasitics) = powerBlock.W_gross - powerBlock.W_net;
  der(E_PB_net) = powerBlock.W_net;
  E_checks = E_solars - E_loss_curtailments - E_loss_availabilitys - E_loss_opticals - E_loss_defocuss - E_loss_rcvs - E_loss_tank - E_PB_loss_thermal_to_electric - E_PB_loss_parasitics - E_PB_net;
//************************************ Equations below exist to close the model
  particleReceiver.raw_input[1] = H_rcv "Particle receiver drop height - input to the static receiver surrogate model";
  particleReceiver.raw_input[2] = ar_rec "Particle receiver aspect ratio (H/W) - input to the static receiver surrogate model";
  particleReceiver.raw_input[3] = T_hot_set "Particle receiver outlet target temperature (K) - input to the static receiver surrogate model";
  particleReceiver.raw_input[4] = particleReceiver.heat.Q_flow "Particle receiver incident heat (W) - input to the static receiver surrogate model";
  particleReceiver.raw_input[5] = particleReceiver.T_in "Particle receiver inlet temp. (K) - input to the static receiver surrogate model";
  particleReceiver.raw_input[6] = particleReceiver.Tamb "Particle receiver ambient temp. (K) - input to the static receiver surrogate model";
  particleReceiver.raw_input[7] = 1 + 0.13929 * particleReceiver.Wspd * Modelica.Math.exp(-1 * ((abs(Modelica.SIunits.Conversions.to_deg(particleReceiver.Wdir) - 180) - 105) / 30) ^ 2) "Particle receiver wind factor (product of wind speed and wind direction) - input to the static receiver surrogate model";
//************************************ Assigning the input to the surrogate model*/
  powerBlock.raw_input[1] = powerBlock.load;
  powerBlock.raw_input[2] = tankHot.medium.T;
  powerBlock.raw_input[3] = data.Tdry;
  der(accumulated_m) = particleReceiver.fluid_a.m_flow;
//************************************ Dual tank dispatch optimiser variables
  E_TES_init = if use_TES then tankHot.L / 100 * TES_capacity else 0 "Initial level of TES [MWhth] -- TES_capacity has been multiplied by upstreammultiplier";
  H2stg_init = tank_H2.L / 100 * H2stg_capacity "Initial level of H2 storage [MWhth]";
  if set_dispatch_optimiser_dual_tank == true then
//******************** If dispatch optimiser dual tank is ON then populate the array
    der(counter_dual_tank) = 1;
    when counter_dual_tank > 0 then
      time_simul_dual_tank = time "Rounding down the time";
      for i in 1:integer(n_horizon / dt_dual_tank) loop
//******************** Obtain the declination and solar hour angle to interpolate optical eff for the next horizon
        t_forecast_dual_tank[i] = if time_simul_dual_tank + i * dt_dual_tank * 3600 < 3.1536e7 then time_simul_dual_tank + i * dt_dual_tank * 3600 else time_simul_dual_tank + i * dt_dual_tank * 3600 - 3.1536e7;
        (dec_horizon_dual_tank[i], hra_horizon_dual_tank[i]) = SolarTherm.Models.Sources.SolarFunctions.PSA_Algorithm(t_forecast_dual_tank[i], lon, lat, t_zone, year);
        eta_opt_horizon_calm_dual_tank[i] = SolarTherm.Utilities.opt_eff_horizon(CV.to_deg(dec_horizon_dual_tank[i]), CV.to_deg(hra_horizon_dual_tank[i]), opt_eff_calm);
        eta_opt_horizon_windy_dual_tank[i] = SolarTherm.Utilities.opt_eff_horizon(CV.to_deg(dec_horizon_dual_tank[i]), CV.to_deg(hra_horizon_dual_tank[i]), opt_eff_windy);
//******************** Obtain the wind speed for the next horizon
        Wspd_horizon_dual_tank[i] = SolarTherm.Utilities.wind_speed_horizon(t_forecast_dual_tank[i], wea_table);
//******************** Obtain the DNI for the next horizon
        dni_dual_tank[i] = SolarTherm.Utilities.DNI_horizon(t_forecast_dual_tank[i], wea_table);
//******************** Obtain the GHI for the next horizon
        dhi_dual_tank[i] = SolarTherm.Utilities.DHI_horizon(t_forecast_dual_tank[i], wea_table);
        albedo_horizon_dual_tank[i] = SolarTherm.Utilities.Albedo_horizon(t_forecast_dual_tank[i], wea_table);
        Tdry_horizon_dual_tank[i] = Modelica.SIunits.Conversions.from_degC(SolarTherm.Utilities.Tdry_horizon(t_forecast_dual_tank[i], wea_table));
        zen_horizon_dual_tank[i] = SolarTherm.Models.Sources.SolarFunctions.solarZenith(dec_horizon_dual_tank[i], if abs(time - 43200) < 1e-3 then 1e-6 else hra_horizon_dual_tank[i], lat);
        if (cos(zen_horizon_dual_tank[i]) * sin(CV.from_deg(lat)) - sin(dec_horizon_dual_tank[i])) / (sin(zen_horizon_dual_tank[i]) * cos(CV.from_deg(lat))) > 0.999 then
          azi_horizon_dual_tank[i] = CN.pi;
        else
          azi_horizon_dual_tank[i] = SolarTherm.Models.Sources.SolarFunctions.solarAzimuth(dec_horizon_dual_tank[i], if abs(time - 43200) < 1e-3 then 1e-6 else hra_horizon_dual_tank[i], lat) + CN.pi;
        end if;
        ele_horizon_dual_tank[i] = SolarTherm.Models.Sources.SolarFunctions.elevationAngle(dec_horizon_dual_tank[i], if abs(time - 43200) < 1e-3 then 1e-6 else hra_horizon_dual_tank[i], lat);
//****************** Call the PV function to populate array of PV net power in
        P_PV_in_z[i] = if PV_fraction_final > 1e-3 then SolarTherm.Utilities.PVArray_func(N_parallel_final_PV, N_series_final_PV, dhi_dual_tank[i], albedo_horizon_dual_tank[i], Wspd_horizon_dual_tank[i], Tdry_horizon_dual_tank[i], PV_Target, azi_s, ele_s, lat, ele_min, dni_dual_tank[i], dec_horizon_dual_tank[i], if abs(time - 43200) < 1e-3 then 1e-6 else hra_horizon_dual_tank[i], zen_horizon_dual_tank[i], ele_horizon_dual_tank[i], azi_horizon_dual_tank[i]) / 1e6 else 0 "Net PV power production over the next forecast horizon in MWe";
//******************** Obtain the optical efficiency for the next horizon
        if set_swaying_optical_eff == true then
//********************* Assuming linear relationship between effective slope error vs. optical efficiency
          slope_error_runtime_dual_tank[i] = slope_error + (slope_error_windy - slope_error) / Wspd_max * max(Wspd_horizon_dual_tank[i], 0);
          if on_CSP then
//************************ Allow for no solar field input by zeroing the optical efficiency
            if slope_error_runtime_dual_tank[i] < slope_error_windy then
              eta_opt_horizon_dual_tank[i] = eta_opt_horizon_windy_dual_tank[i] + (eta_opt_horizon_calm_dual_tank[i] - eta_opt_horizon_windy_dual_tank[i]) / (slope_error_windy - slope_error) * (slope_error_runtime_dual_tank[i] - slope_error);
            else
              eta_opt_horizon_dual_tank[i] = eta_opt_horizon_windy[i];
            end if;
          else
            eta_opt_horizon_dual_tank[i] = 0;
          end if;
        else
//********************* No relationship between effective slope error vs. optical efficiency
          slope_error_runtime_dual_tank[i] = -1;
          if on_CSP then
            eta_opt_horizon_dual_tank[i] = eta_opt_horizon_calm_dual_tank[i];
          else
            eta_opt_horizon_dual_tank[i] = 0;
          end if;
        end if;
//********************* Populate the eta_rcv_horizon_dual_tank_array
        eta_rcv_horizon_dual_tank[i] = SolarTherm.Utilities.eta_rcv_horizon(particleReceiver.session, {H_rcv, ar_rec, T_hot_set, eta_opt_horizon_dual_tank[i] * A_field * dni_dual_tank[i], T_cold_set, Tdry_horizon_dual_tank[i], 1 + 0.13929 * Wspd_horizon_dual_tank[i] * Modelica.Math.exp(-1 * ((abs(Modelica.SIunits.Conversions.to_deg(Wspd_horizon_dual_tank[i]) - 180) - 105) / 30) ^ 2)}, particleReceiver.inputsize, X_max_rcv, X_min_rcv, particleReceiver.out_max, particleReceiver.out_min);
      end for;
      optimalSolution = SolarTherm.Utilities.LinProgFuncDualTank(c0, c1, c2, c3, n_horizon, dt_dual_tank, c_ratio, etaPB, eta_opt_horizon_dual_tank, eta_rcv_horizon_dual_tank, etaAEL, etaBurner, etaResistive, dni_dual_tank, P_PV_in_z, P_AEL_nameplate, A_sf, PB_size, TES_capacity * hot_tnk_full_ub / 100, Q_TES_HX_max, E_TES_init, E_TES_min, Q_H2_reactor_max, LHV, delta_h_H2, H2stg_capacity, H2stg_init, H2stg_min, time_simul_dual_tank);
//******************* Grab the results for simulation purpose
// PV solutions
      fraction_P_PV_AEL = optimalSolution[1];
      fraction_P_PV_TES = optimalSolution[2];
      fraction_P_PV_dumped = optimalSolution[3];
// TES solutions
      fraction_Q_TES_PB = optimalSolution[4];
      fraction_Q_TES_HX = optimalSolution[5];
// H2 stg solutions
      fraction_Q_H2_reactor = optimalSolution[6];
      fraction_Q_H2_burner = optimalSolution[7];
// SF solutions
      fraction_Q_SF_TES = optimalSolution[8];
      fraction_Q_SF_dumped = optimalSolution[9];
//*************************** Aditional equations to close the system level model
      if fraction_P_PV_AEL > (-1) then
// LP solved
        pV_Power_Splitter.fraction_P_PV_AEL = fraction_P_PV_AEL;
        pV_Power_Splitter.fraction_P_PV_TES = fraction_P_PV_TES;
        pV_Power_Splitter.fraction_P_PV_dumped = fraction_P_PV_dumped;
        controlHot.logic.fraction_Q_TES_PB = fraction_Q_TES_PB;
        heatExchanger_ParticleGasIndirectH2.fraction_Q_TES_HX = fraction_Q_TES_HX;
        tank_H2.fraction_Q_H2_reactor = fraction_Q_H2_reactor;
        tank_H2.fraction_Q_TES_HX = fraction_Q_TES_HX;
        tank_H2.fraction_Q_H2_burner = fraction_Q_H2_burner;
        particleReceiver.fraction_Q_SF_dumped = fraction_Q_SF_dumped;
      else
// LP unsolved
        if PVArray.P_net > electrolyser.P_electro_requested then
          pV_Power_Splitter.fraction_P_PV_AEL = electrolyser.P_electro_requested / PVArray.P_net;
          pV_Power_Splitter.fraction_P_PV_TES = (PVArray.P_net - electrolyser.P_electro_requested) / PVArray.P_net;
          controlHot.logic.fraction_Q_TES_PB = 0;
        else
          pV_Power_Splitter.fraction_P_PV_AEL = 1;
          pV_Power_Splitter.fraction_P_PV_TES = 0;
          controlHot.logic.fraction_Q_TES_PB = if PB_fraction_final < 1e-3 then 0 else min(1, (electrolyser.P_electro_requested - PVArray.P_net) / (upstreammultiplier * P_net));
        end if;
        pV_Power_Splitter.fraction_P_PV_dumped = 1 - fraction_P_PV_AEL - fraction_P_PV_TES;
// TES solutions
        if tank_H2.on_discharge then
          if tank_H2.on_discharge_TES then
            tank_H2.fraction_Q_H2_reactor = 1;
            tank_H2.fraction_Q_H2_burner = 0;
            tank_H2.fraction_Q_TES_HX = 1;
            heatExchanger_ParticleGasIndirectH2.fraction_Q_TES_HX = 1;
          else
            tank_H2.fraction_Q_H2_reactor = 1;
            tank_H2.fraction_Q_H2_burner = 1;
            tank_H2.fraction_Q_TES_HX = 0;
            heatExchanger_ParticleGasIndirectH2.fraction_Q_TES_HX = 0;
          end if;
        else
          tank_H2.fraction_Q_H2_reactor = 0;
          tank_H2.fraction_Q_H2_burner = 0;
          tank_H2.fraction_Q_TES_HX = 0;
          heatExchanger_ParticleGasIndirectH2.fraction_Q_TES_HX = 0;
        end if;
        particleReceiver.fraction_Q_SF_dumped = 0;
      end if;
//*************************** End
      reinit(counter_dual_tank, const_t_dual_tank);
    end when;
  else
//******************** If dispatch optimiser dual tank is OFF then populate the array
    der(counter_dual_tank) = 0;
    time_simul_dual_tank = time;
    for i in 1:integer(n_horizon / dt_dual_tank) loop
//******************** Obtain the declination and solar hour angle to interpolate optical eff for the next horizon
      t_forecast_dual_tank[i] = -1;
      dec_horizon_dual_tank[i] = -1;
      hra_horizon_dual_tank[i] = -1;
      eta_opt_horizon_calm_dual_tank[i] = -1;
      eta_opt_horizon_windy_dual_tank[i] = -1;
//******************** Obtain the wind speed for the next horizon
      Wspd_horizon_dual_tank[i] = -1;
//******************** Obtain the DNI for the next horizon
      dni_dual_tank[i] = -1;
//******************** Obtain the GHI for the next horizon
      dhi_dual_tank[i] = -1;
      Tdry_horizon_dual_tank[i] = -1;
      zen_horizon_dual_tank[i] = -1;
      azi_horizon_dual_tank[i] = -1;
      ele_horizon_dual_tank[i] = -1;
      albedo_horizon_dual_tank[i] = -1;
//****************** Call the PV function to populate array of PV net power in
      P_PV_in_z[i] = -1;
//******************** Obtain the optical efficiency for the next horizon
      slope_error_runtime_dual_tank[i] = -1;
      eta_opt_horizon_dual_tank[i] = -1;
//********************* Populate the eta_rcv_horizon_dual_tank_array
      eta_rcv_horizon_dual_tank[i] = -1;
    end for;
    optimalSolution = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
//******************* Grab the results for simulation purpose
// PV solutions
    if PVArray.P_net > P_hybrid_system_final then
      fraction_P_PV_AEL = P_hybrid_system_final / PVArray.P_net;
      fraction_P_PV_TES = if set_no_thermal then 0 else (PVArray.P_net - P_hybrid_system_final) / PVArray.P_net "If no thermal, then no TES";
    else
      fraction_P_PV_AEL = 1;
      fraction_P_PV_TES = 0;
    end if;
    fraction_P_PV_dumped = 1 - fraction_P_PV_AEL - fraction_P_PV_TES;
// TES solutions
    fraction_Q_TES_PB = if set_no_thermal or PB_fraction_final < 1e-3 then 0 else 1;
    fraction_Q_TES_HX = if set_no_thermal then 0 else 1 "If no thermal, then no TES, no indirect particle gas HX";
// H2 stg solutions
    fraction_Q_H2_reactor = -1;
    fraction_Q_H2_burner = -1;
// SF solutions
    fraction_Q_SF_TES = -1;
    fraction_Q_SF_dumped = -1;
//*************************** Aditional equations to close the system level model
    pV_Power_Splitter.fraction_P_PV_AEL = fraction_P_PV_AEL;
//fraction_P_PV_AEL;
    pV_Power_Splitter.fraction_P_PV_TES = fraction_P_PV_TES;
//fraction_P_PV_TES;
    pV_Power_Splitter.fraction_P_PV_dumped = fraction_P_PV_dumped;
//fraction_P_PV_dumped;
    controlHot.logic.fraction_Q_TES_PB = fraction_Q_TES_PB;
    heatExchanger_ParticleGasIndirectH2.fraction_Q_TES_HX = fraction_Q_TES_HX;
    tank_H2.fraction_Q_H2_reactor = 1;
    tank_H2.fraction_Q_TES_HX = 1;
    tank_H2.fraction_Q_H2_burner = 1;
    particleReceiver.fraction_Q_SF_dumped = 0;
//*************************** End
  end if;
//************************************ Dispatch optimiser variables
  SLmax = E_max * 2.77778e-10;
  DEmax = Q_flow_des / 1e6;
  SLinit = tankHot.L / 100 * SLmax;
  if set_dispatch_optimiser == true then
    der(counter) = 1;
  else
    der(counter) = 0;
  end if;
  when counter > 0 then
    time_simul = floor(time) "Rounding down the time";
    if set_dispatch_optimiser then
//******************** If dispatch optimiser is on then populate the array
      for i in 1:integer(horizon / dt) loop
//******************** Obtain the declination and solar hour angle to interpolate optical eff for the next horizon
        (dec_horizon[i], hra_horizon[i]) = SolarTherm.Models.Sources.SolarFunctions.PSA_Algorithm(if time_simul + (i - 1) * dt * 3600 < 3.1536e7 then time_simul + (i - 1) * dt * 3600 else time_simul + (i - 1) * dt * 3600 - 3.1536e7, lon, lat, t_zone, year);
        eta_opt_horizon_calm[i] = SolarTherm.Utilities.opt_eff_horizon(CV.to_deg(dec_horizon[i]), CV.to_deg(hra_horizon[i]), opt_eff_calm);
        eta_opt_horizon_windy[i] = SolarTherm.Utilities.opt_eff_horizon(CV.to_deg(dec_horizon[i]), CV.to_deg(hra_horizon[i]), opt_eff_windy);
//******************** Obtain the wind speed for the next horizon
        Wspd_horizon[i] = SolarTherm.Utilities.wind_speed_horizon(if time_simul + (i - 1) * dt * 3600 < 3.1536e7 then time_simul + (i - 1) * dt * 3600 else time_simul + (i - 1) * dt * 3600 - 3.1536e7, wea_table);
//******************** Obtain the DNI for the next horizon
        DNI_horizon[i] = SolarTherm.Utilities.DNI_horizon(if time_simul + (i - 1) * dt * 3600 < 3.1536e7 then time_simul + (i - 1) * dt * 3600 else time_simul + (i - 1) * dt * 3600 - 3.1536e7, wea_table);
//******************** Obtain the price for the next horizon
        pri_horizon[i] = SolarTherm.Utilities.pri_horizon(if time_simul + (i - 1) * dt * 3600 < 3.1536e7 then time_simul + (i - 1) * dt * 3600 else time_simul + (i - 1) * dt * 3600 - 3.1536e7, pri_table);
//******************** Obtain the optical efficiency for the next horizon
        if set_swaying_optical_eff == true then
          slope_error_runtime[i] = slope_error + (slope_error_windy - slope_error) / Wspd_max * max(Wspd_horizon[i], 0);
          if slope_error_runtime[i] < slope_error_windy then
//********************* Assuming linear relationship between effective slope error vs. optical efficiency
            eta_opt_horizon[i] = (eta_opt_horizon_windy[i] + (eta_opt_horizon_calm[i] - eta_opt_horizon_windy[i]) / (slope_error_windy - slope_error) * (slope_error_runtime[i] - slope_error)) * eta_rcv_assumption;
          else
            eta_opt_horizon[i] = eta_opt_horizon_windy[i] * eta_rcv_assumption;
          end if;
        else
          slope_error_runtime[i] = -1;
          eta_opt_horizon[i] = eta_opt_horizon_calm[i];
        end if;
      end for;
      optimalDispatch = SolarTherm.Utilities.LinProgFunc_particle(horizon, dt, time_simul, eta_opt_horizon, DNI_horizon, pri_horizon, etaG, t_storage, DEmax, SLmax, SLinit, SLminrel, A_field);
    else
//******************** If dispatch optimiser is OFF then populate the array with non-sense values
      optimalDispatch = DEmax "If dispatch optimiser is unused, the dispatched energy is always at maximum whenever possible";
      for i in 1:integer(horizon / dt) loop
        dec_horizon[i] = 0;
        hra_horizon[i] = 0;
        eta_opt_horizon_calm[i] = -1;
        eta_opt_horizon_windy[i] = -1;
        eta_opt_horizon[i] = -1;
        slope_error_runtime[i] = -1;
        Wspd_horizon[i] = -1;
        DNI_horizon[i] = -1;
        pri_horizon[i] = -1;
      end for;
    end if;
    reinit(counter, const_t);
  end when;
//************* Equation here is to close the power block controller models  since the component is lack of 1 equation to close optimalMassFlow variable
  controlHot.logic.optimalMassFlow = if set_dispatch_optimiser == true then min(optimalDispatch / DEmax * m_flow_blk, m_flow_blk) else 0;
  dummyRatio = optimalDispatch / DEmax;
  der(E_resource) = max(sun.dni * n_helios * A_helio, 0.0);
  der(E_losses_optical) = (1 - heliostatsField.nu) * max(heliostatsField.solar.dni * heliostatsField.n_h * heliostatsField.A_h, 0.0);
  der(E_losses_availability) = (1 - he_av_design) * max(heliostatsField.nu * heliostatsField.solar.dni * heliostatsField.n_h * heliostatsField.A_h, 0.0);
  der(E_losses_curtailment) = if heliostatsField.on_hf == true then 0 else (1 - he_av_design) * max(heliostatsField.nu * heliostatsField.solar.dni * heliostatsField.n_h * heliostatsField.A_h, 0.0);
//der(E_losses_defocus) = if heliostatsField.on_internal then if heliostatsField.defocus_internal then abs(heliostatsField.Q_net - heliostatsField.Q_raw) else 0 else 0;
  der(E_losses_defocus) = particleReceiver.heat.Q_flow - particleReceiver.Q_incident;
  der(E_resource_after_optical_eff) = max(sun.dni * n_helios * A_helio, 0.0) * heliostatsField.nu;
  der(E_helio_incident) = if heliostatsField.on_hf then heliostatsField.n_h * heliostatsField.A_h * max(0.0, heliostatsField.solar.dni) else 0.0;
  der(E_helio_raw) = heliostatsField.Q_raw;
  der(E_helio_net) = heliostatsField.Q_net;
  der(E_helio_spillage) = heliostatsField.Q_spil;
  der(E_helio_cosine) = heliostatsField.Q_cosine;
  der(E_recv_incident) = particleReceiver.heat.Q_flow;
  der(E_recv_net) = particleReceiver.Q_rcv;
  der(E_pb_input) = powerBlock.Q_HX;
  der(E_pb_gross) = -powerBlock.W_gross;
  der(E_pb_net) = powerBlock.W_net;
  P_elec = sum1.y;
  der(E_elec) = sum1.y;
  R_spot = market.profit;
  der(TOD_W) = market.price.price * powerBlock.W_net;
//************************************ Connections from data
  connect(DNI_input.y, sun.dni) annotation(
    Line(points = {{-113, 70}, {-102, 70}, {-102, 69.8}, {-82.6, 69.8}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(Wspd_input.y, heliostatsField.Wspd) annotation(
    Line(points = {{-112.7, 30}, {-100, 30}, {-100, 29.54}, {-87.68, 29.54}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(Pres_input.y, tankCold.p_top) annotation(
    Line(points = {{55, 28}, {49.5, 28}, {49.5, -8}, {50, -8}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(Pres_input.y, tankHot.p_top) annotation(
    Line(points = {{55, 28}, {46, 28}, {8, 28}, {8, 78}, {30.5, 78.3}, {32, 74}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(Tamb_input.y, powerBlock.T_amb) annotation(
    Line(points = {{118, 80}, {102.4, 80}, {102.4, 34.4}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(Tamb_input.y, tankHot.T_amb) annotation(
    Line(points = {{118, 80}, {22, 80}, {22, 74}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(tankCold.T_amb, Tamb_input.y) annotation(
    Line(points = {{57, -8}, {57, 20}, {92, 20}, {92, 42}, {118, 42}, {118, 80}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
//************************************ Fluid connections
  connect(liftHX.fluid_b, powerBlock.fluid_a) annotation(
    Line(points = {{70, 44}, {86, 44}, {86, 29.46}, {98.08, 29.46}}, color = {0, 127, 255}));
//************************************ controlHot connections
  connect(tankHot.L, controlHot.L_mea) annotation(
    Line(points = {{38, 68}, {38, 72.5}, {46, 72.5}}, color = {0, 0, 127}));
  connect(controlHot.m_flow, liftHX.m_flow) annotation(
    Line(points = {{59, 69}, {79, 69}, {79, 54}}, color = {0, 0, 127}));
  connect(controlHot.defocus, or1.u1) annotation(
    Line(points = {{52, 77}, {52, 86}, {-106, 86}, {-106, 8}, {-102.8, 8}}, color = {255, 0, 255}, pattern = LinePattern.Dash));
//************************************ Solar field connections i.e. solar.heat port and control
  connect(sun.solar, heliostatsField.solar) annotation(
    Line(points = {{-72, 60}, {-72, 36}}, color = {255, 128, 0}));
  connect(or1.y, heliostatsField.defocus) annotation(
    Line(points = {{-93.6, 8}, {-92, 8}, {-92, 8.8}, {-87.68, 8.8}}, color = {255, 0, 255}, pattern = LinePattern.Dash));
//************************************ PowerBlock connections
  connect(parasities_input.y, powerBlock.parasities) annotation(
    Line(points = {{147, 44}, {147, 40.85}, {109.6, 40.85}, {109.6, 34.4}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(always_on.y, heliostatsField.on_hopper) annotation(
    Line(points = {{-116, -4}, {-112, -4}, {-112, 20}, {-88, 20}, {-88, 20}}, color = {255, 0, 255}));
  connect(controlHot.rampingout, powerBlock.ramping) annotation(
    Line(points = {{59, 72}, {106, 72}, {106, 34}}, color = {255, 0, 255}));
  connect(heliostatsField.heat, particleReceiver.heat) annotation(
    Line(points = {{-56, 28}, {-42, 28}, {-42, 34}, {-42, 34}}, color = {191, 0, 0}));
  connect(heliostatsField.on, particleReceiver.on) annotation(
    Line(points = {{-72, 2}, {-72, 2}, {-72, -14}, {-38, -14}, {-38, 16}, {-28, 16}, {-28, 16}}, color = {255, 0, 255}));
  connect(liftRC.fluid_b, particleReceiver.fluid_a) annotation(
    Line(points = {{-8, -26}, {-22, -26}, {-22, 16}, {-22, 16}}, color = {0, 127, 255}));
  connect(Wind_dir.y, particleReceiver.Wdir) annotation(
    Line(points = {{-116, 52}, {-30, 52}, {-30, 42}, {-30, 42}}, color = {0, 0, 127}));
  connect(Tamb_input.y, particleReceiver.Tamb) annotation(
    Line(points = {{118, 80}, {-22, 80}, {-22, 42}, {-22, 42}}, color = {0, 0, 127}));
  connect(booleanExpression.y, or1.u2) annotation(
    Line(points = {{-116, -22}, {-106, -22}, {-106, 4}, {-102, 4}, {-102, 4}}, color = {255, 0, 255}));
  connect(tankCold.L, particleReceiver.level) annotation(
    Line(points = {{46, -14}, {24, -14}, {24, 24}, {-20, 24}}, color = {0, 0, 127}));
  connect(Wspd_input.y, wspdScaler.WspdGround) annotation(
    Line(points = {{-112, 30}, {-98, 30}, {-98, 50}, {-59, 50}, {-59, 45}}, color = {0, 0, 127}));
  connect(wspdScaler.WspdOut, particleReceiver.Wspd) annotation(
    Line(points = {{-43, 45}, {-26, 45}, {-26, 42}}, color = {0, 0, 127}));
  connect(powerBlock.fluid_b, LiftCold.fluid_a) annotation(
    Line(points = {{96, 14}, {84, 14}, {84, -6}, {118, -6}, {118, -22}, {112, -22}, {112, -22}}, color = {0, 127, 255}));
  connect(LiftCold.fluid_b, tankCold.fluid_a) annotation(
    Line(points = {{100, -22}, {76, -22}, {76, -14}, {62, -14}, {62, -13}}, color = {0, 127, 255}));
  connect(tankHot.fluid_b, liftHX.fluid_a) annotation(
    Line(points = {{38, 57}, {38, 44}, {81, 44}}, color = {0, 127, 255}));
  connect(tankHot.fluid_b2, lift.fluid_a) annotation(
    Line(points = {{38, 54}, {136, 54}, {136, 102}, {184, 102}, {184, 102}}, color = {255, 0, 0}, thickness = 1));
  connect(lift.fluid_b, heatExchanger_ParticleGas.particle_port_in) annotation(
    Line(points = {{196, 102}, {248, 102}, {248, 34}}, color = {255, 0, 0}, thickness = 1));
  connect(heatExchanger_ParticleGas.particle_port_out, tankCold.fluid_a) annotation(
    Line(points = {{248, -2}, {254, -2}, {254, -14}, {62, -14}, {62, -13}}, color = {255, 0, 0}, thickness = 1));
  connect(tankHot.L, heatExchanger_ParticleGas.L_hot_tank) annotation(
    Line(points = {{38, 68}, {106, 68}, {106, 114}, {288, 114}, {288, 10}, {252, 10}}, color = {0, 0, 127}));
  connect(Tamb_input.y, heatExchanger_ParticleGas.T_amb) annotation(
    Line(points = {{118, 80}, {116, 80}, {116, 160}, {320, 160}, {320, 22}, {252, 22}, {252, 22}}, color = {0, 0, 127}));
  connect(electrolyser.H2_mdot_out, SMR.H2_electrolyser) annotation(
    Line(points = {{238, -56}, {248, -56}, {248, -84}, {130, -84}, {130, -108}, {132, -108}}, color = {255, 0, 0}, pattern = LinePattern.Dash));
//************************************ End of The System Model ********************************************\\
  connect(SMR.H2_SMR, sumH2.u[1]) annotation(
    Line(points = {{168, -142}, {240, -142}, {240, -116}, {286, -116}, {286, -116}}, color = {255, 0, 0}, thickness = 1));
  connect(electrolyser.H2_mdot_out, sumH2.u[2]) annotation(
    Line(points = {{238, -56}, {270, -56}, {270, -116}, {286, -116}, {286, -116}}, color = {255, 0, 0}, thickness = 1));
  connect(tankCold.L, simpleElectricalHeater.L) annotation(
    Line(points = {{46, -14}, {24, -14}, {24, -104}, {27, -104}, {27, -169}}, color = {0, 85, 0}, thickness = 1));
  connect(simpleElectricalHeater.particle_port_out, tankHot.fluid_a) annotation(
    Line(points = {{-15, -177}, {-38, -177}, {-38, -60}, {12, -60}, {12, 70}, {16, 70}}, color = {0, 85, 0}, thickness = 1));
  connect(particleReceiver.m_flow_out, sumMdot.u[1]) annotation(
    Line(points = {{-22, 28}, {1, 28}, {1, 31}, {24, 31}}, color = {255, 0, 0}, pattern = LinePattern.Dash, thickness = 0.5));
  connect(simpleElectricalHeater.mdot_heater, sumMdot.u[2]) annotation(
    Line(points = {{17, -169}, {17, -66}, {18, -66}, {18, 10}, {24, 10}, {24, 31}}, color = {255, 0, 0}, pattern = LinePattern.Dash, thickness = 0.5));
  connect(sumMdot.y, controlHot.m_flow_in) annotation(
    Line(points = {{24, 45}, {24, 48}, {46, 48}, {46, 65.5}}, color = {255, 0, 0}, pattern = LinePattern.Dash, thickness = 0.5));
  connect(particleReceiver.fluid_b, tankHot.fluid_a) annotation(
    Line(points = {{-20, 38}, {-8, 38}, {-8, 70}, {16, 70}, {16, 70}}, color = {0, 127, 255}));
  connect(tankCold.fluid_b2, simpleElectricalHeater.particle_port_in) annotation(
    Line(points = {{46, -28}, {40, -28}, {40, -128}, {47, -128}, {47, -177}}, color = {85, 85, 0}, thickness = 1));
  connect(tankCold.fluid_b, liftRC.fluid_a) annotation(
    Line(points = {{46, -24}, {6, -24}, {6, -26}, {6, -26}}, color = {0, 127, 255}));
  connect(controlHot.T_mea, tankHot.T_mea) annotation(
    Line(points = {{49, 62}, {49, 52}, {42, 52}, {42, 64}, {40, 64}}, color = {255, 0, 0}, pattern = LinePattern.Dash, thickness = 0.5));
  connect(controlHot.T_amb_input, Tamb_input.y) annotation(
    Line(points = {{53, 62}, {53, 52}, {64, 52}, {64, 80}, {118, 80}}, color = {255, 0, 0}, pattern = LinePattern.Dash, thickness = 0.5));
  connect(Tamb_input.y, tank_H2.T_amb) annotation(
    Line(points = {{118, 80}, {114, 80}, {114, 196}, {380, 196}, {380, -27}, {409, -27}}, color = {0, 0, 127}));
  connect(electrolyser.W_d_PV, simpleElectricalHeater.W_PV) annotation(
    Line(points = {{238, -52}, {260, -52}, {260, -88}, {50, -88}, {50, -106}, {47, -106}, {47, -169}}, color = {0, 0, 127}));
  connect(electrolyser.W_d_PB, simpleElectricalHeater.W_PB) annotation(
    Line(points = {{238, -46}, {284, -46}, {284, -78}, {36, -78}, {36, -169}}, color = {0, 0, 127}));
  connect(tankHot.fluid_b2, liftH2HX.fluid_a) annotation(
    Line(points = {{38, 54}, {94, 54}, {94, -135}, {81, -135}}, color = {0, 127, 255}));
  connect(liftH2HX.fluid_b, heatExchanger_ParticleGasIndirectH2.particle_port_in) annotation(
    Line(points = {{74, -134}, {74, -189}, {176, -189}, {176, -230}}, color = {0, 127, 255}));
  connect(heatExchanger_ParticleGasIndirectH2.particle_port_out, tankCold.fluid_a) annotation(
    Line(points = {{196, -230}, {196, -12}, {62, -12}}, color = {0, 127, 255}));
  connect(tank_H2.charging, electrolyser.H2_tank_charging) annotation(
    Line(points = {{451, -68}, {444, -68}, {444, 2}, {222, 2}, {222, -36}}, color = {255, 0, 255}));
  connect(tank_H2.charging, controlHot.H2_tank_charging) annotation(
    Line(points = {{451, -68}, {456, -68}, {456, 58}, {56, 58}, {56, 62}}, color = {255, 0, 255}));
  connect(tank_H2.discharging, heatExchanger_ParticleGasIndirectH2.on_discharge) annotation(
    Line(points = {{452, -80}, {528, -80}, {528, -158}, {190, -158}, {190, -232}}, color = {255, 0, 255}));
  connect(h2_splitter.flow_out_2, simpleGasBurner.fuel_in) annotation(
    Line(points = {{338, -272}, {219, -272}, {219, -317}, {122, -317}}, color = {0, 0, 127}));
  connect(tank_H2.discharging, h2_recycle.on) annotation(
    Line(points = {{452, -80}, {464, -80}, {464, -136}, {340, -136}, {340, -173}}, color = {255, 0, 255}));
  connect(controlHot.on_discharge, tank_H2.on_discharge_TES) annotation(
    Line(points = {{58, 74}, {432, 74}, {432, -28}}, color = {255, 0, 255}));
  connect(h2_splitter.flow_out_1, sum.u[1]) annotation(
    Line(points = {{338, -248}, {302, -248}, {302, -192}, {278, -192}, {278, -192}}, color = {0, 0, 127}));
  connect(h2_recycle.out, sum.u[2]) annotation(
    Line(points = {{318, -194}, {280, -194}, {280, -192}, {278, -192}}, color = {0, 0, 127}));
  connect(sum.y, heatExchanger_ParticleGasIndirectH2.gas_in) annotation(
    Line(points = {{256, -192}, {244, -192}, {244, -222}, {196, -222}}, color = {0, 0, 127}));
  connect(tank_H2.H2_out, decompressingUnit.H2_in) annotation(
    Line(points = {{452, -98}, {690, -98}, {690, -246}, {587, -246}}, color = {0, 0, 127}));
  connect(decompressingUnit.H2_out, h2_splitter.flow_in) annotation(
    Line(points = {{462, -247}, {442, -247}, {442, -258}, {400, -258}}, color = {0, 0, 127}));
  connect(Albedo_input.y, PVArray.albedo) annotation(
    Line(points = {{-110, -88}, {-102, -88}, {-102, -78}, {-90, -78}, {-90, -78}}, color = {0, 0, 127}));
  connect(DHI_input.y, PVArray.DHI) annotation(
    Line(points = {{-110, -68}, {-92, -68}, {-92, -68}, {-90, -68}}, color = {0, 0, 127}));
  connect(Tamb_input.y, PVArray.Tdry) annotation(
    Line(points = {{118, 80}, {-106, 80}, {-106, -104}, {-90, -104}, {-90, -104}}, color = {0, 0, 127}));
  connect(sun.solar, PVArray.sun) annotation(
    Line(points = {{-72, 60}, {-68, 60}, {-68, -64}, {-68, -64}}, color = {255, 170, 0}));
  connect(Wspd_input.y, PVArray.Wspd) annotation(
    Line(points = {{-112, 30}, {-98, 30}, {-98, -94}, {-90, -94}, {-90, -94}}, color = {0, 0, 127}));
  connect(PVArray.W_net, pV_Power_Splitter.P_PV_in) annotation(
    Line(points = {{-44, -86}, {-32, -86}, {-32, -84}, {-24, -84}, {-24, -84}}, color = {0, 0, 127}));
  connect(pV_Power_Splitter.P_PV_AEL, sum1.u[1]) annotation(
    Line(points = {{-2, -84}, {88, -84}, {88, -56}, {140, -56}}, color = {255, 0, 255}, thickness = 1));
  connect(pV_Power_Splitter.P_PV_AEL, controlHot.PV_input) annotation(
    Line(points = {{-2, -84}, {40, -84}, {40, 68}, {46, 68}, {46, 70}}, color = {255, 0, 255}, thickness = 1));
  connect(pV_Power_Splitter.P_PV_AEL, electrolyser.W_electrolyser_PV) annotation(
    Line(points = {{-2, -84}, {188, -84}, {188, -60}, {196, -60}, {196, -58}}, color = {0, 255, 0}, thickness = 1));
  connect(sumTESHeater.y, simpleElectricalHeater.W_electric) annotation(
    Line(points = {{-12, -133}, {-12, -132}, {6, -132}, {6, -169}}, color = {0, 85, 0}, thickness = 1));
  connect(electrolyser.W_dumped, sumTESHeater.u[2]) annotation(
    Line(points = {{238, -62}, {244, -62}, {244, -102}, {-12, -102}, {-12, -119}}, color = {0, 85, 0}, thickness = 1));
  connect(tank_H2.emergency_burner_signal, heatExchanger_ParticleGasIndirectH2.emergency_burner) annotation(
    Line(points = {{452, -44}, {476, -44}, {476, -264}, {180, -264}, {180, -232}, {180, -232}}, color = {255, 0, 255}));
  connect(heatExchanger_ParticleGasIndirectH2.gas_out, simpleExchanger.HTF_in) annotation(
    Line(points = {{175, -222}, {150, -222}, {150, -256}, {60, -256}, {60, -254}}, color = {0, 0, 127}));
  connect(simpleGasBurner.Q_flow_burner, simpleExchanger.Q_in) annotation(
    Line(points = {{102, -318}, {-156, -318}, {-156, -176}, {24, -176}, {24, -220}, {26, -220}}, color = {0, 0, 127}));
  connect(controlHot.on_discharge, heatExchanger_ParticleGasIndirectH2.on_discharge_TES) annotation(
    Line(points = {{58, 74}, {222, 74}, {222, -254}, {182, -254}, {182, -232}, {180, -232}}, color = {255, 0, 255}));
  connect(tank_H2.mdot_h2_burner, h2_splitter.mdot_h2_burner) annotation(
    Line(points = {{450, -108}, {510, -108}, {510, -196}, {374, -196}, {374, -232}, {374, -232}}, color = {0, 0, 127}));
  connect(h2_splitter.flow_out_1, h2_recycle.H2_feedstock) annotation(
    Line(points = {{338, -248}, {332, -248}, {332, -216}, {422, -216}, {422, -192}, {364, -192}, {364, -196}}, color = {0, 0, 127}));
  connect(h2_recycle.out, simpleExchanger.HTF_rec) annotation(
    Line(points = {{318, -194}, {294, -194}, {294, -286}, {94, -286}, {94, -238}, {60, -238}, {60, -238}}, color = {0, 0, 127}));
  connect(upstream_multiplier.y, numberup.u1) annotation(
    Line(points = {{100, 120}, {126, 120}, {126, 16}, {130, 16}, {130, 14}}, color = {0, 0, 127}));
  connect(powerBlock.W_net, numberup.u2) annotation(
    Line(points = {{116, 22}, {122, 22}, {122, 2}, {130, 2}, {130, 2}}, color = {0, 0, 127}));
  connect(numberup.y, sum1.u[2]) annotation(
    Line(points = {{154, 8}, {172, 8}, {172, -30}, {124, -30}, {124, -56}, {140, -56}, {140, -56}}, color = {255, 0, 255}, thickness = 1));
  connect(sum1.y, market.W_net) annotation(
    Line(points = {{164, -56}, {178, -56}, {178, 44}, {196, 44}, {196, 44}}, color = {255, 0, 255}, thickness = 1));
  connect(numberup.y, electrolyser.W_electrolyser_PB) annotation(
    Line(points = {{154, 8}, {188, 8}, {188, -46}, {194, -46}, {194, -48}, {196, -48}}, color = {0, 255, 0}, thickness = 1));
  connect(pV_Power_Splitter.P_PV_TES, pv_divider.u2) annotation(
    Line(points = {{-12, -94}, {-12, -94}, {-12, -100}, {-34, -100}, {-34, -112}, {-82, -112}, {-82, -124}, {-82, -124}}, color = {0, 85, 0}, thickness = 1));
  connect(upstream_numberdown.y, pv_divider.u1) annotation(
    Line(points = {{-122, -118}, {-70, -118}, {-70, -124}, {-70, -124}}, color = {0, 85, 0}, thickness = 1));
  connect(pv_divider.y, sumTESHeater.u[1]) annotation(
    Line(points = {{-76, -146}, {-76, -146}, {-76, -156}, {-50, -156}, {-50, -114}, {-12, -114}, {-12, -118}, {-12, -118}}, color = {0, 85, 0}, thickness = 1));
  connect(sumH2.y, compressorTrainH2.H2_in) annotation(
    Line(points = {{310, -116}, {314, -116}, {314, 42}, {344, 42}, {344, 16}, {344, 16}}, color = {255, 0, 0}, thickness = 1));
  connect(compressorTrainH2.H2_out, tank_H2.H2_in) annotation(
    Line(points = {{344, -20}, {344, -50}, {364, -50}}, color = {0, 0, 127}));
protected
  annotation(
    Diagram(coordinateSystem(extent = {{-140, -120}, {160, 140}}, initialScale = 0.1), graphics = {Text(lineColor = {217, 67, 180}, extent = {{4, 92}, {40, 90}}, textString = "defocus strategy", fontSize = 9), Text(origin = {-8, -20}, lineColor = {217, 67, 180}, extent = {{-58, -18}, {-14, -40}}, textString = "on/off strategy", fontSize = 9), Text(origin = {12, 24}, extent = {{-52, 8}, {-4, -12}}, textString = "Receiver", fontSize = 6, fontName = "CMU Serif"), Text(origin = {12, 4}, extent = {{-110, 4}, {-62, -16}}, textString = "Heliostats Field", fontSize = 6, fontName = "CMU Serif"), Text(origin = {4, -8}, extent = {{-80, 86}, {-32, 66}}, textString = "Sun", fontSize = 6, fontName = "CMU Serif"), Text(origin = {-4, 2}, extent = {{0, 58}, {48, 38}}, textString = "Hot Tank", fontSize = 6, fontName = "CMU Serif"), Text(origin = {8, 2}, extent = {{30, -24}, {78, -44}}, textString = "Cold Tank", fontSize = 6, fontName = "CMU Serif"), Text(origin = {4, -2}, extent = {{80, 12}, {128, -8}}, textString = "Power Block", fontSize = 6, fontName = "CMU Serif"), Text(origin = {56, 12}, extent = {{112, 16}, {160, -4}}, textString = "Market", fontSize = 6, fontName = "CMU Serif"), Text(origin = {20, 4}, extent = {{-6, 20}, {42, 0}}, textString = "Receiver Control", fontSize = 6, fontName = "CMU Serif"), Text(origin = {2, 32}, extent = {{30, 62}, {78, 42}}, textString = "Power Block Control", fontSize = 6, fontName = "CMU Serif"), Text(origin = {-6, -14}, extent = {{-146, -26}, {-98, -46}}, textString = "Data Source", fontSize = 7, fontName = "CMU Serif"), Text(origin = {0, -40}, extent = {{-10, 8}, {10, -8}}, textString = "Lift Receiver", fontSize = 6, fontName = "CMU Serif"), Text(origin = {110, -42}, extent = {{-14, 8}, {14, -8}}, textString = "LiftCold", fontSize = 6, fontName = "CMU Serif"), Text(origin = {85, 59}, extent = {{-19, 11}, {19, -11}}, textString = "LiftHX", fontSize = 6, fontName = "CMU Serif")}),
    Icon(coordinateSystem(extent = {{-140, -120}, {160, 140}})),
    experiment(StopTime = 3.1536e+07, StartTime = 0, Tolerance = 1e-06, Interval = 3600),
    __Dymola_experimentSetupOutput,
    Documentation(revisions = "<html>
	<ul>
	<li> Philipe Gunawan Gan (Jan 2022) :<br>Released Hybrid CSP-PV. </li>
	</ul>

	</html>"),
    __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
end Hybrid_CSP_PV_Particle_HotH2_Burner_OptDispatch2;
