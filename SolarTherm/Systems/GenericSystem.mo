within SolarTherm.Systems;

model GenericSystem
  import SI = Modelica.SIunits;
  import nSI = Modelica.SIunits.Conversions.NonSIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import FI = SolarTherm.Models.Analysis.Finances;
  import SolarTherm.Types.Solar_angles;
  // Input Parameters
  // *********************
  parameter Boolean match_sam = false "Configure to match SAM output";
  parameter Boolean fac_fra_wrt_dni = false "true if DNI fraction is to be used to calculate heat losses in receiver
  		and false if inlet thermal power fraction is to be used to calculate heat losses in receiver";
  parameter Boolean fixed_field = true "true if the size of the solar field is fixed";
  // if fixed_field=false, P_gro will be fixed and given as input, otherwise, R_des will be fixed and given as input
  parameter String wea_file "Weather file";
  parameter String opt_file "Optical efficiency file";
  parameter Solar_angles angles = Solar_angles.ele_azi "Angles used in the lookup table file";
  parameter String pri_file = "" "Electricity price file";
  parameter Real wdelay[8] = {0, 0, 0, 0, 0, 0, 0, 0} "Weather file delays";
  parameter Real SM "Solar multiple";
  parameter SI.Power P_gro(fixed = if fixed_field then false else true) "Power block gross rating at design";
  parameter SI.RadiantPower R_des(fixed = if fixed_field then true else false) "Input power for receiver at design";
  parameter SI.Length H_tower = 0 "Tower height";
  parameter SI.Area A_receiver "Receiver area";
  parameter SI.Area A0_receiver "Receiver reference area";
  parameter SI.Efficiency eff_cyc = 0.37 "Efficiency of power cycle at design point";
  parameter SI.Efficiency eff_opt = 1 "Efficiency of optics at design point (max in opt_file)";
  parameter Real t_storage(unit = "h") = 6 "Hours of storage";
  parameter Real ini_frac(min = 0, max = 1) = 0.0 "Initial fraction charged";
  parameter Boolean const_dispatch = true "Constant dispatch of energy";
  parameter nSI.Angle_deg deploy_angle = 0 "Altitude angle to start tracking";
  parameter nSI.Angle_deg stow_angle = 0 "Altitude angle to stop tracking";
  parameter SI.Irradiance dni_start = 0.001 "DNI at which solar field is switched on";
  parameter SI.Irradiance dni_stop = 0 "DNI at which solar field is switched off";
  parameter Real t_blk_heat(unit = "h") = 0 "Hours to heat up power block at blk_heat";
  parameter Real t_blk_cool(unit = "h") = 1 "Hours to cool power block";
  parameter Real blk_disp = 0.2 "Min fraction of power to power block";
  parameter Real blk_heat = 0.1 "Fraction of power to block for heating";
  parameter Real tnk_min_start = 0.1 "Minimum fraction of tank to start dispatch";
  // Parameters at the nominal design point of the system
  parameter SI.Temperature rec_T_amb_des = 298.15 "Ambient temperature at design point";
  parameter SI.Temperature tnk_T_amb_des = 298.15 "Ambient temperature at design point";
  parameter SI.Temperature blk_T_amb_des = 298.15 "Ambient temperature at design point";
  parameter SI.Temperature par_T_amb_des = 298.15 "Ambient temperature at design point";
  parameter SI.Irradiance dni_des = 1000 "DNI at design point";
  parameter Real rec_fr = 0.01 "Receiver loss fraction of radiance at design point";
  parameter Real tnk_fr = 0.01 "Tank loss fraction of tank in one day at design point";
  parameter Real par_fr = 0.01 "Parasitics fraction of power block rating at design point";
  parameter Real par_fix_fr = 0.0 "Fixed parasitics as fraction of net rating";
  // If using SAM values for rec_cf, then convert according to:
  // {c0, c1, c2, c3} -> {0, c0, c1, c2, c3}
  parameter Real rec_cf[:] = {1} "Receiver coefficients";
  parameter Real rec_ca[:] = {1} "Receiver coefficients";
  parameter Real rec_cw[:] = {1} "Receiver coefficients";
  parameter Real tnk_cf[:] = {1} "Tank coefficients";
  parameter Real tnk_ca[:] = {1} "Tank coefficients";
  parameter Real blk_cf[:] = {1} "Power block coefficients";
  parameter Real blk_ca[:] = {1} "Power block coefficients";
  parameter Real par_cf[:] = {1} "Parasitics coefficients";
  parameter Real par_ca[:] = {1} "Parasitics coefficients";
  parameter Real tnk_crit_lb(min = 0, max = 1) = 0.001;
  parameter Real tnk_crit_ub(min = 0, max = 1) = 0.02;
  parameter Real tnk_empty_lb(min = 0, max = 1) = 0.01;
  parameter Real tnk_empty_ub(min = 0, max = 1) = 0.02;
  parameter Real tnk_full_lb(min = 0, max = 1) = 0.98;
  parameter Real tnk_full_ub(min = 0, max = 1) = 0.99;
  // Contingencies should be included
  parameter Real land_mult = 1 "Land area multiplier";
  parameter FI.AreaPrice pri_field = 0 "Field cost per design aperture area";
  parameter FI.AreaPrice pri_site = 0 "Site improvements cost per area";
  parameter FI.AreaPrice pri_land = 0 "Land cost per area";
  parameter FI.Money pri_tower = 0 "Fixed tower cost";
  parameter Real idx_pri_tower = 0 "Tower cost scaling exponent";
  parameter FI.Money pri_receiver = 0 "Receiver reference cost";
  parameter Real idx_pri_receiver = 0 "Receiver cost scaling exponent";
  parameter FI.EnergyPrice pri_storage = 0 "Storage cost per energy capacity";
  parameter FI.PowerPrice pri_block = 0 "Power block cost per gross rated power";
  parameter FI.PowerPrice pri_bop = 0 "Balance of plant cost per gross rated power";
  parameter Real pri_om_name(unit = "$/W/year") = 0 "O&M cost per nameplate per year";
  parameter Real pri_om_prod(unit = "$/J/year") = 0 "O&M cost per production per year";
  parameter Real r_disc = 0.05 "Discount rate";
  parameter Integer t_life(unit = "year") = 25 "Lifetime of plant";
  parameter Integer t_cons(unit = "year") = 0 "Years of construction";
  // Calculated Parameters
  // *********************
  parameter SI.HeatFlowRate Q_flow_des = if fixed_field then if match_sam then R_des / ((1 + rec_fr) * SM) else R_des * (1 - rec_fr) / SM else P_gro / eff_cyc "Heat to power block at design";
  parameter SI.Energy E_max = t_storage * 3600 * Q_flow_des "Maximum tank stored energy";
  parameter Boolean storage = true "Storage component present";
  parameter SI.Area A_field = R_des / eff_opt / dni_des "Field area";
  parameter SI.Area A_land = land_mult * A_field "Land area";
  parameter SI.Power P_net = (1 - par_fr) * P_gro "Power block net rating at design";
  parameter SI.Power P_name = P_net "Nameplate power";
  parameter FI.Money C_field = A_field * pri_field "Field cost";
  parameter FI.Money C_site = A_field * pri_site "Site improvements cost";
  parameter FI.Money C_tower = pri_tower * exp(idx_pri_tower * H_tower) "Tower cost";
  parameter FI.Money C_receiver = pri_receiver * (A_receiver / A0_receiver) ^ idx_pri_receiver "Receiver cost";
  parameter FI.Money C_storage = E_max * pri_storage "Storage cost";
  parameter FI.Money C_block = P_gro * pri_block "Power block cost";
  parameter FI.Money C_bop = P_gro * pri_bop "Balance of plant cost";
  parameter FI.Money C_land = A_land * pri_land "Land cost";
  parameter FI.Money C_cap = C_field + C_site + C_tower + C_receiver + C_storage + C_block + C_bop + C_land "Capital costs";
  parameter FI.MoneyPerYear C_year = P_name * pri_om_name "Cost per year";
  parameter Real C_prod(unit = "$/J/year") = pri_om_prod "Cost per production per year";
  // System components
  // *********************
  SolarTherm.Models.Sources.Weather.WeatherSource wea(file = wea_file, delay = wdelay);
  SolarTherm.Models.CSP.CRS.HeliostatsField.SteeredCL CL(redeclare model OptEff = SolarTherm.Models.CSP.CRS.HeliostatsField.FileOE(angles = angles, file = opt_file), orient_north = wea.orient_north, A = if match_sam then 1.0273 * A_field else A_field, steer_rate = 0.002, target_error = 0.0001, actual_0 = 0.0, dni_start = dni_start, dni_stop = dni_stop);
  // ~8 minutes till fully on sun
  // if large can be large source of missing energy
  Models.CSP.CRS.Receivers.GenericRC RC(match_sam = match_sam, fac_fra_wrt_dni = fac_fra_wrt_dni, Q_flow_loss_des = if match_sam then rec_fr * SM * Q_flow_des else rec_fr * R_des, R_des = R_des, I_des = dni_des, T_amb_des = rec_T_amb_des, cf = rec_cf, ca = rec_ca, cw = rec_cw);
  SolarTherm.Models.Storage.Tank.GenericST ST(E_max = E_max, E_0 = E_max * ini_frac, Q_flow_loss_des = tnk_fr * E_max / (24 * 3600), T_amb_des = tnk_T_amb_des, cf = tnk_cf, ca = tnk_ca) if storage;
  SolarTherm.Models.PowerBlocks.GenericStartPB PB(eff_des = eff_cyc, Q_flow_des = Q_flow_des, T_amb_des = blk_T_amb_des, cf = blk_cf, ca = blk_ca, Q_flow_disp = blk_disp * Q_flow_des, Q_flow_heat = blk_heat * Q_flow_des, t_heat = t_blk_heat * 3600, t_cool = t_blk_cool * 3600);
  SolarTherm.Models.PowerBlocks.GenericParasitics par(P_par_des = par_fr * P_gro, P_gross_des = P_gro, T_amb_des = par_T_amb_des, cf = par_cf, ca = par_ca);
  SolarTherm.Models.Control.TankDispatch dis(full_lb = tnk_full_lb * E_max, full_ub = tnk_full_ub * E_max, empty_lb = tnk_empty_lb * E_max, empty_ub = tnk_empty_ub * E_max, crit_lb = tnk_crit_lb * E_max, crit_ub = tnk_crit_ub * E_max, level_start = tnk_min_start * E_max, disp_min = blk_disp * Q_flow_des, heat_rate = blk_heat * Q_flow_des) if storage;
  // Needs to be configured in instantiation if not const_dispatch
  SolarTherm.Models.Sources.Schedule.Scheduler sch if not const_dispatch;
  Models.Analysis.Performance per(schedule = true, pri_file = pri_file);
  // Variables
  Real sched;
  SI.Power P_elec "Net electrical power out";
  SI.Energy E_elec "Generated electricity";
  FI.Money R_spot "Spot market revenue";
  // Equations

initial equation
  if fixed_field then
    P_gro = Q_flow_des * eff_cyc;
  else
    R_des = if match_sam then SM * Q_flow_des * (1 + rec_fr) else SM * Q_flow_des / (1 - rec_fr);
  end if;
equation
  connect(wea.wbus, CL.wbus);
  connect(wea.wbus, RC.wbus);
  if storage then
    connect(wea.wbus, ST.wbus);
  end if;
  connect(wea.wbus, PB.wbus);
  connect(wea.wbus, par.wbus);
  connect(CL.R_foc, RC.R);
  if storage then
    connect(RC.Q_flow, ST.Q_flow_in);
    connect(ST.Q_flow_out, PB.Q_flow);
    connect(ST.Q_flow_in, dis.flow_in);
    dis.flow_tar = sched * Q_flow_des;
    connect(ST.E, dis.level);
    connect(PB.heated, dis.heated);
    connect(dis.flow_dis, ST.Q_flow_set);
  else
    connect(RC.Q_flow, PB.Q_flow);
  end if;
  connect(PB.P, par.P_gen);
  P_elec = PB.P - par.P_par - par_fix_fr * P_net;
  connect(P_elec, per.P_elec);
  per.P_sch = sched * P_name;
  connect(per.E_elec, E_elec);
  connect(per.R_spot, R_spot);
  if wea.wbus.azi <= 180 and wea.wbus.alt >= deploy_angle or wea.wbus.azi >= 180 and wea.wbus.alt >= stow_angle then
    CL.target = 1;
  else
    CL.target = 0;
  end if;
  if const_dispatch then
    sched = 1;
  else
    connect(wea.wbus, sch.wbus);
    sched = sch.v;
  end if;
end GenericSystem;