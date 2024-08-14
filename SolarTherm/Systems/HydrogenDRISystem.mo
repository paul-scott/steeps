within SolarTherm.Systems;

model HydrogenDRISystem
    // Imports
    import SI = Modelica.SIunits;
    import CN = Modelica.Constants;
    import CV = Modelica.SIunits.Conversions;
    import FI = SolarTherm.Models.Analysis.Finances;
    import nSI = Modelica.SIunits.Conversions.NonSIunits;
    import SolarTherm.Utilities.Tables.STMotab;

    extends Modelica.Icons.Example;

    // Renewable energy input
    parameter String pv_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/pv_gen_Pilbara 3_1.0MWe.motab");
    parameter String wind_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Renewable/wind_gen_Pilbara 3_320.0MWe.motab");
    parameter String schedule_file = Modelica.Utilities.Files.loadResource("modelica://SolarTherm/Data/Schedules/mip_schedule.motab");
    parameter STMotab.STMotab pv_motab = STMotab.STMotab(pv_file);
    parameter STMotab.STMotab wind_motab = STMotab.STMotab(wind_file);
    parameter Modelica.SIunits.Power pv_ref_size = 1e6 "PV farm reference size";
    parameter Modelica.SIunits.Power wind_ref_size = 320e6 "Wind farm reference size";
    parameter Modelica.SIunits.Power P_elec_min = 1e6;
    parameter Modelica.SIunits.Efficiency pv_fraction = 1 "Fraction of pv capacity at design";
    parameter Real renewable_multiple = 39.03883123206496 "Oversizing factor of renewable power to process power demand";
    parameter Modelica.SIunits.Power P_elec_max = renewable_multiple * P_process_des "Maximum Combined PV/Wind electrical output";
    parameter Modelica.SIunits.HeatFlowRate P_process_des = 77500e3 "Process power demand at design";

	Modelica.Blocks.Sources.CombiTimeTable scheduler(
		tableOnFile=true,
		tableName="power",
		smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
		columns=2:12,
		fileName=schedule_file)
		annotation (Placement(transformation(extent={{-90,50},{-70,70}})));	

    SolarTherm.Models.Sources.GridInput renewable_input(
        P_elec_max = P_elec_max, 
        P_elec_min = P_elec_min, 
        P_elec_pv_ref_size = pv_ref_size, 
        P_elec_wind_ref_size = wind_ref_size, 
        pv_file = pv_file, 
        pv_fraction = pv_fraction, 
        wind_file = wind_file, 
        renewable_multiple = renewable_multiple);

    Modelica.Blocks.Sources.BooleanExpression curtail(y=false);
    Modelica.Blocks.Sources.RealExpression P_curtail(y=P_ren_curtail);

    // Parameters
    parameter SI.Efficiency eff_process = 1 "Power block efficiency";
    parameter Real t_storage(unit="h") = 17.408955285 "Hours of storage";
    parameter SI.Energy E_max = 1607612.566834e3*3600 "Max stored energy";

    parameter SI.Energy E_up_u = 0.95*E_max "Upper energy limit";
    parameter SI.Energy E_up_l = 0.93*E_max "Upper energy limit";
    parameter SI.Energy E_low_u = 0.07*E_max "Lower energy limit";
    parameter SI.Energy E_low_l = 0.05*E_max "Lower energy limit";
    parameter SI.Energy E_start = 0.45*E_max "Lower energy limit";
    parameter SI.Efficiency eta_PHES_charge = sqrt(0.79) "Charging efficiency of the PHES";
	parameter SI.Efficiency eta_PHES_discharge = sqrt(0.79) "Discharging efficiency of the PHES";
	parameter SI.Efficiency eta_H2_in = 0.99 "Efficiency of storing hydrogen";
	parameter SI.Efficiency eta_H2_out = 0.99 "Efficiency of discharging hydrogen";
	parameter SI.Energy E_h2_max = 13568613.701532e3*3600;
	parameter SI.Energy E_h2_start = 0.35*E_h2_max;

    SI.HeatFlowRate P_direct;
    SI.HeatFlowRate P_EES_in;
    SI.HeatFlowRate P_EES_out;
    SI.Power P_elec_ely;  
    SI.Power P_elec;
    SI.Power P_ren_curtail;

    SI.Power E_h2_ely_out;  // Declarado aquí
    SI.Efficiency eta_ely = 0.7;  // Eficiencia del electrolizador, se debe ajustar según el diseño
    SI.Power E_h2_direct;  // Declarado aquí
    SI.Power E_h2_stor_in;  // Declarado aquí
    SI.Power E_h2_stor_out;  // Declarado aquí
    SI.Power E_h2;  // Declarado aquí
    SI.Power m_steel;  // Declarado aquí
    SI.Efficiency f_reduction = 0.9;  // Factor de reducción, se debe ajustar según el proceso
    SI.Efficiency f_smelter = 0.1;  // Factor de fundición, se debe ajustar según el proceso

    SI.Energy E_phes(min=0, max=E_max) "Stored energy";
    SI.Energy E_h2_stor;  // Declarado aquí

    Boolean full "True if the storage tank is full";
    Real P_elec_in;

function dispatch
    input Real t;
    output Real P_curtail;
    output Real P_direct;
    output Real P_EES_in;
    output Real P_EES_out;
    output Real E_h2_stor_in;
    output Real E_h2_stor_out;
    external "C" st_mip(t, P_curtail, P_direct, P_EES_in, P_EES_out, E_h2_stor_in, E_h2_stor_out);
    annotation (IncludeDirectory = "modelica://SolarTherm/Systems", Include = "#include \"forecast.c\"");
end dispatch;

initial equation
    E_phes = E_start;
    E_h2_stor = E_h2_start;
    if E_phes > E_up_u then
        full = true;
    elseif E_phes < E_up_l then
        full = false;
    else
        full = true;
    end if;

algorithm
    when E_phes > E_up_u then
        full := true;
    elsewhen E_phes < E_up_l then
        full := false;
    end when;

equation
    //t,P_curt,P_direct,P_EES_in,P_EES_out,H2_st_in,H2_st_out,P_ely,pv_out,wind_out
    P_elec_in = scheduler.y[1] + scheduler.y[2] + scheduler.y[3] + scheduler.y[7];
    connect(renewable_input.curtail,curtail.y);
    connect(renewable_input.P_schedule,P_curtail.y);
    P_ren_curtail = scheduler.y[1];
    P_direct = scheduler.y[2];
    P_EES_in = scheduler.y[3];
	P_EES_out = scheduler.y[4];
	E_h2_stor_in = scheduler.y[5];
	E_h2_stor_out = scheduler.y[6];

	// Pump hydro energy system
    P_elec_ely = max(0,renewable_input.electricity - P_ren_curtail - P_EES_in - P_direct);
    
    der(E_phes) = P_EES_in * eta_PHES_charge - P_EES_out / eta_PHES_discharge;

    P_elec = P_ren_curtail + P_direct + P_elec_ely;

	// Hydrogen subsystem
	E_h2_ely_out = P_elec_ely * eta_ely;
	E_h2_ely_out = E_h2_direct + E_h2_stor_in;
	E_h2 = E_h2_direct + E_h2_stor_out;


	der(E_h2_stor) = E_h2_stor_in * eta_H2_in - E_h2_stor_out / eta_H2_out;

	// Steel Production
	m_steel = E_h2 * f_reduction + P_elec * f_smelter;

    annotation(experiment(StartTime = 0, StopTime = 864000, Interval = 300, Tolerance = 1e-06),
 Diagram);
end HydrogenDRISystem;
