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
		columns=2:10,
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

    Modelica.Blocks.Sources.BooleanExpression curtail(y=full);
    Modelica.Blocks.Sources.RealExpression P_curtail(y=P_ren_curtail);

    // Parameters
    parameter SI.Efficiency eff_process = 1 "Power block efficiency";
    parameter Real t_storage(unit="h") = 8 "Hours of storage";
    parameter SI.Energy E_max = P_process_des*t_storage*3600 "Max stored energy";

    parameter SI.Energy E_up_u = 0.95*E_max "Upper energy limit";
    parameter SI.Energy E_up_l = 0.93*E_max "Upper energy limit";
    parameter SI.Energy E_low_u = 0.07*E_max "Lower energy limit";
    parameter SI.Energy E_low_l = 0.05*E_max "Lower energy limit";
    parameter SI.Energy E_start = 0.3*E_max "Lower energy limit";

    SI.HeatFlowRate P_direct;
    SI.HeatFlowRate P_EES_in;
    SI.HeatFlowRate P_EES_out;
    SI.Power P_elec_ely;  
    SI.Power P_elec;
    SI.Power P_ren_curtail;

    SI.MassFlowRate m_h2_ely_out;  // Declarado aquí
    SI.Efficiency eta_ely = 0.7;  // Eficiencia del electrolizador, se debe ajustar según el diseño
    SI.MassFlowRate m_h2_direct;  // Declarado aquí
    SI.MassFlowRate m_h2_stor_in;  // Declarado aquí
    SI.MassFlowRate m_h2_stor_out;  // Declarado aquí
    SI.MassFlowRate m_h2;  // Declarado aquí
    SI.MassFlowRate m_steel;  // Declarado aquí
    SI.Efficiency f_reduction = 0.9;  // Factor de reducción, se debe ajustar según el proceso
    SI.Efficiency f_smelter = 0.1;  // Factor de fundición, se debe ajustar según el proceso

    SI.Energy E_phes(min=0, max=E_max) "Stored energy";
    SI.Mass m_h2_stor;  // Declarado aquí

    Boolean full "True if the storage tank is full";
    Real P_elec_in;

function dispatch
    input Real t;
    output Real P_curtail;
    output Real P_direct;
    output Real P_EES_in;
    output Real P_EES_out;
    output Real m_h2_stor_in;
    output Real m_h2_stor_out;
    external "C" st_mip(t, P_curtail, P_direct, P_EES_in, P_EES_out, m_h2_stor_in, m_h2_stor_out);
    annotation (IncludeDirectory = "modelica://SolarTherm/Systems", Include = "#include \"forecast.c\"");
end dispatch;

initial equation
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
    (P_ren_curtail, P_direct, P_EES_in, P_EES_out, m_h2_stor_in, m_h2_stor_out) = dispatch(time);

	// Pump hydro energy system
    P_EES_in + P_direct + P_elec_ely = renewable_input.electricity;
    
    der(E_phes) = P_EES_in - P_EES_out;

    P_elec = P_EES_out + P_direct;

	// Hydrogen subsystem
	m_h2_ely_out = P_elec_ely * eta_ely;
	m_h2_ely_out = m_h2_direct + m_h2_stor_in;
	m_h2 = m_h2_direct + m_h2_stor_out;

	der(m_h2_stor) = m_h2_stor_in - m_h2_stor_out;

	// Steel Production
	m_steel = m_h2 * f_reduction + P_elec * f_smelter;

    annotation(experiment(StartTime = 0, StopTime = 864000, Interval = 300, Tolerance = 1e-06),
 Diagram);
end HydrogenDRISystem;
