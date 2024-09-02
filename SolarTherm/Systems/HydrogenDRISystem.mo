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
    parameter Modelica.SIunits.HeatFlowRate P_nom_des = 92344e3 "Process power demand at design";
    parameter Modelica.SIunits.HeatFlowRate h2_nom_des = 347140e3 "H2 Process power demand at design";

	Modelica.Blocks.Sources.CombiTimeTable scheduler(
		tableOnFile=true,
		tableName="power",
		smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
		columns=2:13,
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

    // Electrical energy sub-system
    parameter SI.Energy E_up_u = 0.95*E_max "Upper energy limit";
    parameter SI.Energy E_up_l = 0.93*E_max "Upper energy limit";
    parameter SI.Energy E_low_u = 0.07*E_max "Lower energy limit";
    parameter SI.Energy E_low_l = 0.05*E_max "Lower energy limit";
    parameter SI.Energy E_start = 0.45*E_max "Lower energy limit";
    parameter SI.Efficiency eta_PHES_charge = sqrt(0.79) "Charging efficiency of the PHES";
	parameter SI.Efficiency eta_PHES_discharge = sqrt(0.79) "Discharging efficiency of the PHES";
	
	// Hydrogen sub-system
	parameter SI.Efficiency eta_H2_in = 0.99 "Efficiency of storing hydrogen";
	parameter SI.Efficiency eta_H2_out = 0.99 "Efficiency of discharging hydrogen";
	parameter SI.Energy E_h2_max = 13568613.701532e3*3600;
	parameter SI.Energy E_h2_start = 0.35*E_h2_max;
    parameter SI.Efficiency eta_electroliser = 0.7;
    parameter SI.Efficiency eta_H2_burner = 0.95;
    parameter Real HHV_H2 = 39.2 "kWh/kg higher heating value of hydrogen";
	parameter Real LHV_H2 = 33.33 "kWh/kg lower heating value of hydrogen";


    // Variables
    SI.HeatFlowRate P_direct "Bypassed Electrical Power";
    SI.HeatFlowRate P_EES_in "EES Electrical Power Input";
    SI.HeatFlowRate P_EES_out "EES Electrical Power Output";
    SI.Power P_elec_electroliser "Electrical Power Input to the Electroliser";
    SI.Power P_ren_curtail "Curtailed Electrical Power";

    SI.Power E_h2_electroliser_out "Chemical Energy of H2 from the electroliser";
    SI.Power E_h2_direct "Chemical Energy of H2 bypassed from the electroliser";
    SI.Power E_h2_stor_in "H2 Chemical Energy Storage Input";
    SI.Power E_h2_stor_out "H2 Chemical Energy Storage Output";
    SI.Power E_h2_burner_in "H2 Burner Input";

    SI.Energy E_EES(min=0, max=E_max) "EES Stored Energy";
    SI.Energy E_h2_stor "H2 Stored Chemical Energy";

	SI.Power P_load "Total supplied electrical load";
	SI.Power E_h2_load "Total supplied H2 energy load";
	
	SI.Energy E_load;
	SI.Energy E_h2;
	
	Real CF_elec;
	Real CF_h2;
    Boolean full "True if the storage tank is full";

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
    E_EES = E_start;
    E_h2_stor = E_h2_start;
    if E_EES > E_up_u then
        full = true;
    elseif E_EES < E_up_l then
        full = false;
    else
        full = true;
    end if;

algorithm
    when E_EES > E_up_u then
        full := true;
    elsewhen E_EES < E_up_l then
        full := false;
    end when;

equation
    connect(renewable_input.curtail,curtail.y);
    connect(renewable_input.P_schedule,P_curtail.y);
    P_direct = scheduler.y[2];
    P_EES_in = scheduler.y[3];
	P_EES_out = scheduler.y[4];
	E_h2_stor_in = scheduler.y[5];
	E_h2_stor_out = scheduler.y[6];
    P_elec_electroliser = scheduler.y[7];
	E_h2_burner_in = scheduler.y[12];

	// Pump hydro energy system
    P_ren_curtail = max(0, renewable_input.electricity - P_EES_in - P_direct - P_elec_electroliser);
    
    P_load = P_direct + P_EES_out;
    
    der(E_load) = P_load;
    CF_elec = E_load/(P_nom_des*8760*3600);
    
    der(E_EES) = P_EES_in * eta_PHES_charge - P_EES_out / eta_PHES_discharge;

	// Hydrogen subsystem
	E_h2_electroliser_out = P_elec_electroliser * eta_electroliser * LHV_H2/HHV_H2;
	E_h2_electroliser_out = E_h2_direct + E_h2_stor_in;

	der(E_h2_stor) = E_h2_stor_in * eta_H2_in - E_h2_stor_out / eta_H2_out;
	
	E_h2_load + E_h2_burner_in = E_h2_direct + E_h2_stor_out;
	der(E_h2) = E_h2_load;
	CF_h2 = E_h2/(h2_nom_des*8760*3600);

    annotation(experiment(StartTime = 0, StopTime = 31536000, Interval = 300, Tolerance = 1e-06),
 Diagram);
end HydrogenDRISystem;