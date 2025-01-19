within SolarTherm.Models.Storage.Thermocline.Annular;

model Annular_Storage_Section_SM //Stationary Momentum Version (SM)
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import Tables = Modelica.Blocks.Tables;
  //Initialize Material Packages
  replaceable package Fluid_Package = SolarTherm.Materials.PartialMaterial "Fluid Package";
  replaceable package Filler_Package = SolarTherm.Materials.PartialMaterial "Filler Package";
  replaceable package Tank_Package = SolarTherm.Materials.SS316L "Tank Package (steel shell)";
  replaceable package Encapsulation_Package = Filler_Package "Encapsulation Package, default is the same as Filler package, effectively no encapsulation";
  //Fluid Material States
  Fluid_Package.State fluid_in "Model which calculates properties at inlet of the section";
  Fluid_Package.State fluid_out "Model which calculates properties at outlet of the section";
  //Interfacial heat transfer Settings
  parameter Integer Correlation = 1 "1=WakaoKaguei, 2=MelissariArgyropolus, 3=Conservative, 4=Bellan, 5=Laminar, 6=Laminar+Turbulent, 7 = Nie";
  parameter Real f_hc = 1.0 "Scaling factor for convective heat transfer. Lower this if fouling is suspected, default 1.0";
  //Height offset for plotting purposes
  parameter SI.Length z_offset = 0.0 "Amount of height offset if there is a tank below it";
  //Tank Design parameters
  parameter SI.Energy E_max = 144e9 "Design storage capacity";
  parameter SI.Energy E_unit = rho_f_avg * V_fluid * (h_f_max - h_f_min) + rho_p * V_solid * (h_p_max - h_p_min);
  parameter SI.Energy E_unit_solid = rho_p * V_solid * (h_p_max - h_p_min) "This one assumes only storage in the solid phase";
  //parameter Real ar = 2.0 "Tank Aspect ratio H/D";
  //parameter Real eta = 0.22 "Porosity";
  //parameter Real d_p = 0.015 "Diameter of sphere (particle) (m)";
  //parameter Real t_e = d_p/(2.0*N_p) "Thickness of encapsulation, default is such that it is at a value that preserves equidistant radii discretizations (m)";
  parameter SI.Length D_pipe = 0.10 "Air pipe diameter (m)";
  parameter SI.Length D_solid = 0.20 "Solid cylinder outer diameter (m)";
  parameter SI.Length L_pipe = 10.0 "Length of the pipe in x-axis (m)";
  parameter Real Multiplier = E_max / E_unit;
  //Temperature Bounds
  parameter SI.Temperature T_min = CV.from_degC(100) "Design cold Temperature of everything in the tank (K)";
  parameter SI.Temperature T_max = CV.from_degC(600) "Design hot Temperature of everything in the tank (K)";
  parameter SI.Temperature T_start = T_min "Initial (uniform) temperature of all components (K), defaults to T_min";
  //Calculated Tank Design Parameters
  parameter SI.Length r_pipe = 0.5 * D_pipe "Air pipe radius (m)";
  parameter SI.Length r_solid = 0.5 * D_solid "Solid cylinder outer radius (m)";
  parameter SI.Volume V_solid = CN.pi * (r_solid * r_solid - r_pipe * r_pipe) * L_pipe "Total volume of solid phase (m3)";
  parameter SI.Volume V_fluid = CN.pi * (r_pipe * r_pipe) * L_pipe "Total volume of fluid phase (m3)";
  //parameter SI.Length H_tank = (4 * E_max / (CN.pi * (1 / ar) ^ 2 * (rho_f_avg * (h_f_max - h_f_min) * eta + rho_p * (h_p_max - h_p_min) * (1.0 - eta)))) ^ (1 / 3);
  //parameter SI.Diameter D_tank = H_tank / ar;
  //parameter SI.Area A = CN.pi * D_tank * D_tank / 4.0 "Cross sectional area of tank";
  //Thermal Losses
  //parameter SI.Area A_loss_tank = CN.pi*D_tank*D_tank*0.5 + CN.pi*D_tank*H_tank "Heat loss area (m2)";
  parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.4 "Heat loss coeff of surfaces (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_wall = U_loss_tank "Cylinder wall heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_side = U_loss_tank "Circle heat loss coeff (W/m2K)";
  //parameter SI.CoefficientOfHeatTransfer U_top = U_loss_tank "Top circle heat loss coeff (W/m2K)";
  //parameter SI.CoefficientOfHeatTransfer U_bot = U_loss_tank "Bottom circle heat loss coeff (W/m2K)";
  //Inititalize temperature and enthalpy profile
  parameter SI.Temperature T_f_start[N_f] = fill(T_start, N_f);
  parameter SI.Temperature T_p_start[N_f, N_p] = fill(fill(T_start, N_p), N_f);
  //parameter SI.Temperature T_e_start = T_start;
  parameter SI.SpecificEnthalpy h_f_start[N_f] = fill(Fluid_Package.h_Tf(T_start, 0.0), N_f) "Defaults to uniform";
  parameter SI.SpecificEnthalpy h_p_start[N_f, N_p] = fill(fill(Filler_Package.h_Tf(T_start, 0.0), N_p), N_f) "Defaults to uniform";
  //Property bounds
  //Fluid
  parameter SI.SpecificEnthalpy h_f_min = Fluid_Package.h_Tf(T_min, 0) "Starting enthalpy of the HTF";
  parameter SI.SpecificEnthalpy h_f_max = Fluid_Package.h_Tf(T_max, 0) "Starting enthalpy of the HTF";
  parameter SI.Density rho_f_min = Fluid_Package.rho_Tf(T_min, 0);
  parameter SI.Density rho_f_max = Fluid_Package.rho_Tf(T_max, 0);
  parameter SI.Density rho_f_avg = (rho_f_min + rho_f_max) / 2;
  //parameter SI.DynamicViscosity mu_f_avg = (Fluid_Package.mu_T(T_min)+Fluid_Package.mu_T(T_max))/2.0;
  //Filler
  parameter SI.SpecificEnthalpy h_p_max = Filler_Package.h_Tf(T_max, 1.0);
  parameter SI.SpecificEnthalpy h_p_min = Filler_Package.h_Tf(T_min, 0.0);
  parameter SI.Density rho_p_min = Filler_Package.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p_max = Filler_Package.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p = min(rho_p_min, rho_p_max) "kg/m3";
  //Encapsulation
  //parameter SI.SpecificEnthalpy h_e_max = Encapsulation_Package.h_Tf(T_max, 1.0);
  //parameter SI.SpecificEnthalpy h_e_min = Encapsulation_Package.h_Tf(T_min, 0.0);
  //parameter SI.Density rho_e_min = Encapsulation_Package.rho_Tf(T_min, 0.0);
  //parameter SI.Density rho_e_max = Encapsulation_Package.rho_Tf(T_max, 1.0);
  //parameter SI.Density rho_e = min(rho_e_min,rho_e_max) "kg/m3";
  //Discretization
  parameter Integer N_f = 10 "Number of finite volume elements in fluid";
  parameter Integer N_p = 5 "Number of finite volume elements in filler, including encapsulation";
  parameter SI.Length dx = L_pipe / N_f "discretization in the x-axis";
  parameter SI.Length dr = (r_solid - r_pipe) / N_p "discretization in the r-axis";
  //parameter SI.Length dr[N_p] = cat(1,fill(0.5 * ((d_p - 2.0*t_e) / (N_p - 1)),N_p-1),{t_e}) "radial thickness of each particle discretization, with last one being the encapsulation";
  //Initialise Fluid Array
  parameter SI.Length z_f[N_f] = Z_position(L_pipe, N_f);
  parameter SI.Area A_fx = CN.pi * r_pipe * r_pipe "Cross sectional area of the fluid column";
  parameter SI.Area A_fr = 2.0 * CN.pi * r_pipe * dx "Cylindrical area of the fluid column for one element";
  parameter SI.Volume dV_f = A_fx * dx "Volume of fluid element";
  parameter SI.Mass dm_f = rho_f_avg*dV_f "Mass of fluid element, using average density";
  SI.Temperature T_f[N_f] "(K)";
  SI.SpecificEnthalpy h_f[N_f](start = h_f_start) "J/kg";
  SI.Density rho_f[N_f] "kg/m3";
  //Plotting
  parameter Real ZDH[N_f] = Relative_Tank_Axes(L_pipe, N_f) "Non-dimensional tank vertical axis";
  //Operational Controls
  Integer State(start = 2) "operational state 2=standby, 3=discharge, 1=charge";
  //Inlet and outlet enthalpies and temperatures
  SI.SpecificEnthalpy h_in "Inlet Enthalpy depends on mass flow direction (J/kg)";
  SI.SpecificEnthalpy h_out "Outlet Enthalpy depends on mass flow direction (J/kg)";
  SI.Temperature T_in "Inlet Temperature depends on mass flow direction";
  SI.Temperature T_out "Outlet Temperature depends on mass flow direction";
  SI.Temperature T_amb;
  //Mass flow rates and superficial velocity
  SI.MassFlowRate m_flow(start = 0.0) "kg/s";
  SI.MassFlowRate m_flow_unit(start = 0.0) "kg/s";
  SI.Velocity u_flow[N_f] "Fluid velocity for each element, wrt the x-axis (m/s)";
  //SI.Velocity u_0 "Fluid velocity through empty cross section = u_flow/eta (m/s)";
  //Analytics
  SI.Energy E_stored(start = 0.0) "Make sure the tank starts from T_min for this to be correct";
  Real Level(start = 0.0) "Tank energy charging level (0-1)";
  //SI.HeatFlowRate Q_loss_total "Heat loss from the entire surface area";
  //Initialise Particle
  SI.Temperature T_p[N_f, N_p](start = T_p_start) "Temperature of particle elements";
  //Filler Geometry
  //parameter Real N_spheres_total = (N_f*6*(1-eta)*A*dx/(CN.pi*(d_p^3))) "Total number of spheres in the tank";
  parameter SI.Mass m_pj[N_p] = Annular.Annular_Masses(r_pipe, r_solid, dx, N_p, rho_p);
  parameter SI.Area A_px[N_p] = Annular.Annular_CSA(r_pipe, r_solid, N_p) "Annulus cross sectional area of each radial solid element";
  parameter SI.Length r_pj[N_p] = Annular.Annular_Radii(r_pipe, r_solid, N_p) "midpoint radius of each radial solid element";
  //Initialise Encapsulation
  //SI.Temperature T_e[N_f](start=T_e_start) "Temperature of encapsulation elements";
  //SI.SpecificEnthalpy h_e[N_f](start = h_e_min) "J/kg";
  //SI.ThermalConductivity k_e[N_f] "W/mK";
  //Calculated Pumping Losses
  //parameter Real eff_pump = 1.00 "Pump electricity to work efficiency";
  //SI.Pressure p_drop_total "Sum of all pressure drops";
  //SI.Power W_loss_pump "losses due to pressure drop";
  //Cost breakdown
  //parameter Real C_fluid = max(rho_f_max,rho_f_min)*eta*(CN.pi*D_tank*D_tank*H_tank/4.0)*Fluid_Package.cost;
  //parameter Real C_section = C_fluid + C_filler + C_insulation + C_tank + C_encapsulation;
  //parameter Real C_insulation = if U_loss_tank > 1e-3 then (16.72/U_loss_tank + 0.04269)*A_loss_tank else 0.0;
  //parameter Real C_insulation = if U_loss_tank > 1e-3 then CpA_external_insulation(T_max,U_loss_tank)*A_loss_tank else 0.0;
  //parameter SI.Density rho_bulk = eta*max(rho_f_max,rho_f_min)+(1.0-eta)*rho_p "Bulk density used to calculate hydrostatic pressure";
  //parameter Real C_tank = C_shell(max(rho_f_max,rho_f_min),H_tank,D_tank,Tank_Package.sigma_yield(T_max),Tank_Package.rho_Tf(298.15,0.0),Tank_Package.cost,4.0);
  //parameter Real C_filler = (if abs(Filler_Package.MM - Encapsulation_Package.MM) < 1e-6 then rho_p*(1.0-eta)*(CN.pi*D_tank*D_tank*H_tank/4.0)*Filler_Package.cost else rho_p*(1.0-eta)*(CN.pi*D_tank*D_tank*H_tank/4.0)*Filler_Package.cost*(((d_p-2*t_e)/d_p)^3));
  //parameter Real C_encapsulation = (if abs(Filler_Package.MM - Encapsulation_Package.MM) < 1e-6 then 0.0 else rho_e*(1.0-eta)*(CN.pi*D_tank*D_tank*H_tank/4.0)*Encapsulation_Package.cost*(1-(((d_p-2*t_e)/d_p)^3)));
  //Filler Surface Area Correction
  //parameter Real f_surface = 1.0 "Don't touch this";
  //Initialise Filler surface temperature
  SI.Temperature T_s[N_f](start = T_f_start);
  //parameter SI.Length r_p[N_p] = cat(1,Particle_Radii(d_p-2*t_e,N_p-1),{(d_p/2)-(t_e/2)}) "Radii of each particle element centre";
  //Filler mass-liquid fraction
  //Real f_p[N_f, N_p](start=fill(fill(0.0,N_p),N_f)) "Mass liquid fraction of filler";
  //Measured outlet temperature
  Real T_outlet_degC "Outlet temperature in degrees Celcius";
  //Convection Properties
  Real Pe[N_f] "Peclet Number";
  Real Bi[N_f] "Biot Number";
  Real Re[N_f] "Reynolds";
  Real Pr[N_f] "Prandtl";
  Real Nu[N_f] "Nusselt";
  Real h_c[N_f] "Heat transfer coeff (W/m2K)";
  //Real f_fric[N_f] "Friction factor";
  //Flow-Scheme dependent numbers
  //Real Nu_low[N_f] "Nusselt number at Re=2000";
  //Real Nu_high[N_f] "Nusselt number at Re=4000";
  //Real f_low[N_f] "Friction factor at Re=2000";
  //Real f_high[N_f] "Friction factor at Re=4000";
  //Integer flow_scheme[N_f](start = fill(1,N_f)) "Flow scheme (1=Laminar, 2=Transition, 3=Turbulent)";

  //SI.ThermalConductance U_in[N_f, N_p] "K/W"; //Obsolete
  SI.ThermalConductance U_in[N_f, N_p] "K/W";
  SI.ThermalConductance U_right[N_f, N_p] "K/W";
  //SI.ThermalConductance U_out_e[N_f] "K/W";
  protected
  //Filler Properties
  SI.SpecificEnthalpy h_p[N_f, N_p](start = h_p_start) "J/kg";
  SI.ThermalConductivity k_p[N_f, N_p] "W/mK";
  //parameter SI.Mass m_p[N_f, N_p] = fill(cat(1,Particle_Masses(d_p-2*t_e, N_p-1, rho_p),{rho_e*(1/6)*CN.pi*((d_p^3)-((d_p-2*t_e)^3))}), N_f) "Masses of each particle element";
  //parameter SI.Mass m_e = rho_e*(1/6)*CN.pi*((d_p^3)-((d_p-2*t_e)^3)) "Masses of encapsulation in one particle";
  //Pressure Drop
  //SI.Pressure p_drop[N_f] "Pressure drop across each mesh element";
  //Thermal Losses
  //SI.HeatFlowRate Q_loss_wall[N_f] "Heat loss from the wall";
  //SI.HeatFlowRate Q_loss_top "Heat loss from the top";
  //SI.HeatFlowRate Q_loss_bot "Heat loss from the bottom";
  //Fluid Properties
  SI.ThermalConductivity k_f[N_f] "W/mK";
  //SI.ThermalConductivity k_eff[N_f] "W/mK";
  SI.DynamicViscosity mu_f[N_f] "Pa.s";
  SI.DynamicViscosity mu_f_wall[N_f] "Pa.s (for the wall surface for convection)";
  SI.SpecificHeatCapacity c_pf[N_f] "J/kgK";
  Fluid_Package.State fluid[N_f](each h_start = h_f_min) "Fluid object array";
  Fluid_Package.State fluid_wall[N_f](each h_start = h_f_min) "Fluid object array for the surface of the wall for convection calcs";
  //Try filler state "Remove this if using function-based calculation"
  Filler_Package.State filler[N_f, N_p] "Filler object array";
  //Encapsulation_Package.State encapsulation[N_f] "Encapsulation object array";
  Real der_h_f[N_f] "Rate of change of specific enthalpy of fluid";
  //Real der_h_p[N_f,N_p] "Rate of change of specific enthalpy of solid";
algorithm
//Fluid equations
  if State == 1 then //Fluid is flowing downwards uflow is negative
    der_h_f[1] := 
    ((-2.0 * k_f[1] * k_f[2]) * (T_f[1] - T_f[2]) * A_fx / ((k_f[1] + k_f[2]) * dx) 
    +  m_flow_unit * (h_f[1] - h_f[2]) 
    - h_c[1] * (T_f[1] - T_s[1]) * A_fr ) / dm_f;
    h_out := h_f[1];
    for i in 2:N_f - 1 loop
      der_h_f[i] := 
      (2.0 * k_f[i - 1] * k_f[i] * (T_f[i - 1] - T_f[i]) * A_fx / ((k_f[i - 1] + k_f[i]) * dx) 
      - 2.0 * k_f[i] * k_f[i + 1] * (T_f[i] - T_f[i + 1]) * A_fx / ((k_f[i] + k_f[i + 1]) * dx) 
      + m_flow_unit * (h_f[i] - h_f[i + 1]) 
      - h_c[i] * (T_f[i] - T_s[i]) * A_fr) / dm_f;
    end for;
    der_h_f[N_f] := 
    (2.0 * k_f[N_f - 1] * k_f[N_f] * (T_f[N_f - 1] - T_f[N_f]) * A_fx / ((k_f[N_f - 1] + k_f[N_f]) * dx) 
    + m_flow_unit * (h_f[N_f] - h_in) 
    - h_c[N_f] * (T_f[N_f] - T_s[N_f]) * A_fr) / dm_f;
  else
    der_h_f[1] := 
    ((-2.0 * k_f[1] * k_f[2] * (T_f[1] - T_f[2]) * A_fx / ((k_f[1] + k_f[2]) * dx)) 
    + m_flow_unit * (h_in - h_f[1]) 
    - h_c[1] * (T_f[1] - T_s[1]) * A_fr) / dm_f;
    for i in 2:N_f - 1 loop
      der_h_f[i] := 
      (2.0 * k_f[i - 1] * k_f[i] * (T_f[i - 1] - T_f[i]) * A_fx / ((k_f[i - 1] + k_f[i]) * dx) 
      - 2.0 * k_f[i] * k_f[i + 1] * (T_f[i] - T_f[i + 1]) * A_fx / ((k_f[i] + k_f[i + 1]) * dx) 
      + m_flow_unit * (h_f[i - 1] - h_f[i]) 
      - h_c[i] * (T_f[i] - T_s[i]) * A_fr) / dm_f;
    end for;
    der_h_f[N_f] := 
    (2.0 * k_f[N_f - 1] * k_f[N_f] * (T_f[N_f - 1] - T_f[N_f]) * A_fx / ((k_f[N_f - 1] + k_f[N_f]) * dx) 
    + m_flow_unit * (h_f[N_f - 1] - h_f[N_f]) 
    - h_c[N_f] * (T_f[N_f] - T_s[N_f]) * A_fr) / dm_f;
    h_out := h_f[N_f];
  end if;
  
  /*
  for i in 1:N_f loop
    if Correlation == 1 then //Liquids
      Nu[i] := Nusselt_Liquid(Re[i],Pr[i],T_f[i],T_s[i],mu_f[i],mu_f_wall[i]);
    else //Gas
      Nu[i] := Nusselt_Gas(Re[i],Pr[i],T_f[i],T_s[i]);
    end if;
  end for;
  */
  /*
//Solid Equations  
  for i in 2:N_f-1 loop
//Innermost solid annulus
    der_h_p[i, 1] := (U_in[i, 1] * (T_s[i] - T_p[i, 1]) - U_in[i, 2] * (T_p[i, 1] - T_p[i, 2]) + U_right[i - 1, 1] * (T_p[i - 1, 1] - T_p[i, 1]) - U_right[i, 1] * (T_p[i, 1] - T_p[i + 1, 1])) / m_pj[1];
//Middle solid annuli
    for j in 2:N_p - 1 loop
      der_h_p[i, j] := (U_in[i, j] * (T_p[i, j - 1] - T_p[i, j]) - U_in[i, j + 1] * (T_p[i, j] - T_p[i, j + 1]) + U_right[i - 1, j] * (T_p[i - 1, j] - T_p[i, j]) - U_right[i, j] * (T_p[i, j] - T_p[i + 1, j])) / m_pj[j];
    end for;
//Outermost solid annulus
    der_h_p[i, N_p] := (U_in[i, N_p] * (T_p[i, N_p - 1] - T_p[i, N_p]) + U_right[i - 1, N_p] * (T_p[i - 1, N_p] - T_p[i, N_p]) - U_right[i, N_p] * (T_p[i, N_p] - T_p[i + 1, N_p]) - U_wall * (T_p[i, N_p] - T_amb) * CN.pi * D_solid * dx) / m_pj[N_p];
  end for;
//Left end of the pipe i = 1
//Innermost solid annulus
  der_h_p[1, 1] := (U_in[1, 1] * (T_s[1] - T_p[1, 1]) - U_in[1, 2] * (T_p[1, 1] - T_p[1, 2]) - U_right[1, 1] * (T_p[1, 1] - T_p[2, 1]) - U_side * (T_p[1, 1] - T_amb) * A_px[1]) / m_pj[1];
//Middle solid annuli
  for j in 2:N_p - 1 loop
    der_h_p[1, j] := (U_in[1, j] * (T_p[1, j - 1] - T_p[1, j]) - U_in[1, j + 1] * (T_p[1, j] - T_p[1, j + 1]) - U_right[1, j] * (T_p[1, j] - T_p[2, j]) - U_side * (T_p[1, j] - T_amb) * A_px[j]) / m_pj[j];
  end for;
//Outermost solid annulus
  der_h_p[1, N_p] := (U_in[1, N_p] * (T_p[1, N_p - 1] - T_p[1, N_p]) - U_right[1, N_p] * (T_p[1, N_p] - T_p[2, N_p]) - U_wall * (T_p[1, N_p] - T_amb) * CN.pi * D_solid * dx - U_side * (T_p[1, N_p] - T_amb) * A_px[N_p]) / m_pj[N_p];
//Right end of the pipe i = N_f
//Innermost solid annulus
  der_h_p[N_f, 1] := (U_in[N_f, 1] * (T_s[N_f] - T_p[N_f, 1]) - U_in[N_f, 2] * (T_p[N_f, 1] - T_p[N_f, 2]) + U_right[N_f - 1, 1] * (T_p[N_f - 1, 1] - T_p[N_f, 1]) - U_side * (T_p[N_f, 1] - T_amb) * A_px[1]) / m_pj[1];
//Middle solid annuli
  for j in 2:N_p - 1 loop
    der_h_p[N_f, j] := (U_in[N_f, j] * (T_p[N_f, j - 1] - T_p[N_f, j]) - U_in[N_f, j + 1] * (T_p[N_f, j] - T_p[N_f, j + 1]) + U_right[N_f - 1, j] * (T_p[N_f - 1, j] - T_p[N_f, j]) - U_side * (T_p[N_f, j] - T_amb) * A_px[j]) / m_pj[j];
  end for;
//Outermost solid annulus
  der_h_p[N_f, N_p] := (U_in[N_f, N_p] * (T_p[N_f, N_p - 1] - T_p[N_f, N_p]) + U_right[N_f - 1, N_p] * (T_p[N_f - 1, N_p] - T_p[N_f, N_p]) - U_wall * (T_p[N_f, N_p] - T_amb) * CN.pi * D_solid * dx - U_side * (T_p[N_f, N_p] - T_amb) * A_px[N_p]) / m_pj[N_p];
  */
  
//Charging (Mass flows top to bottom)
//Bottom Charging Fluid Node
//End Bottom Charging Fluid Node
//Middle Charging Fluid Nodes
//End Middle Charging Fluid Nodes
//Top Charging Fluid Node
//End Top Charging Fluid Node
//Discharge (Mass flows bottom to top)
//Bottom Discharge Node
//End Bottom Discharge Node
//Middle Discharge Nodes
//End Middle Discharge Nodes
//Top Discharge Node
//Re, Pr, Nu Calculation
//Convection Equations
/*
  for i in 1:N_f loop
    
    Bi[i] := Nu[i] * k_f[i] / (6.0 * k_p[i, N_p]);
//Use outermost shell conductivity
    Pe[i] := Re[i] * Pr[i];
    Pr[i] := max(1.0e-8,c_pf[i] * mu_f[i] / k_f[i]);
    Re[i] := max(10.0, rho_f[i] * D_pipe * abs(u_flow[i]) / mu_f[i]);
    f_low[i] := 64.0 / 2000.0;
    f_high[i] := (1.82 * log10(4000.0) - 1.64) ^ (-2.0);
    Nu_low[i] := 3.66 * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
    
    
    if Correlation == 1 then
//Liquids
      if T_s[i] > T_f[i] then
//Heating the Fluid
        Nu_high[i] := f_high[i] / 8.0 * 4000.0 * Pr[i] / (1.07 + 12.7 * (f_high[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0)) * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
      else
//Cooling the Fluid
        Nu_high[i] := f_high[i] / 8.0 * 4000.0 * Pr[i] / (1.07 + 12.7 * (f_high[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0)) * (mu_f[i] / mu_f_wall[i]) ^ 0.25;
      end if;
    else
//Air Humble,Lowdermilk, Desmon
      if T_s[i] > T_f[i] + 0.1 then
//Heating the Fluid
        Nu_high[i] := 0.023 * 4000.0 ^ 0.8 * (Pr[i] ^ 0.4) * ((max(1.0e-8,T_s[i] / T_f[i])) ^ (-0.55));
      else
//Cooling the Fluid, exponent is zero
        Nu_high[i] := 0.023 * 4000.0 ^ 0.8 * (Pr[i] ^ 0.4);
      end if;
    end if;
    if Re[i] < 2000.0 then
//Laminar
//if pre(flow_scheme[i]) == 1 then
      f_fric[i] := 64.0 / Re[i];
      Nu[i] := 3.66 * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
    elseif Re[i] > 4000.0 then
//Turbulent
//elseif pre(flow_scheme[i]) == 3 then
      f_fric[i] := (1.82 * log10(Re[i]) - 1.64) ^ (-2.0);
      if Correlation == 1 then
//Liquids
        if T_s[i] > T_f[i] then
//Heating the Fluid
          Nu[i] := f_fric[i] / 8.0 * Re[i] * Pr[i] / (1.07 + 12.7 * (f_fric[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0)) * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
        else
//Cooling the Fluid
          Nu[i] := f_fric[i] / 8.0 * Re[i] * Pr[i] / (1.07 + 12.7 * (f_fric[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0)) * (mu_f[i] / mu_f_wall[i]) ^ 0.25;
        end if;
      else
        if T_s[i] > T_f[i] + 0.1 then
//Heating the Fluid
          Nu[i] := 0.023 * (Re[i] ^ 0.8) * (Pr[i] ^ 0.4) * ((max(1.0e-8,T_s[i] / T_f[i])) ^ (-0.55));
        else
//Cooling the Fluid
          Nu[i] := 0.023 * (Re[i] ^ 0.8) * (Pr[i] ^ 0.4);
        end if;
      end if;
    else
      f_fric[i] := f_low[i] + (f_high[i] - f_low[i]) * SolarTherm.Utilities.Phis((Re[i] - 2000.0) / 2000.0);
      Nu[i] := Nu_low[i] + (Nu_high[i] - Nu_low[i]) * SolarTherm.Utilities.Phis((Re[i] - 2000.0) / 2000.0);
    end if;
    
    h_c[i] := Nu[i] * k_f[i] / D_pipe;
    
    
  end for;
*/

/*
  for i in 1:N_f loop
    Pr[i] := c_pf[i] * mu_f[i] / k_f[i];
    if abs(u_flow) > 1e-12 then
      Re[i] := rho_f_avg * D_pipe * abs(u_flow) / mu_f[i];
    else
      Re[i] := 1.0; //Just avoid zero
    end if;
    //
    if Re[i] <= 2000.0 then //Laminar
      f_fric[i] := 64.0/max(Re[i],1.0); //just avoid divison by zero. f_fric does not affect Nu
      Nu[i] := 4.36;
    elseif Re[i] >= 4000.0 then //turbulent
      f_fric[i] := (1.82*log10(Re[i]) - 1.64)^(-2.0);
      Nu[i] := (f_fric[i]/8.0)*(Re[i] - 1000.0)*Pr[i]/(1.0 + 12.7*((f_fric[i]/8.0)^(0.5))*(((Pr[i])^0.66)-1.0));
    else //transition
      f_fric[i] := 64.0/Re[i] + SolarTherm.Utilities.Phis((Re[i]-2.0e3)/2.0e3)*((1.82*log10(Re[i]) - 1.64)^(-2.0) - 64.0/Re[i]);
      Nu[i] := 4.36 + SolarTherm.Utilities.Phis((Re[i]-2.0e3)/2.0e3)*((f_fric[i]/8.0)*(Re[i] - 1000.0)*Pr[i]/(1.0 + 12.7*((f_fric[i]/8.0)^ (0.5))*(((Pr[i])^0.66) -1.0)) - 4.36);
    end if;
  end for;
  */
initial equation
  for i in 1:N_f loop
    fluid[i].h = h_f_start[i];
    for j in 1:N_p loop
      filler[i, j].h = h_p_start[i, j];
    end for;
//encapsulation[i].h = h_p_start[i,N_p];
  end for;
equation
  m_flow_unit = m_flow/Multiplier;
  for i in 1:N_f loop
    der_h_f[i] = der(h_f[i]);
    u_flow[i] = m_flow_unit / (rho_f[i] * A_fx);
  end for;
  
  //for i in 1:N_f loop
    //for j in 1:N_p loop
      //der_h_p[i,j] = der(h_p[i,j]);
    //end for;
  //end for;
//Determine which operational state: In this version, standby and discharge are lumped.
  if m_flow < 0.0 then
//mass is flowing downwards so charging
    State = 1;
  else
//mass is flowing upwards so discharging
    State = 3;
  end if;
  //u_flow = m_flow / (rho_f_avg * A_fx);
//positive if flowing upwards (discharge)
//u_0 = u_flow*eta; //Velocity through empty cross-section
//Fluid inlet and outlet properties
  fluid_in.h = h_in;
  fluid_out.h = h_out;
  fluid_in.T = T_in;
  fluid_out.T = T_out;
//Wall surface properties for convection
  for i in 1:N_f loop
    fluid_wall[i].T = T_s[i];
    //mu_f_wall[i] = max(1.0e-8,fluid_wall[i].mu);
  end for;
//Fluid Property evaluation SolarSalt
  for i in 1:N_f loop
    h_f[i] = fluid[i].h;
    T_f[i] = fluid[i].T;
    c_pf[i] = fluid[i].cp;
    k_f[i] = fluid[i].k;
    rho_f[i] = fluid[i].rho; //added this for stationary momentum (SM) assumption
//k_eff[i] = eta*fluid[i].k; //Effective thermal conductivity of fluid (weighted by porosity)
    mu_f[i] = max(1.0e-8,fluid[i].mu);
    mu_f_wall[i] = max(1.0e-8,fluid_wall[i].mu);
  end for;
//Particle Property evaluation quartzite and sand
  for i in 1:N_f loop
    for j in 1:N_p loop
/*
      //Function-based calculation of filler state
      f_p[i,j] = Filler_Package.f_h(h_p[i,j]);
      T_p[i,j] = Filler_Package.T_h(h_p[i,j]);
      k_p[i,j] = Filler_Package.k_Tf(T_p[i,j],f_p[i,j]); 
      */
      filler[i, j].h = h_p[i, j];
      T_p[i, j] = filler[i, j].T;
//f_p[i,j] = filler[i,j].f;
      k_p[i, j] = filler[i, j].k;
    end for;
//encapsulation property evaluation
//encapsulation[i].h = h_p[i,N_p];
//T_p[i,N_p] = encapsulation[i].T;
//k_p[i,N_p] = encapsulation[i].k;
//f_p[i,N_p] = encapsulation[i].f;
  end for;

  
//Convection Equations
  for i in 1:N_f loop
    
    Bi[i] = Nu[i] * k_f[i] / (6.0 * k_p[i, N_p]);
//Use outermost shell conductivity
    Pe[i] = Re[i] * Pr[i];
    Pr[i] = max(1.0e-8,c_pf[i] * mu_f[i] / k_f[i]);
    Re[i] = max(10.0, rho_f[i] * D_pipe * abs(u_flow[i]) / mu_f[i]);
    
    
    if Correlation == 1 then //Liquids
      Nu[i] = SolarTherm.Utilities.Nusselt.Internal_Flow.Nusselt_Liquid(Re[i],Pr[i],T_f[i],T_s[i],mu_f[i],mu_f_wall[i]);
    else //Gas
      Nu[i] = SolarTherm.Utilities.Nusselt.Internal_Flow.Nusselt_Gas(Re[i],Pr[i],T_f[i],T_s[i]);
    end if;
    
    /*
    f_low[i] = 64.0 / 2000.0;
    f_high[i] = (1.82 * log10(4000.0) - 1.64) ^ (-2.0);
    Nu_low[i] = 3.66 * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
    
    
    if Correlation == 1 then
//Liquids
      if T_s[i] > T_f[i] then
//Heating the Fluid
        Nu_high[i] = (((f_high[i] / 8.0) * 4000.0 * Pr[i] )/ (1.07 + 12.7 * (f_high[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0))) * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
      else
//Cooling the Fluid
        Nu_high[i] = (((f_high[i] / 8.0) * 4000.0 * Pr[i] )/ (1.07 + 12.7 * (f_high[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0))) * (mu_f[i] / mu_f_wall[i]) ^ 0.25;
      end if;
    else
//Air Humble,Lowdermilk, Desmon
      if T_s[i] > T_f[i] + 0.01 then
//Heating the Fluid
        Nu_high[i] = 0.023 * (4000.0 ^ 0.8) * (Pr[i] ^ 0.4) * ((max(1.0e-8,T_s[i] / T_f[i])) ^ (-0.55));
      else
//Cooling the Fluid, exponent is zero
        Nu_high[i] = 0.023 * (4000.0 ^ 0.8) * (Pr[i] ^ 0.4);
      end if;
    end if;
    if Re[i] < 2000.0 then
//Laminar
//if pre(flow_scheme[i]) == 1 then
      f_fric[i] = 64.0 / Re[i];
      Nu[i] = 3.66 * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
    elseif Re[i] > 4000.0 then
//Turbulent
//elseif pre(flow_scheme[i]) == 3 then
      f_fric[i] = (1.82 * log10(Re[i]) - 1.64) ^ (-2.0);
      if Correlation == 1 then
//Liquids
        if T_s[i] > T_f[i] then
//Heating the Fluid
          Nu[i] = (((f_fric[i] / 8.0) * Re[i] * Pr[i]) / (1.07 + 12.7 * (f_fric[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0))) * (mu_f[i] / mu_f_wall[i]) ^ 0.11;
        else
//Cooling the Fluid
          Nu[i] = (((f_fric[i] / 8.0) * Re[i] * Pr[i]) / (1.07 + 12.7 * (f_fric[i] / 8.0) ^ 0.5 * (Pr[i] ^ 0.667 - 1.0))) * (mu_f[i] / mu_f_wall[i]) ^ 0.25;
        end if;
      else
        if T_s[i] > T_f[i] + 0.01 then
//Heating the Fluid
          Nu[i] = 0.023 * (Re[i] ^ 0.8) * (Pr[i] ^ 0.4) * ((max(1.0e-8,T_s[i] / T_f[i])) ^ (-0.55));
        else
//Cooling the Fluid
          Nu[i] = 0.023 * (Re[i] ^ 0.8) * (Pr[i] ^ 0.4);
        end if;
      end if;
    else
      f_fric[i] = f_low[i] + (f_high[i] - f_low[i]) * SolarTherm.Utilities.Phis((Re[i] - 2000.0) / 2000.0);
      Nu[i] = Nu_low[i] + (Nu_high[i] - Nu_low[i]) * SolarTherm.Utilities.Phis((Re[i] - 2000.0) / 2000.0);
    end if;
    */

    h_c[i] = f_hc * Nu[i] * k_f[i] / D_pipe;
   
    
  end for;

  //Solid Equations
  //Radial Conductance
  for i in 1:N_f loop
    U_in[i, 1] = 2.0 * CN.pi * dx / (log(r_pj[1] / (r_pj[1] - 0.5 * dr)) / k_p[i, 1]);
    for j in 2:N_p loop
      U_in[i, j] = 2.0 * CN.pi * dx / (log((r_pj[j] - 0.5 * dr) / (r_pj[j] - dr)) / k_p[i, j - 1] + log(r_pj[j] / (r_pj[j] - 0.5 * dr)) / k_p[i, j]);
    end for;
//
    U_in[i, 1] * (T_s[i] - T_p[i, 1]) = h_c[i] * A_fr * (T_f[i] - T_s[i]); //links fluid and filler equations
    //reconciles fluid and filler heat transfer
  end for;
//Axial Conductance
  for j in 1:N_p loop
    for i in 1:N_f - 1 loop
      U_right[i, j] = 2.0 * A_px[j] * k_p[i, j] * k_p[i + 1, j] / (dx * (k_p[i, j] + k_p[i + 1, j]));
    end for;
    U_right[N_f, j] = 2.0 * k_p[N_f, j] * A_px[j] / dx;
  end for;
  
  for i in 2:N_f - 1 loop
//not the ends of the pipe
//Innermost solid annulus
    der(h_p[i, 1]) = (U_in[i, 1] * (T_s[i] - T_p[i, 1]) - U_in[i, 2] * (T_p[i, 1] - T_p[i, 2]) + U_right[i - 1, 1] * (T_p[i - 1, 1] - T_p[i, 1]) - U_right[i, 1] * (T_p[i, 1] - T_p[i + 1, 1])) / m_pj[1];
//Middle solid annuli
    for j in 2:N_p - 1 loop
      der(h_p[i, j]) = (U_in[i, j] * (T_p[i, j - 1] - T_p[i, j]) - U_in[i, j + 1] * (T_p[i, j] - T_p[i, j + 1]) + U_right[i - 1, j] * (T_p[i - 1, j] - T_p[i, j]) - U_right[i, j] * (T_p[i, j] - T_p[i + 1, j])) / m_pj[j];
    end for;
//Outermost solid annulus
    der(h_p[i, N_p]) = (U_in[i, N_p] * (T_p[i, N_p - 1] - T_p[i, N_p]) + U_right[i - 1, N_p] * (T_p[i - 1, N_p] - T_p[i, N_p]) - U_right[i, N_p] * (T_p[i, N_p] - T_p[i + 1, N_p]) - U_wall * (T_p[i, N_p] - T_amb) * CN.pi * D_solid * dx) / m_pj[N_p];
  end for;
//Left end of the pipe i = 1
//Innermost solid annulus
  der(h_p[1, 1]) = (U_in[1, 1] * (T_s[1] - T_p[1, 1]) - U_in[1, 2] * (T_p[1, 1] - T_p[1, 2]) - U_right[1, 1] * (T_p[1, 1] - T_p[2, 1]) - U_side * (T_p[1, 1] - T_amb) * A_px[1]) / m_pj[1];
//Middle solid annuli
  for j in 2:N_p - 1 loop
    der(h_p[1, j]) = (U_in[1, j] * (T_p[1, j - 1] - T_p[1, j]) - U_in[1, j + 1] * (T_p[1, j] - T_p[1, j + 1]) - U_right[1, j] * (T_p[1, j] - T_p[2, j]) - U_side * (T_p[1, j] - T_amb) * A_px[j]) / m_pj[j];
  end for;
//Outermost solid annulus
  der(h_p[1, N_p]) = (U_in[1, N_p] * (T_p[1, N_p - 1] - T_p[1, N_p]) - U_right[1, N_p] * (T_p[1, N_p] - T_p[2, N_p]) - U_wall * (T_p[1, N_p] - T_amb) * CN.pi * D_solid * dx - U_side * (T_p[1, N_p] - T_amb) * A_px[N_p]) / m_pj[N_p];
//Right end of the pipe i = N_f
//Innermost solid annulus
  der(h_p[N_f, 1]) = (U_in[N_f, 1] * (T_s[N_f] - T_p[N_f, 1]) - U_in[N_f, 2] * (T_p[N_f, 1] - T_p[N_f, 2]) + U_right[N_f - 1, 1] * (T_p[N_f - 1, 1] - T_p[N_f, 1]) - U_side * (T_p[N_f, 1] - T_amb) * A_px[1]) / m_pj[1];
//Middle solid annuli
  for j in 2:N_p - 1 loop
    der(h_p[N_f, j]) = (U_in[N_f, j] * (T_p[N_f, j - 1] - T_p[N_f, j]) - U_in[N_f, j + 1] * (T_p[N_f, j] - T_p[N_f, j + 1]) + U_right[N_f - 1, j] * (T_p[N_f - 1, j] - T_p[N_f, j]) - U_side * (T_p[N_f, j] - T_amb) * A_px[j]) / m_pj[j];
  end for;
//Outermost solid annulus
  der(h_p[N_f, N_p]) = (U_in[N_f, N_p] * (T_p[N_f, N_p - 1] - T_p[N_f, N_p]) + U_right[N_f - 1, N_p] * (T_p[N_f - 1, N_p] - T_p[N_f, N_p]) - U_wall * (T_p[N_f, N_p] - T_amb) * CN.pi * D_solid * dx - U_side * (T_p[N_f, N_p] - T_amb) * A_px[N_p]) / m_pj[N_p];
  
//Particle Surface equations energy balance
//Heat loss calculations, different form than the equations above as they were in terms of rho*dh/dt not m*dh/dt
//Q_loss_top = U_top*CN.pi*D_tank*D_tank*0.25*(T_f[N_f]-T_amb);
//Q_loss_bot = U_bot*CN.pi*D_tank*D_tank*0.25*(T_f[1]-T_amb);
//for i in 2:N_f-1 loop
//Q_loss_wall[i] = U_wall*CN.pi*D_tank*(T_f[i]-T_amb)*dx;
//end for;
//Q_loss_wall[1] = U_wall*CN.pi*D_tank*(T_f[1]-T_amb)*dx;
//Q_loss_wall[N_f] = U_wall*CN.pi*D_tank*(T_f[N_f]-T_amb)*dx;
//Q_loss_total = Q_loss_top + sum(Q_loss_wall) + Q_loss_bot;
//End heat loss calculations
//Calculated Pumping losses
//for i in 1:N_f loop
//p_drop[i] = dx*(((600*((1-eta)^2)*mu_f[i]*abs(m_flow))/((eta^3)*(d_p^2)*rho_f_avg*CN.pi*(D_tank^2)))+((28*(1-eta)*(m_flow^2))/((eta^3)*d_p*rho_f_avg*CN.pi*CN.pi*(D_tank^4))));
//end for;
//p_drop_total = sum(p_drop);
//W_loss_pump = (abs(m_flow)/rho_f_avg)*p_drop_total/eff_pump;
//Analyics
  der(E_stored) = Multiplier*abs(m_flow_unit) * (h_in - h_out);
  Level = E_stored / (Multiplier*E_unit);
  if m_flow > 1.0e-3 then
//Discharging, outlet is the top
    T_outlet_degC = T_f[N_f] - 273.15;
  elseif m_flow < (-1.0e-3) then
//Charging, outlet is the bottom
    T_outlet_degC = T_f[1] - 273.15;
  else
//No flow, output reference temperature
    T_outlet_degC = 25.0;
  end if;
  annotation(
    Documentation(revisions = "<html>
		<p>By Zebedee Kee on 03/12/2020</p>
		</html>", info = "<html>
		<p>This model contains the heat-transfer calculations of a thermocline packed bed storage tank with spherical filler geometry. This model does not contain any fluid connectors, for the CSP component with connectors, see Thermocline_Spheres_SingleTank. Variables fluid_top and fluid_bot provides the enthalpy-temperature relationship of the fluid material. Depending on whether m_flow is positive (discharging, fluid flowing upwards) or negative (charging, fluid flowing downwards), the charging/discharging equations are applied. In this iteration of the model, discharging and standby are lumped into one state.</p>
		</html>"));
end Annular_Storage_Section_SM;
