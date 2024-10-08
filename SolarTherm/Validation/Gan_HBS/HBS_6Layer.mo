within SolarTherm.Validation.Gan_HBS;

model HBS_6Layer
  //Thermocline tank with spherical fillers - Heat Transfer Model
  //Lumped Capacitance Model, suitable for air fluid - rock filler systems
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  import Tables = Modelica.Blocks.Tables;
  //Initialize Material Packages
  replaceable package Fluid_Package = SolarTherm.Materials.PartialMaterial "Fluid Package";
  replaceable package Filler_Package_1 = SolarTherm.Materials.PartialMaterial "Filler Package layer 1";
  replaceable package Filler_Package_2 = SolarTherm.Materials.PartialMaterial "Filler Package layer 1";
  replaceable package Filler_Package_3 = SolarTherm.Materials.PartialMaterial "Filler Package layer 1";
  replaceable package Filler_Package_4 = SolarTherm.Materials.PartialMaterial "Filler Package layer 1";
  replaceable package Filler_Package_5 = SolarTherm.Materials.PartialMaterial "Filler Package layer 1";
  replaceable package Filler_Package_6 = SolarTherm.Materials.PartialMaterial "Filler Package layer 1";
  replaceable package Tank_Package = SolarTherm.Materials.SS316L "Tank Package (steel shell)";
  //replaceable package Encapsulation_Package = Filler_Package "Encapsulation Package, default is the same as Filler package, effectively no encapsulation";
  //Fluid Material States
  Fluid_Package.State fluid_in "Model which calculates properties at inlet of the section";
  Fluid_Package.State fluid_out "Model which calculates properties at outlet of the section";
  //Interfacial heat transfer Settings
  parameter Integer Correlation = 1 "Internal Flow, Gas";
  //Height offset for plotting purposes
  parameter SI.Length z_offset = 0.0 "Amount of height offset if there is a tank below it";
  //Tank Design parameters
  parameter SI.Energy E_max = 144e9 "Design storage capacity";
  parameter Real ar = 2.0 "Tank Aspect ratio H/D";
  parameter Real eta = 0.40 "Structural porosity of the firebrick (m)";
  parameter Real d_p = 0.02 "Diameter of the pores (m)";
  //parameter Real s_p = 0.04 "Separation of the pores (m)";
  //parameter Real t_e = d_p / (2.0 * N_p) "Thickness of encapsulation, default is such that it is at a value that preserves equidistant radii discretizations (m)";
  //Temperature Bounds
  parameter SI.Temperature T_min = CV.from_degC(620) "Design cold Temperature of everything in the tank (K)";
  parameter SI.Temperature T_max = CV.from_degC(820) "Design hot Temperature of everything in the tank (K)";
  parameter SI.Temperature T_start = T_min "Initial (uniform) temperature of all components (K), defaults to T_min";
  //Calculated Tank Design Parameters
  parameter SI.Length H_tank = (4 * E_max / (CN.pi * (1 / ar) ^ 2 * (rho_f_avg * (h_f_max - h_f_min) * eta + rho_p * (h_p_max - h_p_min) * (1.0 - eta)))) ^ (1 / 3);
  parameter SI.Diameter D_tank = H_tank / ar;
  parameter SI.Area A = CN.pi * D_tank * D_tank / 4.0 "Cross sectional area of tank";
  //Thermal Losses
  SI.Temperature T_amb;
  parameter SI.Area A_loss_tank = CN.pi * D_tank * D_tank * 0.5 + CN.pi * D_tank * H_tank "Heat loss area (m2)";
  parameter SI.CoefficientOfHeatTransfer U_loss_tank = 0.1 "Heat loss coeff of surfaces (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_wall = U_loss_tank "Cylinder wall heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_top = U_loss_tank "Top circle heat loss coeff (W/m2K)";
  parameter SI.CoefficientOfHeatTransfer U_bot = U_loss_tank "Bottom circle heat loss coeff (W/m2K)";
  //Inititalize temperature and enthalpy profile
  parameter SI.Temperature T_f_start[N_f] = fill(T_start, N_f);
  parameter SI.Temperature T_p_start[N_f] = fill(T_start, N_f);
  //parameter SI.Temperature T_e_start = T_start;
  parameter SI.SpecificEnthalpy h_f_start[N_f] = fill(Fluid_Package.h_Tf(T_start, 0.0), N_f) "Defaults to uniform";
  parameter SI.SpecificEnthalpy h_p_start[N_f];// = fill(Filler_Package.h_Tf(T_start, 0.0), N_f) "Defaults to uniform";
  //Property bounds
  //Fluid
  parameter SI.SpecificEnthalpy h_f_min = Fluid_Package.h_Tf(T_min, 0) "Starting enthalpy of the HTF";
  parameter SI.SpecificEnthalpy h_f_max = Fluid_Package.h_Tf(T_max, 0) "Starting enthalpy of the HTF";
  parameter SI.Density rho_f_min = Fluid_Package.rho_Tf(T_min, 0);
  parameter SI.Density rho_f_max = Fluid_Package.rho_Tf(T_max, 0);
  parameter SI.Density rho_f_avg = (rho_f_min + rho_f_max) / 2;
  //Filler
  parameter SI.SpecificEnthalpy h_p_max = Filler_Package_1.h_Tf(T_max, 1.0);
  parameter SI.SpecificEnthalpy h_p_min = Filler_Package_1.h_Tf(T_min, 0.0);
  
  parameter SI.Density rho_p1_min = Filler_Package_1.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p1_max = Filler_Package_1.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p1 = min(rho_p1_min, rho_p1_max) "kg/m3";
  parameter SI.Density rho_p2_min = Filler_Package_2.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p2_max = Filler_Package_2.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p2 = min(rho_p2_min, rho_p2_max) "kg/m3";
  parameter SI.Density rho_p3_min = Filler_Package_3.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p3_max = Filler_Package_3.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p3 = min(rho_p3_min, rho_p3_max) "kg/m3";
  parameter SI.Density rho_p4_min = Filler_Package_4.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p4_max = Filler_Package_4.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p4 = min(rho_p4_min, rho_p4_max) "kg/m3";
  parameter SI.Density rho_p5_min = Filler_Package_5.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p5_max = Filler_Package_5.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p5 = min(rho_p5_min, rho_p5_max) "kg/m3";
  parameter SI.Density rho_p6_min = Filler_Package_6.rho_Tf(T_min, 0.0);
  parameter SI.Density rho_p6_max = Filler_Package_6.rho_Tf(T_max, 1.0);
  parameter SI.Density rho_p6 = min(rho_p6_min, rho_p6_max) "kg/m3";
  //Encapsulation
  //parameter SI.SpecificEnthalpy h_e_max = Encapsulation_Package.h_Tf(T_max, 1.0);
  //parameter SI.SpecificEnthalpy h_e_min = Encapsulation_Package.h_Tf(T_min, 0.0);
  //parameter SI.Density rho_e_min = Encapsulation_Package.rho_Tf(T_min, 0.0);
  //parameter SI.Density rho_e_max = Encapsulation_Package.rho_Tf(T_max, 1.0);
  //parameter SI.Density rho_e = min(rho_e_min, rho_e_max) "kg/m3";
  //Discretization
  parameter SI.Length dz = H_tank / N_f "discretization vertical length of fluid";
  //parameter SI.Length dr[N_p] = cat(1,fill(0.5 * ((d_p - 2.0*t_e) / (N_p - 1)),N_p-1),{t_e}) "radial thickness of each particle discretization, with last one being the encapsulation";
  parameter Integer N_f = 100 "Number of finite volume elements in fluid";
  parameter Integer N_p = 1 "Only One Element";
  //Initialise Fluid Array
  parameter SI.Length z_f[N_f] = SolarTherm.Models.Storage.Thermocline.Z_position(H_tank, N_f) .+ z_offset;
  SI.Temperature T_f[N_f] "(K)";
  SI.SpecificEnthalpy h_f[N_f];
  //(start = h_f_start) "J/kg";
  //Plotting
  parameter Real ZDH[N_f] = SolarTherm.Models.Storage.Thermocline.Relative_Tank_Axes(H_tank, N_f) "Non-dimensional tank vertical axis";
  //Operational Controls
  Integer State(start = 2) "operational state 2=standby, 3=discharge, 1=charge";
  //Inlet and outlet enthalpies and temperatures
  SI.SpecificEnthalpy h_in "Inlet Enthalpy depends on mass flow direction (J/kg)";
  SI.SpecificEnthalpy h_out "Outlet Enthalpy depends on mass flow direction (J/kg)";
  SI.Temperature T_in "Inlet Temperature depends on mass flow direction";
  SI.Temperature T_out "Outlet Temperature depends on mass flow direction";
  //Mass flow rates and superficial velocity
  SI.MassFlowRate m_flow(start = 0.0) "kg/s";
  SI.Velocity u_flow "m/s";
  //SI.Velocity u_0 "Fluid velocity through empty cross section = u_flow/eta (m/s)";
  //Analytics
  SI.Energy E_stored(start = 0.0) "Make sure the tank starts from T_min for this to be correct";
  Real Level(start = 0.0) "Tank energy charging level (0-1)";
  SI.HeatFlowRate Q_loss_total "Heat loss from the entire surface area";
  
  //Initialise Encapsulation
  //SI.Temperature T_e[N_f](start=T_e_start) "Temperature of encapsulation elements";
  //SI.SpecificEnthalpy h_e[N_f](start = h_e_min) "J/kg";
  //SI.ThermalConductivity k_e[N_f] "W/mK";
  //Calculated Pumping Losses
  parameter Real eff_pump = 1.00 "Pump electricity to work efficiency";
  SI.Pressure p_drop_total "Sum of all pressure drops";
  SI.Power W_loss_pump "losses due to pressure drop";
  //Cost breakdown
  parameter Real C_fluid=0.0;// = max(rho_f_max, rho_f_min) * eta * (CN.pi * D_tank * D_tank * H_tank / 4.0) * Fluid_Package.cost;
  parameter Real C_section = C_fluid + C_filler + C_insulation + C_tank + C_encapsulation;
  //parameter Real C_insulation = if U_loss_tank > 1e-3 then (16.72/U_loss_tank + 0.04269)*A_loss_tank else 0.0;
  parameter Real C_insulation =0.0;//= if U_loss_tank > 1e-3 then CpA_external_insulation(T_max, U_loss_tank) * A_loss_tank else 0.0;
  parameter Real C_tank = 0.0;// C_shell(max(rho_f_max, rho_f_min), H_tank, D_tank, Tank_Package.sigma_yield(T_max), Tank_Package.rho_Tf(298.15, 0.0), 4.0);
  parameter Real C_filler= 0.0;// = rho_p * (1.0 - eta) * (CN.pi * D_tank * D_tank * H_tank / 4.0) * Filler_Package.cost;
  parameter Real C_encapsulation = 0.0;
  //Filler Surface Area Correction
  parameter Real f_surface = 1.0 "Don't touch this";
  //Initialise Filler surface temperature
  SI.Temperature T_s[N_f](start = T_f_start);
  //parameter SI.Length r_p[N_p] = cat(1,Particle_Radii(d_p-2*t_e,N_p-1),{(d_p/2)-(t_e/2)}) "Radii of each particle element centre";
  parameter SI.Volume V_p[N_f] = fill((1.0-eta)*CN.pi*D_tank*D_tank*H_tank/(4.0*N_f), N_f) "Solid volume of each filler slice (m3)";
  parameter SI.Density rho_p[N_f] = cat(1,fill(rho_p1,5), fill(rho_p2,37),fill(rho_p3,16),fill(rho_p4,8), fill(rho_p5,6), fill(rho_p6,28));
  parameter SI.Mass m_p[N_f] = V_p.*rho_p;
  
  parameter SI.Area A_loss_top = 0.25*CN.pi*D_tank*D_tank "Surface area of the top and bottom of the tank for heat loss calculations (m2)";
  parameter SI.Area A_loss_wall_i = CN.pi*D_tank*dz "Surface area of the side of each element for heat loss calculations (m2)";
  //parameter SI.Mass m_p[N_f] = fill((1.0-eta)*CN.pi*D_tank*D_tank*H_tank*rho_p/(4.0*N_f), N_f) "Masses of each particle";
  SI.SpecificEnthalpy h_p[N_f] "J/kg";
  SI.Temperature T_p[N_f](start = T_p_start) "Temperature of particle elements";

  //Initialise Particle
  
  //Filler mass-liquid fraction
  Real f_p[N_f] "Mass liquid fraction of filler";
  Real Bi[N_f] "Biot Number";
  //Convection Properties
  Real Re[N_f] "Reynolds";
  Real Pr[N_f] "Prandtl";
  Real Nu[N_f] "Nusselt";
  //Real h_v[N_f] "Volumetric heat transfer coeff (W/m3K)";
  SI.CoefficientOfHeatTransfer h_c[N_f] "Convective heat transfer coefficient of internal flow (W/m2K)";
  //SI.ThermalConductance U_in[N_f, N_p] "K/W"; //Obsolete
  //SI.ThermalConductance U_out[N_f, N_p] "K/W";
  //SI.ThermalConductance U_out_e[N_f] "K/W";
  SI.ThermalConductance U_up[N_f] "Thermal conductance between a filler element and the one on top of it (K/W)";
  //Filler Properties
  
  SI.ThermalConductivity k_p[N_f] "W/mK";
  //Filler Geometry
  //parameter Real N_spheres_total = N_f * 6 * (1 - eta) * A * dz / (CN.pi * d_p ^ 3) "Total number of spheres in the tank";
  
  //parameter SI.Mass m_e = rho_e*(1/6)*CN.pi*((d_p^3)-((d_p-2*t_e)^3)) "Masses of encapsulation in one particle";
  //Pressure Drop
  SI.Pressure p_drop[N_f] "Pressure drop across each mesh element";
  //Thermal Losses
  SI.HeatFlowRate Q_loss_wall[N_f] "Heat loss from the wall";
  SI.HeatFlowRate Q_loss_top "Heat loss from the top";
  SI.HeatFlowRate Q_loss_bot "Heat loss from the bottom";
  //Fluid Properties
  SI.ThermalConductivity k_f[N_f] "W/mK";
  //SI.ThermalConductivity k_eff[N_f] "W/mK";
  SI.DynamicViscosity mu_f[N_f] "Pa.s";
  SI.SpecificHeatCapacity c_pf[N_f] "J/kgK";
protected
  Fluid_Package.State fluid[N_f](each h_start = h_f_min) "Fluid object array";
  //Try filler state "Remove this if using function-based calculation"
  Filler_Package_1.State filler1[5] "Filler object array";
  Filler_Package_2.State filler2[37] "Filler object array";
  Filler_Package_3.State filler3[16] "Filler object array";
  Filler_Package_4.State filler4[8] "Filler object array";
  Filler_Package_5.State filler5[6] "Filler object array";
  Filler_Package_6.State filler6[28] "Filler object array";
  //.State encapsulation[N_f] "Encapsulation object array";
  Real der_h_f[N_f] "rate of change of h_f calculated explicitly";
algorithm
//Fluid Equations
  if State == 1 then
    der_h_f[1] := ((-2.0 * k_f[1] * k_f[2]) * (T_f[1] - T_f[2]) / ((k_f[1] + k_f[2]) * dz * dz) + rho_f_avg * u_flow * (h_f[1] - h_f[2]) / dz - 4.0 * h_c[1] * (T_f[1] - T_s[1]) / (eta * d_p) - U_bot * (T_f[1] - T_amb) / (eta * dz) - U_wall * CN.pi * D_tank * (T_f[1] - T_amb) / (eta * A)) / rho_f_avg;
    h_out := h_f[1];
    for i in 2:N_f - 1 loop
      der_h_f[i] := (2.0 * k_f[i - 1] * k_f[i] * (T_f[i - 1] - T_f[i]) / ((k_f[i - 1] + k_f[i]) * dz * dz) - 2.0 * k_f[i] * k_f[i + 1] * (T_f[i] - T_f[i + 1]) / ((k_f[i] + k_f[i + 1]) * dz * dz) + rho_f_avg * u_flow * (h_f[i] - h_f[i + 1]) / dz - 4.0 * h_c[i] * (T_f[i] - T_s[i]) / (eta * d_p) - U_wall * CN.pi * D_tank * (T_f[i] - T_amb) / (eta * A)) / rho_f_avg;
    end for;
    der_h_f[N_f] := (2.0 * k_f[N_f - 1] * k_f[N_f] * (T_f[N_f - 1] - T_f[N_f]) / ((k_f[N_f - 1] + k_f[N_f]) * dz * dz) + rho_f_avg * u_flow * (h_f[N_f] - h_in) / dz - 4.0 * h_c[N_f] * (T_f[N_f] - T_s[N_f]) / (eta * d_p) - U_wall * CN.pi * D_tank * (T_f[N_f] - T_amb) / (eta * A) - U_top * (T_f[N_f] - T_amb) / (eta * dz)) / rho_f_avg;
  else
    der_h_f[1] := ((-2.0 * k_f[1] * k_f[2] * (T_f[1] - T_f[2]) / ((k_f[1] + k_f[2]) * dz * dz)) + rho_f_avg * u_flow * (h_in - h_f[1]) / dz - 4.0 * h_c[1] * (T_f[1] - T_s[1]) / (eta * d_p) - U_bot * (T_f[1] - T_amb) / (eta * dz) - U_wall * CN.pi * D_tank * (T_f[1] - T_amb) / (eta * A)) / rho_f_avg;
    for i in 2:N_f - 1 loop
      der_h_f[i] := (2.0 * k_f[i - 1] * k_f[i] * (T_f[i - 1] - T_f[i]) / ((k_f[i - 1] + k_f[i]) * dz * dz) - 2.0 * k_f[i] * k_f[i + 1] * (T_f[i] - T_f[i + 1]) / ((k_f[i] + k_f[i + 1]) * dz * dz) + rho_f_avg * u_flow * (h_f[i - 1] - h_f[i]) / dz - 4.0 * h_c[i] * (T_f[i] - T_s[i]) / (eta * d_p) - U_wall * CN.pi * D_tank * (T_f[i] - T_amb) / (eta * A)) / rho_f_avg;
    end for;
    der_h_f[N_f] := (2.0 * k_f[N_f - 1] * k_f[N_f] * (T_f[N_f - 1] - T_f[N_f]) / ((k_f[N_f - 1] + k_f[N_f]) * dz * dz) + rho_f_avg * u_flow * (h_f[N_f - 1] - h_f[N_f]) / dz - 4.0 * h_c[N_f] * (T_f[N_f] - T_s[N_f]) / (eta * d_p) - U_wall * CN.pi * D_tank * (T_f[N_f] - T_amb) / (eta * A) - U_top * (T_f[N_f] - T_amb) / (eta * dz)) / rho_f_avg;
    h_out := h_f[N_f];
  end if;
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
initial equation
  for i in 1:N_f loop
    fluid[i].h = h_f_start[i];
//for j in 1:N_p-1 loop
    //filler[i].h = h_p_start[i];
//end for;
//encapsulation[i].h = h_p_start[i,N_p];
  end for;
  for i in 1:5 loop
    filler1[i].h = h_p_start[i];
  end for;
  for i in 6:42 loop
    filler2[i-5].h = h_p_start[i];
  end for;
  for i in 43:58 loop
    filler3[i-42].h = h_p_start[i];
  end for;
  for i in 59:66 loop
    filler4[i-58].h = h_p_start[i];
  end for;
  for i in 67:72 loop
    filler5[i-66].h = h_p_start[i];
  end for;
  for i in 73:100 loop
    filler6[i-72].h = h_p_start[i];
  end for;
  
  
  
  
equation
  for i in 1:N_f loop
    T_s[i] = T_p[i];
  end for;
//Fluid Equations
  for i in 1:N_f loop
    der(h_f[i]) = der_h_f[i];
  end for;
//Determine which operational state: In this version, standby and discharge are lumped.
  if m_flow < 0.0 then
//mass is flowing downwards so charging
    State = 1;
  else
//mass is flowing upwards so discharging
    State = 3;
  end if;
  u_flow = m_flow / (eta * rho_f_avg * A);
//positive if flowing upwards (discharge)
  //u_0 = u_flow * eta;
//Velocity through empty cross-section
//Fluid inlet and outlet properties
  fluid_in.h = h_in;
  fluid_out.h = h_out;
  fluid_in.T = T_in;
  fluid_out.T = T_out;
/*
  //Fluid Equations
  if State == 1 then
  //Charging (Mass flows top to bottom)
  //Bottom Charging Fluid Node
    rho_f_avg * der(h_f[1]) =
    (-2.0*k_f[1]*k_f[2])*(T_f[1]-T_f[2])/((k_f[1]+k_f[2])*dz*dz)
    + (rho_f_avg*u_flow)*(h_f[1]-h_f[2])/dz
    - h_v[1]*(T_f[1] - T_s[1])/eta
    - U_bot*(T_f[1]-T_amb)/(eta*dz) 
    - U_wall*CN.pi*D_tank*(T_f[1]-T_amb)/(eta*A);
    
    h_out = h_f[1];
  //End Bottom Charging Fluid Node
  //Middle Charging Fluid Nodes
    for i in 2:N_f - 1 loop
      rho_f_avg*der(h_f[i]) = 
      2.0*k_f[i - 1]*k_f[i]*(T_f[i-1]-T_f[i])/((k_f[i-1]+k_f[i])*dz*dz)
      - 2.0*k_f[i]*k_f[i+1]*(T_f[i]-T_f[i + 1])/((k_f[i]+k_f[i+1])*dz*dz)
      + (rho_f_avg*u_flow)*(h_f[i]-h_f[i+1])/dz
      - h_v[i]*(T_f[i]-T_s[i])/eta
      - U_wall*CN.pi*D_tank*(T_f[i]-T_amb)/(eta*A);
    end for;
  //End Middle Charging Fluid Nodes
  //Top Charging Fluid Node
    rho_f_avg*der(h_f[N_f]) = 
    2.0*k_f[N_f-1]*k_f[N_f]*(T_f[N_f-1]-T_f[N_f])/((k_f[N_f-1]+k_f[N_f])*dz*dz)
    + (rho_f_avg*u_flow)*(h_f[N_f]-h_in)/dz
    - h_v[N_f]*(T_f[N_f]-T_s[N_f])/eta
    - U_wall*CN.pi*D_tank*(T_f[N_f]-T_amb)/(eta*A)
    - U_top*(T_f[N_f]-T_amb)/(eta*dz);
  //End Top Charging Fluid Node
  else
  //Discharge (Mass flows bottom to top)
  //Bottom Discharge Node
    rho_f_avg*der(h_f[1]) =
    -2.0*k_f[1]*k_f[2]*(T_f[1]-T_f[2])/((k_f[1]+k_f[2])*dz*dz)
    + (rho_f_avg*u_flow)*(h_in-h_f[1])/dz
    - h_v[1]*(T_f[1]-T_s[1])/eta
    - U_bot*(T_f[1]-T_amb)/(eta*dz)
    - U_wall*CN.pi*D_tank*(T_f[1]-T_amb)/(eta*A);
  //End Bottom Discharge Node
  //Middle Discharge Nodes
    for i in 2:N_f - 1 loop
      rho_f_avg*der(h_f[i]) =
      2.0*k_f[i-1]*k_f[i]*(T_f[i-1]-T_f[i])/((k_f[i-1]+k_f[i])*dz*dz)
      - 2.0*k_f[i]*k_f[i + 1]*(T_f[i]-T_f[i+1])/((k_f[i]+k_f[i+1])*dz*dz)
      + (rho_f_avg*u_flow)*(h_f[i-1]-h_f[i])/dz
      - h_v[i]*(T_f[i]-T_s[i])/eta
      - U_wall*CN.pi*D_tank*(T_f[i]-T_amb)/(eta*A);
    end for;
  //End Middle Discharge Nodes
  //Top Discharge Node
    rho_f_avg*der(h_f[N_f]) =
    2.0*k_f[N_f-1]*k_f[N_f]*(T_f[N_f-1]-T_f[N_f])/((k_f[N_f-1]+k_f[N_f])*dz*dz)
    + (rho_f_avg*u_flow)*(h_f[N_f-1]-h_f[N_f])/dz
    - h_v[N_f]*(T_f[N_f]-T_s[N_f])/eta
    - U_wall*CN.pi*D_tank*(T_f[N_f]-T_amb)/(eta*A)
    - U_top*(T_f[N_f]-T_amb)/(eta*dz);
    
    h_out = h_f[N_f];
  end if;
  */
//Fluid Property evaluation SolarSalt
  for i in 1:N_f loop
    fluid[i].h = h_f[i];
    T_f[i] = fluid[i].T;
    c_pf[i] = fluid[i].cp;
    k_f[i] = fluid[i].k;
//k_eff[i] = eta*fluid[i].k; //Effective thermal conductivity of fluid (weighted by porosity)
    mu_f[i] = fluid[i].mu;
  end for;
//Particle Property evaluation
  for i in 1:5 loop
    filler1[i].h = h_p[i];
    T_p[i] = filler1[i].T;
    f_p[i] = filler1[i].f;
    k_p[i] = filler1[i].k;
  end for;
  for i in 6:42 loop
    filler2[i-5].h = h_p[i];
    T_p[i] = filler2[i-5].T;
    f_p[i] = filler2[i-5].f;
    k_p[i] = filler2[i-5].k;
  end for;
  for i in 43:58 loop
    filler3[i-42].h = h_p[i];
    T_p[i] = filler3[i-42].T;
    f_p[i] = filler3[i-42].f;
    k_p[i] = filler3[i-42].k;
  end for;
  for i in 59:66 loop
    filler4[i-58].h = h_p[i];
    T_p[i] = filler4[i-58].T;
    f_p[i] = filler4[i-58].f;
    k_p[i] = filler4[i-58].k;
  end for;
  for i in 67:72 loop
    filler5[i-66].h = h_p[i];
    T_p[i] = filler5[i-66].T;
    f_p[i] = filler5[i-66].f;
    k_p[i] = filler5[i-66].k;
  end for;
  for i in 73:100 loop
    filler6[i-72].h = h_p[i];
    T_p[i] = filler6[i-72].T;
    f_p[i] = filler6[i-72].f;
    k_p[i] = filler6[i-72].k;
  end for;
  
//Convection Equations
  for i in 1:N_f loop
    Bi[i] = Nu[i] * k_f[i] / (6.0 * k_p[i]);
    if abs(u_flow) > 1e-12 then
//There is actually mass flowing
      Re[i] = rho_f_avg * d_p * abs(u_flow) / mu_f[i];
      Pr[i] = c_pf[i] * mu_f[i] / k_f[i];
      Nu[i] = SolarTherm.Validation.Gan_HBS.Nusselt_Gas_Gan(Re[i],Pr[i]);//,T_f[i],T_s[i]);
    else
      Re[i] = 0;
      Pr[i] = 0;
      Nu[i] = 3.66;
    end if;
    //h_v[i] = f_surface * 6.0 * (1.0 - eta) * Nu[i] * k_f[i] / (d_p * d_p);
    h_c[i] = Nu[i]*k_f[i]/d_p;
    
//Note that filler surface area correction factor is applied elsewhere.
  end for;
//Particle Equations
  //Vertical Thermal Conductance
  for i in 1:N_f-1 loop
    U_up[i] = 0.5*k_p[i]*k_p[i+1]*(1.0-eta)*CN.pi*D_tank*D_tank/((k_p[i]+k_p[i+1])*dz);
  end for;
  U_up[N_f] = 0.0; //There is nothing on top of it 


  m_p[1] * der(h_p[1]) = -1.0*U_up[1]*(T_p[1]-T_p[2]) +
    h_c[1] * D_tank * D_tank * CN.pi * dz * (T_f[1] - T_s[1]) / d_p - U_wall*A_loss_wall_i*(T_p[1]-T_amb) - U_bot*A_loss_top*(T_p[1]-T_amb);
  for i in 2:N_f-1 loop
    //m_p[i] * der(h_p[i]) = CN.pi * d_p ^ 3 * h_v[i] / (6.0 * (1.0 - eta)) * (T_f[i] - T_s[i]);
    m_p[i] * der(h_p[i]) = U_up[i-1]*(T_p[i-1]-T_p[i]) - U_up[i]*(T_p[i]-T_p[i+1]) + h_c[i] * D_tank * D_tank * CN.pi * dz * (T_f[i] - T_s[i]) / d_p;
  end for;
  m_p[N_f] * der(h_p[N_f]) = U_up[N_f-1]*(T_p[N_f-1]-T_p[N_f]) + h_c[N_f] * D_tank * D_tank * CN.pi * dz * (T_f[N_f] - T_s[N_f]) / d_p - U_wall*A_loss_wall_i*(T_p[N_f]-T_amb) - U_top*A_loss_top*(T_p[N_f]-T_amb);
//Heat loss calculations, different form than the equations above as they were in terms of rho*dh/dt not m*dh/dt
  Q_loss_top = U_top * CN.pi * D_tank * D_tank * 0.25 * (T_f[N_f] - T_amb);
  Q_loss_bot = U_bot * CN.pi * D_tank * D_tank * 0.25 * (T_f[1] - T_amb);
  for i in 1:N_f loop
    Q_loss_wall[i] = U_wall * CN.pi * D_tank * (T_f[i] - T_amb) * dz;
  end for;
  Q_loss_total = Q_loss_top + sum(Q_loss_wall) + Q_loss_bot;
//End heat loss calculations
//Analyics
  der(E_stored) = abs(m_flow) * (h_in - h_out) - Q_loss_total;
  Level = E_stored / E_max;
//Calculated Pumping losses
  for i in 1:N_f loop
    p_drop[i] = dz * (600 * (1 - eta) ^ 2 * mu_f[i] * abs(m_flow) / (eta ^ 3 * d_p ^ 2 * rho_f_avg * CN.pi * D_tank ^ 2) + 28 * (1 - eta) * m_flow ^ 2 / (eta ^ 3 * d_p * rho_f_avg * CN.pi * CN.pi * D_tank ^ 4));
  end for;
  p_drop_total = sum(p_drop);
  W_loss_pump = abs(m_flow) / rho_f_avg * p_drop_total / eff_pump;
  annotation(
    Documentation(revisions = "<html>
		<p>By Zebedee Kee on 03/12/2020</p>
		</html>", info = "<html>
		<p>This model contains the heat-transfer calculations of a thermocline packed bed storage tank with spherical filler geometry. This model does not contain any fluid connectors, for the CSP component with connectors, see Thermocline_Spheres_SingleTank. Variables fluid_top and fluid_bot provides the enthalpy-temperature relationship of the fluid material. Depending on whether m_flow is positive (discharging, fluid flowing upwards) or negative (charging, fluid flowing downwards), the charging/discharging equations are applied. In this iteration of the model, discharging and standby are lumped into one state.</p>
		</html>"));

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end HBS_6Layer;