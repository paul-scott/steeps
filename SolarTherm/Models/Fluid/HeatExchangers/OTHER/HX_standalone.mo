within SolarTherm.Models.Fluid.HeatExchangers;
model HX_standalone
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import MA = Modelica.Math;
  import SolarTherm.{Models,Media};
  import Modelica.Math.Vectors;
  replaceable package Medium1 = Media.Sodium.Sodium_pT "Medium props for Sodium";
  replaceable package Medium2 = Media.ChlorideSalt.ChlorideSalt_pT "Medium props for Molten Salt";
  
  //Design Parameters
  parameter SI.HeatFlowRate Q_d_des = 50e6 "Design Heat Flow Rate";
  parameter SI.Temperature T_Na1_des = 740 + 273.15 "Desing Sodium Hot Fluid Temperature";
  parameter SI.Temperature T_MS1_des = 500 + 273.15 "Desing Molten Salt Cold Fluid Temperature";
  parameter SI.Temperature T_MS2_des = 720 + 273.15 "Desing Molten Salt Hot Fluid Temperature";
  parameter SI.Pressure p_Na1_des = 101325 "Design Sodium Inlet Pressure";
  parameter SI.Pressure p_MS1_des = 101325 "Design Molten Salt Inlet Pressure";
  
  //Sweep Parameters
  //parameter SI.Length d_o[19] = {6.35e-3,9.53e-3,12.70e-3,15.88e-3,19.05e-3,22.23e-3,25.40e-3,28.58e-3,31.75e-3,34.93e-3,38.10e-3,41.28e-3,44.45e-3,47.63e-3,50.80e-3,53.98e-3,57.15e-3,60.33e-3,63.50e-3} "Outer Tube Diameter";
  //parameter SI.Length L[10]={2,4,6,8,10,12,14,16,18,20} "Tube Length";
  //parameter Integer N_p[2]={4,8} "Tube passes number";
  //parameter Integer layout[2]={1,2} "Tube layout";
  //parameter SI.Temperature T_Na2[11] = {630+273.15,640+273.15,650+273.15,660+273.15,670+273.15,680+273.15,690+273.15,700+273.15,710+273.15,720+273.15,730+273.15} "Sodium Cold Fluid Temperature";
  parameter SI.Length d_o[7]={28.58e-3,31.75e-3,34.93e-3,38.10e-3,41.28e-3,44.45e-3,50.80e-3} "Outer Tube Diameter";
  parameter SI.Length L[6]={4,6,8,10,12,14} "Tube Length";
  parameter Integer N_p[1]={4} "Tube passes number";
  parameter Integer layout[1]={2} "Tube layout";
  parameter SI.Temperature T_Na2_des[6] = {650+273.15,670+273.15,680+273.15,690+273.15,700+273.15,710+273.15} "Sodium Cold Fluid Temperature";

  //Auxiliary parameters
  parameter Integer num_dim=5;
  parameter Integer dim1=size(d_o,1);
  parameter Integer dim2=size(L,1);
  parameter Integer dim3=size(N_p,1);
  parameter Integer dim4=size(layout,1);
  parameter Integer dim5=size(T_Na2_des,1);
  parameter Integer dim_tot=dim1*dim2*dim3*dim4*dim5;
  parameter Real vec[dim_tot,num_dim](each fixed=false);
  parameter Integer iter(fixed=false);
  parameter Boolean optimize_and_run(fixed=false);
  parameter Real m_flow_min_frac=0.25;
  parameter Real m_flow_max_frac=1.10;
  
  //Design Parameters
  parameter SI.MassFlowRate m_flow_Na_des[dim_tot](each fixed=false) "Sodium mass flow rate";
  parameter SI.MassFlowRate m_flow_MS_des[dim_tot](each fixed=false) "Molten-Salt mass flow rate";
  parameter Real F_des[dim_tot](each fixed=false) "Temperature correction factor";  
  parameter SI.ThermalConductance UA_des[dim_tot](each fixed=false) "UA";
  parameter SI.Area A_tot[dim_tot](each fixed=false) "Exchange Area";
  parameter SI.CoefficientOfHeatTransfer U_des[dim_tot](each fixed=false) "Heat tranfer coefficient";
  parameter Integer N_t[dim_tot](each fixed=false) "Number of tubes";
  parameter SI.Pressure Dp_tube_des[dim_tot](each fixed=false) "Tube-side pressure drop";
  parameter SI.Pressure Dp_shell_des[dim_tot](each fixed=false) "Shell-side pressure drop";
  parameter SI.CoefficientOfHeatTransfer h_s_des[dim_tot](each fixed=false) "Shell-side Heat tranfer coefficient";
  parameter SI.CoefficientOfHeatTransfer h_t_des[dim_tot](each fixed=false) "Tube-side Heat tranfer coefficient";
  parameter SI.Length D_s[dim_tot](each fixed=false) "Shell Diameter";
  parameter SI.Velocity v_Na_des[dim_tot](each fixed=false) "Sodium velocity in tubes";
  parameter SI.Velocity v_max_MS_des[dim_tot](each fixed=false) "Molten Salt velocity in shell";
  parameter SI.Volume V_HX_des[dim_tot](each fixed=false) "Heat-Exchanger Total Volume";
  parameter SI.Mass m_HX_des[dim_tot](each fixed=false) "Heat-Exchanger Total Mass";
  parameter Real C_BEC[dim_tot](unit = fill("€",dim_tot),each fixed=false) "Bare cost @2018";
  parameter Real C_pump_des[dim_tot](unit = fill("€/year",dim_tot),each fixed=false) "Annual pumping cost";
  parameter Real TAC_des[dim_tot](unit = fill("€/year",dim_tot),each fixed=false) "Total Annualized Cost";
  parameter Real ex_eff_des[dim_tot](each fixed=false) "HX Exergetic Efficiency";
  parameter Real en_eff_des[dim_tot](each fixed=false) "HX Energetic Efficiency";
  
  //Optimal Parameter Values
  parameter Integer result(fixed=false);
  parameter Real TAC_min(unit = "€/year",fixed=false) "Minimum Total Annualized Cost";
  parameter SI.Area A_tot_opt(fixed=false) "Optimal Exchange Area";
  parameter SI.CoefficientOfHeatTransfer U_opt_des(fixed=false) "Optimal Heat tranfer coefficient";
  parameter Integer N_t_opt(fixed=false) "Optimal Number of tubes";
  parameter SI.Pressure Dp_tube_opt_des(fixed=false) "Optimal Tube-side pressure drop";
  parameter SI.Pressure Dp_shell_opt_des(fixed=false) "Optimal Shell-side pressure drop";
  parameter SI.CoefficientOfHeatTransfer h_s_opt_des(fixed=false) "Optimal Shell-side Heat tranfer coefficient";
  parameter SI.CoefficientOfHeatTransfer h_t_opt_des(fixed=false) "Optimal Tube-side Heat tranfer coefficient";
  parameter SI.Length D_s_opt(fixed=false) "Optimal Shell Diameter";
  parameter SI.Velocity v_Na_opt_des(fixed=false) "Optimal Sodium velocity in tubes";
  parameter SI.Velocity v_max_MS_opt_des(fixed=false) "Optimal Molten Salt velocity in shell";
  parameter SI.Volume V_HX_opt_des(fixed=false) "Optimal Heat-Exchanger Total Volume";
  parameter SI.Mass m_HX_opt_des(fixed=false) "Optimal Heat-Exchanger Total Mass";
  parameter Real C_BEC_opt(unit = "€",fixed=false) "Optimal Bare cost @2018";
  parameter Real C_pump_opt_des(unit = "€/year",fixed=false) "Optimal Annual pumping cost";
  parameter SI.Length d_o_opt(fixed=false) "Optimal Outer Tube Diameter";
  parameter SI.Length L_opt(fixed=false) "Optimal Tube Length";
  parameter Integer N_p_opt(fixed=false) "Optimal Tube passes number";
  parameter Integer layout_opt(fixed=false) "Optimal Tube Layout";
  parameter SI.Temperature T_Na2_opt_des(fixed=false) "Optimal outlet sodium temperature";
  parameter SI.MassFlowRate m_flow_Na_opt_des(fixed=false) "Optimal Sodium mass flow rate";
  parameter SI.MassFlowRate m_flow_MS_opt_des(fixed=false) "Optimal Molten-Salt mass flow rate";
  parameter Real F_opt_des(fixed=false) "Optimal Temperature correction factor";  
  parameter SI.ThermalConductance UA_opt_des(fixed=false) "Optimal UA";
  parameter Real ex_eff_opt_des(fixed=false) "Optimal HX Exergetic Efficiency";
  parameter Real en_eff_opt_des(fixed=false) "Optimal HX Energetic Efficiency";
  parameter SI.MassFlowRate m_flow_min_Na(fixed=false) "Minimum Sodium mass flow rate";
  parameter SI.MassFlowRate m_flow_max_Na(fixed=false) "Maximum Sodium mass flow rate";
  
  //Variables
  SI.MassFlowRate m_flow_Na "Sodium mass flow rate";
  SI.MassFlowRate m_flow_MS "Molten Salt mass flow rate";
  SI.Temperature T_Na1(start=740+273.15,nominal=740+273.15) "Sodium Hot Fluid Temperature";
  SI.Temperature T_MS1 "Molten Salt Cold Fluid Temperature";
  SI.Temperature T_MS2 "Molten Salt Hot Fluid Temperature";
  SI.Temperature T_Na2(start=690+273.15,nominal=690+273.15) "Sodium Cold Fluid Temperature";
  SI.Pressure p_Na1 "Sodium Inlet Pressure";
  SI.Pressure p_MS1 "Molten Salt Inlet Pressure";
  SI.CoefficientOfHeatTransfer U "Heat tranfer coefficient";
  SI.CoefficientOfHeatTransfer h_s "Shell-side Heat tranfer coefficient";
  SI.CoefficientOfHeatTransfer h_t "Tube-side Heat tranfer coefficient";
  SI.HeatFlowRate Q "Design Heat Flow Rate";
  Real F "Temperature correction factor";
  SI.TemperatureDifference DT1 "Sodium-Molten Salt temperature difference 1";
  SI.TemperatureDifference DT2 "Sodium-Molten Salt temperature difference 2";
  SI.TemperatureDifference LMTD "Logarithmic mean temperature difference";
  SI.Pressure Dp_tube "Tube-side pressure drop";
  SI.Pressure Dp_shell "Shell-side pressure drop";
  SI.Velocity v_Na "Sodium velocity in tubes";
  SI.Velocity v_max_MS "Molten Salt velocity in shell";
  
  //Fluid Properties
  SI.Temperature Tm_Na(start=((690+740)/2)+273.15,nominal=((690+740)/2)+273.15) "Mean Sodium Fluid Temperature";
  SI.Temperature Tm_MS "Mean Molten Salts Fluid Temperature";
  SI.ThermalConductivity k_Na "Sodium Conductivity @mean temperature";
  SI.ThermalConductivity k_MS "Molten Salts Conductivity @mean temperature";
  SI.Density rho_Na "Sodium density @mean temperature";
  SI.Density rho_MS "Molten Salts density @mean temperature";
  SI.DynamicViscosity mu_Na "Sodium dynamic viscosity @mean temperature";
  SI.DynamicViscosity mu_MS "Molten Salts  dynamic viscosity @mean temperature";
  SI.DynamicViscosity mu_Na_wall "Sodium dynamic viscosity @wall temperature";
  SI.DynamicViscosity mu_MS_wall "Molten salts dynamic viscosity @wall temperature";
  SI.SpecificHeatCapacity cp_Na "Sodium specific heat capacity @mean temperature";
  SI.SpecificHeatCapacity cp_MS "Molten Salts specific heat capacity @mean temperature";
  Medium1.ThermodynamicState state_mean_Na;
  Medium1.ThermodynamicState state_input_Na;
  Medium1.ThermodynamicState state_output_Na;
  Medium2.ThermodynamicState state_mean_MS;
  Medium2.ThermodynamicState state_wall_MS;
  Medium2.ThermodynamicState state_input_MS;
  Medium2.ThermodynamicState state_output_MS;
  
initial algorithm
optimize_and_run:=true;
if optimize_and_run==true then
  iter:=1;
  for ll in 1:dim5 loop
    for ww in 1:dim4 loop
      for ii in 1:dim3 loop
        for jj in 1:dim2 loop
          for kk in 1:dim1 loop
          vec[iter,1]:=d_o[kk];
          vec[iter,2]:=L[jj];
          vec[iter,3]:=N_p[ii];
          vec[iter,4]:=layout[ww];
          vec[iter,5]:=T_Na2_des[ll];
          (m_flow_Na_des[iter], m_flow_MS_des[iter], F_des[iter], UA_des[iter], N_t[iter], U_des[iter], A_tot[iter], Dp_tube_des[iter], Dp_shell_des[iter], TAC_des[iter], h_s_des[iter], h_t_des[iter], D_s[iter], v_Na_des[iter], v_max_MS_des[iter], V_HX_des[iter], m_HX_des[iter], C_BEC[iter], C_pump_des[iter], ex_eff_des[iter], en_eff_des[iter]):= Design_HX(Q_d=Q_d_des, T_Na1=T_Na1_des, T_MS1=T_MS1_des, T_MS2=T_MS2_des, d_o=d_o[kk], L=L[jj], N_p=N_p[ii], layout=layout[ww], T_Na2=T_Na2_des[ll], p_MS1 = p_MS1_des, p_Na1 = p_Na1_des, c_e = 0.13, r = 0.05, H_y = 4500);
          iter:=iter+1;
          end for;
        end for;
      end for;
    end for;
  end for; 
  TAC_min:=min(TAC_des);
  result:=Vectors.find(TAC_min,TAC_des);
  A_tot_opt:=A_tot[result];
  U_opt_des:=U_des[result];
  N_t_opt:=N_t[result];
  Dp_tube_opt_des:=Dp_tube_des[result];
  Dp_shell_opt_des:=Dp_shell_des[result];
  h_s_opt_des:=h_s_des[result];
  h_t_opt_des:=h_t_des[result];
  D_s_opt:=D_s[result];
  v_Na_opt_des:=v_Na_des[result];
  v_max_MS_opt_des:=v_max_MS_des[result];
  V_HX_opt_des:=V_HX_des[result];
  m_HX_opt_des:=m_HX_des[result];
  C_BEC_opt:=C_BEC[result];
  C_pump_opt_des:=C_pump_des[result];
  d_o_opt:=vec[result,1];
  L_opt:=vec[result,2];
  N_p_opt:=integer(vec[result,3]);
  layout_opt:=integer(vec[result,4]);
  T_Na2_opt_des:=vec[result,5];
  m_flow_Na_opt_des:=m_flow_Na_des[result];
  m_flow_MS_opt_des:=m_flow_MS_des[result];
  m_flow_min_Na:=m_flow_min_frac*m_flow_Na_opt_des;
  m_flow_max_Na:=m_flow_max_frac*m_flow_Na_opt_des;
  F_opt_des:=F_des[result];
  UA_opt_des:=UA_des[result];
  ex_eff_opt_des:=ex_eff_des[result];
  en_eff_opt_des:=en_eff_des[result];
  else 
  TAC_min:=1.4439e+04;
  A_tot_opt:=19.3982;
  U_opt_des:=1.1953e+03;
  N_t_opt:=24;
  Dp_tube_opt_des:=3.0651e+04;
  Dp_shell_opt_des:=277.3451;
  h_s_opt_des:=1.6392e+03;
  h_t_opt_des:=4.1181e+04;
  D_s_opt:=0.4520;
  v_Na_opt_des:=2.0690;
  v_max_MS_opt_des:=0.1565;
  V_HX_opt_des:=1.0056;
  m_HX_opt_des:=2.0715e+03;
  C_BEC_opt:=1.7387e+05;
  C_pump_opt_des:=487.8192;
  d_o_opt:=0.0508;
  L_opt:=6;
  N_p_opt:=4;
  layout_opt:=2;
  T_Na2_opt_des:=963.1500;
  m_flow_Na_opt_des:=15.9719;
  m_flow_MS_opt_des:=4.2386;
  F_opt_des:=0.5711;
  UA_opt_des:=2.3186e+04;
  ex_eff_opt_des:=0.9461;
  en_eff_opt_des:=0.9167;
  m_flow_min_Na:=m_flow_min_frac*m_flow_Na_opt_des;
  m_flow_max_Na:=m_flow_max_frac*m_flow_Na_opt_des;
  end if;
  
equation
  p_MS1=101325; //Substitute with inlet pressure
  p_Na1=101325; //Substitute with inlet pressure
  T_Na1=740+273.15; //Substitute with inlet temperature
  T_MS1=500+273.15; //Substitute with inlet temperature
  m_flow_Na=((m_flow_max_Na-m_flow_min_Na)*time/10)+m_flow_max_Na; //Substitute with inlet flow rate
  
  T_MS2=720+273.15; //Impose value
  
  //Molten Salt properties
  Tm_MS=(T_MS1+T_MS2)/2;
  state_mean_MS=Medium2.setState_pTX(p_MS1, Tm_MS);
  state_input_MS=Medium2.setState_pTX(p_MS1, T_MS1);
  state_output_MS=Medium2.setState_pTX(p_MS1, T_MS2);
  state_wall_MS=Medium2.setState_pTX(p_MS1, Tm_Na);
  rho_MS=Medium2.density(state_mean_MS);
  cp_MS=Medium2.specificHeatCapacityCp(state_mean_MS);
  mu_MS=Medium2.dynamicViscosity(state_mean_MS);
  k_MS=Medium2.thermalConductivity(state_mean_MS);
  mu_MS_wall=Medium2.dynamicViscosity(state_wall_MS);
  //Sodium properties
  Tm_Na=(T_Na1+T_Na2)/2;
  state_mean_Na=Medium1.setState_pTX(p_Na1, Tm_Na);
  state_input_Na=Medium1.setState_pTX(p_Na1, T_Na1);
  state_output_Na=Medium1.setState_pTX(p_Na1, T_Na2);
  rho_Na=Medium1.density(state_mean_Na);
  cp_Na=Medium1.specificHeatCapacityCp(state_mean_Na);
  mu_Na=Medium1.dynamicViscosity(state_mean_Na);
  mu_Na_wall=mu_Na;
  k_Na=Medium1.thermalConductivity(state_mean_Na);
  //Problem  
  DT1=T_Na1-T_MS2;
  DT2=T_Na2-T_MS1;
  LMTD=(DT1-DT2)/MA.log(DT1 / DT2);
  F=TempCorrFactor(T_Na1=T_Na1, T_Na2=T_Na2, T_MS1=T_MS1, T_MS2=T_MS2);
  (U, h_s, h_t)=HTCs(d_o=d_o_opt, N_p=N_p_opt, layout=layout_opt, N_t=N_t_opt, state_mean_Na=state_mean_Na, state_mean_MS=state_mean_MS, state_wall_MS=state_wall_MS, m_flow_Na=m_flow_Na, m_flow_MS=m_flow_MS);
  Q=m_flow_Na*cp_Na*(T_Na1-T_Na2);
  Q=m_flow_MS*cp_MS*(T_MS2-T_MS1);
  Q=U*A_tot_opt*F*LMTD;
  (Dp_tube, Dp_shell, v_Na, v_max_MS)=Dp_losses(d_o=d_o_opt, N_p=N_p_opt, layout=layout_opt, N_t=N_t_opt, L=L_opt, state_mean_Na=state_mean_Na, state_mean_MS=state_mean_MS, state_wall_MS=state_wall_MS, m_flow_Na=m_flow_Na, m_flow_MS=m_flow_MS);
    
end HX_standalone;
