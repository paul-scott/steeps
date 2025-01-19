within SolarTherm.Calculators;

model HX_Temperature_Profiles_v2 "Checks the temperature profiles within a heat exchanger, the hot stream contains a mixture of (A = H2) and (B = H2O)"
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import MA = Modelica.Blocks.Math;

  replaceable package Medium_H_A = Modelica.Media.IdealGases.SingleGases.H2;
  replaceable package Medium_H_B = Modelica.Media.Water.WaterIF97_ph;
  replaceable package Medium_C = Modelica.Media.IdealGases.SingleGases.H2;

  parameter SI.Pressure p_H_A = 0.7657e5 "Pressure of hot stream component A (Pa)";
  parameter SI.Pressure p_H_B = 0.2343e5 "Pressure of hot stream component B (Pa)";
  parameter SI.Pressure p_C = 1.0e5 "Pressure of cold stream (Pa)";
  
  parameter SI.Temperature T_H_in = 600.0 + 273.15 "Hot inlet temperature (K)";
  parameter SI.Temperature T_C_in = 54.0 + 273.15 "Cold inlet temperature (K)";
  
  parameter SI.MassFlowRate m_flow_H_A = 1.87 "Mass flow rate of hot stream component A (kg/s)";
  parameter SI.MassFlowRate m_flow_H_B = 5.10 "Mass flow rate of hot stream component B (kg/s)";
  parameter SI.MassFlowRate m_flow_C = 6.3 "Mass flow rate of cold stream (kg/s)";
  
  parameter Real eff_HX = 0.95 "Heat exchanger effectiveness (-)";
  
  parameter SI.SpecificHeatCapacity cp_H_A = (Medium_H_A.specificEnthalpy_pTX(p_H_A,T_H_in,{1.0}) - Medium_H_A.specificEnthalpy_pTX(p_H_A,T_C_in,{1.0}))/(T_H_in-T_C_in)"Approximate c_p of the hot stream component A (J/kgK)";
  parameter SI.SpecificHeatCapacity cp_H_B = (Medium_H_B.specificEnthalpy_pTX(p_H_B,T_H_in,{1.0}) - Medium_H_B.specificEnthalpy_pTX(p_H_B,max(T_C_in,Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p_H_B)+1.0),{1.0}))/(T_H_in-T_C_in)"Approximate c_p of the hot stream component B (J/kgK)";
  
  parameter SI.Temperature T_sat = max(T_C_in,Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.tsat(p_H_B));

  parameter SI.SpecificHeatCapacity cp_C = (Medium_C.specificEnthalpy_pTX(p_C,T_H_in,{1.0}) - Medium_C.specificEnthalpy_pTX(p_C,T_C_in,{1.0}))/(T_H_in-T_C_in)"Approximate c_p of the cold stream (J/kgK)";
  
  parameter SI.ThermalConductance C_H = m_flow_H_A*cp_H_A + m_flow_H_B*cp_H_B "Heat capacity rate of hot stream (W/K)";
  parameter SI.ThermalConductance C_C = m_flow_C*cp_C "Heat capacity rate of cold stream (W/K)";
  
  parameter SI.ThermalConductance C_max = max(C_H,C_C) "Heat capacity rate of the larger stream (W/K)";
  parameter SI.ThermalConductance C_min = min(C_H,C_C) "Heat capacity rate of the smaller stream (W/K)";
  parameter Real C_r = C_min/C_max "Heat capacity rate ratio of the heat exchanger (-)";
  
  parameter SI.HeatFlowRate Q_flow_max = C_min*(T_H_in-T_C_in) "Maximum possible heat flow rate across the HX (J/s)";
  
  parameter SI.HeatFlowRate Q_flow_HX = eff_HX*Q_flow_max "Expected heat flow rate across the HX (J/s)";
  
  //Medium_H.BaseProperties State_H[100] "Property state of the hot media";
  //Medium_C.BaseProperties State_C[100] "Property state of the cold media";
  
  parameter SI.HeatFlowRate Delta_Q = Q_flow_HX/100.0 "Incremental heat flow rate across 1/100th of the HX length (J/s)";
  SI.SpecificEnthalpy h_H_A[101];
  SI.SpecificEnthalpy h_H_B[101];
  SI.SpecificEnthalpy h_C[101];
  
  SI.Temperature T_H[101];
  SI.Temperature T_C[101];
  
  SI.TemperatureDifference Delta_T[101];
  SI.TemperatureDifference Delta_T_avg;
  SI.TemperatureDifference Delta_T_min;
  
  parameter SI.SpecificEnthalpy h_H_A_in = Medium_H_A.specificEnthalpy_pTX(p_H_A,T_H_in,{1.0});
  parameter SI.SpecificEnthalpy h_H_B_in = Medium_H_B.specificEnthalpy_pTX(p_H_B,T_H_in,{1.0});
  parameter SI.SpecificEnthalpy h_C_in = Medium_C.specificEnthalpy_pTX(p_C,T_C_in,{1.0});
  
  parameter SI.HeatFlowRate H_flow_H_tot[101] = linspace(H_flow_H_tot_in - Q_flow_HX, H_flow_H_tot_in, 101);// = m_flow_H_A*h_H_A_in + m_flow_H_B*h_H_B_in;
  
  parameter SI.HeatFlowRate H_flow_H_tot_in = m_flow_H_A*h_H_A_in + m_flow_H_B*h_H_B_in;
  
  //parameter SI.SpecificEnthalpy 
  
  //SI.HeatFlowRate  H_flow_H_tot[11];
//initial algorithm
  //for i in 1:11 loop
    //H_flow_H_tot[i] := H_flow_H_tot_in - (11-i)*Delta_Q;
  //end for;
  
//algorithm
  //h_H[100] := h_H_in;
  //h_C[1] := h_C_in;
  //for i in 1:11 loop
    //h_C[i] := h_C_in + (i-1)*Delta_Q/m_flow_C;
    //h_H[i] := h_H_in - (101-i)*Delta_Q/m_flow_H;
    //T_C[i] := Medium_C.temperature_phX(p_C,h_C[i],{1.0});
    //T_H[i] := Medium_H_A.temperature_phX(p_H_A,h_H_A[i],{1.0});
    
    //Delta_T[i] := T_H[i]-T_C[i];
  //end for;
  
  //Delta_T_avg := sum(Delta_T)/11.0;
  //Delta_T_min := min(Delta_T);
//initial equation
  //h_H_A[11] = Medium_H_A.specificEnthalpy_pTX(p_H_A,T_H_in,{1.0});
  //h_H_B[11] = Medium_H_B.specificEnthalpy_pTX(p_H_B,T_H_in,{1.0});


equation
  
  for i in 1:101 loop
    h_H_A[i] = Medium_H_A.specificEnthalpy_pTX(p_H_A,T_H[i],{1.0});
    h_H_B[i] = Medium_H_B.specificEnthalpy_pTX(p_H_B,T_H[i],{1.0});
   // T_H[i] = Medium_H_A.temperature_phX(p_H_A,h_H_A[i],{1.0});
   // T_H[i] = Medium_H_B.temperature_phX(p_H_B,h_H_B[i],{1.0});
    
    h_C[i] = h_C_in + (i-1)*Delta_Q/m_flow_C;
    //H_flow_H_tot[i] = H_H_in_tot - (11-i)*Delta_Q;
    
    m_flow_H_A*h_H_A[i] + m_flow_H_B*h_H_B[i] = H_flow_H_tot[i];
    T_C[i] = Medium_C.temperature_phX(p_C,h_C[i],{1.0});
    //T_H[i] = Medium_H_A.temperature_phX(p_H_A,h_H_A[i],{1.0});
    Delta_T[i] = T_H[i]-T_C[i];

  //T_H[100] := T_H_in;
  //T_C[1] := T_C_in;
  
  //for i in 1:99 loop
   // h_C[i+1] = h_C[i] + Delta_Q/m_flow_C;
    //h_H[100-i] = h_H[101-i] - Delta_Q/m_flow_H;
   // T_H[i+1] = Medium_H.temperature_phX(p_H,h_H[i+1],{1.0});
   // T_C[i+1] = Medium_C.temperature_phX(p_C,h_C[i+1],{1.0});
  end for;
  Delta_T_avg = sum(Delta_T)/101.0;
  Delta_T_min = min(Delta_T);

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),experiment(StopTime = 120, StartTime = 0, Tolerance = 1.0e-5, Interval = 1, maxStepSize = 1, initialStepSize = 1));
end HX_Temperature_Profiles_v2;