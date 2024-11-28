within SolarTherm.Calculators;

model HX_Temperature_Profiles_v1 "Checks the temperature profiles within a heat exchanger"
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import MA = Modelica.Blocks.Math;

  replaceable package Medium_H = SolarTherm.Media.SolidParticles.IOE_Dehydroxylated_ph;
  replaceable package Medium_C = Modelica.Media.IdealGases.SingleGases.H2;

  
  parameter SI.Pressure p_H = 1.0e5 "Pressure of hot stream (Pa)";
  parameter SI.Pressure p_C = 1.0e5 "Pressure of cold stream (Pa)";
  
  parameter SI.Temperature T_H_in = 797.0 + 273.15 "Hot inlet temperature (K)";
  parameter SI.Temperature T_C_in = 236.0 + 273.15 "Cold inlet temperature (K)";
  
  parameter SI.MassFlowRate m_flow_H = 93.86 "Mass flow rate of hot stream (kg/s)";
  parameter SI.MassFlowRate m_flow_C = 6.3 "Mass flow rate of cold stream (kg/s)";
  
  parameter Real eff_HX = 0.80 "Heat exchanger effectiveness (-)";
  
  parameter SI.SpecificHeatCapacity cp_H = (Medium_H.specificEnthalpy_pTX(p_H,T_H_in,{1.0}) - Medium_H.specificEnthalpy_pTX(p_H,T_C_in,{1.0}))/(T_H_in-T_C_in)"Approximate c_p of the hot stream (J/kgK)";

  parameter SI.SpecificHeatCapacity cp_C = (Medium_C.specificEnthalpy_pTX(p_C,T_H_in,{1.0}) - Medium_C.specificEnthalpy_pTX(p_C,T_C_in,{1.0}))/(T_H_in-T_C_in)"Approximate c_p of the cold stream (J/kgK)";
  
  parameter SI.ThermalConductance C_H = m_flow_H*cp_H "Heat capacity rate of hot stream (W/K)";
  parameter SI.ThermalConductance C_C = m_flow_C*cp_C "Heat capacity rate of cold stream (W/K)";
  
  parameter SI.ThermalConductance C_max = max(C_H,C_C) "Heat capacity rate of the larger stream (W/K)";
  parameter SI.ThermalConductance C_min = min(C_H,C_C) "Heat capacity rate of the smaller stream (W/K)";
  parameter Real C_r = C_min/C_max "Heat capacity rate ratio of the heat exchanger (-)";
  
  parameter SI.HeatFlowRate Q_flow_max = C_min*(T_H_in-T_C_in) "Maximum possible heat flow rate across the HX (J/s)";
  
  parameter SI.HeatFlowRate Q_flow_HX = eff_HX*Q_flow_max "Expected heat flow rate across the HX (J/s)";
  
  //Medium_H.BaseProperties State_H[100] "Property state of the hot media";
  //Medium_C.BaseProperties State_C[100] "Property state of the cold media";
  
  parameter SI.HeatFlowRate Delta_Q = Q_flow_HX/100.0 "Incremental heat flow rate across 1/100th of the HX length (J/s)";
  SI.SpecificEnthalpy h_H[101];
  SI.SpecificEnthalpy h_C[101];
  
  SI.Temperature T_H[101];
  SI.Temperature T_C[101];
  
  SI.TemperatureDifference Delta_T[101];
  SI.TemperatureDifference Delta_T_avg;
  SI.TemperatureDifference Delta_T_min;
  
  parameter SI.SpecificEnthalpy h_H_in = Medium_H.specificEnthalpy_pTX(p_H,T_H_in,{1.0});
  parameter SI.SpecificEnthalpy h_C_in = Medium_C.specificEnthalpy_pTX(p_C,T_C_in,{1.0});
  
algorithm
  //h_H[100] := h_H_in;
  //h_C[1] := h_C_in;
  for i in 1:101 loop
    h_C[i] := h_C_in + (i-1)*Delta_Q/m_flow_C;
    h_H[i] := h_H_in - (101-i)*Delta_Q/m_flow_H;
    
    T_C[i] := Medium_C.temperature_phX(p_C,h_C[i],{1.0});
    T_H[i] := Medium_H.temperature_phX(p_H,h_H[i],{1.0});
    
    Delta_T[i] := T_H[i]-T_C[i];
  end for;
  
  Delta_T_avg := sum(Delta_T)/101.0;
  Delta_T_min := min(Delta_T);
  
  //T_H[100] := T_H_in;
  //T_C[1] := T_C_in;
  
  //for i in 1:99 loop
   // h_C[i+1] = h_C[i] + Delta_Q/m_flow_C;
    //h_H[100-i] = h_H[101-i] - Delta_Q/m_flow_H;
   // T_H[i+1] = Medium_H.temperature_phX(p_H,h_H[i+1],{1.0});
   // T_C[i+1] = Medium_C.temperature_phX(p_C,h_C[i+1],{1.0});
  //end for;
    

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),experiment(StopTime = 120, StartTime = 0, Tolerance = 1.0e-5, Interval = 1, maxStepSize = 1, initialStepSize = 1));
end HX_Temperature_Profiles_v1;