within SolarTherm.Calculators;

model Interpolate_ThermalConductivity
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import MA = Modelica.Blocks.Math;

  //parameter SI.Temperature T_data[13] = {298.15,373.15,473.15,588.15,673.15,773.15,873.15,1073.15,1273.15,1473.15,1673.15,1873.15,2073.15};Al2O3
  //parameter SI.ThermalConductivity k_data[13] = {25.104, 22.8028, 19.6648, 15.4808, 10.6692, 14.644, 8.9956, 6.6944, 6.276, 5.8576, 5.4392, 5.8576, 7.1128}; Al2O3
  parameter SI.Temperature T_data[5] = {473.15,673.15,1073.15,1473.15,1873.15}; //SiO2
  parameter SI.ThermalConductivity k_data[5] = {1.046,1.2552,1.6736,2.092,2.5104}; //SiO2
  
  SI.Temperature T;
  SI.ThermalConductivity k;
  
equation
  T = 300.0+time*10.0;
  k = Modelica.Math.Vectors.interpolate(T_data,k_data,T);
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)),experiment(StopTime = 120, StartTime = 0, Tolerance = 1.0e-5, Interval = 1, maxStepSize = 1, initialStepSize = 1));
end Interpolate_ThermalConductivity;