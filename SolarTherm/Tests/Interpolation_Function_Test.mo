within SolarTherm.Tests;

model Interpolation_Function_Test
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  import CV = Modelica.SIunits.Conversions;
  extends Modelica.Icons.Example;
  
  SI.Temperature T;
  SI.ThermalConductivity k;
  SI.Temperature T_array[8] = {373.15,473.15,673.15,873.15,1073.15,1273.15,1473.15,1673.15};
  SI.ThermalConductivity k_array[8] = {4.8534,4.3514,3.6819,3.3472,3.1798,3.0125,3.0125,3.0125};
equation
  //T = 300.0 + 0.1*time;
  T = 900.0 + 600.0*sin(time/100);
  //k = SolarTherm.Utilities.Interpolation.Interpolate1D(T_array,k_array,T);
  k = Modelica.Math.Vectors.interpolate(T_array,k_array,T)
  

annotation(experiment(StopTime = 12000, StartTime = 0, Tolerance = 1e-3, Interval = 10.0),
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Interpolation_Function_Test;