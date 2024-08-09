within SolarTherm.Validation.Gan_HBS;

function Initialise_h_f_start_discharging
  input Real z_start_data[:];
  input Real T_data[:];
  input Real z_array[100];
  output Real h_array[100];
algorithm
  for i in 1:100 loop
    h_array[i] := SolarTherm.Validation.Gan_HBS.Materials.Air_533kPa.h_Tf( Modelica.Math.Vectors.interpolate(z_start_data,T_data,z_array[i]),0.0);
  end for;

annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Initialise_h_f_start_discharging;