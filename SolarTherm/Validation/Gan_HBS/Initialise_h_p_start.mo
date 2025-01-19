within SolarTherm.Validation.Gan_HBS;

function Initialise_h_p_start
  input Real z_start_data[:];
  input Real T_data[:];
  input Real z_array[100];
  output Real h_array[100];
algorithm
  for i in 1:5 loop
    h_array[i] := SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_1.h_Tf( Modelica.Math.Vectors.interpolate(z_start_data,T_data,z_array[i]),0.0);
  end for;
  for i in 6:42 loop
    h_array[i] := SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_2.h_Tf( Modelica.Math.Vectors.interpolate(z_start_data,T_data,z_array[i]),0.0);
  end for;
  for i in 43:58 loop
    h_array[i] := SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_3.h_Tf( Modelica.Math.Vectors.interpolate(z_start_data,T_data,z_array[i]),0.0);
  end for;
  for i in 59:66 loop
    h_array[i] := SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_4.h_Tf( Modelica.Math.Vectors.interpolate(z_start_data,T_data,z_array[i]),0.0);
  end for;
  for i in 67:72 loop
    h_array[i] := SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_5.h_Tf( Modelica.Math.Vectors.interpolate(z_start_data,T_data,z_array[i]),0.0);
  end for;
  for i in 73:100 loop
    h_array[i] := SolarTherm.Validation.Gan_HBS.Materials.Checkerbrick_6.h_Tf( Modelica.Math.Vectors.interpolate(z_start_data,T_data,z_array[i]),0.0);
  end for;
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Initialise_h_p_start;