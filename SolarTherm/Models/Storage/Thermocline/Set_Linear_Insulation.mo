within SolarTherm.Models.Storage.Thermocline;

function Set_Linear_Insulation
  input Real U_loss_bot;
  input Real U_loss_top;
  input Real H_tank;
  input Real z_f[N_f];
  input Integer N_f;
  output Real U_loss_array[N_f];
algorithm
  for i in 1:N_f loop
    U_loss_array[i] := Modelica.Math.Vectors.interpolate({0.0,H_tank},{U_loss_bot,U_loss_top},z_f[i]);
  end for;
annotation(
    Diagram(coordinateSystem(preserveAspectRatio = false)));
end Set_Linear_Insulation;