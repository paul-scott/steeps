within SolarTherm.Models.Storage.Thermocline;
//This is to non-dimensionalise a vertical tank axis
function Relative_Tank_Axes "Input tank height, and number of fluid CVs, output array of non-dimensional z_positions"
  input SI.Length H_tank;
  input Integer N_f;
  output SI.Length ZdH[N_f];
protected
  SI.Length dz "Height step";
algorithm
  dz := H_tank / N_f;
  for i in 1:N_f loop
    ZdH[i] := (i - 0.5) * dz / H_tank;
  end for;
end Relative_Tank_Axes;