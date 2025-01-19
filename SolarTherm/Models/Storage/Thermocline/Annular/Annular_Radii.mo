within SolarTherm.Models.Storage.Thermocline.Annular;
function Annular_Radii
  import SI = Modelica.SIunits;
  import CN = Modelica.Constants;
  input SI.Length r_pipe "Fluid pipe radius";
  input SI.Length r_solid "Solid phase outer radius";
  input Integer N_p "Number of mesh elements";
  output SI.Length r[N_p] "Midpoint radii";
protected
  SI.Length dr "Radius step";
algorithm
  dr := (r_solid-r_pipe)/N_p;
  for j in 1:N_p loop
    r[j] := r_pipe + ((j - 0.5) * dr);
  end for;
end Annular_Radii;