within SolarTherm.Models.Storage.Thermocline;
//This is for a thermocline tank with horizontal cylindrical filler.
function Cylinder_Radii "Input outer filler particle diameter and number of filler CVs, output array of midpoint radii"
  input SI.Length d_p "Outer diameter";
  input Integer N_p "Number of mesh elements";
  output SI.Length r[N_p] "Midpoint radii";
protected
  SI.Mass m[N_p];
  SI.Length dr "Radius step";
algorithm
  dr := d_p / (2.0 * N_p);
  for i in 1:N_p loop
    r[i] := (i - 0.5) * dr;
  end for;
end Cylinder_Radii;