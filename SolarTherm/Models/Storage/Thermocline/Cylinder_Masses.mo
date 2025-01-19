within SolarTherm.Models.Storage.Thermocline;
//This is for a thermocline tank with horizontal cylindrical filler.
function Cylinder_Masses "Input outer diameter, number of filler CVs, density, outputs array of cylinder CV masses"
  input SI.Length d_p "Outer diameter";
  input Integer N_p "Number of mesh elements";
  input SI.Density rho "Density";
  input SI.Length L "Cylinder Length";
  output SI.Mass m[N_p];
protected
  SI.Length dr "Radius step";
  SI.Length r[N_p] "Midpoint radii";
algorithm
  dr := d_p / (2.0 * N_p);
  for i in 1:N_p loop
    r[i] := (i - 0.5) * dr;
    m[i] := 2.0*CN.pi*rho*L*r[i]*dr;
  end for;
end Cylinder_Masses;