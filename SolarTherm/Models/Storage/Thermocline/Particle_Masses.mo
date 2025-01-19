within SolarTherm.Models.Storage.Thermocline;
//This is for a thermocline with spherical filler.
function Particle_Masses "Input outer diameter, number of filler CVs, density, outputs array of particle CV masses"
  input SI.Length d_p "Outer diameter";
  input Integer N_p "Number of mesh elements";
  input SI.Density rho "Density";
  output SI.Mass m[N_p];
protected
  SI.Length dr "Radius step";
  SI.Length r[N_p] "Midpoint radii";
algorithm
  dr := d_p / (2.0 * N_p);
  for i in 1:N_p loop
    r[i] := (i - 0.5) * dr;
    m[i] := rho * 4.0 * CN.pi * (r[i] ^ 2 * dr + dr ^ 3 / 12.0);
  end for;
end Particle_Masses;