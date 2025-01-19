within SolarTherm.Utilities.DNI_Models;
//DNI as a function of zenith angle of the sun. A. B. Meinel & M. P. Meinel (1976), Applied Solar Energy: An Introduction, Addison Wesley, USA, 41, 180.
function Meinel
  input Real zenith_deg "Zenith angle (deg)";
  output Real DNI "Direct normal irradiance W/m2";
protected
  Real zenith_rad = zenith_deg*Modelica.Constants.pi/180.0 "Zenith angle (rad)";
  Real AM = 1.0/Modelica.Math.cos(zenith_rad);
algorithm
  DNI := 1367.0*0.7^(AM^0.678);
end Meinel;